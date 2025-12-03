#include <cuda_runtime.h>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include "../mlipper_util.h"
#include "partial_likelihood.cuh"

__device__ inline void scale_double_pow2(double &x, int shift) {
    long long bits = __double_as_longlong(x);
    bits += ((long long)shift << 52);   // exponent += shift
    x = __longlong_as_double(bits);
}

__device__ __forceinline__ size_t per_node_span(const DeviceTree& D) {
    return (size_t)D.sites * (size_t)D.rate_cats * (size_t)D.states;
}

// ===== Common helpers =====
__device__ __forceinline__ unsigned int* site_scaler_ptr_base(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site,
    unsigned int rate_cats)
{
    return D.d_site_scaler
        ? D.d_site_scaler + (size_t)site * (size_t)rate_cats
        : nullptr;
}

// ===== Per-site computations =====
template<int RATE_CATS>
__device__ __forceinline__ void compute_tip_tip_site_ratecat(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site)
{
    // lookup-based helper (uses precomputed D.d_lookup)
    const size_t span     = (size_t)4 * RATE_CATS;
    const size_t per_node = per_node_span(D);

    const unsigned char* left_tip  = D.d_tipchars + (size_t)op.left_tip_index  * D.sites;
    const unsigned char* right_tip = D.d_tipchars + (size_t)op.right_tip_index * D.sites;

    const unsigned int j = (unsigned int)left_tip[site];
    const unsigned int k = (unsigned int)right_tip[site];

    const size_t lookup_idx = (((size_t)j << op.log2_stride) + k) * span;
    const size_t parent_off = (size_t)op.parent_id * per_node + (size_t)site * span;

    const double4* __restrict__ src_vec =
        reinterpret_cast<const double4*>(D.d_lookup + lookup_idx);
    double4* __restrict__ dst_vec =
        reinterpret_cast<double4*>(D.d_clv_pool + parent_off);

    #pragma unroll
    for (int r = 0; r < RATE_CATS; ++r) {
        dst_vec[r] = src_vec[r];
    }

    if (D.d_site_scaler) {
        unsigned int* site_scaler_ptr =
            site_scaler_ptr_base(D, op, site, RATE_CATS);

        #pragma unroll
        for (int r = 0; r < RATE_CATS; ++r) {
            double* pout = reinterpret_cast<double*>(dst_vec + r);

            double maxv = fmax(fmax(pout[0], pout[1]),
                               fmax(pout[2], pout[3]));
            int expv;
            frexp(maxv, &expv);
            if (expv < SCALE_THRESHOLD_EXPONENT) {
                unsigned int shift = SCALE_THRESHOLD_EXPONENT - expv;
                site_scaler_ptr[r] += shift;

                #pragma unroll
                for (int s = 0; s < 4; ++s) {
                    scale_double_pow2(pout[s], shift);
                }
            }
        }
    }
}

template<int RATE_CATS>
__device__ __forceinline__ void compute_tip_tip_site_ratecat_nolookup(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site)
{
    // on-the-fly helper (no lookup table)
    const size_t span     = (size_t)4 * RATE_CATS;
    const size_t per_node = per_node_span(D);

    const unsigned char* left_tip  = D.d_tipchars + (size_t)op.left_tip_index  * D.sites;
    const unsigned char* right_tip = D.d_tipchars + (size_t)op.right_tip_index * D.sites;

    const unsigned int j = (unsigned int)left_tip[site];
    const unsigned int k = (unsigned int)right_tip[site];

    const unsigned int jmask_base = D.d_tipmap[j];
    const unsigned int kmask_base = D.d_tipmap[k];

    const double* __restrict__ jmat_base =
        D.d_pmat + (size_t)op.left_id  * RATE_CATS * 4 * 4;
    const double* __restrict__ kmat_base =
        D.d_pmat + (size_t)op.right_id * RATE_CATS * 4 * 4;

    const size_t parent_off = (size_t)op.parent_id * per_node + (size_t)site * span;
    double* __restrict__ dst = D.d_clv_pool + parent_off;

    unsigned int* site_scaler_ptr =
        site_scaler_ptr_base(D, op, site, RATE_CATS);

    #pragma unroll
    for (int r = 0; r < RATE_CATS; ++r) {
        const double* __restrict__ jmat = jmat_base + (size_t)r * 4 * 4;
        const double* __restrict__ kmat = kmat_base + (size_t)r * 4 * 4;
        double* __restrict__ Pout = dst + (size_t)r * 4;

        double col_scale_max_val = 0.0;

        const double* Lrow = jmat;
        const double* Rrow = kmat;
        #pragma unroll
        for (int i = 0; i < 4; ++i) {
            double termj = 0.0;
            double termk = 0.0;
            unsigned int jmask = jmask_base;
            unsigned int kmask = kmask_base;

            #pragma unroll
            for (int m = 0; m < 4; ++m) {
                if (jmask & 1u) termj += Lrow[m];
                if (kmask & 1u) termk += Rrow[m];
                jmask >>= 1;
                kmask >>= 1;
            }

            Pout[i] = termj * termk;
            if (Pout[i] > col_scale_max_val) col_scale_max_val = Pout[i];

            Lrow += 4;
            Rrow += 4;
        }

        if (site_scaler_ptr) {
            double* pout = Pout;
            double maxv = fmax(fmax(pout[0], pout[1]),
                               fmax(pout[2], pout[3]));
            int expv;
            frexp(maxv, &expv);
            if (expv < SCALE_THRESHOLD_EXPONENT) {
                unsigned int shift = SCALE_THRESHOLD_EXPONENT - expv;
                site_scaler_ptr[r] += shift;

                #pragma unroll
                for (int s = 0; s < 4; ++s)
                    scale_double_pow2(pout[s], shift);
            }
        }
    }
}

__device__ __forceinline__ void compute_tip_tip_site_4_generic(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site)
{
    const size_t span     = (size_t)4 * (size_t)D.rate_cats;
    const size_t per_node = per_node_span(D);

    const unsigned char* left_tip  = D.d_tipchars + (size_t)op.left_tip_index  * D.sites;
    const unsigned char* right_tip = D.d_tipchars + (size_t)op.right_tip_index * D.sites;

    const unsigned int j = (unsigned int)left_tip[site];
    const unsigned int k = (unsigned int)right_tip[site];

    const unsigned int jmask_base = D.d_tipmap[j];
    const unsigned int kmask_base = D.d_tipmap[k];

    const double* __restrict__ jmat_base =
        D.d_pmat + (size_t)op.left_id  * D.rate_cats * 4 * 4;
    const double* __restrict__ kmat_base =
        D.d_pmat + (size_t)op.right_id * D.rate_cats * 4 * 4;

    const size_t parent_off = (size_t)op.parent_id * per_node + (size_t)site * span;
    double* __restrict__ dst = D.d_clv_pool + parent_off;

    unsigned int* site_scaler_ptr =
        site_scaler_ptr_base(D, op, site, (unsigned int)D.rate_cats);

    for (int r = 0; r < D.rate_cats; ++r) {
        const double* __restrict__ jmat = jmat_base + (size_t)r * 4 * 4;
        const double* __restrict__ kmat = kmat_base + (size_t)r * 4 * 4;
        double* __restrict__ Pout = dst + (size_t)r * 4;

        double col_scale_max_val = 0.0;

        const double* Lrow = jmat;
        const double* Rrow = kmat;
        for (int i = 0; i < 4; ++i) {
            double termj = 0.0;
            double termk = 0.0;
            unsigned int jmask = jmask_base;
            unsigned int kmask = kmask_base;

            for (int m = 0; m < 4; ++m) {
                if (jmask & 1u) termj += Lrow[m];
                if (kmask & 1u) termk += Rrow[m];
                jmask >>= 1;
                kmask >>= 1;
            }

            Pout[i] = termj * termk;
            if (Pout[i] > col_scale_max_val) col_scale_max_val = Pout[i];

            Lrow += 4;
            Rrow += 4;
        }

        if (site_scaler_ptr) {
            double maxv = fmax(fmax(Pout[0], Pout[1]),
                               fmax(Pout[2], Pout[3]));
            int expv;
            frexp(maxv, &expv);
            if (expv < SCALE_THRESHOLD_EXPONENT) {
                unsigned int shift = SCALE_THRESHOLD_EXPONENT - expv;
                site_scaler_ptr[r] += shift;
                for (int s = 0; s < 4; ++s) {
                    scale_double_pow2(Pout[s], shift);
                }
            }
        }
    }
}

__device__ __forceinline__ void compute_tip_tip_site_generic(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site)
{
    const unsigned int states = (unsigned int)D.states;
    const unsigned int rate_cats = (unsigned int)D.rate_cats;
    const size_t span = (size_t)states * rate_cats;
    const size_t per_node = per_node_span(D);

    const unsigned char* left_tip  = D.d_tipchars + (size_t)op.left_tip_index  * D.sites;
    const unsigned char* right_tip = D.d_tipchars + (size_t)op.right_tip_index * D.sites;

    const unsigned int j = (unsigned int)left_tip[site];
    const unsigned int k = (unsigned int)right_tip[site];

    const size_t src_off = (((size_t)j << op.log2_stride) + k) * span;
    const size_t dst_off = (size_t)op.parent_id * per_node + (size_t)site * span;

    for (size_t t = 0; t < span; ++t) {
        D.d_clv_pool[dst_off + t] = D.d_lookup[src_off + t];
    }

    if (D.d_site_scaler) {
        unsigned int* site_scaler_ptr =
            site_scaler_ptr_base(D, op, site, rate_cats);

        for (unsigned int r = 0; r < rate_cats; ++r) {
            double* pout = D.d_clv_pool + dst_off + (size_t)r * states;
            double maxv = 0.0;
            for (unsigned int s = 0; s < states; ++s) {
                double v = pout[s];
                if (v > maxv) maxv = v;
            }
            int expv;
            frexp(maxv, &expv);
            if (expv < SCALE_THRESHOLD_EXPONENT) {
                unsigned int shift = SCALE_THRESHOLD_EXPONENT - expv;
                site_scaler_ptr[r] += shift;
                for (unsigned int s = 0; s < states; ++s) {
                    scale_double_pow2(pout[s], shift);
                }
            }
        }
    }
}

template<int RATE_CATS>
__device__ __forceinline__ void compute_tip_inner_site_ratecat(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site)
{
    const size_t span     = (size_t)4 * RATE_CATS;
    const size_t per_node = per_node_span(D);

    const bool tip_on_left = op.left_tip_index >= 0;
    const int  tip_index   = tip_on_left ? op.left_tip_index  : op.right_tip_index;
    const int  inner_id    = tip_on_left ? op.right_id : op.left_id;
    const int  tip_node_id = tip_on_left ? op.left_id  : op.right_id;

    const unsigned char* d_left_tip = D.d_tipchars + (size_t)tip_index * D.sites;
    const double* d_right_clv = D.d_clv_pool + (size_t)inner_id * per_node;

    const double* d_Lmat = D.d_pmat + (size_t)tip_node_id * RATE_CATS * 4 * 4;
    const double* d_Rmat = D.d_pmat + (size_t)inner_id * RATE_CATS * 4 * 4;

    const size_t site_off = (size_t)site * span;
    const unsigned int tmask = D.d_tipmap[d_left_tip[site]];

    unsigned int* site_scaler_ptr =
        site_scaler_ptr_base(D, op, site, RATE_CATS);

    for (int r = 0; r < RATE_CATS; ++r) {
        const double* Lmat = d_Lmat + (size_t)r * 4 * 4;
        const double* Rmat = d_Rmat + (size_t)r * 4 * 4;
        const double* Rclv = d_right_clv + site_off + (size_t)r * 4;
        double* Pout = D.d_clv_pool + (size_t)op.parent_id * per_node + site_off + (size_t)r * 4;
        double col_scale_max_val = 0.0;

        const double r0 = Rclv[0];
        const double r1 = Rclv[1];
        const double r2 = Rclv[2];
        const double r3 = Rclv[3];

        const double* Lrow = Lmat;
        const double* Rrow = Rmat;
        for (int i = 0; i < 4; ++i) {
            double lefterm = 0.0;
            unsigned int lstate = tmask;
            if (lstate & 1u) lefterm += Lrow[0];
            if (lstate & 2u) lefterm += Lrow[1];
            if (lstate & 4u) lefterm += Lrow[2];
            if (lstate & 8u) lefterm += Lrow[3];

            double righterm = Rrow[0] * r0 + Rrow[1] * r1 + Rrow[2] * r2 + Rrow[3] * r3;
            Pout[i] = lefterm * righterm;
            if (Pout[i] > col_scale_max_val) col_scale_max_val = Pout[i];
            Lrow += 4;
            Rrow += 4;
        }

        if (site_scaler_ptr) {
            int scaling_exponent;
            frexp(col_scale_max_val, &scaling_exponent);
            if (scaling_exponent < SCALE_THRESHOLD_EXPONENT) {
                site_scaler_ptr[r] += SCALE_THRESHOLD_EXPONENT - scaling_exponent;
                #pragma unroll
                for (int i = 0; i < 4; ++i) {
                    scale_double_pow2(Pout[i], SCALE_THRESHOLD_EXPONENT - scaling_exponent);
                }
            }
        }
    }
}

__device__ __forceinline__ void compute_tip_inner_site_generic(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site)
{
    const unsigned int states = (unsigned int)D.states;
    const unsigned int rate_cats = (unsigned int)D.rate_cats;
    const size_t span     = (size_t)states * rate_cats;
    const size_t per_node = per_node_span(D);

    const bool tip_on_left = op.left_tip_index >= 0;
    const int  tip_index   = tip_on_left ? op.left_tip_index  : op.right_tip_index;
    const int  inner_id    = tip_on_left ? op.right_id : op.left_id;
    const int  tip_node_id = tip_on_left ? op.left_id  : op.right_id;

    const unsigned char* d_left_tip = D.d_tipchars + (size_t)tip_index * D.sites;
    const double* d_right_clv = D.d_clv_pool + (size_t)inner_id * per_node;

    const double* d_Lmat = D.d_pmat + (size_t)tip_node_id * D.rate_cats * states * states;
    const double* d_Rmat = D.d_pmat + (size_t)inner_id * D.rate_cats * states * states;

    const size_t site_off = (size_t)site * span;
    const unsigned int tmask = D.d_tipmap[d_left_tip[site]];
    unsigned int* site_scaler_ptr =
        site_scaler_ptr_base(D, op, site, rate_cats);

    for (unsigned int r = 0; r < rate_cats; ++r) {
        double col_scale_max_val = 0.0;
        int scaling_exponent;
        const double* Lmat = d_Lmat + (size_t)r * states * states;
        const double* Rmat = d_Rmat + (size_t)r * states * states;
        const double* Rclv = d_right_clv + site_off + (size_t)r * states;
        double* Pout = D.d_clv_pool + (size_t)op.parent_id * per_node + site_off + (size_t)r * states;

        const double* Lrow = Lmat;
        const double* Rrow = Rmat;
        for (unsigned int i = 0; i < states; ++i) {
            double lefterm = 0.0, righterm = 0.0;
            unsigned int lstate = tmask;
            for (unsigned int j = 0; j < states; ++j) {
                if (lstate & 1u) lefterm += Lrow[j];
                righterm += Rrow[j] * Rclv[j];
                lstate >>= 1;
            }
            Pout[i] = lefterm * righterm;
            if (Pout[i] > col_scale_max_val) col_scale_max_val = Pout[i];
            Lrow += states;
            Rrow += states;
        }
        if (site_scaler_ptr) {
            frexp(col_scale_max_val, &scaling_exponent);
            if (scaling_exponent < SCALE_THRESHOLD_EXPONENT) {
                site_scaler_ptr[r] += SCALE_THRESHOLD_EXPONENT - scaling_exponent;
                for (unsigned int i = 0; i < states; ++i) {
                    scale_double_pow2(Pout[i], SCALE_THRESHOLD_EXPONENT - scaling_exponent);
                }
            }
        }
    }
}

template<int RATE_CATS>
__device__ __forceinline__ void compute_inner_inner_site_ratecat(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site)
{
    const size_t span     = (size_t)RATE_CATS * 4;
    const size_t per_node = per_node_span(D);
    const size_t site_off = (size_t)site * span;

    const double* d_left_clv  = D.d_clv_pool + (size_t)op.left_id  * per_node;
    const double* d_right_clv = D.d_clv_pool + (size_t)op.right_id * per_node;
    const double* d_left_mat  = D.d_pmat + (size_t)op.left_id  * RATE_CATS * 4 * 4;
    const double* d_right_mat = D.d_pmat + (size_t)op.right_id * RATE_CATS * 4 * 4;

    unsigned int* site_scaler_ptr =
        site_scaler_ptr_base(D, op, site, RATE_CATS);

    for (int r = 0; r < RATE_CATS; ++r) {
        const double* Lclv = d_left_clv  + site_off + (size_t)r * 4;
        const double* Rclv = d_right_clv + site_off + (size_t)r * 4;

        const double* Lmat = d_left_mat  + (size_t)r * 4 * 4;
        const double* Rmat = d_right_mat + (size_t)r * 4 * 4;

        double* Pout = D.d_clv_pool + (size_t)op.parent_id * per_node + site_off + (size_t)r * 4;
        double col_scale_max_val = 0.0;

        const double l0 = Lclv[0];
        const double l1 = Lclv[1];
        const double l2 = Lclv[2];
        const double l3 = Lclv[3];
        const double r0 = Rclv[0];
        const double r1 = Rclv[1];
        const double r2 = Rclv[2];
        const double r3 = Rclv[3];

        const double lt0 = Lmat[0]*l0 + Lmat[1]*l1 + Lmat[2]*l2 + Lmat[3]*l3;
        const double lt1 = Lmat[4]*l0 + Lmat[5]*l1 + Lmat[6]*l2 + Lmat[7]*l3;
        const double lt2 = Lmat[8]*l0 + Lmat[9]*l1 + Lmat[10]*l2 + Lmat[11]*l3;
        const double lt3 = Lmat[12]*l0 + Lmat[13]*l1 + Lmat[14]*l2 + Lmat[15]*l3;

        const double rt0 = Rmat[0]*r0 + Rmat[1]*r1 + Rmat[2]*r2 + Rmat[3]*r3;
        const double rt1 = Rmat[4]*r0 + Rmat[5]*r1 + Rmat[6]*r2 + Rmat[7]*r3;
        const double rt2 = Rmat[8]*r0 + Rmat[9]*r1 + Rmat[10]*r2 + Rmat[11]*r3;
        const double rt3 = Rmat[12]*r0 + Rmat[13]*r1 + Rmat[14]*r2 + Rmat[15]*r3;

        Pout[0] = lt0 * rt0;
        Pout[1] = lt1 * rt1;
        Pout[2] = lt2 * rt2;
        Pout[3] = lt3 * rt3;

        col_scale_max_val = fmax(fmax(Pout[0], Pout[1]), fmax(Pout[2], Pout[3]));

        if (site_scaler_ptr) {
            int scaling_exponent;
            frexp(col_scale_max_val, &scaling_exponent);
            if (scaling_exponent < SCALE_THRESHOLD_EXPONENT) {
                site_scaler_ptr[r] += SCALE_THRESHOLD_EXPONENT - scaling_exponent;
                #pragma unroll
                for (int j = 0; j < 4; ++j) {
                    scale_double_pow2(Pout[j], SCALE_THRESHOLD_EXPONENT - scaling_exponent);
                }
            }
        }
    }
}

__device__ __forceinline__ void compute_inner_inner_site_generic(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site)
{
    const unsigned int states = (unsigned int)D.states;
    const unsigned int rate_cats = (unsigned int)D.rate_cats;
    const size_t span = (size_t)states * (size_t)rate_cats;
    const size_t per_node = per_node_span(D);
    const size_t site_off = (size_t)site * span;

    const double* d_left_clv  = D.d_clv_pool + (size_t)op.left_id  * per_node;
    const double* d_right_clv = D.d_clv_pool + (size_t)op.right_id * per_node;
    const double* d_left_mat  = D.d_pmat + (size_t)op.left_id  * D.rate_cats * states * states;
    const double* d_right_mat = D.d_pmat + (size_t)op.right_id * D.rate_cats * states * states;

    unsigned int* site_scaler_ptr =
        site_scaler_ptr_base(D, op, site, rate_cats);

    for (unsigned int r = 0; r < rate_cats; ++r) {
        const double* Lclv = d_left_clv  + site_off + (size_t)r * states;
        const double* Rclv = d_right_clv + site_off + (size_t)r * states;

        const double* Lmat = d_left_mat  + (size_t)r * states * states;
        const double* Rmat = d_right_mat + (size_t)r * states * states;

        double* Pout = D.d_clv_pool + (size_t)op.parent_id * per_node + site_off + (size_t)r * states;
        double col_scale_max_val = 0.0;

        const double* Lrow = Lmat;
        const double* Rrow = Rmat;
        for (unsigned int j = 0; j < states; ++j) {
            double lt = 0.0, rt = 0.0;
            #pragma unroll
            for (unsigned int k = 0; k < states; ++k) {
                lt += Lrow[k] * Lclv[k];
                rt += Rrow[k] * Rclv[k];
            }
            Pout[j] = lt * rt;
            if (Pout[j] > col_scale_max_val) col_scale_max_val = Pout[j];
            Lrow += states;
            Rrow += states;
        }

        if (site_scaler_ptr) {
            int scaling_exponent;
            frexp(col_scale_max_val, &scaling_exponent);
            if (scaling_exponent < SCALE_THRESHOLD_EXPONENT) {
                site_scaler_ptr[r] += SCALE_THRESHOLD_EXPONENT - scaling_exponent;
                for (unsigned int j = 0; j < states; ++j) {
                    scale_double_pow2(Pout[j], SCALE_THRESHOLD_EXPONENT - scaling_exponent);
                }
            }
        }
    }
}

// ===== Kernels =====
__global__ void UpdatePartialTipTipKernel(
    const DeviceTree D,
    const NodeOpInfo* __restrict__ ops,
    int num_ops)
{
    unsigned int tid  = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int step = blockDim.x * gridDim.x;

    for (unsigned int site = tid; site < D.sites; site += step) {
        for (int i = 0; i < num_ops; ++i) {
            compute_tip_tip_site_generic(D, ops[i], site);
        }
    }
}

template<int RATE_CATS>
__global__ void UpdatePartialTipTipKernel_states_4_ratecat(
    const DeviceTree D,
    const NodeOpInfo* __restrict__ ops,
    int num_ops)
{
    unsigned int tid  = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int step = blockDim.x * gridDim.x;

    for (unsigned int site = tid; site < D.sites; site += step) {
        for (int i = 0; i < num_ops; ++i) {
            compute_tip_tip_site_ratecat_nolookup<RATE_CATS>(D, ops[i], site);
        }
    }
}

__global__ void UpdatePartialTipTipKernel_states_4_generic(
    const DeviceTree D,
    const NodeOpInfo* __restrict__ ops,
    int num_ops)
{
    unsigned int tid  = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int step = blockDim.x * gridDim.x;

    for (unsigned int site = tid; site < D.sites; site += step) {
        for (int i = 0; i < num_ops; ++i) {
            compute_tip_tip_site_4_generic(D, ops[i], site);
        }
    }
}

__global__ void UpdatePartialTipInnerKernel(
    const DeviceTree D,
    const NodeOpInfo* __restrict__ ops,
    int num_ops)
{
    unsigned int tid  = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int step = blockDim.x * gridDim.x;

    for (unsigned int site = tid; site < D.sites; site += step) {
        for (int i = 0; i < num_ops; ++i) {
            compute_tip_inner_site_generic(D, ops[i], site);
        }
    }
}

template<int RATE_CATS>
__global__ void UpdatePartialTipInnerKernel_states_4_ratecat(
    const DeviceTree D,
    const NodeOpInfo* __restrict__ ops,
    int num_ops)
{
    unsigned int tid  = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int step = blockDim.x * gridDim.x;

    for (unsigned int site = tid; site < D.sites; site += step) {
        for (int i = 0; i < num_ops; ++i) {
            compute_tip_inner_site_ratecat<RATE_CATS>(D, ops[i], site);
        }
    }
}

__global__ void UpdatePartialInnerInnerKernel(
    const DeviceTree D,
    const NodeOpInfo* __restrict__ ops,
    int num_ops)
{
    unsigned int tid  = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int step = blockDim.x * gridDim.x;

    for (unsigned int site = tid; site < D.sites; site += step) {
        for (int i = 0; i < num_ops; ++i) {
            compute_inner_inner_site_generic(D, ops[i], site);
        }
    }
}

template<int RATE_CATS>
__global__ void  UpdatePartialInnerInnerKernel_states_4_ratecat(
    const DeviceTree D,
    const NodeOpInfo* __restrict__ ops,
    int num_ops)
{
    unsigned int tid  = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int step = blockDim.x * gridDim.x;

    for (unsigned int site = tid; site < D.sites; site += step) {
        for (int i = 0; i < num_ops; ++i) {
            compute_inner_inner_site_ratecat<RATE_CATS>(D, ops[i], site);
        }
    }
}

__global__ void create_lookup_kernel(
    const DeviceTree D,
    const NodeOpInfo* __restrict__ ops,
    unsigned int maxstates,
    unsigned int log2_maxstates)
{
    const unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned int total_pairs = maxstates * maxstates;
    if (tid >= total_pairs) return;

    const NodeOpInfo& op = ops[0];
    const unsigned int j = tid / maxstates;
    const unsigned int k = tid % maxstates;

    const size_t span = (size_t)D.states * (size_t)D.rate_cats;
    double* lh_statepair =
        D.d_lookup + (((size_t)j << log2_maxstates) + (size_t)k) * span;

    const double* jmat = D.d_pmat + (size_t)op.left_id  * D.rate_cats * D.states * D.states;
    const double* kmat = D.d_pmat + (size_t)op.right_id * D.rate_cats * D.states * D.states;

    unsigned int index = 0;

    for (unsigned int n = 0; n < (unsigned int)D.rate_cats; ++n) {
        for (unsigned int i = 0; i < (unsigned int)D.states; ++i) {
            double termj = 0.0;
            double termk = 0.0;

            unsigned int jstate = D.d_tipmap[j];
            unsigned int kstate = D.d_tipmap[k];

            for (unsigned int m = 0; m < (unsigned int)D.states; ++m) {
                if (jstate & 1u) termj += jmat[m];
                if (kstate & 1u) termk += kmat[m];
                jstate >>= 1;
                kstate >>= 1;
            }

            jmat += D.states;
            kmat += D.states;

            lh_statepair[index++] = termj * termk;
        }
    }
}

// ===== Host launchers =====
static NodeOpInfo make_node_op(
    int parent_id,
    int left_id,
    int right_id,
    int left_tip_index,
    int right_tip_index,
    NodeOpType type,
    int log2_stride = 0)
{
    NodeOpInfo op{};
    op.parent_id = parent_id;
    op.left_id = left_id;
    op.right_id = right_id;
    op.left_tip_index = left_tip_index;
    op.right_tip_index = right_tip_index;
    op.op_type = static_cast<int>(type);
    op.log2_stride = log2_stride;
    return op;
}

static NodeOpInfo* upload_single_op(const NodeOpInfo& op, cudaStream_t stream) {
    NodeOpInfo* d_op = nullptr;
    CUDA_CHECK(cudaMalloc(&d_op, sizeof(NodeOpInfo)));
    CUDA_CHECK(cudaMemcpyAsync(
        d_op,
        &op,
        sizeof(NodeOpInfo),
        cudaMemcpyHostToDevice,
        stream));
    return d_op;
}

void Launch_Update_Partial_InnerInner(
    const DeviceTree& D,
    int parent_id,
    int left_id,
    int right_id,
    cudaStream_t stream)
{
    validate_states_rate(D.states, D.rate_cats, 64, 8);

    NodeOpInfo op = make_node_op(
        parent_id,
        left_id,
        right_id,
        -1,
        -1,
        OP_INNER_INNER);
    NodeOpInfo* d_op = upload_single_op(op, stream);

    int block = 256;
    int grid = (D.sites + block - 1) / block;

    if (D.states == 4) {
        switch (D.rate_cats) {
            case 1:
                UpdatePartialInnerInnerKernel_states_4_ratecat<1><<<grid, block, 0, stream>>>(
                    D,
                    d_op,
                    1);
                break;
            case 4:
                UpdatePartialInnerInnerKernel_states_4_ratecat<4><<<grid, block, 0, stream>>>(
                    D,
                    d_op,
                    1);
                break;
            case 8:
                UpdatePartialInnerInnerKernel_states_4_ratecat<8><<<grid, block, 0, stream>>>(
                    D,
                    d_op,
                    1);
                break;
            default:
                UpdatePartialInnerInnerKernel<<<grid, block, 0, stream>>>(
                    D,
                    d_op,
                    1);
                break;
        }
    } else {
        UpdatePartialInnerInnerKernel<<<grid, block, 0, stream>>>(
            D,
            d_op,
            1);
    }
    CUDA_CHECK(cudaFree(d_op));
}

void partial_likelihood::compute_inner_inner(
    const DeviceTree& D,
    int parent_id,
    int left_id,
    int right_id,
    cudaStream_t stream)
{
    Launch_Update_Partial_InnerInner(
        D,
        parent_id,
        left_id,
        right_id,
        stream);
}

void Launch_Update_Partial_TipTip(
    const DeviceTree& D,
    int parent_id,
    int left_node_id,
    int right_node_id,
    int left_tip_index,
    int right_tip_index,
    cudaStream_t stream)
{
    validate_states_rate(D.states, D.rate_cats, 64, 8);

    const unsigned int tipmap_size = D.states + 1;
    const int log2_stride = ceil_log2_u32(tipmap_size);

    NodeOpInfo op = make_node_op(
        parent_id,
        left_node_id,
        right_node_id,
        left_tip_index,
        right_tip_index,
        OP_TIP_TIP,
        log2_stride);
    NodeOpInfo* d_op = upload_single_op(op, stream);

    int block = 4;
    int grid  = (tipmap_size * tipmap_size + block - 1) / block;
    create_lookup_kernel<<<grid, block, 0, stream>>>(
        D,
        d_op,
        tipmap_size,
        log2_stride);

    block = 256;
    grid  = (D.sites + block - 1) / block;
    if (D.states == 4) {
        switch (D.rate_cats) {
            case 1:
                UpdatePartialTipTipKernel_states_4_ratecat<1><<<grid, block, 0, stream>>>(
                    D, d_op, 1);
                break;
            case 4:
                UpdatePartialTipTipKernel_states_4_ratecat<4><<<grid, block, 0, stream>>>(
                    D, d_op, 1);
                break;
            case 8:
                UpdatePartialTipTipKernel_states_4_ratecat<8><<<grid, block, 0, stream>>>(
                    D, d_op, 1);
                break;
            default:
                UpdatePartialTipTipKernel_states_4_generic<<<grid, block, 0, stream>>>(
                    D,
                    d_op,
                    1);
                break;
        }
    } else {
        UpdatePartialTipTipKernel<<<grid, block, 0, stream>>>(
            D,
            d_op,
            1);
    }

    CUDA_CHECK(cudaFree(d_op));
}

void partial_likelihood::compute_tip_tip(
    const DeviceTree& D,
    int parent_id,
    int left_node_id,
    int right_node_id,
    int left_tip_index,
    int right_tip_index,
    cudaStream_t stream)
{
    Launch_Update_Partial_TipTip(
        D,
        parent_id,
        left_node_id,
        right_node_id,
        left_tip_index,
        right_tip_index,
        stream
    );
}

void Launch_Update_Partial_TipInner(
    const DeviceTree& D,
    int parent_id,
    int tip_node_id,
    int inner_node_id,
    int tip_index,
    cudaStream_t stream)
{
    validate_states_rate(D.states, D.rate_cats, 64, 8);

    const int block = 256;
    const int grid  = (D.sites + block - 1) / block;

    NodeOpInfo op = make_node_op(
        parent_id,
        tip_node_id,
        inner_node_id,
        tip_index,
        -1,
        OP_TIP_INNER);
    NodeOpInfo* d_op = upload_single_op(op, stream);

    if (D.states == 4) {
        switch (D.rate_cats) {
            case 1:
                UpdatePartialTipInnerKernel_states_4_ratecat<1><<<grid, block, 0, stream>>>(
                    D,
                    d_op,
                    1);
                break;
            case 4:
                UpdatePartialTipInnerKernel_states_4_ratecat<4><<<grid, block, 0, stream>>>(
                    D,
                    d_op,
                    1);
                break;
            case 8:
                UpdatePartialTipInnerKernel_states_4_ratecat<8><<<grid, block, 0, stream>>>(
                    D,
                    d_op,
                    1);
                break;
            default:
                UpdatePartialTipInnerKernel<<<grid, block, 0, stream>>>(
                    D,
                    d_op,
                    1);
                break;
        }
    } else {
        UpdatePartialTipInnerKernel<<<grid, block, 0, stream>>>(
            D,
            d_op,
            1);
    }
    CUDA_CHECK(cudaFree(d_op));
}

void partial_likelihood::compute_tip_inner(
    const DeviceTree& D,
    int parent_id,
    int tip_node_id,
    int inner_node_id,
    int tip_index,
    cudaStream_t stream)
{
    Launch_Update_Partial_TipInner(
        D,
        parent_id,
        tip_node_id,
        inner_node_id,
        tip_index,
        stream
    );
}

void partial_likelihood::compute_tip_inner_swap(
    const DeviceTree& D,
    int parent_id,
    int tip_node_id,
    int inner_node_id,
    int tip_index,
    cudaStream_t stream)
{
    Launch_Update_Partial_TipInner(
        D,
        parent_id,
        tip_node_id,
        inner_node_id,
        tip_index,
        stream
    );
}

__global__ void Rtree_Likelihood_Site_Parallel_Kernel(
    const DeviceTree D,
    const NodeOpInfo* ops,
    int num_ops
) {
    unsigned int tid  = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int step = blockDim.x * gridDim.x;

    for (unsigned int site = tid; site < D.sites; site += step) {
        for (int i = 0; i < num_ops; ++i) {
            const NodeOpInfo& op = ops[i];
            switch (op.op_type) {
                case OP_TIP_TIP:
                    if (D.states == 4) {
                        switch (D.rate_cats) {
                            case 1:
                                compute_tip_tip_site_ratecat_nolookup<1>(D, op, site);
                                break;
                            case 4:
                                compute_tip_tip_site_ratecat_nolookup<4>(D, op, site);
                                break;
                            case 8:
                                compute_tip_tip_site_ratecat_nolookup<8>(D, op, site);
                                break;
                            default:
                                compute_tip_tip_site_4_generic(D, op, site);
                                break;
                        }
                    } else {
                        compute_tip_tip_site_generic(D, op, site);
                    }
                    break;
                case OP_TIP_INNER:
                    if (D.states == 4) {
                        switch (D.rate_cats) {
                            case 1:
                                compute_tip_inner_site_ratecat<1>(D, op, site);
                                break;
                            case 4:
                                compute_tip_inner_site_ratecat<4>(D, op, site);
                                break;
                            case 8:
                                compute_tip_inner_site_ratecat<8>(D, op, site);
                                break;
                            default:
                                compute_tip_inner_site_generic(D, op, site);
                                break;
                        }
                    } else {
                        compute_tip_inner_site_generic(D, op, site);
                    }
                    break;
                case OP_INNER_INNER:
                    if (D.states == 4) {
                        switch (D.rate_cats) {
                            case 1:
                                compute_inner_inner_site_ratecat<1>(D, op, site);
                                break;
                            case 4:
                                compute_inner_inner_site_ratecat<4>(D, op, site);
                                break;
                            case 8:
                                compute_inner_inner_site_ratecat<8>(D, op, site);
                                break;
                            default:
                                compute_inner_inner_site_generic(D, op, site);
                                break;
                        }
                    } else {
                        compute_inner_inner_site_generic(D, op, site);
                    }
                    break;
                default:
                    break;
            }
        }
    }
}
