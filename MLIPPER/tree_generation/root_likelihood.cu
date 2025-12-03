#include "root_likelihood.cuh"

#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>
#include <cuda_runtime.h>
#include "../mlipper_util.h"

namespace root_likelihood {
constexpr int MAX_STATES = 64;
constexpr int MAX_RATECATS = 8;

__constant__ double c_frequencies[MAX_STATES];
__constant__ double c_rate_weights[MAX_RATECATS];

template<int RC, int SITES_PER_THREAD=2>
__global__ __launch_bounds__(256, 2)
void RootLikelihoodCalculation_states_4_1_4_8(
    std::size_t sites,
    const double* __restrict__ d_root_clv,
    const unsigned* __restrict__ d_pattern_w,
    const unsigned* __restrict__ d_site_scaler,
    const int* __restrict__ d_invar_indices,
    double invar_proportion,
    double* __restrict__ d_site_loglk_out)
{
    const int tid0 = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = gridDim.x * blockDim.x;

    double pi0 = c_frequencies[0];
    double pi1 = c_frequencies[1];
    double pi2 = c_frequencies[2];
    double pi3 = c_frequencies[3];

    double wr[RC];
    #pragma unroll
    for (int r=0; r<RC; ++r) wr[r] = c_rate_weights[r];

    for (int base = tid0 * SITES_PER_THREAD; base < (int)sites; base += stride * SITES_PER_THREAD) {
        #pragma unroll
        for (int t=0; t<SITES_PER_THREAD; ++t) {
            int i = base + t;
            if (i >= (int)sites) break;

            const double* __restrict__ clv_site = d_root_clv + (size_t)i * RC * 4;
            double sum_rate = 0.0;

            #pragma unroll
            for (int r=0; r<RC; ++r) {
                const double4* v = reinterpret_cast<const double4*>(clv_site + r*4);
                const double4 a = v[0];
                double val =
                    fma(a.x, pi0,
                    fma(a.y, pi1,
                    fma(a.z, pi2,
                        a.w * pi3)));
                if (d_site_scaler) {
                    unsigned int shift = d_site_scaler[(size_t)i * RC + r];
                    if (shift) val = ldexp(val, -static_cast<int>(shift));
                }
                sum_rate = fma(wr[r], val, sum_rate);
            }

            double site_sum = (1.0 - invar_proportion) * sum_rate;
            if (d_invar_indices) {
                int inv_idx = d_invar_indices[i];
                if (inv_idx >= 0) {
                    const double piv[4] = {pi0, pi1, pi2, pi3};
                    site_sum += invar_proportion * piv[inv_idx];
                }
            }
            const double eps = 1e-300;
            double loglk = log(site_sum > eps ? site_sum : eps);

            if (d_pattern_w) loglk *= static_cast<double>(d_pattern_w[i]);
            d_site_loglk_out[i] = loglk;
        }
    }
}

template<int RC, int SITES_PER_THREAD=2>
__global__ __launch_bounds__(256, 2)
void RootLikelihoodCalculation_states_5_1_4_8(
    std::size_t sites,
    const double* __restrict__ d_root_clv,
    const unsigned* __restrict__ d_pattern_w,
    const unsigned* __restrict__ d_site_scaler,
    const int* __restrict__ d_invar_indices,
    double invar_proportion,
    double* __restrict__ d_site_loglk_out)
{
    const int tid0   = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = gridDim.x  * blockDim.x;

    const double pi0 = c_frequencies[0];
    const double pi1 = c_frequencies[1];
    const double pi2 = c_frequencies[2];
    const double pi3 = c_frequencies[3];
    const double pi4 = c_frequencies[4];

    double wr[RC];
    #pragma unroll
    for (int r=0; r<RC; ++r) wr[r] = c_rate_weights[r];

    for (std::size_t base = std::size_t(tid0) * SITES_PER_THREAD;
         base < sites;
         base += std::size_t(stride) * SITES_PER_THREAD)
    {
        #pragma unroll
        for (int t=0; t<SITES_PER_THREAD; ++t) {
            const std::size_t i = base + t;
            if (i >= sites) break;

            const double* __restrict__ clv_site = d_root_clv + i * (std::size_t)RC * 5;
            double sum_rate = 0.0;

            #pragma unroll
            for (int r=0; r<RC; ++r) {
                const double* __restrict__ cr = clv_site + r*5;
                const double c0 = cr[0];
                const double c1 = cr[1];
                const double c2 = cr[2];
                const double c3 = cr[3];
                const double c4 = cr[4];
                double val = c0*pi0 + c1*pi1 + c2*pi2 + c3*pi3 + c4*pi4;
                if (d_site_scaler) {
                    unsigned int shift = d_site_scaler[i * (std::size_t)RC + r];
                    if (shift) val = ldexp(val, -static_cast<int>(shift));
                }
                sum_rate = fma(wr[r], val, sum_rate);
            }

            double site_sum = (1.0 - invar_proportion) * sum_rate;
            const int inv_idx = d_invar_indices ? d_invar_indices[i] : -1;
            if (inv_idx >= 0) site_sum += invar_proportion * c_frequencies[inv_idx];

            const double eps = 1e-300;
            double loglk = log(site_sum > eps ? site_sum : eps);

            if (d_pattern_w) loglk *= static_cast<double>(d_pattern_w[i]);
            d_site_loglk_out[i] = loglk;
        }
    }
}

template<int RC, int SITES_PER_THREAD=1>
__global__ __launch_bounds__(256, 2)
void RootLikelihoodCalculation_states_20_1_4_8(
    std::size_t sites,
    const double* __restrict__ d_root_clv,
    const unsigned* __restrict__ d_pattern_w,
    const unsigned* __restrict__ d_site_scaler,
    const int* __restrict__ d_invar_indices,
    double invar_proportion,
    double* __restrict__ d_site_loglk_out)
{
    const int tid0   = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = gridDim.x  * blockDim.x;

    const double4 pi0 = *reinterpret_cast<const double4*>(&c_frequencies[0]);
    const double4 pi1 = *reinterpret_cast<const double4*>(&c_frequencies[4]);
    const double4 pi2 = *reinterpret_cast<const double4*>(&c_frequencies[8]);
    const double4 pi3 = *reinterpret_cast<const double4*>(&c_frequencies[12]);
    const double4 pi4 = *reinterpret_cast<const double4*>(&c_frequencies[16]);

    double wr[RC];
    #pragma unroll
    for (int r=0; r<RC; ++r) wr[r] = c_rate_weights[r];

    for (int base = tid0 * SITES_PER_THREAD; base < (int)sites; base += stride * SITES_PER_THREAD) {
        #pragma unroll
        for (int t=0; t<SITES_PER_THREAD; ++t) {
            const int i = base + t;
            if (i >= (int)sites) break;

            const double* __restrict__ clv_site = d_root_clv + (size_t)i * RC * 20;
            double sum_rate = 0.0;

            #pragma unroll
            for (int r=0; r<RC; ++r) {
                const double* __restrict__ cr = clv_site + r*20;
                const double4 v0 = *reinterpret_cast<const double4*>(&cr[ 0]);
                const double4 v1 = *reinterpret_cast<const double4*>(&cr[ 4]);
                const double4 v2 = *reinterpret_cast<const double4*>(&cr[ 8]);
                const double4 v3 = *reinterpret_cast<const double4*>(&cr[12]);
                const double4 v4 = *reinterpret_cast<const double4*>(&cr[16]);
                double s =
                    fma(v0.x, pi0.x, fma(v0.y, pi0.y, fma(v0.z, pi0.z, v0.w * pi0.w))) +
                    fma(v1.x, pi1.x, fma(v1.y, pi1.y, fma(v1.z, pi1.z, v1.w * pi1.w))) +
                    fma(v2.x, pi2.x, fma(v2.y, pi2.y, fma(v2.z, pi2.z, v2.w * pi2.w))) +
                    fma(v3.x, pi3.x, fma(v3.y, pi3.y, fma(v3.z, pi3.z, v3.w * pi3.w))) +
                    fma(v4.x, pi4.x, fma(v4.y, pi4.y, fma(v4.z, pi4.z, v4.w * pi4.w)));
                if (d_site_scaler) {
                    unsigned int shift = d_site_scaler[(size_t)i * RC + r];
                    if (shift) s = ldexp(s, -static_cast<int>(shift));
                }
                sum_rate = fma(wr[r], s, sum_rate);
            }

            double site_sum = (1.0 - invar_proportion) * sum_rate;
            const int inv_idx = d_invar_indices ? d_invar_indices[i] : -1;
            if (inv_idx >= 0) site_sum += invar_proportion * c_frequencies[inv_idx];

            const double eps = 1e-300;
            double loglk = log(site_sum > eps ? site_sum : eps);

            if (d_pattern_w) loglk *= static_cast<double>(d_pattern_w[i]);
            d_site_loglk_out[i] = loglk;
        }
    }
}

template<int STATES>
// Generic kernel for arbitrary state/rate counts; one thread processes one site.
__global__ void RootLikelihoodCalculation_general(
    std::size_t sites,
    int states,
    int rate_cats,
    const double* __restrict__ d_root_clv,
    const unsigned* __restrict__ d_pattern_w,
    const unsigned* __restrict__ d_site_scaler,
    const int* __restrict__ d_invar_indices,
    double invar_proportion,
    double* __restrict__ d_site_loglk_out)
{
    const unsigned int tid0 = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned int stride = gridDim.x * blockDim.x;

    for (std::size_t site = tid0; site < sites; site += stride) {
        const double* clv_site = d_root_clv + site * rate_cats * (size_t)states;
        double sum_rate = 0.0;

        for (int r = 0; r < rate_cats; ++r) {
            const double* cr = clv_site + (size_t)r * states;
            double val = 0.0;
#pragma unroll
            for (int s = 0; s < STATES; ++s) {
                if (s >= states) break;
                val = fma(cr[s], c_frequencies[s], val);
            }
            if (d_site_scaler) {
                unsigned int shift = d_site_scaler[(size_t)site * rate_cats + r];
                if (shift) val = ldexp(val, -static_cast<int>(shift));
            }
            sum_rate = fma(c_rate_weights[r], val, sum_rate);
        }

        double site_sum = (1.0 - invar_proportion) * sum_rate;
        if (d_invar_indices) {
            int inv_idx = d_invar_indices[site];
            if (inv_idx >= 0) site_sum += invar_proportion * c_frequencies[inv_idx];
        }
        double safe = site_sum > 1e-300 ? site_sum : 1e-300;
        double loglk = log(safe);
        if (d_pattern_w) loglk *= static_cast<double>(d_pattern_w[site]);
        d_site_loglk_out[site] = loglk;
    }
}

// Host wrapper: copy model constants, launch specialized/general kernels, and sum log-likelihood.
double ComputeRootLogLikelihood(
    const double* d_root_clv,
    std::size_t sites,
    int states,
    int rate_cats,
    const double* h_frequencies,
    const double* h_rate_weights,
    const unsigned int* d_site_scaler,
    cudaStream_t stream,
    const unsigned* d_pattern_w,
    const int* d_invar_indices,
    double invar_proportion,
    RootLikelihoodWorkspace* workspace)
{
    if (!d_root_clv || !h_frequencies || !h_rate_weights) {
        throw std::runtime_error("Null pointer in root likelihood inputs.");
    }
    validate_states_rate(states, rate_cats, MAX_STATES, MAX_RATECATS);

    CHECK_CUDA(cudaMemcpyToSymbol(c_frequencies, h_frequencies, sizeof(double) * states));
    CHECK_CUDA(cudaMemcpyToSymbol(c_rate_weights, h_rate_weights, sizeof(double) * rate_cats));

    DeviceBuffer d_per_site_owner;
    double* d_per_site = nullptr;
    if (workspace) {
        workspace->ensure(sites);
        d_per_site = workspace->data();
    } else {
        d_per_site_owner = DeviceBuffer(sites);
        d_per_site = d_per_site_owner.get();
    }
    CHECK_CUDA(cudaMemsetAsync(d_per_site, 0, sizeof(double) * sites, stream));

    int block = 256;
    int grid = (sites + block - 1) / block;

    switch (states) {
        case 4:
            switch (rate_cats) {
                case 1:
                    RootLikelihoodCalculation_states_4_1_4_8<1,4><<<grid, block, 0, stream>>>(
                        sites, d_root_clv, d_pattern_w, d_site_scaler, d_invar_indices, invar_proportion, d_per_site);
                    break;
                case 4:
                    RootLikelihoodCalculation_states_4_1_4_8<4,4><<<grid, block, 0, stream>>>(
                        sites, d_root_clv, d_pattern_w, d_site_scaler, d_invar_indices, invar_proportion, d_per_site);
                    break;
                case 8:
                    RootLikelihoodCalculation_states_4_1_4_8<8,4><<<grid, block, 0, stream>>>(
                        sites, d_root_clv, d_pattern_w, d_site_scaler, d_invar_indices, invar_proportion, d_per_site);
                    break;
                default:
                    goto GENERAL;
            }
            break;
        case 5:
            switch (rate_cats) {
                case 1:
                    RootLikelihoodCalculation_states_5_1_4_8<1,4><<<grid, block, 0, stream>>>(
                        sites, d_root_clv, d_pattern_w, d_site_scaler, d_invar_indices, invar_proportion, d_per_site);
                    break;
                case 4:
                    RootLikelihoodCalculation_states_5_1_4_8<4,4><<<grid, block, 0, stream>>>(
                        sites, d_root_clv, d_pattern_w, d_site_scaler, d_invar_indices, invar_proportion, d_per_site);
                    break;
                case 8:
                    RootLikelihoodCalculation_states_5_1_4_8<8,4><<<grid, block, 0, stream>>>(
                        sites, d_root_clv, d_pattern_w, d_site_scaler, d_invar_indices, invar_proportion, d_per_site);
                    break;
                default:
                    goto GENERAL;
            }
            break;
        case 20:
            switch (rate_cats) {
                case 1:
                    RootLikelihoodCalculation_states_20_1_4_8<1,4><<<grid, block, 0, stream>>>(
                        sites, d_root_clv, d_pattern_w, d_site_scaler, d_invar_indices, invar_proportion, d_per_site);
                    break;
                case 4:
                    RootLikelihoodCalculation_states_20_1_4_8<4,4><<<grid, block, 0, stream>>>(
                        sites, d_root_clv, d_pattern_w, d_site_scaler, d_invar_indices, invar_proportion, d_per_site);
                    break;
                case 8:
                    RootLikelihoodCalculation_states_20_1_4_8<8,4><<<grid, block, 0, stream>>>(
                        sites, d_root_clv, d_pattern_w, d_site_scaler, d_invar_indices, invar_proportion, d_per_site);
                    break;
                default:
                    goto GENERAL;
            }
            break;
        default:
        GENERAL:
            RootLikelihoodCalculation_general<MAX_STATES><<<grid, block, 0, stream>>>(
                sites, states, rate_cats, d_root_clv, d_pattern_w, d_site_scaler, d_invar_indices, invar_proportion, d_per_site);
            break;
    }

    CHECK_CUDA(cudaGetLastError());

    std::vector<double> host_lk(sites);
    CHECK_CUDA(cudaMemcpyAsync(host_lk.data(), d_per_site, sizeof(double) * sites,
                               cudaMemcpyDeviceToHost, stream));
    CHECK_CUDA(cudaStreamSynchronize(stream));

    double total = 0.0;
    for (double v : host_lk) total += v;

    return total;
}

} // namespace root_likelihood
