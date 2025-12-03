// core_likelihood.cu
#include "core_likelihood.cuh" 
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <cuda_runtime.h>
#include <cub/cub.cuh>
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/discard_iterator.h>


namespace core_likelihood {
    __constant__ double c_frequencies[MAX_STATES];
    __constant__ double c_rate_weights[MAX_RATECATS];
    __constant__ double c_scale_minlh[SCALE_MAX_DIFF];
    __constant__ unsigned int c_tipmap[256];

}


void device_reduce_sum(const double* d_per_site, size_t sites, double* d_out_sum) {
    
    cudaError_t err;
    void*  d_temp_storage = nullptr;
    size_t temp_bytes     = 0;

    err = cub::DeviceReduce::Sum(d_temp_storage, temp_bytes, d_per_site, d_out_sum, sites);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "[Gpu_ERROR] DeviceReduce::Sum (query temp_bytes) failed: %s\n",cudaGetErrorString(err));
        exit(1);
    }

    err = cudaMalloc((void**)&d_temp_storage, temp_bytes);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "[Gpu_ERROR] cudaMalloc for temporary storage failed (requested %zu bytes): %s\n", temp_bytes, cudaGetErrorString(err));
        exit(1);
    }

    err = cub::DeviceReduce::Sum(d_temp_storage, temp_bytes, d_per_site, d_out_sum, sites);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "[Gpu_ERROR] DeviceReduce::Sum (execution) failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    cudaDeviceSynchronize();

    err = cudaFree(d_temp_storage);
    if (err != cudaSuccess)
    {
        fprintf(stderr,
                "[Gpu_ERROR] cudaFree (temporary storage) failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }
}

__device__ __forceinline__ double dot_states4_vec(const double* __restrict__ clv_r,
                                                  const double* __restrict__ freqs) {
    const double4* v = reinterpret_cast<const double4*>(clv_r);
    const double4 a = v[0];
    return fma(a.x, freqs[0],
           fma(a.y, freqs[1],
           fma(a.z, freqs[2],
               a.w * freqs[3])));
}

template<int STATES>
__device__ __forceinline__ double dot_states_scalar(const double* __restrict__ clv_r,
                                                    const double* __restrict__ freqs) {
    double s = 0.0;
#pragma unroll
    for (int j = 0; j < STATES; ++j) s = fma(clv_r[j], freqs[j], s);
    return s;
}

template<int RC, int SITES_PER_THREAD=2>
__global__ __launch_bounds__(256, 2) void RootLikelihoodCalculation_states_4_1_4_8(
    std::size_t sites,
    const double* __restrict__ d_root_clv,      // [sites * RC * 4]
    const unsigned* __restrict__ d_pattern_w,   // [sites] or nullptr
    const unsigned* __restrict__ d_site_scaler, // [sites] or nullptr
    const int* __restrict__ d_invar_indices,    // [sites] or nullptr; -1 if not
    double invar_proportion,
    double* __restrict__ d_site_loglk_out)
{
    const int tid0 = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = gridDim.x * blockDim.x;

    double pi0 = core_likelihood::c_frequencies[0];
    double pi1 = core_likelihood::c_frequencies[1];
    double pi2 = core_likelihood::c_frequencies[2];
    double pi3 = core_likelihood::c_frequencies[3];

    double wr[RC];
#pragma unroll
    for (int r=0; r<RC; ++r) wr[r] = core_likelihood::c_rate_weights[r];

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
                sum_rate = fma(wr[r], val, sum_rate);
            }

            double site_sum = (1.0 - invar_proportion) * sum_rate;
            if (d_invar_indices) {
                int inv_idx = d_invar_indices[i];
                if (inv_idx >= 0) {
                    const double piv[4] = {pi0,pi1,pi2,pi3};
                    site_sum += invar_proportion * piv[inv_idx];
                }
            }

            const double eps = 1e-300;
            double loglk = log(site_sum > eps ? site_sum : eps);

            if (d_pattern_w) loglk *= (double)d_pattern_w[i];
            if (d_site_scaler) loglk += (double)d_site_scaler[i] * LOG_SCALE_THRESHOLD;

            d_site_loglk_out[i] = loglk;
        }
    }
}

template<int RC, int SITES_PER_THREAD=2>
__global__ __launch_bounds__(256, 2)
void RootLikelihoodCalculation_states_5_1_4_8(
    std::size_t sites,
    const double* __restrict__ d_root_clv,      // [sites * RC * 5]
    const unsigned* __restrict__ d_pattern_w,   // [sites] or nullptr
    const unsigned* __restrict__ d_site_scaler, // [sites] or nullptr
    const int* __restrict__ d_invar_indices,    // [sites] or nullptr; -1 if not invar
    double invar_proportion,
    double* __restrict__ d_site_loglk_out)
{
    const int tid0   = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = gridDim.x  * blockDim.x;

    // 常數記憶體：用標量讀，避免 vector 對齊要求
    const double pi0 = core_likelihood::c_frequencies[0];
    const double pi1 = core_likelihood::c_frequencies[1];
    const double pi2 = core_likelihood::c_frequencies[2];
    const double pi3 = core_likelihood::c_frequencies[3];
    const double pi4 = core_likelihood::c_frequencies[4];

    double wr[RC];
    #pragma unroll
    for (int r=0; r<RC; ++r) wr[r] = core_likelihood::c_rate_weights[r];

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
                sum_rate   = fma(wr[r], val, sum_rate);
            }

            double site_sum = (1.0 - invar_proportion) * sum_rate;

            if (d_invar_indices) {
                const int inv_idx = d_invar_indices[i];
                if (inv_idx >= 0) site_sum += invar_proportion * core_likelihood::c_frequencies[inv_idx];
            }

            const double eps = 1e-300;
            double loglk = log(site_sum > eps ? site_sum : eps);

            if (d_pattern_w)  loglk *= (double)d_pattern_w[i];
            if (d_site_scaler) loglk += (double)d_site_scaler[i] * LOG_SCALE_THRESHOLD;

            d_site_loglk_out[i] = loglk;
        }
    }
}


template<int RC, int SITES_PER_THREAD=1>
__global__ void RootLikelihoodCalculation_states_20_1_4_8(
    std::size_t sites,
    const double* __restrict__ d_root_clv,      // [sites * RC * 20]
    const unsigned* __restrict__ d_pattern_w,   // [sites] or nullptr
    const unsigned* __restrict__ d_site_scaler, // [sites] or nullptr
    const int* __restrict__ d_invar_indices,    // [sites] or nullptr; -1 if not invar
    double invar_proportion,
    double* __restrict__ d_site_loglk_out)
{
    const int tid0   = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = gridDim.x  * blockDim.x;

    const double4 pi0 = *reinterpret_cast<const double4*>(&core_likelihood::c_frequencies[ 0]); // 0..3
    const double4 pi1 = *reinterpret_cast<const double4*>(&core_likelihood::c_frequencies[ 4]); // 4..7
    const double4 pi2 = *reinterpret_cast<const double4*>(&core_likelihood::c_frequencies[ 8]); // 8..11
    const double4 pi3 = *reinterpret_cast<const double4*>(&core_likelihood::c_frequencies[12]); // 12..15
    const double4 pi4 = *reinterpret_cast<const double4*>(&core_likelihood::c_frequencies[16]); // 16..19

    double wr[RC];
    #pragma unroll
    for (int r=0; r<RC; ++r) wr[r] = core_likelihood::c_rate_weights[r];

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

                const double4 v0 = *reinterpret_cast<const double4*>(&cr[ 0]); // 0..3
                const double4 v1 = *reinterpret_cast<const double4*>(&cr[ 4]); // 4..7
                const double4 v2 = *reinterpret_cast<const double4*>(&cr[ 8]); // 8..11
                const double4 v3 = *reinterpret_cast<const double4*>(&cr[12]); // 12..15
                const double4 v4 = *reinterpret_cast<const double4*>(&cr[16]); // 16..19

                double s =
                    fma(v0.x, pi0.x, fma(v0.y, pi0.y, fma(v0.z, pi0.z, v0.w * pi0.w))) +
                    fma(v1.x, pi1.x, fma(v1.y, pi1.y, fma(v1.z, pi1.z, v1.w * pi1.w))) +
                    fma(v2.x, pi2.x, fma(v2.y, pi2.y, fma(v2.z, pi2.z, v2.w * pi2.w))) +
                    fma(v3.x, pi3.x, fma(v3.y, pi3.y, fma(v3.z, pi3.z, v3.w * pi3.w))) +
                    fma(v4.x, pi4.x, fma(v4.y, pi4.y, fma(v4.z, pi4.z, v4.w * pi4.w)));

                sum_rate = fma(wr[r], s, sum_rate);
            }

            double site_sum = (1.0 - invar_proportion) * sum_rate;
            if (d_invar_indices) {
                const int inv_idx = d_invar_indices[i]; // 0..19 或 -1
                if (inv_idx >= 0) {
                    site_sum += invar_proportion * core_likelihood::c_frequencies[inv_idx];
                }
            }

            const double eps = 1e-300;
            double loglk = log(site_sum > eps ? site_sum : eps);

            if (d_pattern_w) loglk *= (double)d_pattern_w[i];
            if (d_site_scaler) loglk += (double)d_site_scaler[i] * LOG_SCALE_THRESHOLD;

            d_site_loglk_out[i] = loglk;
        }
    }
}



// template<int STATES, int RC, bool HAS_PW, bool HAS_SCALER, bool HAS_INVAR, bool VEC4=false>
// __global__ void RootLikelihoodCalculation_states_4(
//     std::size_t sites,
//     const double* __restrict__ d_root_clv,
//     const unsigned* __restrict__ d_pattern_weights,
//     const unsigned* __restrict__ d_site_scaler,
//     const int* __restrict__ d_invar_indices,
//     double invar_proportion,
//     double* __restrict__ d_site_loglk_out)
// {
//     const std::size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
//     const int stride = gridDim.x * blockDim.x;
    
//     if (tid >= sites) return;

//     const double* __restrict__ clv_site = d_root_clv + tid * (RC * STATES);

//     double sum_rate = 0.0;
// #pragma unroll
//     for (int r = 0; r < RC; ++r) {
//         const double* __restrict__ clv_r = clv_site + r * STATES;

//         double val = 0.0;
//         if constexpr (VEC4 && STATES == 4) {
//             val = dot_states4_vec(clv_r, core_likelihood::c_frequencies);
//         } else {
//             val = dot_states_scalar<STATES>(clv_r, core_likelihood::c_frequencies);
//         }
//         sum_rate = fma(core_likelihood::c_rate_weights[r], val, sum_rate);
//     }

//     double site_sum = (1.0 - invar_proportion) * sum_rate;
//     if constexpr (HAS_INVAR) {
//         if (d_invar_indices) {
//             const int inv_idx = d_invar_indices[tid];
//             if (inv_idx >= 0) {
//                 site_sum += invar_proportion * core_likelihood::c_frequencies[inv_idx];
//             }
//         }
//     }

//     const double eps = 1e-300;
//     double loglk = log(site_sum > eps ? site_sum : eps);

//     if constexpr (HAS_PW) {
//         if (d_pattern_weights) loglk *= static_cast<double>(d_pattern_weights[tid]);
//     }

//     if constexpr (HAS_SCALER) {
//         if (d_site_scaler) loglk += static_cast<double>(d_site_scaler[tid]) * LOG_SCALE_THRESHOLD;
//     }

//     d_site_loglk_out[tid] = loglk;
// }


__global__ void RootLikelihoodCalculation(
    std::size_t sites,
    int states,
    int rate_cats,
    const double* d_root_clv,
    const unsigned* d_pattern_weights,
    const unsigned* d_site_scaler,
    const int*  d_invar_indices,
    double invar_proportion,
    double* d_site_loglk_out)
{
    const unsigned lane   = threadIdx.x & 31; // 0..31
    const unsigned warpId = (blockIdx.x * blockDim.x + threadIdx.x) >> 5;
    if (warpId >= sites) return;

    const double* clv_site = d_root_clv + warpId * rate_cats * states;

    double sum_rate = 0.0;
    const unsigned fullmask = 0xFFFFFFFFu;
    const unsigned valid = __ballot_sync(fullmask, lane < states);

    for (unsigned r = 0; r < rate_cats; ++r) {
        double val = 0.0;
        if (lane < states) {
            const double* clv_r = clv_site + r * states;
            val = clv_r[lane] * core_likelihood::c_frequencies[lane];
        }
        for (int offset = 16; offset > 0; offset >>= 1) {
            double other = __shfl_down_sync(valid, val, offset);
            if (lane + offset < states) val += other;
        }

        if (lane == 0) {
            sum_rate += core_likelihood::c_rate_weights[r] * val;
        }
    }

    if (lane == 0) {
        double site_sum = (1.0 - invar_proportion) * sum_rate;
        if (d_invar_indices) {
            const int inv_idx = d_invar_indices[warpId];
            if (inv_idx >= 0) site_sum += invar_proportion * core_likelihood::c_frequencies[inv_idx];
        }

        const double eps = 1e-300;
        double loglk = log(site_sum > eps ? site_sum : eps);

        if (d_pattern_weights) loglk *= static_cast<double>(d_pattern_weights[warpId]);
        if (d_site_scaler) loglk += static_cast<double>(d_site_scaler[warpId]) * LOG_SCALE_THRESHOLD;
        d_site_loglk_out[warpId] = loglk;
    }
}

void core_likelihood::Likelihood_Root::Initialize(const Param& p)
{
    cudaError_t err;
    err = cudaMalloc((void**) &d_per_site, sizeof(double) * p.sites);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_per_site allocation failed: %s\n",cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMemset(d_per_site, 0, sizeof(double) * p.sites);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_per_site reset failed: %s\n",cudaGetErrorString(err));
        exit(1);
    }
   
    err = cudaMalloc((void**)&d_root_clv, sizeof(double) * p.sites * p.rate_cats * p.states);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_root_clv allocation failed: %s\n",cudaGetErrorString(err));
        exit(1);      
    }
    err = cudaMemset(d_root_clv, 0, sizeof(double) * p.sites * p.rate_cats  * p.states);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_root_clv reset failed: %s\n",
                cudaGetErrorString(err));
        exit(1);
    }
    
    err = cudaMalloc((void**)&d_pattern_weights, sizeof(unsigned) * p.sites);
        if (err != cudaSuccess){
            fprintf(stderr, "[Gpu_ERROR] d_pattern_weights allocation failed: %s\n",cudaGetErrorString(err));
            exit(1);
        }
        err = cudaMemset(d_pattern_weights, 0, sizeof(unsigned) * p.sites);
        if (err != cudaSuccess){
            fprintf(stderr, "[Gpu_ERROR] d_pattern_weights reset failed: %s\n",cudaGetErrorString(err));
            exit(1);
        }
    
    err = cudaMalloc((void**)&d_site_scaler, sizeof(unsigned) * p.sites);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_site_scaler allocation failed: %s\n",cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMemset(d_site_scaler, 0, sizeof(unsigned) * p.sites);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_site_scaler reset failed: %s\n",cudaGetErrorString(err));
        exit(1);
    }

    err = cudaMalloc((void**)&d_invar_indices, sizeof(int) * p.sites);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_invar_indices allocation failed: %s\n",cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMemset(d_invar_indices, 0, sizeof(int) * p.sites);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_invar_indices reset failed: %s\n",cudaGetErrorString(err));
        exit(1);
    }

    err = cudaMalloc((void**)&d_out_sum, sizeof(double));
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_out_sum allocation failed: %s\n",cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMemset(d_out_sum, 0, sizeof(double));
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_out_sum reset failed: %s\n",cudaGetErrorString(err));
        exit(1);
    }

}

void core_likelihood::Likelihood_Root::ConstructionOnGpu(
    const double* h_root_clv,
    const double* h_frequencies,
    const double* h_rate_weights,
    const unsigned int*    h_pattern_weights,
    const unsigned int*    h_site_scaler,
    const int*    h_invar_indices,
    const Param& p)

{
    cudaError_t err;

    err = cudaMemcpy(d_root_clv, h_root_clv, sizeof(double) * p.sites * p.rate_cats * p.states, cudaMemcpyHostToDevice);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_root_clv copy failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }

    err = cudaMemcpyToSymbol(core_likelihood::c_frequencies, h_frequencies, sizeof(double) * p.states);
    if (err != cudaSuccess) { std::fprintf(stderr, "cpy c_frequencies: %s\n", cudaGetErrorString(err)); std::abort(); }

    err = cudaMemcpyToSymbol(core_likelihood::c_rate_weights, h_rate_weights, sizeof(double) * p.rate_cats);
    if (err != cudaSuccess) { std::fprintf(stderr, "cpy c_rate_weights: %s\n", cudaGetErrorString(err)); std::abort(); }

    err = cudaMemcpy(d_pattern_weights, h_pattern_weights, sizeof(int) * p.sites, cudaMemcpyHostToDevice);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_pattern_weights copy failed: %s\n",cudaGetErrorString(err));
        exit(1);
    }
    if(h_site_scaler){
        err = cudaMemcpy(d_site_scaler, h_site_scaler, sizeof(unsigned) * p.sites, cudaMemcpyHostToDevice);
        if (err != cudaSuccess){
            fprintf(stderr, "[Gpu_ERROR] d_site_scaler copy failed: %s\n",cudaGetErrorString(err));
            exit(1);
        }
    }
    
    err = cudaMemcpy(d_invar_indices, h_invar_indices, sizeof(int) * p.sites, cudaMemcpyHostToDevice);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_invar_indices copy failed: %s\n",cudaGetErrorString(err));
        exit(1); 
    }


    cudaHostAlloc((void**)&total_loglik, sizeof(double), cudaHostAllocMapped);
    cudaHostGetDevicePointer((void**)&d_out_sum, total_loglik, 0);
}   


void core_likelihood::Likelihood_Root::ComputeLikelihood(
    const Param& p)
{
    cudaError_t err;
    cudaStream_t stream = 0;
    dim3 block(256), grid((p.sites + block.x*4 - 1)/(block.x*4)); 
        
    switch(p.states) {
        case 4:
            switch (p.rate_cats) {
                case 1:
                    RootLikelihoodCalculation_states_4_1_4_8<1,4><<<grid,block,0,stream>>>(
                                p.sites, 
                                d_root_clv, 
                                d_pattern_weights, 
                                d_site_scaler, 
                                d_invar_indices, 
                                p.invar_proportion, 
                                d_per_site);
                    break;
                case 4:
                    RootLikelihoodCalculation_states_4_1_4_8<4,4><<<grid,block,0,stream>>>(
                                p.sites, 
                                d_root_clv, 
                                d_pattern_weights, 
                                d_site_scaler, 
                                d_invar_indices, 
                                p.invar_proportion, 
                                d_per_site);
                    break;
                case 8:
                    RootLikelihoodCalculation_states_4_1_4_8<8,4><<<grid,block,0,stream>>>(
                                p.sites, 
                                d_root_clv, 
                                d_pattern_weights, 
                                d_site_scaler, 
                                d_invar_indices, 
                                p.invar_proportion, 
                                d_per_site);
                    break;
                default:
                    goto GENERAL;
            }
            break;
        case 5:
            switch (p.rate_cats) {
                case 1:
                    RootLikelihoodCalculation_states_5_1_4_8<1,4><<<grid, block, 0, stream>>>(
                                p.sites,
                                d_root_clv,
                                d_pattern_weights,
                                d_site_scaler,
                                d_invar_indices,
                                p.invar_proportion,
                                d_per_site);
                    break;
                case 4:
                    RootLikelihoodCalculation_states_5_1_4_8<4,4><<<grid, block, 0, stream>>>(
                            p.sites,
                            d_root_clv,
                            d_pattern_weights,
                            d_site_scaler,
                            d_invar_indices,
                            p.invar_proportion,
                            d_per_site);
                    break;
                case 8:
                    RootLikelihoodCalculation_states_5_1_4_8<8,4><<<grid, block, 0, stream>>>(
                            p.sites,
                            d_root_clv,
                            d_pattern_weights,
                            d_site_scaler,
                            d_invar_indices,
                            p.invar_proportion,
                            d_per_site);
                    break;
                default:
                    goto GENERAL; 
            }
            break;
        case 20:
            switch (p.rate_cats) {
                case 1:
                    RootLikelihoodCalculation_states_20_1_4_8<1,4><<<grid, block, 0, stream>>>(
                                p.sites,
                                d_root_clv,
                                d_pattern_weights,
                                d_site_scaler,
                                d_invar_indices,
                                p.invar_proportion,
                                d_per_site);
                    break;
                case 4:
                    RootLikelihoodCalculation_states_20_1_4_8<4,4><<<grid, block, 0, stream>>>(
                                p.sites,
                                d_root_clv,
                                d_pattern_weights,
                                d_site_scaler,
                                d_invar_indices,
                                p.invar_proportion,
                                d_per_site);
                    break;
                case 8:
                    RootLikelihoodCalculation_states_20_1_4_8<8,4><<<grid, block, 0, stream>>>(
                                p.sites,
                                d_root_clv,
                                d_pattern_weights,
                                d_site_scaler,
                                d_invar_indices,
                                p.invar_proportion,
                                d_per_site);
                    break; 
                default:
                    goto GENERAL;
            }
            break;
        default:
            GENERAL:
                RootLikelihoodCalculation<<<grid, block, 0, stream>>>(
                    p.sites,
                    p.states,
                    p.rate_cats,
                    d_root_clv,
                    d_pattern_weights,
                    d_site_scaler,
                    d_invar_indices,
                    p.invar_proportion,
                    d_per_site);
                break;
        }

    device_reduce_sum(d_per_site, p.sites, total_loglik);

    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "[Gpu_ERROR] CUDA_Kernel failed: %s\n",cudaGetErrorString(err));
        exit(1);
    }
    cudaDeviceSynchronize();
}

void core_likelihood::Likelihood_Root::PrintLikelihood() const
{
    printf("[Likelihood_Root] Total log-likelihood = %.12f \n", total_loglik);
}

void core_likelihood::Likelihood_Root::CleanUp()
{
    cudaFree(d_per_site);
    cudaFree(d_out_sum);
    cudaFree(d_root_clv);
    cudaFree(d_pattern_weights);
    cudaFree(d_site_scaler);
    cudaFree(d_invar_indices);
}


void core_likelihood::Likelihood_Tip_Inner::Initialize(const Param& p)
{
    cudaError_t err;

    err = cudaMalloc((void**) &d_per_site, sizeof(double) * p.sites);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_per_site allocation failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMemset(d_per_site, 0, sizeof(double) * p.sites);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_per_site reset failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }

    err = cudaMalloc((void**) &d_out_sum, sizeof(double));
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_out_sum allocation failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMemset(d_out_sum, 0, sizeof(double));
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_out_sum reset failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }

    err = cudaMalloc((void**) &d_parent_clv, sizeof(double) * p.sites * p.rate_cats * p.states);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_parent_clv allocation failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMemset(d_parent_clv, 0, sizeof(double) * p.sites * p.rate_cats * p.states);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_parent_clv reset failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }

    if(p.per_rate_scaling)  {
        err = cudaMalloc((void**) &d_parent_scaler, sizeof(unsigned int) * p.sites * p.rate_cats);
        if (err != cudaSuccess){
            fprintf(stderr, "[Gpu_ERROR] d_parent_scaler allocation failed: %s\n", cudaGetErrorString(err));
            exit(1);
        }
        err = cudaMemset(d_parent_scaler, 0, sizeof(unsigned int) * p.sites * p.rate_cats);
        if (err != cudaSuccess){
            fprintf(stderr, "[Gpu_ERROR] d_parent_scaler reset failed: %s\n", cudaGetErrorString(err));
            exit(1);
        }
    } else {
        err = cudaMalloc((void**) &d_parent_scaler, sizeof(unsigned int) * p.sites);
        if (err != cudaSuccess){
            fprintf(stderr, "[Gpu_ERROR] d_parent_scaler allocation failed: %s\n", cudaGetErrorString(err));
            exit(1);
        }
        err = cudaMemset(d_parent_scaler, 0, sizeof(unsigned int) * p.sites);
        if (err != cudaSuccess){
            fprintf(stderr, "[Gpu_ERROR] d_parent_scaler reset failed: %s\n", cudaGetErrorString(err));
            exit(1);
        }
    }
    

    err = cudaMalloc((void**) &d_tipchars, sizeof(unsigned char) * p.sites);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_tipchars allocation failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMemset(d_tipchars, 0, sizeof(unsigned char) * p.sites);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_tipchars reset failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }

    err = cudaMalloc((void**) &d_p_matrix, sizeof(double) * p.rate_cats * p.states * p.states);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_p_matrix allocation failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMemset(d_p_matrix, 0, sizeof(double) * p.rate_cats * p.states * p.states);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_p_matrix reset failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }

    err = cudaMalloc((void**) &d_pattern_weights, sizeof(unsigned) * p.sites);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_pattern_weights allocation failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMemset(d_pattern_weights, 0, sizeof(unsigned int) * p.sites);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_pattern_weights reset failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }

    err = cudaMalloc((void**) &d_invar_indices, sizeof(int) * p.sites);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_invar_indices allocation failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMemset(d_invar_indices, 0, sizeof(int) * p.sites);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_invar_indices reset failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }
}



void core_likelihood::Likelihood_Tip_Inner::ConstructionOnGpu(
    const double* h_parent_clv,
    const unsigned int* h_parent_scaler,
    const unsigned char* h_tipchars,
    const unsigned int* h_tipmap,
    const double* h_matrix,
    const double* h_frequencies,
    const double* h_rate_weights,
    const unsigned int* h_pattern_weights,
    const int*    h_invar_indices,
    const Param& p
    )
{
    cudaError_t err;

    err = cudaMemcpy(d_parent_clv, h_parent_clv, sizeof(double) * p.sites * p.rate_cats * p.states, cudaMemcpyHostToDevice);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_parent_clv copy failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }

    if(h_parent_scaler)
    {
        const size_t n = p.per_rate_scaling ? (p.sites * p.rate_cats) : p.sites;
        err = cudaMemcpy(d_parent_scaler, h_parent_scaler, sizeof(unsigned int) * n, cudaMemcpyHostToDevice);
        if (err != cudaSuccess){
            fprintf(stderr, "[Gpu_ERROR] d_parent_scaler copy failed: %s\n", cudaGetErrorString(err));
            exit(1);
        }
    }

    err = cudaMemcpy(d_tipchars, h_tipchars, sizeof(unsigned char) * p.sites, cudaMemcpyHostToDevice);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_tipchars copy failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }

    err = cudaMemcpyToSymbol(core_likelihood::c_tipmap, h_tipmap, sizeof(unsigned int) * 256);
    if (err != cudaSuccess) { std::fprintf(stderr, "cpy c_tipmap: %s\n", cudaGetErrorString(err)); exit(1); }

    err = cudaMemcpy(d_p_matrix, h_matrix, sizeof(double) * p.rate_cats * p.states * p.states, cudaMemcpyHostToDevice);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_p_matrix copy failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }

    err = cudaMemcpyToSymbol(core_likelihood::c_frequencies, h_frequencies, sizeof(double) * p.states);
    if (err != cudaSuccess) {
        std::fprintf(stderr, "cpy c_frequencies: %s\n", cudaGetErrorString(err)); std::abort(); 
    }

    err = cudaMemcpyToSymbol(core_likelihood::c_rate_weights, h_rate_weights, sizeof(double) * p.rate_cats);

    err = cudaMemcpy(d_pattern_weights, h_pattern_weights, sizeof(unsigned int) * p.sites, cudaMemcpyHostToDevice);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_pattern_weights copy failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMemcpy(d_invar_indices, h_invar_indices, sizeof(int) * p.sites, cudaMemcpyHostToDevice);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_invar_indices copy failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }

    std::vector<double> scale_minlh(SCALE_MAX_DIFF);
    double scale_factor = 1.0;
    
    if(p.per_rate_scaling){
        for (int i = 0; i < SCALE_MAX_DIFF; ++i) {
            scale_factor *= SCALE_THRESHOLD;
            scale_minlh[i] = scale_factor;
        }
        
        err = cudaMemcpyToSymbol(core_likelihood::c_scale_minlh, scale_minlh.data(), SCALE_MAX_DIFF * sizeof(double), cudaMemcpyHostToDevice);
        if (err != cudaSuccess){
            fprintf(stderr, "[Gpu_ERROR] c_scale_minlh copy failed: %s\n", cudaGetErrorString(err));
            exit(1);
        }
    
        err = cudaMalloc((void**)&d_site_scaling_min, sizeof(unsigned) * p.sites);
        if(err != cudaSuccess){
            fprintf(stderr, "[Gpu_ERROR] d_site_scaling_min allocation failed: %s\n", cudaGetErrorString(err));
            exit(1);
        }
        size_t N = p.sites * p.rate_cats;
        thrust::device_ptr<const unsigned> temp_reduce_data(d_parent_scaler);

        thrust::counting_iterator<unsigned> rate_scaling_temp_ptr(0);
        auto keys_begin = thrust::make_transform_iterator(rate_scaling_temp_ptr, KeyBySite{static_cast<unsigned>(p.rate_cats)});

        thrust::reduce_by_key(
            keys_begin, keys_begin + N,
            temp_reduce_data,
            thrust::make_discard_iterator(),
            thrust::device_pointer_cast(d_site_scaling_min),
            thrust::equal_to<unsigned>(),
            thrust::minimum<unsigned>());

        err = cudaGetLastError();
        if (err != cudaSuccess)
        {   
            fprintf(stderr, "[Gpu_ERROR] thrust::reduce_by_key failed: %s\n", cudaGetErrorString(err));
            exit(1);
        }
    }
    else{
        d_site_scaling_min = nullptr;
    }

}

template<int RC, int  SITES_PER_THREAD = 2>
__global__ __launch_bounds__(256, 2)
void TipInnerLikelihoodCalculation_states_4_1_4_8(
    std::size_t sites,
    const double* __restrict__ d_parent_clv,        // [sites * RC * 4]
    const unsigned* __restrict__ d_parent_scaler,   // [sites * RC] or nullptr
    const unsigned char* __restrict__ d_tipchars,   // [sites]
    const double* __restrict__ d_p_matrix,          // [RC * 4 * 4]
    const unsigned* __restrict__ d_pattern_weights, // [sites] or nullptr
    double invar_proportion,
    const int* __restrict__ d_invar_indices,        // [sites] or nullptr
    const unsigned* __restrict__ d_site_scaling_min,// [sites] or nullptr
    double* __restrict__ d_site_loglk_out,
    bool per_rate_scaling)
{
    const int tid0   = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = gridDim.x  * blockDim.x;

    const double pi0 = core_likelihood::c_frequencies[0];
    const double pi1 = core_likelihood::c_frequencies[1];
    const double pi2 = core_likelihood::c_frequencies[2];
    const double pi3 = core_likelihood::c_frequencies[3];

    double wr[RC];
    #pragma unroll
    for (int r = 0; r < RC; ++r) wr[r] = core_likelihood::c_rate_weights[r];

    for (std::size_t base = std::size_t(tid0) * SITES_PER_THREAD;
         base < sites;
         base += std::size_t(stride) * SITES_PER_THREAD)
    {
        #pragma unroll
        for (int t = 0; t < SITES_PER_THREAD; ++t) {
            const std::size_t i = base + t;
            if (i >= sites) break;

            const unsigned tip_state = static_cast<unsigned>(d_tipchars[i]);
            const unsigned mask      = core_likelihood::c_tipmap[tip_state];

            const double* __restrict__ clv_site = d_parent_clv + i * (std::size_t)RC * 4;

            double site_sum = 0.0;

            #pragma unroll
            for (int r = 0; r < RC; ++r) {
                const double* __restrict__ pmat_r = d_p_matrix + r * 16;
                const double* __restrict__ cr     = clv_site + r * 4;

                double terma_r = 0.0;

                // j = 0..3
                #pragma unroll
                for (int j = 0; j < 4; ++j) {
                    // sum over tip-allowed k
                    double termb = 0.0;

                    // k=0..3；用 mask 判斷 tip 可接受狀態
                    // (pmat_r row-major: j*4 + k)
                    if (mask & (1u << 0)) termb += pmat_r[j*4 + 0];
                    if (mask & (1u << 1)) termb += pmat_r[j*4 + 1];
                    if (mask & (1u << 2)) termb += pmat_r[j*4 + 2];
                    if (mask & (1u << 3)) termb += pmat_r[j*4 + 3];

                    // clv * pi * termb
                    const double pi_j = (j==0)?pi0 : (j==1)?pi1 : (j==2)?pi2 : pi3;
                    terma_r = fma(cr[j] * pi_j, termb, terma_r);
                }

                if (per_rate_scaling && d_parent_scaler && d_site_scaling_min) {
                    const int diff =
                        min(int(d_parent_scaler[i * RC + r]) - int(d_site_scaling_min[i]), SCALE_MAX_DIFF);
                    if (diff > 0) terma_r *= core_likelihood::c_scale_minlh[diff - 1];
                }

                if (invar_proportion > 0.0 && d_invar_indices) {
                    const int inv_idx = d_invar_indices[i]; // -1：非 invar；0..3：狀態
                    const double inv_lk = (inv_idx >= 0)
                        ? ((inv_idx==0)?pi0 : (inv_idx==1)?pi1 : (inv_idx==2)?pi2 : pi3)
                        : 0.0;
                    site_sum += wr[r] * ((1.0 - invar_proportion) * terma_r + invar_proportion * inv_lk);
                    continue;
                }
                site_sum += wr[r] * terma_r;
            }

            // loglk
            const double eps = 1e-300;
            double loglk = log(site_sum > eps ? site_sum : eps);

            if (d_site_scaling_min) loglk += static_cast<double>(d_site_scaling_min[i]) * LOG_SCALE_THRESHOLD;

            if (d_pattern_weights) loglk *= double(d_pattern_weights[i]);

            d_site_loglk_out[i] = loglk;
        }
    }
}

template<int RC, int SITES_PER_THREAD = 2>
__global__ __launch_bounds__(256, 2)
void TipInnerLikelihoodCalculation_states_5_1_4_8(
    std::size_t sites,
    const double* __restrict__ d_parent_clv,        // [sites * RC * 5]
    const unsigned* __restrict__ d_parent_scaler,   // [sites * RC] or nullptr (per-rate)
    const unsigned char* __restrict__ d_tipchars,   // [sites]
    const double* __restrict__ d_p_matrix,          // [RC * 5 * 5], row-major (j*5 + k)
    const unsigned* __restrict__ d_pattern_w,       // [sites] or nullptr
    double invar_proportion,
    const int* __restrict__ d_invar_indices,        // [sites] or nullptr; -1 if not invar
    const unsigned* __restrict__ d_site_scaling_min,// [sites] or nullptr (global min per-site scaler)
    double* __restrict__ d_site_loglk_out,
    bool per_rate_scaling                           // ★ 新增：是否啟用 per-rate scaling
){
    const int tid0   = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = gridDim.x  * blockDim.x;

    // 常數頻率：標量讀，與你 root 5-state 版本一致
    const double pi0 = core_likelihood::c_frequencies[0];
    const double pi1 = core_likelihood::c_frequencies[1];
    const double pi2 = core_likelihood::c_frequencies[2];
    const double pi3 = core_likelihood::c_frequencies[3];
    const double pi4 = core_likelihood::c_frequencies[4];

    // rate-weights 放暫存器
    double wr[RC];
    #pragma unroll
    for (int r = 0; r < RC; ++r) wr[r] = core_likelihood::c_rate_weights[r];

    for (std::size_t base = std::size_t(tid0) * SITES_PER_THREAD;
         base < sites;
         base += std::size_t(stride) * SITES_PER_THREAD)
    {
        #pragma unroll
        for (int t = 0; t < SITES_PER_THREAD; ++t) {
            const std::size_t i = base + t;
            if (i >= sites) break;

            const unsigned tip_state = static_cast<unsigned>(d_tipchars[i]);
            const unsigned mask      = core_likelihood::c_tipmap[tip_state]; // 5-state bitmask

            const double* __restrict__ clv_site = d_parent_clv + i * (std::size_t)RC * 5;

            double sum_rate = 0.0;

            #pragma unroll
            for (int r = 0; r < RC; ++r) {
                const double* __restrict__ cr     = clv_site + r * 5;      // CLV[j]
                const double* __restrict__ pmat_r = d_p_matrix + r * 25;   // 5*5

                // 計算：terma_r = Σ_j [ clv[j] * pi[j] * Σ_{k∈tipmask} P(j->k) ]
                double terma_r = 0.0;

                // j=0..4（標量版，避免任何 16B 對齊要求）
                {
                    // j=0
                    double termb = 0.0;
                    if (mask & (1u << 0)) termb += pmat_r[0*5 + 0];
                    if (mask & (1u << 1)) termb += pmat_r[0*5 + 1];
                    if (mask & (1u << 2)) termb += pmat_r[0*5 + 2];
                    if (mask & (1u << 3)) termb += pmat_r[0*5 + 3];
                    if (mask & (1u << 4)) termb += pmat_r[0*5 + 4];
                    terma_r = fma(cr[0] * pi0, termb, terma_r);
                }
                {
                    // j=1
                    double termb = 0.0;
                    if (mask & (1u << 0)) termb += pmat_r[1*5 + 0];
                    if (mask & (1u << 1)) termb += pmat_r[1*5 + 1];
                    if (mask & (1u << 2)) termb += pmat_r[1*5 + 2];
                    if (mask & (1u << 3)) termb += pmat_r[1*5 + 3];
                    if (mask & (1u << 4)) termb += pmat_r[1*5 + 4];
                    terma_r = fma(cr[1] * pi1, termb, terma_r);
                }
                {
                    // j=2
                    double termb = 0.0;
                    if (mask & (1u << 0)) termb += pmat_r[2*5 + 0];
                    if (mask & (1u << 1)) termb += pmat_r[2*5 + 1];
                    if (mask & (1u << 2)) termb += pmat_r[2*5 + 2];
                    if (mask & (1u << 3)) termb += pmat_r[2*5 + 3];
                    if (mask & (1u << 4)) termb += pmat_r[2*5 + 4];
                    terma_r = fma(cr[2] * pi2, termb, terma_r);
                }
                {
                    // j=3
                    double termb = 0.0;
                    if (mask & (1u << 0)) termb += pmat_r[3*5 + 0];
                    if (mask & (1u << 1)) termb += pmat_r[3*5 + 1];
                    if (mask & (1u << 2)) termb += pmat_r[3*5 + 2];
                    if (mask & (1u << 3)) termb += pmat_r[3*5 + 3];
                    if (mask & (1u << 4)) termb += pmat_r[3*5 + 4];
                    terma_r = fma(cr[3] * pi3, termb, terma_r);
                }
                {
                    // j=4
                    double termb = 0.0;
                    if (mask & (1u << 0)) termb += pmat_r[4*5 + 0];
                    if (mask & (1u << 1)) termb += pmat_r[4*5 + 1];
                    if (mask & (1u << 2)) termb += pmat_r[4*5 + 2];
                    if (mask & (1u << 3)) termb += pmat_r[4*5 + 3];
                    if (mask & (1u << 4)) termb += pmat_r[4*5 + 4];
                    terma_r = fma(cr[4] * pi4, termb, terma_r);
                }

                // per-rate scaling（由參數 per_rate_scaling 控制）
                if (per_rate_scaling && d_parent_scaler && d_site_scaling_min) {
                    const int diff =
                        min(int(d_parent_scaler[i * RC + r]) - int(d_site_scaling_min[i]), SCALE_MAX_DIFF);
                    if (diff > 0) terma_r *= core_likelihood::c_scale_minlh[diff - 1];
                }

                // 混入 invariant
                if (invar_proportion > 0.0 && d_invar_indices) {
                    const int inv_idx = d_invar_indices[i]; // -1 或 0..4
                    const double inv_lk = (inv_idx >= 0)
                                            ? ((inv_idx==0)?pi0:(inv_idx==1)?pi1:(inv_idx==2)?pi2:(inv_idx==3)?pi3:pi4)
                                            : 0.0;
                    sum_rate += wr[r] * ((1.0 - invar_proportion) * terma_r + invar_proportion * inv_lk);
                } else {
                    sum_rate += wr[r] * terma_r;
                }
            }

            // 與 root 版一致：site_sum -> loglk -> pattern weight / site scaler
            double site_sum = sum_rate; // 已包含 (1-invar)*... + invar*... 的情況
            const double eps = 1e-300;
            double loglk = log(site_sum > eps ? site_sum : eps);
            
            if (d_site_scaling_min) loglk += static_cast<double>(d_site_scaling_min[i]) * LOG_SCALE_THRESHOLD;
            if (d_pattern_w)  loglk *= double(d_pattern_w[i]);
            
            d_site_loglk_out[i] = loglk;
        }
    }
}



template<int RC, int SITES_PER_THREAD = 1>
__global__ __launch_bounds__(256, 2)
void TipInnerLikelihoodCalculation_states_20_1_4_8(
    std::size_t sites,
    const double* __restrict__ d_parent_clv,        // [sites * RC * 20]
    const unsigned* __restrict__ d_parent_scaler,   // [sites * RC] or nullptr
    const unsigned char* __restrict__ d_tipchars,   // [sites]
    const double* __restrict__ d_p_matrix,          // [RC * 20 * 20], row-major (j*20 + k)
    const unsigned* __restrict__ d_pattern_weights, // [sites] or nullptr
    double invar_proportion,
    const int* __restrict__ d_invar_indices,        // [sites] or nullptr; -1 if not invar
    const unsigned* __restrict__ d_site_scaling_min,// [sites] or nullptr
    double* __restrict__ d_site_loglk_out,
    bool per_rate_scaling)
{
    const int tid0   = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = gridDim.x  * blockDim.x;

    // ---- 頻率：double4 分塊（需要 __constant__ 16B 對齊） ----
    const double4 pi0 = *reinterpret_cast<const double4*>(&core_likelihood::c_frequencies[ 0]); //  0.. 3
    const double4 pi1 = *reinterpret_cast<const double4*>(&core_likelihood::c_frequencies[ 4]); //  4.. 7
    const double4 pi2 = *reinterpret_cast<const double4*>(&core_likelihood::c_frequencies[ 8]); //  8..11
    const double4 pi3 = *reinterpret_cast<const double4*>(&core_likelihood::c_frequencies[12]); // 12..15
    const double4 pi4 = *reinterpret_cast<const double4*>(&core_likelihood::c_frequencies[16]); // 16..19

    // ---- rate-weights 放暫存器 ----
    double wr[RC];
    #pragma unroll
    for (int r = 0; r < RC; ++r) wr[r] = core_likelihood::c_rate_weights[r];

    for (int base = tid0 * SITES_PER_THREAD; base < (int)sites; base += stride * SITES_PER_THREAD) {
        #pragma unroll
        for (int t = 0; t < SITES_PER_THREAD; ++t) {
            const int i = base + t;
            if (i >= (int)sites) break;

            const unsigned tip_state = static_cast<unsigned>(d_tipchars[i]);
            const unsigned mask      = core_likelihood::c_tipmap[tip_state];   // 20-state bitmask

            const double* __restrict__ clv_site = d_parent_clv + (size_t)i * RC * 20;

            double site_sum = 0.0;

            #pragma unroll
            for (int r = 0; r < RC; ++r) {
                const double* __restrict__ cr     = clv_site + r * 20;
                const double* __restrict__ pmat_r = d_p_matrix + r * 400;     // 20*20

                // 先把每個 parent 狀態 j 的 termb[j] 算好（依 tip mask 累加此列的 P）
                double termb[20];
                #pragma unroll 1
                for (int j = 0; j < 20; ++j) {
                    double s = 0.0;
                    const double* row = pmat_r + j * 20;
                    #pragma unroll 1
                    for (int k = 0; k < 20; ++k) {
                        if ((mask >> k) & 1u) s += row[k];
                    }
                    termb[j] = s;
                }

                // 用 double4 分塊載入 clv，並與 pi 分塊做乘加（每個元素再乘對應 termb[j]）
                const double4 v0 = *reinterpret_cast<const double4*>(&cr[ 0]); //  0.. 3
                const double4 v1 = *reinterpret_cast<const double4*>(&cr[ 4]); //  4.. 7
                const double4 v2 = *reinterpret_cast<const double4*>(&cr[ 8]); //  8..11
                const double4 v3 = *reinterpret_cast<const double4*>(&cr[12]); // 12..15
                const double4 v4 = *reinterpret_cast<const double4*>(&cr[16]); // 16..19

                double terma_r = 0.0;

                // j = 0..3
                terma_r = fma(v0.x * pi0.x, termb[ 0], terma_r);
                terma_r = fma(v0.y * pi0.y, termb[ 1], terma_r);
                terma_r = fma(v0.z * pi0.z, termb[ 2], terma_r);
                terma_r = fma(v0.w * pi0.w, termb[ 3], terma_r);
                // j = 4..7
                terma_r = fma(v1.x * pi1.x, termb[ 4], terma_r);
                terma_r = fma(v1.y * pi1.y, termb[ 5], terma_r);
                terma_r = fma(v1.z * pi1.z, termb[ 6], terma_r);
                terma_r = fma(v1.w * pi1.w, termb[ 7], terma_r);
                // j = 8..11
                terma_r = fma(v2.x * pi2.x, termb[ 8], terma_r);
                terma_r = fma(v2.y * pi2.y, termb[ 9], terma_r);
                terma_r = fma(v2.z * pi2.z, termb[10], terma_r);
                terma_r = fma(v2.w * pi2.w, termb[11], terma_r);
                // j = 12..15
                terma_r = fma(v3.x * pi3.x, termb[12], terma_r);
                terma_r = fma(v3.y * pi3.y, termb[13], terma_r);
                terma_r = fma(v3.z * pi3.z, termb[14], terma_r);
                terma_r = fma(v3.w * pi3.w, termb[15], terma_r);
                // j = 16..19
                terma_r = fma(v4.x * pi4.x, termb[16], terma_r);
                terma_r = fma(v4.y * pi4.y, termb[17], terma_r);
                terma_r = fma(v4.z * pi4.z, termb[18], terma_r);
                terma_r = fma(v4.w * pi4.w, termb[19], terma_r);

                // per-rate scaling（與 root 相同的查表）
                if (per_rate_scaling && d_parent_scaler && d_site_scaling_min) {
                    const int diff = min(int(d_parent_scaler[i * RC + r]) - int(d_site_scaling_min[i]),
                                         SCALE_MAX_DIFF);
                    if (diff > 0) terma_r *= core_likelihood::c_scale_minlh[diff - 1];
                }

                // 混入 invariant
                if (invar_proportion > 0.0 && d_invar_indices) {
                    const int inv_idx = d_invar_indices[i]; // -1 或 0..19
                    if (inv_idx >= 0) {
                        terma_r = (1.0 - invar_proportion) * terma_r
                                  + invar_proportion * core_likelihood::c_frequencies[inv_idx];
                    }
                }

                site_sum = fma(wr[r], terma_r, site_sum);
            }

            // per-site log-likelihood
            const double eps = 1e-300;
            double loglk = log(site_sum > eps ? site_sum : eps);

            
            if (d_site_scaling_min) loglk += static_cast<double>(d_site_scaling_min[i]) * LOG_SCALE_THRESHOLD;
            if (d_pattern_weights) loglk *= double(d_pattern_weights[i]);
            d_site_loglk_out[i] = loglk;
        }
    }
}


__global__ void TipInnerLikelihoodCalculation(
    std::size_t sites,
    int states,
    int rate_cats,
    const double* d_parent_clv,
    const unsigned* d_parent_scaler,      
    const unsigned char* d_tipchars,
    const double*  d_p_matrix,
    const unsigned* d_pattern_weights,
    double invar_proportion,
    const int*  d_invar_indices,
    const unsigned* d_site_scaling_min,
    double* d_site_loglk_out,
    bool per_rate_scaling
)
{

    const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    const double* clv_site = d_parent_clv + i * rate_cats * states;  
    double site_loglk;
     
    if (i >= sites) return;
    
    double site_sum = 0.0;
   
    for (int r = 0; r < rate_cats; ++r) {
        
        const double* pmat_r   = d_p_matrix + r * states * states;

        double terma_r = 0.0;
        for (int j = 0; j < states; ++j) { // Loop over parent states
            double termb = 0.0;
            unsigned tip_state = static_cast<unsigned>(d_tipchars[i]);
            unsigned mask = core_likelihood::c_tipmap[tip_state];
            for (int k = 0; k < states; ++k)
            {
                if ((mask >> k) & 1u)
                    termb += pmat_r[j*states + k];
            }
            
            terma_r += clv_site[j] * core_likelihood::c_frequencies[j] * termb;
        }

        int site_rate_scaling_diff = (d_site_scaling_min)? min(d_parent_scaler[i * rate_cats + r] - d_site_scaling_min[i], SCALE_MAX_DIFF) : 0;
        if(site_rate_scaling_diff > 0){
            terma_r *= core_likelihood::c_scale_minlh[site_rate_scaling_diff -1];
        }

        if (invar_proportion > 0.0 && d_invar_indices) {
            const int inv_idx = d_invar_indices[i];
            const double inv_site_lk = (inv_idx >= 0) ? core_likelihood::c_frequencies[inv_idx] : 0.0;
            site_sum += core_likelihood::c_rate_weights[r] * ((1.0 - invar_proportion) * terma_r + invar_proportion * inv_site_lk);
        } else {
            site_sum += core_likelihood::c_rate_weights[r] * terma_r;
        }
        clv_site += states;
    }

    site_loglk = log(site_sum);

    if (d_site_scaling_min) {
        site_loglk += static_cast<double>(d_site_scaling_min[i]) * LOG_SCALE_THRESHOLD;
    }

    if (d_pattern_weights) {
        site_loglk *= static_cast<double>(d_pattern_weights[i]);
    }

    d_site_loglk_out[i] = site_loglk;
}

void core_likelihood::Likelihood_Tip_Inner::ComputeLikelihood(
    const Param& p)
{
    cudaError_t err;
    cudaStream_t stream = 0;
    
    dim3 block(256), grid((p.sites + block.x*4 - 1)/(block.x*4)); 
        
    switch(p.states) {
        case 4:
            switch (p.rate_cats) {
                case 1:
                    TipInnerLikelihoodCalculation_states_4_1_4_8<1,4><<<grid, block, 0, stream>>>(
                        p.sites,
                        d_parent_clv,
                        d_parent_scaler,
                        d_tipchars,
                        d_p_matrix,
                        d_pattern_weights,
                        p.invar_proportion,
                        d_invar_indices,
                        d_site_scaling_min,
                        d_per_site,
                        p.per_rate_scaling);
                    break;
                case 4:
                    TipInnerLikelihoodCalculation_states_4_1_4_8<4,4><<<grid, block, 0, stream>>>(
                        p.sites,
                        d_parent_clv,
                        d_parent_scaler,
                        d_tipchars,
                        d_p_matrix,
                        d_pattern_weights,
                        p.invar_proportion,
                        d_invar_indices,
                        d_site_scaling_min,
                        d_per_site,
                        p.per_rate_scaling);
                    break;
                case 8:
                    TipInnerLikelihoodCalculation_states_4_1_4_8<8,4><<<grid, block, 0, stream>>>(
                        p.sites,
                        d_parent_clv,
                        d_parent_scaler,
                        d_tipchars,
                        d_p_matrix,
                        d_pattern_weights,
                        p.invar_proportion,
                        d_invar_indices,
                        d_site_scaling_min,
                        d_per_site,
                        p.per_rate_scaling);
                    break;
                default:
                    goto GENERAL;
            }
            break;
        case 5:
            switch (p.rate_cats) {
                case 1:
                    TipInnerLikelihoodCalculation_states_5_1_4_8<1,4><<<grid, block, 0, stream>>>(
                        p.sites,
                        d_parent_clv,
                        d_parent_scaler,
                        d_tipchars,
                        d_p_matrix,
                        d_pattern_weights,
                        p.invar_proportion,
                        d_invar_indices,
                        d_site_scaling_min,
                        d_per_site,
                        p.per_rate_scaling);
                    break;
                case 4:
                    TipInnerLikelihoodCalculation_states_5_1_4_8<4,4><<<grid, block, 0, stream>>>(
                        p.sites,
                        d_parent_clv,
                        d_parent_scaler,
                        d_tipchars,
                        d_p_matrix,
                        d_pattern_weights,
                        p.invar_proportion,
                        d_invar_indices,
                        d_site_scaling_min,
                        d_per_site,
                        p.per_rate_scaling);
                    break;
                case 8:
                    TipInnerLikelihoodCalculation_states_5_1_4_8<8,4><<<grid, block, 0, stream>>>(
                        p.sites,
                        d_parent_clv,
                        d_parent_scaler,
                        d_tipchars,
                        d_p_matrix,
                        d_pattern_weights,
                        p.invar_proportion,
                        d_invar_indices,
                        d_site_scaling_min,
                        d_per_site,
                        p.per_rate_scaling);
                    break;
                default:
                    goto GENERAL; 
            }
            break;
        case 20:
            switch (p.rate_cats) {
                case 1:
                    TipInnerLikelihoodCalculation_states_20_1_4_8<1,4><<<grid, block, 0, stream>>>(
                        p.sites,
                        d_parent_clv,
                        d_parent_scaler,
                        d_tipchars,
                        d_p_matrix,
                        d_pattern_weights,
                        p.invar_proportion,
                        d_invar_indices,
                        d_site_scaling_min,
                        d_per_site,
                        p.per_rate_scaling);
                    break;
                case 4:
                    TipInnerLikelihoodCalculation_states_20_1_4_8<4,4><<<grid, block, 0, stream>>>(
                        p.sites,
                        d_parent_clv,
                        d_parent_scaler,
                        d_tipchars,
                        d_p_matrix,
                        d_pattern_weights,
                        p.invar_proportion,
                        d_invar_indices,
                        d_site_scaling_min,
                        d_per_site,
                        p.per_rate_scaling);
                    break;
                case 8:
                    TipInnerLikelihoodCalculation_states_20_1_4_8<8,4><<<grid, block, 0, stream>>>(
                        p.sites,
                        d_parent_clv,
                        d_parent_scaler,
                        d_tipchars,
                        d_p_matrix,
                        d_pattern_weights,
                        p.invar_proportion,
                        d_invar_indices,
                        d_site_scaling_min,
                        d_per_site,
                        p.per_rate_scaling);
                    break;
                    break;
                default:
                    goto GENERAL;
            }
            break;
        default:
            GENERAL:
                TipInnerLikelihoodCalculation<<<grid, block, 0, stream>>>(
                    p.sites,
                    p.states,
                    p.rate_cats,
                    d_parent_clv,
                    d_parent_scaler,
                    d_tipchars,
                    d_p_matrix,
                    d_pattern_weights,
                    p.invar_proportion,
                    d_invar_indices,
                    d_site_scaling_min,
                    d_per_site,
                    p.per_rate_scaling);
                break;
    }
    
    device_reduce_sum(d_per_site, p.sites, d_out_sum);

    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "[Gpu_ERROR] device_reduce_sum failed: %s\n",cudaGetErrorString(err));
        exit(1);
    }

    err = cudaMemcpy(&total_loglik, d_out_sum, sizeof(double), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "[Gpu_ERROR] d_out_sum copy failed: %s\n",cudaGetErrorString(err));
        exit(1);
    }
}

void core_likelihood::Likelihood_Tip_Inner::PrintLikelihood() const
{
    printf("[Likelihood_Tip_Inner] Total log-likelihood = %.12f \n", total_loglik);
}
void core_likelihood::Likelihood_Tip_Inner::CleanUp()
{
    cudaFree(d_per_site);
    cudaFree(d_out_sum);
    cudaFree(d_parent_clv);
    cudaFree(d_parent_scaler);
    cudaFree(d_tipchars);
    cudaFree(d_p_matrix);
    cudaFree(d_pattern_weights);
    cudaFree(d_invar_indices);
    if (d_site_scaling_min) cudaFree(d_site_scaling_min);
}

// ==========================================================
// Likelihood_Inner_Inner
// ==========================================================

void core_likelihood::Likelihood_Inner_Inner::Initialize(Param& p)
{
    cudaError_t err;

    err = cudaMalloc((void**) &d_per_site, sizeof(double) * p.sites);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_per_site allocation failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMemset(d_per_site, 0, sizeof(double) * p.sites);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_per_site reset failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }

    err = cudaMalloc((void**) &d_out_sum, sizeof(double));
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_out_sum allocation failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMemset(d_out_sum, 0, sizeof(double));
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_out_sum reset failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }

    err = cudaMalloc((void**) &d_parent_clv, sizeof(double) * p.sites * p.rate_cats * p.states);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_parent_clv allocation failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMemset(d_parent_clv, 0, sizeof(double) * p.sites * p.rate_cats * p.states);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_parent_clv reset failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }

    const int scaler_size = p.per_rate_scaling ? (p.sites * p.rate_cats) : p.sites;

    err = cudaMalloc((void**) &d_parent_scaler, sizeof(unsigned int) * scaler_size);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_parent_scaler allocation failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMemset(d_parent_scaler, 0, sizeof(unsigned int) * scaler_size);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_parent_scaler reset failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }

    err = cudaMalloc((void**) &d_child_clv, sizeof(double) * p.sites * p.rate_cats * p.states);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_child_clv allocation failed: %s\n",  cudaGetErrorString(err));
        exit(1);
    }  

    err = cudaMemset(d_child_clv, 0, sizeof(double) * p.sites * p.rate_cats * p.states);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_child_clv reset failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMalloc((void**) &d_child_scaler, sizeof(unsigned int) * scaler_size);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_child_scaler allocation failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMemset(d_child_scaler, 0, sizeof(unsigned int) * scaler_size);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_child_scaler reset failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }

    err = cudaMalloc((void**) &d_site_scaler, sizeof(unsigned int) * p.sites * p.rate_cats);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_site_scaler allocation failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMemset(d_site_scaler, 0, sizeof(unsigned int) * p.sites * p.rate_cats);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_site_scaler reset failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }

    err = cudaMalloc((void**) &d_p_matrix, sizeof(double) * p.rate_cats * p.states * p.states);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_p_matrix allocation failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMemset(d_p_matrix, 0, sizeof(double) * p.rate_cats * p.states * p.states);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_p_matrix reset failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }

    err = cudaMalloc((void**) &d_pattern_weights, sizeof(unsigned int) * p.sites);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_pattern_weights allocation failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMemset(d_pattern_weights, 0, sizeof(unsigned int) * p.sites);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_pattern_weights reset failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }

    err = cudaMalloc((void**) &d_invar_indices, sizeof(int) * p.sites);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_invar_indices allocation failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMemset(d_invar_indices, 0, sizeof(int) * p.sites);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_invar_indices reset failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }

}

void core_likelihood::Likelihood_Inner_Inner::ConstructionOnGpu(
    const double * h_parent_clv,
    const unsigned int * h_parent_scaler,
    const double * h_child_clv,
    const unsigned int * h_child_scaler,
    const double * h_p_matrix,
    double * h_frequencies,
    const double * h_rate_weights,
    const unsigned int * h_pattern_weights,
    const int * h_invar_indices,
    Param& p)
{
    cudaError_t err;

    err = cudaMemcpy(d_parent_clv, h_parent_clv, sizeof(double) * p.sites * p.rate_cats * p.states, cudaMemcpyHostToDevice);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_parent_clv copy failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }

    const size_t scaler_size = p.per_rate_scaling ? (p.sites * p.rate_cats) : p.sites;
    if(h_parent_scaler)
    {
        err = cudaMemcpy(d_parent_scaler, h_parent_scaler, sizeof(unsigned int) * scaler_size , cudaMemcpyHostToDevice);
        if (err != cudaSuccess){
            fprintf(stderr, "[Gpu_ERROR] d_parent_caler copy failed: %s\n", cudaGetErrorString(err));
            exit(1);
        }
    }

    err = cudaMemcpy(d_child_clv, h_child_clv, sizeof(double) * p.sites * p.rate_cats * p.states, cudaMemcpyHostToDevice);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_child_clv copy failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }

    if(h_child_scaler)
    {
        const size_t scaler_size = p.per_rate_scaling ? p.sites * p.rate_cats : p.sites;
        err = cudaMemcpy(d_child_scaler, h_child_scaler, sizeof(unsigned int) * scaler_size, cudaMemcpyHostToDevice);
        if (err != cudaSuccess){
            fprintf(stderr, "[Gpu_ERROR] d_child_scaler copy failed: %s\n", cudaGetErrorString(err));
            exit(1);
        }
    }

    err = cudaMemcpy(d_p_matrix, h_p_matrix, sizeof(double) * p.rate_cats * p.states * p.states, cudaMemcpyHostToDevice);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_p_matrix copy failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMemcpyToSymbol(core_likelihood::c_frequencies, h_frequencies, sizeof(double) * p.states);
    if (err != cudaSuccess) {
        std::fprintf(stderr, "cpy c_frequencies: %s\n", cudaGetErrorString(err));
        exit(1); 
    }

    err = cudaMemcpyToSymbol(core_likelihood::c_rate_weights, h_rate_weights, sizeof(double) * p.rate_cats);
    if (err != cudaSuccess) {
        std::fprintf(stderr, "cpy c_rate_weights: %s\n", cudaGetErrorString(err));
        exit(1); 
    }
    
    err = cudaMemcpy(d_pattern_weights, h_pattern_weights, sizeof(unsigned int) * p.sites, cudaMemcpyHostToDevice);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_pattern_weights copy failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMemcpy(d_invar_indices, h_invar_indices, sizeof(int) * p.sites, cudaMemcpyHostToDevice);
    if (err != cudaSuccess){
        fprintf(stderr, "[Gpu_ERROR] d_invar_indices copy failed: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    std::vector<double> scale_minlh(SCALE_MAX_DIFF);
    double scale_factor = 1.0;
    if(p.per_rate_scaling){
        
        for (int i = 0; i < SCALE_MAX_DIFF; ++i) {
            scale_factor *= SCALE_THRESHOLD;
            scale_minlh[i] = scale_factor;
        }
    
        size_t sym_size = 0;
        err = cudaGetSymbolSize(&sym_size, c_scale_minlh);
        if (err != cudaSuccess) {
            fprintf(stderr, "[Gpu_ERROR] cudaGetSymbolSize(c_scale_minlh) failed: %s\n",
                    cudaGetErrorString(err));
            std::exit(1);
        }
        const size_t want = sizeof(double) * SCALE_MAX_DIFF;
        if (sym_size != want) {
            fprintf(stderr, "[Gpu_ERROR] c_scale_minlh size mismatch (device=%zu, host=%zu)\n",
                    sym_size, want);
            std::exit(1);
        }

        err = cudaMemcpyToSymbol(core_likelihood::c_scale_minlh, scale_minlh.data(), SCALE_MAX_DIFF * sizeof(double), cudaMemcpyHostToDevice);
        if (err != cudaSuccess){
            fprintf(stderr, "[Gpu_ERROR] c_scale_minlh copy failed: %s\n", cudaGetErrorString(err));
            exit(1);
        }
        
        err = cudaMalloc((void**)&d_site_scaling_min, sizeof(unsigned) * p.sites);
        if(err != cudaSuccess){
            fprintf(stderr, "[Gpu_ERROR] d_site_scaling_min allocation failed: %s\n", cudaGetErrorString(err));
            exit(1);
        }
        
        const size_t N = static_cast<size_t>(p.sites) * static_cast<size_t>(p.rate_cats);

        thrust::counting_iterator<unsigned> idx_begin(0);
        auto keys_begin = thrust::make_transform_iterator(
            idx_begin, KeyBySite{static_cast<unsigned>(p.rate_cats)}
        );

        auto values_begin = thrust::make_transform_iterator(
            idx_begin, SumAndStore{ d_parent_scaler, d_child_scaler, d_site_scaler }
        );
        cudaFree(d_parent_scaler);
        cudaFree(d_child_scaler);

        thrust::reduce_by_key(
            keys_begin, keys_begin + N,
            values_begin,
            thrust::make_discard_iterator(),
            thrust::device_pointer_cast(d_site_scaling_min),
            thrust::equal_to<unsigned>(),
            thrust::minimum<unsigned>()
        );

        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            fprintf(stderr, "[Gpu_ERROR] thrust::reduce_by_key (min per-site) failed: %s\n", cudaGetErrorString(err));
            exit(1);
        }
    }
    else{
        d_site_scaling_min = nullptr;
        d_site_scaler = nullptr;
    }

}

template<int RC, int SITES_PER_THREAD = 2>
__global__ __launch_bounds__(256, 2)
void InnerInnerLikelihoodCalculation_states_4_1_4_8(
    std::size_t sites,
    const double* __restrict__ d_parent_clv,        // [sites * RC * 4]
    const double* __restrict__ d_child_clv,         // [sites * RC * 4]
    const unsigned* __restrict__ d_site_scaler,     // [sites * RC] or nullptr
    const double* __restrict__ d_p_matrix,          // [RC * 4 * 4], row-major
    const unsigned* __restrict__ d_pattern_weights, // [sites] or nullptr
    double invar_proportion,
    const int* __restrict__ d_invar_indices,        // [sites] or nullptr; -1 if not invar
    const unsigned* __restrict__ d_site_scaling_min,// [sites] or nullptr
    double* __restrict__ d_site_loglk_out,
    bool per_rate_scaling)
{
    const int tid0   = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = gridDim.x  * blockDim.x;

    const double pi0 = core_likelihood::c_frequencies[0];
    const double pi1 = core_likelihood::c_frequencies[1];
    const double pi2 = core_likelihood::c_frequencies[2];
    const double pi3 = core_likelihood::c_frequencies[3];

    double wr[RC];
    #pragma unroll
    for (int r = 0; r < RC; ++r) wr[r] = core_likelihood::c_rate_weights[r];

    for (std::size_t base = std::size_t(tid0) * SITES_PER_THREAD;
         base < sites;
         base += std::size_t(stride) * SITES_PER_THREAD)
    {
        #pragma unroll
        for (int t = 0; t < SITES_PER_THREAD; ++t) {
            const std::size_t i = base + t;
            if (i >= sites) break;

            const double* __restrict__ clv_p_site = d_parent_clv + i * (std::size_t)RC * 4;
            const double* __restrict__ clv_c_site = d_child_clv  + i * (std::size_t)RC * 4;

            double site_sum = 0.0;

            #pragma unroll
            for (int r = 0; r < RC; ++r) {
                const double* __restrict__ pmat_r = d_p_matrix + r * 16;   // 4*4
                const double* __restrict__ pr    = clv_p_site + r * 4;
                const double* __restrict__ cr    = clv_c_site + r * 4;

                // terma_r = Σ_j pr[j] * pi[j] * ( Σ_k P_r[j,k] * cr[k] )
                double terma_r = 0.0;

                // j = 0
                {
                    const double termb = fma(pmat_r[0],  cr[0],
                                        fma(pmat_r[1],  cr[1],
                                        fma(pmat_r[2],  cr[2],
                                            pmat_r[3] * cr[3])));
                    terma_r = fma(pr[0] * pi0, termb, terma_r);
                }
                // j = 1
                {
                    const double termb = fma(pmat_r[4],  cr[0],
                                        fma(pmat_r[5],  cr[1],
                                        fma(pmat_r[6],  cr[2],
                                            pmat_r[7] * cr[3])));
                    terma_r = fma(pr[1] * pi1, termb, terma_r);
                }
                // j = 2
                {
                    const double termb = fma(pmat_r[8],  cr[0],
                                        fma(pmat_r[9],  cr[1],
                                        fma(pmat_r[10], cr[2],
                                            pmat_r[11] * cr[3])));
                    terma_r = fma(pr[2] * pi2, termb, terma_r);
                }
                // j = 3
                {
                    const double termb = fma(pmat_r[12], cr[0],
                                        fma(pmat_r[13], cr[1],
                                        fma(pmat_r[14], cr[2],
                                            pmat_r[15] * cr[3])));
                    terma_r = fma(pr[3] * pi3, termb, terma_r);
                }

                // per-rate scaling（min-lh）補償
                if (per_rate_scaling && d_site_scaler && d_site_scaling_min) {
                    const int diff = min(int(d_site_scaler[i * RC + r]) - int(d_site_scaling_min[i]),
                                         (int)SCALE_MAX_DIFF);
                    if (diff > 0) terma_r *= core_likelihood::c_scale_minlh[diff - 1];
                }

                // invariant fastpath（和 tip-inner 相同語義）
                if (invar_proportion > 0.0 && d_invar_indices) {
                    const int inv_idx = d_invar_indices[i]; // -1 非 invar；0..3 為狀態
                    const double inv_lk = (inv_idx >= 0)
                        ? ((inv_idx==0)?pi0 : (inv_idx==1)?pi1 : (inv_idx==2)?pi2 : pi3)
                        : 0.0;
                    site_sum += wr[r] * ((1.0 - invar_proportion) * terma_r + invar_proportion * inv_lk);
                } else {
                    site_sum += wr[r] * terma_r;
                }
            }

            // loglk with 防 underflow
            const double eps = 1e-300;
            double loglk = log(site_sum > eps ? site_sum : eps);

            // site_min（跨 rate 的共同 scaling）補償
            if (d_site_scaling_min) loglk += double(d_site_scaling_min[i]) * LOG_SCALE_THRESHOLD;

            // pattern weights
            if (d_pattern_weights) loglk *= double(d_pattern_weights[i]);

            d_site_loglk_out[i] = loglk;
        }
    }
}

template<int RC, int SITES_PER_THREAD = 2>
__global__ __launch_bounds__(256, 2)
void InnerInnerLikelihoodCalculation_states_5_1_4_8(
    std::size_t sites,
    const double* __restrict__ d_parent_clv,        // [sites * RC * 5]
    const double* __restrict__ d_child_clv,         // [sites * RC * 5]
    const unsigned* __restrict__ d_site_scaler,     // [sites * RC] or nullptr
    const double* __restrict__ d_p_matrix,          // [RC * 5 * 5], row-major
    const unsigned* __restrict__ d_pattern_weights, // [sites] or nullptr
    double invar_proportion,
    const int* __restrict__ d_invar_indices,        // [sites] or nullptr; -1 if not invar
    const unsigned* __restrict__ d_site_scaling_min,// [sites] or nullptr
    double* __restrict__ d_site_loglk_out,
    bool per_rate_scaling)
{
    const int tid0   = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = gridDim.x  * blockDim.x;

    // 常數頻率（5態）
    const double pi0 = core_likelihood::c_frequencies[0];
    const double pi1 = core_likelihood::c_frequencies[1];
    const double pi2 = core_likelihood::c_frequencies[2];
    const double pi3 = core_likelihood::c_frequencies[3];
    const double pi4 = core_likelihood::c_frequencies[4];

    // 每個 rate 的權重
    double wr[RC];
    #pragma unroll
    for (int r = 0; r < RC; ++r) wr[r] = core_likelihood::c_rate_weights[r];

    for (std::size_t base = std::size_t(tid0) * SITES_PER_THREAD;
         base < sites;
         base += std::size_t(stride) * SITES_PER_THREAD)
    {
        #pragma unroll
        for (int t = 0; t < SITES_PER_THREAD; ++t) {
            const std::size_t i = base + t;
            if (i >= sites) break;

            const double* __restrict__ clv_p_site = d_parent_clv + i * (std::size_t)RC * 5;
            const double* __restrict__ clv_c_site = d_child_clv  + i * (std::size_t)RC * 5;

            double site_sum = 0.0;

            #pragma unroll
            for (int r = 0; r < RC; ++r) {
                const double* __restrict__ pmat_r = d_p_matrix + r * 25; // 5*5
                const double* __restrict__ pr    = clv_p_site + r * 5;
                const double* __restrict__ cr    = clv_c_site + r * 5;

                double terma_r = 0.0;

                // j = 0
                {
                    const double termb =
                        fma(pmat_r[0],  cr[0],
                        fma(pmat_r[1],  cr[1],
                        fma(pmat_r[2],  cr[2],
                        fma(pmat_r[3],  cr[3],
                            pmat_r[4] *  cr[4]))));
                    terma_r = fma(pr[0] * pi0, termb, terma_r);
                }
                // j = 1
                {
                    const double termb =
                        fma(pmat_r[5],  cr[0],
                        fma(pmat_r[6],  cr[1],
                        fma(pmat_r[7],  cr[2],
                        fma(pmat_r[8],  cr[3],
                            pmat_r[9] *  cr[4]))));
                    terma_r = fma(pr[1] * pi1, termb, terma_r);
                }
                // j = 2
                {
                    const double termb =
                        fma(pmat_r[10], cr[0],
                        fma(pmat_r[11], cr[1],
                        fma(pmat_r[12], cr[2],
                        fma(pmat_r[13], cr[3],
                            pmat_r[14] * cr[4]))));
                    terma_r = fma(pr[2] * pi2, termb, terma_r);
                }
                // j = 3
                {
                    const double termb =
                        fma(pmat_r[15], cr[0],
                        fma(pmat_r[16], cr[1],
                        fma(pmat_r[17], cr[2],
                        fma(pmat_r[18], cr[3],
                            pmat_r[19] * cr[4]))));
                    terma_r = fma(pr[3] * pi3, termb, terma_r);
                }
                // j = 4
                {
                    const double termb =
                        fma(pmat_r[20], cr[0],
                        fma(pmat_r[21], cr[1],
                        fma(pmat_r[22], cr[2],
                        fma(pmat_r[23], cr[3],
                            pmat_r[24] * cr[4]))));
                    terma_r = fma(pr[4] * pi4, termb, terma_r);
                }

                // per-rate scaling（min-lh）補償
                if (per_rate_scaling && d_site_scaler && d_site_scaling_min) {
                    const int diff = min(int(d_site_scaler[i * RC + r]) - int(d_site_scaling_min[i]),
                                         (int)SCALE_MAX_DIFF);
                    if (diff > 0) terma_r *= core_likelihood::c_scale_minlh[diff - 1];
                }

                // invariant fastpath（內點：用 pi[inv_idx]）
                if (invar_proportion > 0.0 && d_invar_indices) {
                    const int inv_idx = d_invar_indices[i]; // -1: 非 invar；0..4: 狀態
                    double inv_lk = 0.0;
                    if (inv_idx >= 0) {
                        inv_lk = (inv_idx==0)?pi0 :
                                 (inv_idx==1)?pi1 :
                                 (inv_idx==2)?pi2 :
                                 (inv_idx==3)?pi3 : pi4;
                    }
                    site_sum += wr[r] * ((1.0 - invar_proportion) * terma_r + invar_proportion * inv_lk);
                } else {
                    site_sum += wr[r] * terma_r;
                }
            }

            // log-likelihood（防 underflow）
            const double eps = 1e-300;
            double loglk = log(site_sum > eps ? site_sum : eps);

            // site_min（跨 rate 的共同 scaling）補償
            if (d_site_scaling_min) loglk += double(d_site_scaling_min[i]) * LOG_SCALE_THRESHOLD;

            // pattern weights
            if (d_pattern_weights) loglk *= double(d_pattern_weights[i]);

            d_site_loglk_out[i] = loglk;
        }
    }
}

// ===== Protein (20 states) Inner-Inner kernel =====
template<int RC, int SITES_PER_THREAD = 2>
__global__ __launch_bounds__(256, 2)
void InnerInnerLikelihoodCalculation_states_20_1_4_8(
    std::size_t sites,
    const double* __restrict__ d_parent_clv,        // [sites * RC * 20]
    const double* __restrict__ d_child_clv,         // [sites * RC * 20]
    const unsigned* __restrict__ d_site_scaler,     // [sites * RC] or nullptr
    const double* __restrict__ d_p_matrix,          // [RC * 20 * 20], row-major
    const unsigned* __restrict__ d_pattern_weights, // [sites] or nullptr
    double invar_proportion,
    const int* __restrict__ d_invar_indices,        // [sites] or nullptr; -1 if not invar
    const unsigned* __restrict__ d_site_scaling_min,// [sites] or nullptr
    double* __restrict__ d_site_loglk_out,
    bool per_rate_scaling)
{
    const int tid0   = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = gridDim.x  * blockDim.x;

    // 頻率暫存（20 states）
    double pi[20];
    #pragma unroll
    for (int s = 0; s < 20; ++s)
        pi[s] = core_likelihood::c_frequencies[s];

    // rate 權重暫存
    double wr[RC];
    #pragma unroll
    for (int r = 0; r < RC; ++r)
        wr[r] = core_likelihood::c_rate_weights[r];

    for (std::size_t base = std::size_t(tid0) * SITES_PER_THREAD;
         base < sites;
         base += std::size_t(stride) * SITES_PER_THREAD)
    {
        #pragma unroll
        for (int t = 0; t < SITES_PER_THREAD; ++t) {
            const std::size_t i = base + t;
            if (i >= sites) break;

            const double* __restrict__ clv_p_site = d_parent_clv + i * (std::size_t)RC * 20;
            const double* __restrict__ clv_c_site = d_child_clv  + i * (std::size_t)RC * 20;

            double site_sum = 0.0;

            #pragma unroll
            for (int r = 0; r < RC; ++r) {
                const double* __restrict__ pmat_r = d_p_matrix + r * 400; // 20×20
                const double* __restrict__ pr     = clv_p_site + r * 20;
                const double* __restrict__ cr     = clv_c_site + r * 20;

                // ---- compute terma_r = Σ_j pr[j]*pi[j]*Σ_k P[j,k]*cr[k] ----
                double terma_r = 0.0;

                #pragma unroll
                for (int j = 0; j < 20; ++j) {
                    const double* row = pmat_r + j * 20;
                    // dot product row · cr (20 項)
                    double termb =
                        fma(row[0],  cr[0],
                        fma(row[1],  cr[1],
                        fma(row[2],  cr[2],
                        fma(row[3],  cr[3],
                        fma(row[4],  cr[4],
                        fma(row[5],  cr[5],
                        fma(row[6],  cr[6],
                        fma(row[7],  cr[7],
                        fma(row[8],  cr[8],
                        fma(row[9],  cr[9],
                        fma(row[10], cr[10],
                        fma(row[11], cr[11],
                        fma(row[12], cr[12],
                        fma(row[13], cr[13],
                        fma(row[14], cr[14],
                        fma(row[15], cr[15],
                        fma(row[16], cr[16],
                        fma(row[17], cr[17],
                        fma(row[18], cr[18],
                            row[19] * cr[19])))))))))))))))))));
                    terma_r = fma(pr[j] * pi[j], termb, terma_r);
                }

                // per-rate scaling 補償
                if (per_rate_scaling && d_site_scaler && d_site_scaling_min ) {
                    const int diff = min(int(d_site_scaler[i * RC + r]) - int(d_site_scaling_min[i]),
                                         (int)SCALE_MAX_DIFF);
                    if (diff > 0) terma_r *= core_likelihood::c_scale_minlh[diff - 1];
                }

                // invariant fastpath
                if (invar_proportion > 0.0 && d_invar_indices) {
                    const int inv_idx = d_invar_indices[i]; // -1 非 invar；0..19 有效狀態
                    const double inv_lk = (inv_idx >= 0) ? pi[inv_idx] : 0.0;
                    site_sum += wr[r] * ((1.0 - invar_proportion) * terma_r + invar_proportion * inv_lk);
                } else {
                    site_sum += wr[r] * terma_r;
                }
            }

            // log-likelihood with underflow guard
            const double eps = 1e-300;
            double loglk = log(site_sum > eps ? site_sum : eps);

            if (d_site_scaling_min)
                loglk += double(d_site_scaling_min[i]) * LOG_SCALE_THRESHOLD;
            if (d_pattern_weights)
                loglk *= double(d_pattern_weights[i]);

            d_site_loglk_out[i] = loglk;
        }
    }
}


__global__ void InnerInnerLikelihoodCalculation(
    const int sites,
    const int states,
    const int rate_cats,
    const double *d_parent_clv,
    const double *d_child_clv,
    const unsigned int * d_site_scaler,
    const double *d_p_matrix,
    const unsigned int *d_pattern_weights,
    const int *d_invar_indices,
    const double invar_proportion,
    const unsigned* d_site_scaling_min,
    double *d_site_loglk_out,
    const bool per_rate_scaling 
){
    const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    const double* clv_site_p = d_parent_clv + i * rate_cats * states;
    const double* clv_site_c = d_child_clv + i * rate_cats * states;  
    double site_loglk;
    
    if(i >= sites) return;
    double site_sum = 0.0;
    for (int r = 0; r < rate_cats; r++) {
        const double* pmat_r = d_p_matrix + r * states * states;
        double terma_r = 0.0;
        for (int j = 0; j < states; j++) {
            double termb = 0.0;
            for (int k = 0; k < states; k++) {
                termb += pmat_r[j * states + k] * clv_site_c[k];
            }
            terma_r += clv_site_p[j] * core_likelihood::c_frequencies[j] * termb;
        }
        int site_rate_scaling_diff = (d_site_scaling_min)? min(d_site_scaler[i * rate_cats + r] - d_site_scaling_min[i], SCALE_MAX_DIFF) : 0;
        if(site_rate_scaling_diff > 0){
            terma_r *= core_likelihood::c_scale_minlh[site_rate_scaling_diff -1];
        }
        
        if (invar_proportion > 0.0 && d_invar_indices) {
            const int inv_idx = d_invar_indices[i];
            const double inv_site_lk = (inv_idx >= 0) ? core_likelihood::c_frequencies[inv_idx] : 0.0;
            site_sum += core_likelihood::c_rate_weights[r] * ((1.0 - invar_proportion) * terma_r + invar_proportion * inv_site_lk);
        } else {
            site_sum += core_likelihood::c_rate_weights[r] * terma_r;
        }
        
        clv_site_p += states;
        clv_site_c += states;
    }
    
    site_loglk = log(site_sum);

    if (d_site_scaling_min) {
        site_loglk += static_cast<double>(d_site_scaling_min[i]) * LOG_SCALE_THRESHOLD;
    }

    if (d_pattern_weights) {
        site_loglk *= static_cast<double>(d_pattern_weights[i]);
    }

    d_site_loglk_out[i] = site_loglk;
}

void core_likelihood::Likelihood_Inner_Inner::ComputeLikelihood(
    Param& p)
{
    cudaError_t err;
    
    cudaStream_t stream = 0;
    dim3 block(256), grid((p.sites + block.x*4 - 1)/(block.x*4)); 
        
    switch(p.states) {
        case 4:
            switch (p.rate_cats) {
                case 1:
                    InnerInnerLikelihoodCalculation_states_4_1_4_8<1,4><<<grid,block,0,stream>>>(
                                p.sites,
                                d_parent_clv,
                                d_child_clv,
                                d_site_scaler,
                                d_p_matrix,
                                d_pattern_weights,
                                p.invar_proportion,
                                d_invar_indices,
                                d_site_scaling_min,
                                d_per_site,
                                p.per_rate_scaling);
                    break;
                case 4:
                    InnerInnerLikelihoodCalculation_states_4_1_4_8<4,4><<<grid,block,0,stream>>>(
                                p.sites,
                                d_parent_clv,
                                d_child_clv,
                                d_site_scaler,
                                d_p_matrix,
                                d_pattern_weights,
                                p.invar_proportion,
                                d_invar_indices,
                                d_site_scaling_min,
                                d_per_site,
                                p.per_rate_scaling);
                    break;
                case 8:
                    InnerInnerLikelihoodCalculation_states_4_1_4_8<8,4><<<grid,block,0,stream>>>(
                                p.sites,
                                d_parent_clv,
                                d_child_clv,
                                d_site_scaler,
                                d_p_matrix,
                                d_pattern_weights,
                                p.invar_proportion,
                                d_invar_indices,
                                d_site_scaling_min,
                                d_per_site,
                                p.per_rate_scaling);
                    break;
                default:
                    goto GENERAL;
            }
            break;
        case 5:
            switch (p.rate_cats) {
                case 1:
                    InnerInnerLikelihoodCalculation_states_5_1_4_8<1,4><<<grid,block,0,stream>>>(
                                p.sites,
                                d_parent_clv,
                                d_child_clv,
                                d_site_scaler,
                                d_p_matrix,
                                d_pattern_weights,
                                p.invar_proportion,
                                d_invar_indices,
                                d_site_scaling_min,
                                d_per_site,
                                p.per_rate_scaling);
                    break;
                case 4:
                    InnerInnerLikelihoodCalculation_states_5_1_4_8<4,4><<<grid,block,0,stream>>>(
                                p.sites,
                                d_parent_clv,
                                d_child_clv,
                                d_site_scaler,
                                d_p_matrix,
                                d_pattern_weights,
                                p.invar_proportion,
                                d_invar_indices,
                                d_site_scaling_min,
                                d_per_site,
                                p.per_rate_scaling);
                    break;
                case 8:
                    InnerInnerLikelihoodCalculation_states_5_1_4_8<8,4><<<grid,block,0,stream>>>(
                                p.sites,
                                d_parent_clv,
                                d_child_clv,
                                d_site_scaler,
                                d_p_matrix,
                                d_pattern_weights,
                                p.invar_proportion,
                                d_invar_indices,
                                d_site_scaling_min,
                                d_per_site,
                                p.per_rate_scaling);
                    break;
                default:
                    goto GENERAL; 
            }
            break;
        case 20:
            switch (p.rate_cats) {
                case 1:
                    InnerInnerLikelihoodCalculation_states_20_1_4_8<1,4><<<grid,block,0,stream>>>(
                                p.sites,
                                d_parent_clv,
                                d_child_clv,
                                d_site_scaler,
                                d_p_matrix,
                                d_pattern_weights,
                                p.invar_proportion,
                                d_invar_indices,
                                d_site_scaling_min,
                                d_per_site,
                                p.per_rate_scaling);
                    break;
                case 4:
                    InnerInnerLikelihoodCalculation_states_20_1_4_8<4,4><<<grid,block,0,stream>>>(
                                p.sites,
                                d_parent_clv,
                                d_child_clv,
                                d_site_scaler,
                                d_p_matrix,
                                d_pattern_weights,
                                p.invar_proportion,
                                d_invar_indices,
                                d_site_scaling_min,
                                d_per_site,
                                p.per_rate_scaling);
                    break;
                case 8:
                    InnerInnerLikelihoodCalculation_states_20_1_4_8<8,4><<<grid,block,0,stream>>>(
                                p.sites,
                                d_parent_clv,
                                d_child_clv,
                                d_site_scaler,
                                d_p_matrix,
                                d_pattern_weights,
                                p.invar_proportion,
                                d_invar_indices,
                                d_site_scaling_min,
                                d_per_site,
                                p.per_rate_scaling);
                    break;
                default:
                    goto GENERAL; 
            }
            break;
        default:
            GENERAL:
                InnerInnerLikelihoodCalculation<<<grid, block, 0, stream>>>(
                    p.sites,
                    p.states,
                    p.rate_cats,
                    d_parent_clv,
                    d_child_clv,
                    d_site_scaler,
                    d_p_matrix,
                    d_pattern_weights,
                    d_invar_indices,
                    p.invar_proportion,
                    d_site_scaling_min,
                    d_per_site,
                    p.per_rate_scaling);
                break;
    }
    err = cudaGetLastError();
    if (err != cudaSuccess)
    {   
        fprintf(stderr, "[Gpu_ERROR] device_reduce_sum failed: %s\n",cudaGetErrorString(err));
        exit(1);
    }
    device_reduce_sum(d_per_site, p.sites, d_out_sum);
    err = cudaGetLastError();
    if (err != cudaSuccess)
    {   
        fprintf(stderr, "[Gpu_ERROR] device_reduce_sum failed: %s\n",cudaGetErrorString(err));
        exit(1);
    }

    err = cudaMemcpy(&total_loglik, d_out_sum, sizeof(double), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "[Gpu_ERROR] d_out_sum copy failed: %s\n",cudaGetErrorString(err));
        exit(1);
    }
}

void core_likelihood::Likelihood_Inner_Inner::PrintLikelihood() const
{
    printf("[Likelihood_Inner_Inner]  Total log-likelihood = %.12f \n", total_loglik);
}

void core_likelihood::Likelihood_Inner_Inner::CleanUp()
{
    cudaFree(d_per_site);
    cudaFree(d_out_sum);
    cudaFree(d_parent_clv);
    cudaFree(d_parent_scaler);
    cudaFree(d_child_clv);
    cudaFree(d_child_scaler);
    cudaFree(d_p_matrix);
    cudaFree(d_pattern_weights);
    cudaFree(d_invar_indices);
}

double core_likelihood::ComputeRootLogLikelihoodFromDevice(
    const double* d_root_clv,
    std::size_t sites,
    int states,
    int rate_cats,
    double invar_proportion,
    const double* h_frequencies,
    const double* h_rate_weights,
    const unsigned int* d_pattern_weights,
    const unsigned int* d_site_scaler,
    const int* d_invar_indices,
    cudaStream_t stream)
{
    if (!d_root_clv) {
        throw std::runtime_error("Null root CLV pointer.");
    }
    if (!h_frequencies || !h_rate_weights) {
        throw std::runtime_error("Frequencies and rate weights must be provided.");
    }
    if (states <= 0 || states > MAX_STATES) {
        throw std::runtime_error("States out of supported range.");
    }
    if (rate_cats <= 0 || rate_cats > MAX_RATECATS) {
        throw std::runtime_error("Rate categories out of supported range.");
    }

    cudaError_t err = cudaMemcpyToSymbol(core_likelihood::c_frequencies, h_frequencies,
                                         sizeof(double) * states);
    if (err != cudaSuccess) {
        throw std::runtime_error(std::string("Failed to copy frequencies to constant memory: ") + cudaGetErrorString(err));
    }
    err = cudaMemcpyToSymbol(core_likelihood::c_rate_weights, h_rate_weights,
                             sizeof(double) * rate_cats);
    if (err != cudaSuccess) {
        throw std::runtime_error(std::string("Failed to copy rate weights to constant memory: ") + cudaGetErrorString(err));
    }

    double* d_per_site = nullptr;
    err = cudaMalloc((void**)&d_per_site, sizeof(double) * sites);
    if (err != cudaSuccess) {
        throw std::runtime_error(std::string("d_per_site allocation failed: ") + cudaGetErrorString(err));
    }
    err = cudaMemset(d_per_site, 0, sizeof(double) * sites);
    if (err != cudaSuccess) {
        cudaFree(d_per_site);
        throw std::runtime_error(std::string("d_per_site memset failed: ") + cudaGetErrorString(err));
    }

    double* d_out_sum = nullptr;
    err = cudaMalloc((void**)&d_out_sum, sizeof(double));
    if (err != cudaSuccess) {
        cudaFree(d_per_site);
        throw std::runtime_error(std::string("d_out_sum allocation failed: ") + cudaGetErrorString(err));
    }
    err = cudaMemset(d_out_sum, 0, sizeof(double));
    if (err != cudaSuccess) {
        cudaFree(d_per_site);
        cudaFree(d_out_sum);
        throw std::runtime_error(std::string("d_out_sum memset failed: ") + cudaGetErrorString(err));
    }

    dim3 block(256);
    dim3 grid((sites + block.x*4 - 1)/(block.x*4));

    switch(states) {
        case 4:
            switch (rate_cats) {
                case 1:
                    RootLikelihoodCalculation_states_4_1_4_8<1,4><<<grid,block,0,stream>>>(
                                sites,
                                d_root_clv,
                                d_pattern_weights,
                                d_site_scaler,
                                d_invar_indices,
                                invar_proportion,
                                d_per_site);
                    break;
                case 4:
                    RootLikelihoodCalculation_states_4_1_4_8<4,4><<<grid,block,0,stream>>>(
                                sites,
                                d_root_clv,
                                d_pattern_weights,
                                d_site_scaler,
                                d_invar_indices,
                                invar_proportion,
                                d_per_site);
                    break;
                case 8:
                    RootLikelihoodCalculation_states_4_1_4_8<8,4><<<grid,block,0,stream>>>(
                                sites,
                                d_root_clv,
                                d_pattern_weights,
                                d_site_scaler,
                                d_invar_indices,
                                invar_proportion,
                                d_per_site);
                    break;
                default:
                    goto GENERAL_ROOT;
            }
            break;
        case 5:
            switch (rate_cats) {
                case 1:
                    RootLikelihoodCalculation_states_5_1_4_8<1,4><<<grid,block,0,stream>>>(
                                sites,
                                d_root_clv,
                                d_pattern_weights,
                                d_site_scaler,
                                d_invar_indices,
                                invar_proportion,
                                d_per_site);
                    break;
                case 4:
                    RootLikelihoodCalculation_states_5_1_4_8<4,4><<<grid,block,0,stream>>>(
                                sites,
                                d_root_clv,
                                d_pattern_weights,
                                d_site_scaler,
                                d_invar_indices,
                                invar_proportion,
                                d_per_site);
                    break;
                case 8:
                    RootLikelihoodCalculation_states_5_1_4_8<8,4><<<grid,block,0,stream>>>(
                                sites,
                                d_root_clv,
                                d_pattern_weights,
                                d_site_scaler,
                                d_invar_indices,
                                invar_proportion,
                                d_per_site);
                    break;
                default:
                    goto GENERAL_ROOT;
            }
            break;
        case 20:
            switch (rate_cats) {
                case 1:
                    RootLikelihoodCalculation_states_20_1_4_8<1,4><<<grid,block,0,stream>>>(
                                sites,
                                d_root_clv,
                                d_pattern_weights,
                                d_site_scaler,
                                d_invar_indices,
                                invar_proportion,
                                d_per_site);
                    break;
                case 4:
                    RootLikelihoodCalculation_states_20_1_4_8<4,4><<<grid,block,0,stream>>>(
                                sites,
                                d_root_clv,
                                d_pattern_weights,
                                d_site_scaler,
                                d_invar_indices,
                                invar_proportion,
                                d_per_site);
                    break;
                case 8:
                    RootLikelihoodCalculation_states_20_1_4_8<8,4><<<grid,block,0,stream>>>(
                                sites,
                                d_root_clv,
                                d_pattern_weights,
                                d_site_scaler,
                                d_invar_indices,
                                invar_proportion,
                                d_per_site);
                    break;
                default:
                    goto GENERAL_ROOT;
            }
            break;
        default:
        GENERAL_ROOT:
            RootLikelihoodCalculation<<<grid,block,0,stream>>>(
                                sites,
                                states,
                                rate_cats,
                                d_root_clv,
                                d_pattern_weights,
                                d_site_scaler,
                                d_invar_indices,
                                invar_proportion,
                                d_per_site);
            break;
    }

    err = cudaGetLastError();
    if (err != cudaSuccess) {
        cudaFree(d_per_site);
        cudaFree(d_out_sum);
        throw std::runtime_error(std::string("Root likelihood kernel launch failed: ") + cudaGetErrorString(err));
    }

    device_reduce_sum(d_per_site, sites, d_out_sum);

    double total = 0.0;
    err = cudaMemcpy(&total, d_out_sum, sizeof(double), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        cudaFree(d_per_site);
        cudaFree(d_out_sum);
        throw std::runtime_error(std::string("Failed to copy total log-likelihood to host: ") + cudaGetErrorString(err));
    }

    cudaFree(d_per_site);
    cudaFree(d_out_sum);
    return total;
}
