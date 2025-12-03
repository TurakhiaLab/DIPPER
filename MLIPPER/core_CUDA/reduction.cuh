#pragma once
#include <cuda_runtime.h>




static __inline__ __device__
double warp_reduce_sum(double val) {
    for (int offset = 16; offset > 0; offset >>= 1) {
        val += __shfl_down_sync(0xffffffff, val, offset);
    }
    return val;
}

static __inline__ __device__
double block_reduce_sum(double val) {
    
    static __shared__ double shared[8]; // Assumes max 1024 threads per block
    int lane = threadIdx.x & 31;
    int wid  = threadIdx.x >> 5;

    val = warp_reduce_sum(val);
    if (lane == 0) shared[wid] = val;
    __syncthreads();

    val = (threadIdx.x < (blockDim.x >> 5)) ? shared[lane] : 0.0;
    if (wid == 0) {
        val = warp_reduce_sum(val);
    }
    return val;
}
