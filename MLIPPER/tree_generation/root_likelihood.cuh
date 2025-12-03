#pragma once
#ifndef TREE_GENERATION_ROOT_LIKELIHOOD_CUH
#define TREE_GENERATION_ROOT_LIKELIHOOD_CUH

#include <cstddef>
#include <string>
#include <cuda_runtime.h>
#include <vector>
#include "../mlipper_util.h"

namespace root_likelihood {

struct RootLikelihoodWorkspace {
    double* d_per_site{nullptr};
    std::size_t capacity{0};

    ~RootLikelihoodWorkspace() {
        if (d_per_site) cudaFree(d_per_site);
    }

    void ensure(std::size_t count) {
        if (count <= capacity) return;
        if (d_per_site) cudaFree(d_per_site);
        CHECK_CUDA(cudaMalloc(&d_per_site, sizeof(double) * count));
        capacity = count;
    }

    double* data() const { return d_per_site; }
};

    double ComputeRootLogLikelihood(
        const double* d_root_clv,
        std::size_t sites,
        int states,
        int rate_cats,
        const double* h_frequencies,
        const double* h_rate_weights,
        const unsigned int* d_site_scaler,
        cudaStream_t stream = 0,
        const unsigned* d_pattern_w = nullptr,
        const int* d_invar_indices = nullptr,
        double invar_proportion = 0.0,
        RootLikelihoodWorkspace* workspace = nullptr);
}

#endif // TREE_GENERATION_ROOT_LIKELIHOOD_CUH
