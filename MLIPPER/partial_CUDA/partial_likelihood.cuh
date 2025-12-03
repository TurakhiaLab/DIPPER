#pragma once
#ifndef PARTIAL_LIKELIHOOD_CUH
#define PARTIAL_LIKELIHOOD_CUH
#include "tree.hpp"
#include <cuda_runtime.h>
#include <cstddef>
#include <cstdint>
#include "../mlipper_util.h"
#define SCALE_THRESHOLD_EXPONENT -64

enum NodeOpType : int {
    OP_TIP_TIP = 0,
    OP_TIP_INNER = 1,
    OP_INNER_INNER = 2
};

struct NodeOpInfo {
    int parent_id = -1;
    int left_id = -1;
    int right_id = -1;
    int left_tip_index = -1;
    int right_tip_index = -1;
    int op_type = OP_TIP_TIP;
    int log2_stride = 0;
};



namespace partial_likelihood {

    struct Param {
        std::size_t sites;
        int states;
        int rate_cats;
        bool rate_scaling;
        bool per_rate_scaling;
        
        Param() : sites(0), states(0), rate_cats(0), rate_scaling(false), per_rate_scaling(false){};

        Param(std::size_t t_sites, int t_states, int t_rate_cats, bool t_rate_scaling, bool t_per_rate_scaling){
            sites = t_sites;
            states = t_states;
            rate_cats = t_rate_cats;
            rate_scaling = t_rate_scaling;
            per_rate_scaling = t_per_rate_scaling;
        };
    };

    struct Partial_Likelihood_Tip_Tip {

        double*        d_parent_clv      = nullptr;  // [sites * states * rate_cats]
        unsigned int*  d_parent_scaler   = nullptr;  // [sites] or [sites*rate_cats] (optional; zero-initialized if provided)

        const unsigned char* d_left_tipchars   = nullptr; // [sites]  tip 字元編碼（與你現有 pipeline 對齊）
        const unsigned char* d_right_tipchars  = nullptr; // [sites]

        const double*  d_lookup         = nullptr;   // [tipmap_size * tipmap_size * states * rate_cats]
        
        unsigned int tipmap_size = 0;
        unsigned int log2_stride = 0;
 
        bool owns_parent_clv    = false;
        bool owns_parent_scaler = false;
        bool owns_left_chars    = false;
        bool owns_right_chars   = false;
        bool owns_lookup        = false;

        void ConstructionOnGpu(
            const Param& p,
            const unsigned char* h_left_tipchars,   // [sites]
            const unsigned char* h_right_tipchars,  // [sites]
            const double*        h_lookup,          // lookup buffer on host
            std::size_t          lookup_count_doubles,
            cudaStream_t         stream
        );
        void UpdatePartialLikelihood(const Param& p, cudaStream_t stream);
        void CleanUp();
    };

    
    struct Partial_Likelihood_Tip_Inner {

        // device buffers (owned or external)
        const unsigned char* d_left_tipchars   = nullptr; // [sites]
        const unsigned int*  d_tipmap          = nullptr; // [tipmap_size] bitmasks
        const double*        d_right_clv       = nullptr; // [sites * rate_cats * states]
        const double*        d_left_matrix     = nullptr; // [rate_cats * states * states] (row-major: from k -> j)
        const double*        d_right_matrix    = nullptr; // [rate_cats * states * states] (row-major: from j -> k)
        double*              d_parent_clv      = nullptr; // [sites * rate_cats * states]
        unsigned int*        d_parent_scaler   = nullptr; // optional (per-rate)

        // addressing
        unsigned int tipmap_size = 0;


        // ownership flags
        bool owns_left_chars   = false;
        bool owns_tipmap       = false;
        bool owns_right_clv    = false;
        bool owns_left_matrix  = false;
        bool owns_right_matrix = false;
        bool owns_parent_clv   = false;
        bool owns_parent_scaler= false;

        void ConstructionOnGpu(
            const Param& p,
            const unsigned char* h_left_tipchars,   // [sites] or nullptr
            const unsigned int*  h_tipmap,          // [tipmap_size] or nullptr
            const double*        h_right_clv,       // [sites*rate_cats*states] or nullptr
            const double*        h_left_matrix,     // [rate_cats*states*states] or nullptr
            const double*        h_right_matrix,    // [rate_cats*states*states] or nullptr
            cudaStream_t         stream = 0);

        void UpdatePartialLikelihood(const Param& P, cudaStream_t stream = 0) const;

        void CleanUp();
    };

    struct Partial_Likelihood_Inner_Inner {

        // device buffers
        const double* d_left_clv  = nullptr;  // [sites*rate*states]
        const double* d_right_clv = nullptr;  // [sites*rate*states]
        const double* d_left_matrix  = nullptr; // [rate*states*states]
        const double* d_right_matrix = nullptr; // [rate*states*states]
        double*       d_parent_clv   = nullptr; // [sites*rate*states]
        unsigned int* d_parent_scaler = nullptr; // optional

        // ownership
        bool owns_left_clv = false, owns_right_clv = false;
        bool owns_left_matrix = false, owns_right_matrix = false;
        bool owns_parent_clv = false, owns_parent_scaler = false;

        void ConstructionOnGpu(
            const Param& p,
            const double* h_left_clv,      // 可為 nullptr：表示外部已在 GPU
            const double* h_right_clv,
            const double* h_left_matrix,
            const double* h_right_matrix,
            cudaStream_t  stream = 0);

        void UpdatePartialLikelihood(const Param& P, cudaStream_t stream = 0) const;

        void CleanUp();
    };

    void compute_inner_inner(const DeviceTree &D, int parent_id, int left_id, int right_id, cudaStream_t stream);
    void compute_tip_inner_swap(const DeviceTree& D, int parent_id, int tip_node_id, int inner_node_id, int tip_index, cudaStream_t stream);
    void compute_tip_inner(const DeviceTree& D, int parent_id, int tip_node_id, int inner_node_id, int tip_index, cudaStream_t stream);
    void compute_tip_tip(const DeviceTree& D, int parent_id, int left_node_id, int right_node_id, int left_tip_index, int right_tip_index, cudaStream_t stream);
}
__global__ void UpdatePartialTipTipKernel(
    const DeviceTree D,
    const NodeOpInfo* __restrict__ ops,
    int num_ops);

template<int RATE_CATS>
__global__ void UpdatePartialTipTipKernel_states_4_ratecat(
    const DeviceTree D,
    const NodeOpInfo* __restrict__ ops,
    int num_ops);

__global__ void UpdatePartialTipTipKernel_states_4_generic(
    const DeviceTree D,
    const NodeOpInfo* __restrict__ ops,
    int num_ops);

__global__ void UpdatePartialTipInnerKernel(
    const DeviceTree D,
    const NodeOpInfo* __restrict__ ops,
    int num_ops);

template<int RATE_CATS>
__global__ void UpdatePartialTipInnerKernel_states_4_ratecat(
    const DeviceTree D,
    const NodeOpInfo* __restrict__ ops,
    int num_ops);

__global__ void UpdatePartialInnerInnerKernel(
    const DeviceTree D,
    const NodeOpInfo* __restrict__ ops,
    int num_ops);

template<int RATE_CATS>
__global__ void  UpdatePartialInnerInnerKernel_states_4_ratecat(
    const DeviceTree D,
    const NodeOpInfo* __restrict__ ops,
    int num_ops);

__global__ void Rtree_Likelihood_Site_Parallel_Kernel(
    const DeviceTree D,
    const NodeOpInfo* ops,
    int num_ops);
#endif // PARTIAL_LIKELIHOOD_CUH
