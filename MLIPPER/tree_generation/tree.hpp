#pragma once
#include <cuda_runtime.h>
#include <vector>
#include <string>
#include <unordered_map>
#include <stdexcept>
#include <cstddef>
#include <libpll/pll.h>

#include "../pmatrix/pmat.h"


// ===== CUDA error helper =====
#ifndef CUDA_CHECK
#define CUDA_CHECK(expr) do { \
    cudaError_t _e = (expr);  \
    if (_e != cudaSuccess) {  \
      throw std::runtime_error(std::string("[CUDA] ") + cudaGetErrorString(_e)); \
    }                         \
} while(0)
#endif

inline uint8_t encode_state_DNA5(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't':
        case 'U': case 'u': return 3;
        case '-':           return 4; // gap
        default:            return 4; // 先當 gap；若要 IUPAC 可改成 bitmask 方案
    }
}


// 你的 GPU 佈局假設：每個 node 的 CLV 大小一樣（sites * rate_cats * states）
struct TreeNode {
    int   id = -1;
    bool  is_tip = false;
    int   left = -1;    // child id
    int   right = -1;   // child id
    int   parent = -1;  // parent id (root = -1)
    double branch_length_to_parent = 0.0; // Newick 冒號後的數字
    std::string name;   // tip 才會有
    // GPU 偏移：你可以直接用這兩個偏移去 d_clv_pool / d_scaler_pool
    size_t clv_offset = 0;     // elements-based offset (not bytes)
    size_t scaler_offset = 0;  // elements-based offset (per-site 或 per-rate 視你需求)
};

struct TreeBuildResult {
    std::vector<TreeNode> nodes;     // 0..N-1
    int root_id = -1;
    std::vector<int> postorder;      // 後序排列（children -> parent）
    // 也提供 tip name -> node id 對應（方便把 MSA tip CLV 放上對應節點）
    std::unordered_map<std::string,int> tip_node_by_name;
};

TreeBuildResult build_tree_from_newick_with_pll(
    const std::vector<std::string>& msa_tip_names,
    const std::string& newick_text,
    size_t sites,
    int states,
    int rate_cats,
    bool per_rate_scaling);

struct DeviceTree {
    // sizes
    int     N = 0;            // nodes = tips + inners
    int     tips = 0;
    int     inners = 0;
    size_t  sites = 0;
    int     states = 0;
    int     rate_cats = 0;
    bool    per_rate_scaling = false;

    double *d_lambdas  = nullptr;  // [states]
    double *d_V        = nullptr;  // [states*states]  row-major
    double *d_Vinv     = nullptr;  // [states*states]  row-major
    double *d_U        = nullptr;  // [states*states]  (可選；若內核需要 U 直接使用)
    double *d_rate_w   = nullptr;  // [rate_cats]      離散Gamma或其他 rate 類別

    // topology (device)
    int    *d_postorder = nullptr; // [N]
    int    *d_parent    = nullptr; // [N]
    int    *d_left      = nullptr; // [N]
    int    *d_right     = nullptr; // [N]
    uint8_t*d_is_tip    = nullptr; // [N]
    double *d_blen      = nullptr; // [N] branch length to parent

    // tips
    // tip index 0..tips-1 的順序 = 依 TreeBuildResult.nodes 的後序遍歷中遇到的 tips 順序
    // d_tip_node_ids: tipIndex -> nodeId（方便 kernel 需要 nodeId 時查）
    int     *d_tip_node_ids = nullptr;          // [tips]
    uint8_t *d_tipchars     = nullptr;          // [tips * sites], DNA5: 0..4 (含 gap)

    // CLV pool for ALL nodes (internal will be produced; tips usually unused)
    // layout: contiguous by node: node i at [i * per_node_elems ..)
    double  *d_clv_pool     = nullptr;          // [N * sites * rate_cats * states]
    uint64_t*d_clv_offsets  = nullptr;          // [N], elements-based offset (uint64 for device)

    // scaler pool (site or site*rate)
    unsigned *d_site_scaler = nullptr;          // [sites] or [sites*rate_cats]

    double* d_pmat = nullptr;
    double* d_lookup = nullptr;
    unsigned int* d_tipmap = nullptr;

    // convenience
    size_t per_node_elems() const {
        return sites * (size_t)rate_cats * (size_t)states;
    }
    size_t clv_pool_elems() const {
        return (size_t)N * per_node_elems();
    }
    size_t scaler_elems() const {
        return per_rate_scaling ? (sites * (size_t)rate_cats) : sites;
    }
    size_t pmat_per_node_elems() const {
        return (size_t)rate_cats * (size_t)states * (size_t)states;
    }
};

struct HostPacking {
    // 給 cudaMemcpy 的 host 端暫存
    std::vector<int>     postorder, parent, left, right;
    std::vector<uint8_t> is_tip;
    std::vector<double>  blen;
    std::vector<int>     tip_node_ids;      // size = tips
    std::vector<uint8_t> tipchars;          // size = tips * sites

    std::vector<uint64_t> clv_offsets;      // size = N
    std::vector<unsigned> site_scaler;      // size = sites or sites*rate
    std::vector<double>   pmats;
};


HostPacking pack_host_arrays_from_tree_and_msa(
    const TreeBuildResult& T,
    const std::vector<std::string>& msa_tip_names,  // len = tips
    const std::vector<std::string>& msa_rows,       // len = tips, each row length = sites
    size_t sites,
    int states
);
void fill_pmats_in_host_packing(
    const TreeBuildResult&       T,
    HostPacking&                 H,
    const EigResult&             er,
    const std::vector<double>&   pi,           // len = states
    const std::vector<double>&   rate_rates,   // len = rate_cats (per-category rate multipliers)
    int states,
    int rate_cats
);

DeviceTree upload_to_gpu(
    const TreeBuildResult& T,
    const HostPacking& H,
    const EigResult& er,
    const std::vector<double>& rate_weights,
    size_t sites, int states, int rate_cats, bool per_rate_scaling
);

struct BuildToGpuResult{
    DeviceTree dev;
    TreeBuildResult tree;
    HostPacking hostPack;
};

BuildToGpuResult BuildAllToGPU(
    const std::vector<std::string>& msa_tip_names,
    const std::vector<std::string>& msa_rows,
    const std::string& newick_text,
    const std::vector<double>& Q_rowmajor,   // size = states*states
    const std::vector<double>& pi,           // size = states
    const std::vector<double>& rate_rates,   // size = rate_cats
    const std::vector<double>& rate_weights, // size = rate_cats
    size_t sites, int states, int rate_cats, bool per_rate_scaling);

void free_device_tree(DeviceTree& D);
static void throw_if(bool cond, const char* msg);

void launch_init_tip_clv(const DeviceTree& D);

double eval_root_loglikelihood(
    const DeviceTree& D,
    int root_id,
    const std::vector<double>& pi,
    const std::vector<double>& rate_weights,
    cudaStream_t stream
);

double EvaluateTreeLogLikelihood(
    const DeviceTree&      D,
    const TreeBuildResult& T,
    const HostPacking&     H,
    const std::vector<double>& pi,
    const std::vector<double>& rate_weights,
    cudaStream_t stream
);
