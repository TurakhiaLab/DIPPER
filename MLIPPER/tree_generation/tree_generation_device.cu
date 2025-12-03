#pragma once
#include <cuda_runtime.h>
#include <vector>
#include <string>
#include <unordered_map>
#include <stdexcept>
#include <cstdint>
#include "../mlipper_util.h"
#include "tree.hpp" 
#include "pmat.h"
#include "root_likelihood.cuh"
#include "partial_likelihood.cuh"

// forward decl for device likelihood wrapper
double EvaluateTreeLogLikelihood_device(
    const DeviceTree&      D,
    const TreeBuildResult& T,
    const HostPacking&     H,
    const std::vector<double>& pi,
    const std::vector<double>& rate_weights,
    cudaStream_t stream = 0);


// Build transition probability matrices for each node/rate category on host.
void fill_pmats_in_host_packing(
    const TreeBuildResult&       T,
    HostPacking&                 H,
    const EigResult&             er,
    const std::vector<double>&   pi,           // len = states
    const std::vector<double>&   rate_rates,   // len = rate_cats (per-category rate multipliers)
    int states,
    int rate_cats)
{
    const int N = (int)T.nodes.size();
    const size_t per_node = (size_t)rate_cats * states * states;

    H.pmats.assign((size_t)N * per_node, 0.0);

    for (int nid = 0; nid < N; ++nid) {
        const TreeNode& nd = T.nodes[nid];
        if (nd.parent < 0) {
            // root 通常不需要往上的 P 矩陣
            continue;
        }

        double* base = H.pmats.data() + (size_t)nid * per_node;
        double blen  = nd.branch_length_to_parent;

        for (int rc = 0; rc < rate_cats; ++rc) {
            double r = rate_rates[rc];        // rate category multiplier
            double t = blen;                  // branch length
            double p = 0;                   // 先假設沒有 invar 調整

            double* P = base + (size_t)rc * states * states;

            pmatrix_from_triple(
                er.Vinv.data(), er.V.data(), er.lambdas.data(),
                            r, t, p, P, states);
        }
    }
}

// Pack topology and tip encodings from tree/MSA into HostPacking.
HostPacking pack_host_arrays_from_tree_and_msa(
        const TreeBuildResult& T,
        const std::vector<std::string>& msa_tip_names,  // len = tips
        const std::vector<std::string>& msa_rows,       // len = tips, each row length = sites
        size_t sites,
        int states)
{
    if (msa_rows.size() != msa_tip_names.size())
        throw std::runtime_error("MSA rows/names size mismatch.");
    if (msa_rows.empty()) throw std::runtime_error("Empty MSA.");
    if (msa_rows[0].size() != sites)
        throw std::runtime_error("Sites mismatch.");

    const int N = (int)T.nodes.size();

    // 拓樸
    HostPacking H;
    H.postorder = T.postorder;
    H.parent.resize(N, -1);
    H.left.resize(N, -1);
    H.right.resize(N, -1);
    H.is_tip.resize(N, 0);
    H.blen.resize(N, 0.0);
    H.clv_offsets.resize(N, 0);

    for (int i = 0; i < N; ++i) {
        const auto& nd = T.nodes[i];
        H.parent[i] = nd.parent;
        H.left[i]   = nd.left;
        H.right[i]  = nd.right;
        H.is_tip[i] = nd.is_tip ? 1 : 0;
        H.blen[i]   = nd.branch_length_to_parent;
        H.clv_offsets[i] = (uint64_t)nd.clv_offset; // elements-based
    }
    
    // 建立 tip 名稱 -> MSA row 的查表
    std::unordered_map<std::string,int> name2row;
    name2row.reserve(msa_tip_names.size()*2);
    for (int r = 0; r < (int)msa_tip_names.size(); ++r) name2row[msa_tip_names[r]] = r;

    // 收集 tip 在後序中的順序，以及對應 node id
    std::vector<int> tip_node_ids_host;
    tip_node_ids_host.reserve(msa_tip_names.size());
    for (int id : T.postorder) {
        if (T.nodes[id].is_tip) tip_node_ids_host.push_back(id);
    }
    const int tips = (int)tip_node_ids_host.size();

    if (tips != (int)msa_tip_names.size())
        throw std::runtime_error("Tip count in tree != MSA names.");

    H.tip_node_ids = tip_node_ids_host;

    // 依「tipIndex」(0..tips-1) 的固定順序，把對應 MSA row 的序列編碼成 tipchars
    H.tipchars.resize((size_t)tips * sites);

    if (states == 4 || states == 5) {
        for (int t = 0; t < tips; ++t) {
            const int node_id = tip_node_ids_host[t];
            const auto& name  = T.nodes[node_id].name;
            auto it = name2row.find(name);
            if (it == name2row.end())
                throw std::runtime_error("Tip not found in MSA: " + name);
            const std::string& row = msa_rows[it->second];
            for (size_t s = 0; s < sites; ++s)
                H.tipchars[(size_t)t * sites + s] = encode_state_DNA5(row[s]);
        }
    } 
    // else if (states == 20) {
    //     auto encAA = [](char c)->uint8_t {
    //         // TODO: 依你的 protein 映射補齊；暫時把未知當 19
    //         switch (c) {
    //             case 'A': case 'a': return 0;
    //             // ... 其餘 18 種
    //             default: return 19;
    //         }
    //     };
    //     for (int t = 0; t < tips; ++t) {
    //         const int node_id = tip_node_ids_host[t];
    //         const auto& name  = T.nodes[node_id].name;
    //         auto it = name2row.find(name);
    //         if (it == name2row.end())
    //             throw std::runtime_error("Tip not found in MSA: " + name);
    //         const std::string& row = msa_rows[it->second];
    //         for (size_t s = 0; s < sites; ++s)
    //             H.tipchars[(size_t)t * sites + s] = encAA(row[s]);
    //     }
    // } else {
    //     throw std::runtime_error("Unsupported states; expect 4/5/20.");
    // }

    return H;
}

// ===== 把 Host 打包拷到 GPU，建立 DeviceTree =====
// Upload host packing and model parameters to GPU, constructing DeviceTree.
DeviceTree upload_to_gpu(
        const TreeBuildResult& T,
        const HostPacking& H,
        const EigResult& er,
        const std::vector<double>& rate_weights,
        size_t sites, int states, int rate_cats, bool per_rate_scaling)
{
    DeviceTree D;
    D.N = (int)T.nodes.size();
    D.tips = (int)H.tip_node_ids.size();
    D.inners = D.N - D.tips;
    D.sites = sites;
    D.states = states;
    D.rate_cats = rate_cats;
    D.per_rate_scaling = per_rate_scaling;
    
    // --- alloc topology ---
    CUDA_CHECK(cudaMalloc(&D.d_lambdas, sizeof(double) * D.states));
    CUDA_CHECK(cudaMalloc(&D.d_V,       sizeof(double) * D.states * D.states));
    CUDA_CHECK(cudaMalloc(&D.d_Vinv,    sizeof(double) * D.states * D.states));

    CUDA_CHECK(cudaMalloc(&D.d_U,       sizeof(double) * D.states * D.states));
    CUDA_CHECK(cudaMalloc(&D.d_rate_w,  sizeof(double) * D.rate_cats));

    CUDA_CHECK(cudaMemcpy(D.d_lambdas, er.lambdas.data(), sizeof(double) * D.states, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(D.d_V,       er.V.data(),       sizeof(double) * D.states * D.states, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(D.d_Vinv,    er.Vinv.data(),    sizeof(double) * D.states * D.states, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(D.d_U,       er.U.data(),       sizeof(double) * D.states * D.states, cudaMemcpyHostToDevice));
    if ((int)rate_weights.size() != rate_cats) {
        throw std::runtime_error("rate_weights size mismatch.");
    }
    CUDA_CHECK(cudaMemcpy(D.d_rate_w,  rate_weights.data(), sizeof(double) * D.rate_cats, cudaMemcpyHostToDevice));

    CUDA_CHECK(cudaMalloc(&D.d_postorder, sizeof(int) * D.N));
    CUDA_CHECK(cudaMalloc(&D.d_parent,    sizeof(int) * D.N));
    CUDA_CHECK(cudaMalloc(&D.d_left,      sizeof(int) * D.N));
    CUDA_CHECK(cudaMalloc(&D.d_right,     sizeof(int) * D.N));
    CUDA_CHECK(cudaMalloc(&D.d_is_tip,    sizeof(uint8_t) * D.N));
    CUDA_CHECK(cudaMalloc(&D.d_blen,      sizeof(double) * D.N));
    

    CUDA_CHECK(cudaMemcpy(D.d_postorder, H.postorder.data(), sizeof(int) * D.N, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(D.d_parent,    H.parent.data(),    sizeof(int) * D.N, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(D.d_left,      H.left.data(),      sizeof(int) * D.N, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(D.d_right,     H.right.data(),     sizeof(int) * D.N, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(D.d_is_tip,    H.is_tip.data(),    sizeof(uint8_t) * D.N, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(D.d_blen,      H.blen.data(),      sizeof(double) * D.N, cudaMemcpyHostToDevice));

    // --- tips ---
    CUDA_CHECK(cudaMalloc(&D.d_tip_node_ids, sizeof(int) * D.tips));
    CUDA_CHECK(cudaMemcpy(D.d_tip_node_ids, H.tip_node_ids.data(), sizeof(int) * D.tips, cudaMemcpyHostToDevice));

    CUDA_CHECK(cudaMalloc(&D.d_tipchars, sizeof(uint8_t) * (size_t)D.tips * sites));
    CUDA_CHECK(cudaMemcpy(D.d_tipchars, H.tipchars.data(), sizeof(uint8_t) * (size_t)D.tips * sites, cudaMemcpyHostToDevice));

    const unsigned int tipmap_size = D.states + 1;
    // Lookup table for tip-tip was used previously; the current kernels compute
    // tip-tip on the fly, so we skip allocating d_lookup to save memory/time.
    D.d_lookup = nullptr;

    // --- CLV & Offsets ---
    const size_t per_node = (size_t)sites * (size_t)rate_cats * (size_t)states;
    const size_t clv_elems = (size_t)D.N * per_node;
    CUDA_CHECK(cudaMalloc(&D.d_clv_pool,  sizeof(double) * clv_elems));
    CUDA_CHECK(cudaMalloc(&D.d_clv_offsets, sizeof(uint64_t) * D.N));
    CUDA_CHECK(cudaMemcpy(D.d_clv_offsets, H.clv_offsets.data(), sizeof(uint64_t) * D.N, cudaMemcpyHostToDevice));
    
    // 建議把 CLV 池清零（或留給 kernel 覆蓋）
    CUDA_CHECK(cudaMemset(D.d_clv_pool, 0, sizeof(double) * clv_elems));

    CUDA_CHECK(cudaMalloc(&D.d_pmat, sizeof(double) * D.N * D.rate_cats * D.states * D.states));
    CUDA_CHECK(cudaMemcpy(D.d_pmat, H.pmats.data(),sizeof(double) * D.N * D.rate_cats * D.states * D.states, cudaMemcpyHostToDevice));
    // --- scaler ---
    const size_t scaler_len = per_rate_scaling ? (sites * (size_t)rate_cats) : sites;
    CUDA_CHECK(cudaMalloc(&D.d_site_scaler, sizeof(unsigned) * scaler_len));
    // 初始化 0
    CUDA_CHECK(cudaMemset(D.d_site_scaler, 0, sizeof(unsigned) * scaler_len));

    std::vector<unsigned int> tipmap(tipmap_size);
    for (unsigned int j = 0; j < tipmap_size; ++j) {
        if(j == D.states){
            tipmap[j] = 15;  // 0 -> 0001, 1 -> 0010, 2 -> 0100, ...
        }
        else{
            tipmap[j] = 1u << j;  // 0 -> 0001, 1 -> 0010, 2 -> 0100, ...
        }
    }
    CUDA_CHECK(cudaMalloc(&D.d_tipmap, sizeof(unsigned) * tipmap_size));
    CUDA_CHECK(cudaMemcpy(D.d_tipmap, tipmap.data(),  tipmap_size * sizeof(unsigned int), cudaMemcpyHostToDevice));

    return D;
}

// ===== 釋放 GPU 資源 =====
// Release all device buffers held by DeviceTree.
void free_device_tree(DeviceTree& D) {
    auto F = [](void* p){ if(p) cudaFree(p); };
    F(D.d_postorder); F(D.d_parent); F(D.d_left); F(D.d_right);
    F(D.d_is_tip);    F(D.d_blen);
    F(D.d_tip_node_ids); F(D.d_tipchars);
    F(D.d_clv_pool);  F(D.d_clv_offsets);
    F(D.d_site_scaler);
    F(D.d_lambdas); F(D.d_V); F(D.d_Vinv); F(D.d_U); F(D.d_rate_w);
    F(D.d_pmat); F(D.d_lookup); F(D.d_tipmap);
    D = DeviceTree{};
}



// Initialize tip conditional likelihood vectors on device from encoded tipchars.
__global__ void InitTipCLVKernel(DeviceTree D)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total = D.tips * D.sites;
    if (idx >= total) return;
    int tip_idx = idx / D.sites;
    int site    = idx % D.sites;

    int node_id = D.d_tip_node_ids[tip_idx];

    uint8_t encoded = D.d_tipchars[(size_t)tip_idx * D.sites + site];

    double* clv = D.d_clv_pool + D.d_clv_offsets[node_id];

    // 每個 site 的起點
    size_t stride_site  = (size_t)D.rate_cats * D.states;
    size_t base_site    = (size_t)site * stride_site;

    for (int rc = 0; rc < D.rate_cats; ++rc) {
        size_t base_rc = base_site + (size_t)rc * D.states;

        if (encoded < 4) {
            // 明確 A/C/G/T: one-hot
            for (int s = 0; s < D.states; ++s)
                clv[base_rc + s] = (s == encoded) ? 1.0 : 0.0;
        } else {
            // 第 5 個狀態 (N, gap...)：暫時給均勻
            double val = 1.0; // 或者用 pi[s]
            for (int s = 0; s < D.states; ++s)
                clv[base_rc + s] = val;
        }
    }
}

// Launch tip CLV initialization kernel.
void launch_init_tip_clv(const DeviceTree& D)
{
    dim3 block(128);
    dim3 grid( (D.tips * D.sites + block.x - 1) / block.x );

    InitTipCLVKernel<<<grid, block>>>(D);

    CUDA_CHECK(cudaGetLastError());
}

// Compute root log-likelihood on device given a root id and model vectors.
double eval_root_loglikelihood(
    const DeviceTree& D,
    int root_id,
    const std::vector<double>& pi,
    const std::vector<double>& rate_weights,
    cudaStream_t stream)
{
    if (root_id < 0 || root_id >= D.N) {
        throw std::runtime_error("Invalid root id.");
    }
    if ((int)pi.size() != D.states) {
        throw std::runtime_error("pi vector size mismatch.");
    }
    if ((int)rate_weights.size() != D.rate_cats) {
        throw std::runtime_error("rate weights size mismatch.");
    }
    const size_t per_node = D.per_node_elems();
    const double* d_root_clv = D.d_clv_pool + (size_t)root_id * per_node;
    return root_likelihood::ComputeRootLogLikelihood(
        d_root_clv,
        D.sites,
        D.states,
        D.rate_cats,
        pi.data(),
        rate_weights.data(),
        D.d_site_scaler,
        stream,
        nullptr,          // d_pattern_w
        nullptr,          // d_invar_indices
        0.0               // invar_proportion
    );
}

// Assemble a NodeOpInfo record for device execution.
static inline NodeOpInfo make_node_op_device(
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

// Evaluate full tree log-likelihood on device using prepared DeviceTree.
double EvaluateTreeLogLikelihood_device(
    const DeviceTree&      D,
    const TreeBuildResult& T,
    const HostPacking&     H,
    const std::vector<double>& pi,
    const std::vector<double>& rate_weights,
    cudaStream_t stream)
{
    const int N = (int)T.nodes.size();
    std::vector<int> node2tip(N, -1);

    for (int tip_idx = 0; tip_idx < (int)H.tip_node_ids.size(); ++tip_idx) {
        int nid = H.tip_node_ids[tip_idx];
        node2tip[nid] = tip_idx;
    }

    std::vector<NodeOpInfo> ops_host;
    ops_host.reserve(N);

    const unsigned int tipmap_size = (unsigned int)(D.states + 1);
    const int log2_stride = ceil_log2_u32(tipmap_size);

    for (int nid : T.postorder) {
        const TreeNode& nd = T.nodes[nid];
        if (nd.is_tip) continue;

        const int L = nd.left;
        const int R = nd.right;
        const bool tipL = T.nodes[L].is_tip;
        const bool tipR = T.nodes[R].is_tip;

        if (tipL && tipR) {
            int left_tip_idx  = node2tip[L];
            int right_tip_idx = node2tip[R];
            NodeOpInfo op = make_node_op_device(
                nid,
                L,
                R,
                left_tip_idx,
                right_tip_idx,
                OP_TIP_TIP,
                log2_stride);
            ops_host.push_back(op);
        } else if (tipL && !tipR) {
            int left_tip_idx = node2tip[L];
            NodeOpInfo op = make_node_op_device(
                nid,
                L,
                R,
                left_tip_idx,
                -1,
                OP_TIP_INNER);
            ops_host.push_back(op);
        } else if (!tipL && tipR) {
            int right_tip_idx = node2tip[R];
            NodeOpInfo op = make_node_op_device(
                nid,
                L,
                R,
                -1,
                right_tip_idx,
                OP_TIP_INNER);
            ops_host.push_back(op);
        } else {
            NodeOpInfo op = make_node_op_device(
                nid,
                L,
                R,
                -1,
                -1,
                OP_INNER_INNER);
            ops_host.push_back(op);
        }
    }

    NodeOpInfo* d_ops = nullptr;
    const int num_ops = (int)ops_host.size();
    if (num_ops > 0) {
        CUDA_CHECK(cudaMalloc(&d_ops, sizeof(NodeOpInfo) * (size_t)num_ops));
        CUDA_CHECK(cudaMemcpyAsync(
            d_ops,
            ops_host.data(),
            sizeof(NodeOpInfo) * (size_t)num_ops,
            cudaMemcpyHostToDevice,
            stream));
    }
    
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    int num_sms = prop.multiProcessorCount;

    const int block = 512;
    
    int blocks_per_sm_target = 2; 
    int max_blocks = num_sms * blocks_per_sm_target;
    const int grid = (D.sites + block - 1) / block;
    if(grid > max_blocks){
        const int grid = max_blocks;
    }

    if (num_ops > 0) {
        Rtree_Likelihood_Site_Parallel_Kernel<<<grid, block, 0, stream>>>(
            D,
            d_ops,
            num_ops);
        CUDA_CHECK(cudaGetLastError());
    }

    if (d_ops) CUDA_CHECK(cudaFree(d_ops));

    return eval_root_loglikelihood(D, T.root_id, pi, rate_weights, stream);
}
