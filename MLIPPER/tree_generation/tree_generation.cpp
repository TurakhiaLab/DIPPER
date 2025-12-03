#pragma once
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <unordered_map>
#include <vector>
#include <string>

#include "../mlipper_util.h"
#include "tree.hpp"
#include "pmat.h"
#include "core_likelihood.cuh"
#include "partial_likelihood.cuh"

double EvaluateTreeLogLikelihood_device(
    const DeviceTree&      D,
    const TreeBuildResult& T,
    const HostPacking&     H,
    const std::vector<double>& pi,
    const std::vector<double>& rate_weights,
    cudaStream_t stream = 0);

static void throw_if(bool cond, const char* msg) {
    if (cond) throw std::runtime_error(msg);
}

// Parse Newick text and build a rooted TreeBuildResult with topology and offsets.
TreeBuildResult build_tree_from_newick_with_pll(
    const std::vector<std::string>& msa_tip_names,
    const std::string& newick_text,
    size_t sites,
    int states,
    int rate_cats,
    bool per_rate_scaling)
{
    TreeBuildResult out;
    // 1) 用 libpll 解析 Newick（rooted tree）
    const char* filename = "./tree.nwk";
    {
        std::ofstream ofs(filename);
        ofs << newick_text;
    }

    pll_rtree_t* rtree = pll_rtree_parse_newick(filename);

    throw_if(!rtree, "pll_rtree_parse_newick failed (check Newick syntax).");

    // 2) 收集所有節點（libpll 會提供 nodes 陣列與 tips/inner 數量）
    const unsigned int num_tips   = rtree->tip_count;
    const unsigned int num_inners = rtree->inner_count;
    const unsigned int num_nodes  = num_tips + num_inners; // rooted：含 root 在內

    out.nodes.resize(num_nodes);

    // 3) 名稱對齊：MSA tip name -> index
    std::unordered_map<std::string,int> msa_idx;
    msa_idx.reserve(msa_tip_names.size()*2);
    for (int i = 0; i < (int)msa_tip_names.size(); ++i) {
        msa_idx[msa_tip_names[i]] = i;
    }

    // 4) 建 node_id 對應：libpll rooted tree 通常 tips 在前面，inner 在後面，這裡保守點用遍歷拿順序
    //    我們做一次 "postorder" 遍歷，順便決定 child/parent 關係與 id。
    std::vector<pll_rnode_t*> postorder_nodes(num_nodes, nullptr);
    unsigned int count = 0;

    // 定義 callback：只接受一個 node
    auto cb = [](pll_rnode_t *node) -> int {
        return PLL_SUCCESS;
    };

    // 直接讓 libpll 把節點依序填進 outbuffer
    int rc = pll_rtree_traverse(rtree->root,
                                PLL_TREE_TRAVERSE_POSTORDER,
                                cb,                              // callback
                                postorder_nodes.data(),           // outbuffer
                                &count);
    throw_if(rc != PLL_SUCCESS, "pll_rtree_traverse (POSTORDER) failed.");
    throw_if(count != (int)num_nodes, "postorder node count mismatch.");

    // 建立一個 rnode* → id 的對應表
    std::unordered_map<pll_rnode_t*, int> id_of;
    id_of.reserve(num_nodes*2);

    // 先依 postorder 指定 id：0..N-1
    for (int i = 0; i < num_nodes; ++i) {
        id_of[ postorder_nodes[i] ] = i;
        out.nodes[i].id = i;
    }

    // 5) 填 left/right/parent/is_tip/name/branch_length_to_parent
    //    注意：root 的 node->length 是「它連到 parent 的長度」，root 沒 parent → 我們放 0
    for (int i = 0; i < num_nodes; ++i) {
        pll_rnode_t* nd = postorder_nodes[i];
        TreeNode &dst = out.nodes[i];

        // 判斷 tip (leaf)：libpll tip 沒有 child，內部有左右 child
        const bool is_tip = (nd->left == nullptr && nd->right == nullptr);
        dst.is_tip = is_tip;

        // name（tip 才有）
        if (is_tip && nd->label) dst.name = std::string(nd->label);
        else dst.name.clear();

        // parent 與 branch length
        if (nd->parent) {
            dst.parent = id_of[ nd->parent ];
            dst.branch_length_to_parent = nd->length; // 即 Newick 的冒號數字
        } else {
            dst.parent = -1; // root
            dst.branch_length_to_parent = 0.0;
            out.root_id = i;
        }

        // children（內部節點才有）
        if (!is_tip) {
            throw_if(!nd->left || !nd->right, "Internal node missing children (non-binary?).");
            dst.left  = id_of[ nd->left ];
            dst.right = id_of[ nd->right ];
        }
    }

    // 6) 對齊 MSA 名稱：建立 tip name -> node id；同時檢查每個 tip 在 MSA 中都找得到
    out.tip_node_by_name.reserve(num_tips*2);
    unsigned int found_tips = 0;
    for (const TreeNode& tn : out.nodes) {
        if (!tn.is_tip) continue;
        if (tn.name.empty()) {
            throw std::runtime_error("Encountered a tip with empty name.");
        }
        out.tip_node_by_name[tn.name] = tn.id;
        if (msa_idx.find(tn.name) != msa_idx.end()) ++found_tips;
    }
    throw_if(found_tips != num_tips, "Some tips in Newick not found in MSA names.");

    // 7) 產生「整棵樹的 postorder（以 id 表示）」：children 先於 parent
    out.postorder.resize(num_nodes);
    for (int i = 0; i < num_nodes; ++i) out.postorder[i] = out.nodes[i].id; // 因為我們 id 就是 postorder 順序

    // 8) 設定每個 node 的 CLV / Scaler 偏移（以 element 為單位，非 byte）
    const size_t clv_span  = (size_t) rate_cats * (size_t) states;
    const size_t clv_count_per_node = sites * clv_span;
    const size_t scaler_count_per_node = per_rate_scaling
        ? (sites * (size_t)rate_cats)
        : (sites);

    for (auto &tn : out.nodes) {
        tn.clv_offset    = (size_t)tn.id * clv_count_per_node;     // 你可以直接 d_clv_pool + tn.clv_offset 使用
        tn.scaler_offset = (size_t)tn.id * scaler_count_per_node;  // d_scaler_pool + tn.scaler_offset
    }

    // 9) 清理 libpll 結構
    pll_rtree_destroy(rtree, nullptr);

    return out;
}

// End-to-end pipeline: build tree, pack host arrays, compute matrices, upload to GPU.
BuildToGpuResult BuildAllToGPU(
    const std::vector<std::string>& msa_tip_names,
    const std::vector<std::string>& msa_rows,
    const std::string& newick_text,
    const std::vector<double>& Q_rowmajor,   // size = states*states
    const std::vector<double>& pi,           // size = states
    const std::vector<double>& rate_rates,   // size = rate_cats
    const std::vector<double>& rate_weights, // size = rate_cats
    size_t sites, int states, int rate_cats, bool per_rate_scaling)
{
    throw_if(states > 64, "states exceeds MAX_STATES (64).");
    throw_if(rate_cats > 8, "rate_cats exceeds MAX_RATECATS (8).");
    if (Q_rowmajor.size() != (size_t)states*(size_t)states)
        throw std::runtime_error("Q size mismatch.");
    if (pi.size() != (size_t)states)
        throw std::runtime_error("pi size mismatch.");

    // 1) Newick → TreeBuildResult（postorder、拓樸、offset）
    TreeBuildResult T = build_tree_from_newick_with_pll(
        msa_tip_names, newick_text, sites, states, rate_cats, per_rate_scaling);

    // 2) 對齊 MSA → HostPacking（topology arrays, tipchars, offsets, scaler）
    HostPacking H = pack_host_arrays_from_tree_and_msa(
        T, msa_tip_names, msa_rows, sites, states);

    // 3) Host 做 GTR 分解（V/Vinv/λ）
    EigResult Eigen = gtr_eigendecomp_cpu(Q_rowmajor.data(), pi.data(), states);

    fill_pmats_in_host_packing(T, H, Eigen, pi, rate_rates, states, rate_cats);

    // 4) 全部上傳到 GPU（含模型常數）
    DeviceTree D = upload_to_gpu(T, H, Eigen, rate_weights, sites, states, rate_cats, per_rate_scaling);

    return BuildToGpuResult{ std::move(D), std::move(T), std::move(H) };
}

// Evaluate tree log-likelihood using prepared DeviceTree/HostPacking and model vectors.
double EvaluateTreeLogLikelihood(
    const DeviceTree&      D,
    const TreeBuildResult& T,
    const HostPacking&     H,
    const std::vector<double>& pi,
    const std::vector<double>& rate_weights,
    cudaStream_t stream = 0)
{
    return EvaluateTreeLogLikelihood_device(D, T, H, pi, rate_weights, stream);
}
