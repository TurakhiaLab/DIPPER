#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cctype>
#include <filesystem>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <libpll/pll.h>
#if defined(__has_include)
#  if __has_include(<libpll/pllmod_tree.h>)
#    include <libpll/pllmod_tree.h>
#  endif
#endif

#include "tree.hpp"
#include "parse_file.hpp"


static const TreeBuildResult* g_gpu_tree_for_verify = nullptr;
static HostPacking* g_gpu_host_for_verify = nullptr;
static const std::vector<double>* g_gpu_pi_for_verify = nullptr;
static const std::vector<double>* g_gpu_rate_w_for_verify = nullptr;

struct SimpleNode {
    std::string label;
    std::string length_str;
    bool has_length = false;
    std::vector<std::unique_ptr<SimpleNode>> children;
};

static void skip_ws(const std::string& text, size_t& pos) {
    while (pos < text.size() && std::isspace(static_cast<unsigned char>(text[pos]))) ++pos;
}

static std::unique_ptr<SimpleNode> parse_newick_node(const std::string& text, size_t& pos);

static std::unique_ptr<SimpleNode> parse_newick_node(const std::string& text, size_t& pos) {
    skip_ws(text, pos);
    auto node = std::make_unique<SimpleNode>();
    if (pos < text.size() && text[pos] == '(') {
        ++pos;
        while (true) {
            node->children.push_back(parse_newick_node(text, pos));
            skip_ws(text, pos);
            if (pos >= text.size())
                throw std::runtime_error("Unexpected end of Newick while parsing children.");
            if (text[pos] == ',') {
                ++pos;
                continue;
            }
            if (text[pos] == ')') {
                ++pos;
                break;
            }
            throw std::runtime_error("Unexpected character in Newick children list.");
        }
    }

    skip_ws(text, pos);
    size_t label_start = pos;
    while (pos < text.size()) {
        char c = text[pos];
        if (c == ':' || c == ',' || c == ')' || c == ';' || std::isspace(static_cast<unsigned char>(c)))
            break;
        ++pos;
    }
    if (label_start < pos)
        node->label = text.substr(label_start, pos - label_start);

    skip_ws(text, pos);
    if (pos < text.size() && text[pos] == ':') {
        ++pos;
        skip_ws(text, pos);
        size_t len_start = pos;
        while (pos < text.size()) {
            char c = text[pos];
            if (c == ',' || c == ')' || c == ';' || std::isspace(static_cast<unsigned char>(c)))
                break;
            ++pos;
        }
        if (len_start < pos) {
            node->length_str = text.substr(len_start, pos - len_start);
            node->has_length = true;
        }
    }
    return node;
}

static void resolve_polytomies(SimpleNode* node) {
    for (auto& child : node->children) resolve_polytomies(child.get());
    if (node->children.size() <= 2) return;
    std::vector<std::unique_ptr<SimpleNode>> kids = std::move(node->children);
    std::unique_ptr<SimpleNode> combined = std::move(kids[0]);
    for (size_t i = 1; i + 1 < kids.size(); ++i) {
        auto merged = std::make_unique<SimpleNode>();
        merged->children.reserve(2);
        merged->children.push_back(std::move(combined));
        merged->children.push_back(std::move(kids[i]));
        merged->has_length = true;
        merged->length_str = "0";
        combined = std::move(merged);
    }
    node->children.clear();
    node->children.push_back(std::move(combined));
    node->children.push_back(std::move(kids.back()));
}

static void append_newick(const SimpleNode& node, std::string& out) {
    if (!node.children.empty()) {
        out.push_back('(');
        for (size_t i = 0; i < node.children.size(); ++i) {
            append_newick(*node.children[i], out);
            if (i + 1 < node.children.size()) out.push_back(',');
        }
        out.push_back(')');
    }
    if (!node.label.empty()) out += node.label;
    if (node.has_length) {
        out.push_back(':');
        out += node.length_str;
    }
}

static std::string normalize_newick(const std::string& raw) {
    try {
        size_t pos = 0;
        auto root = parse_newick_node(raw, pos);
        skip_ws(raw, pos);
        if (pos < raw.size() && raw[pos] == ';') ++pos;
        skip_ws(raw, pos);
        if (pos != raw.size())
            throw std::runtime_error("Extra characters after Newick tree.");
        resolve_polytomies(root.get());
        std::string normalized;
        append_newick(*root, normalized);
        normalized.push_back(';');
        return normalized;
    } catch (...) {
        return raw;
    }
}

static int cb_full_traversal(pll_unode_t* node)
{
    (void)node;
    return 1;
}

static void set_missing_branch_length(pll_utree_t* tree, float length)
{
    const unsigned int total = tree->tip_count + tree->inner_count;
    for (unsigned int i = 0; i < total; ++i) {
        if (!tree->nodes[i]) continue;
        if (tree->nodes[i]->length == 0.0f) tree->nodes[i]->length = length;
    }
}

static bool has_flag(int argc, char** argv, const char* flag)
{
    if (!flag) return false;
    for (int i = 1; i < argc; ++i) {
        if (!argv[i]) continue;
        if (std::string(argv[i]) == flag) return true;
    }
    return false;
}

static std::string read_config_path_from_args(int argc, char** argv)
{
    std::string path = "data/config.yaml";
    for (int i = 1; i < argc; ++i) {
        if (!argv[i]) continue;
        std::string arg = argv[i];
        if (arg == "--config" && i + 1 < argc && argv[i + 1]) {
            path = argv[++i];
        }
    }
    namespace fs = std::filesystem;
    if (!fs::exists(path)) {
        fs::path alt = "tree_generation/data/config.yaml";
        if (fs::exists(alt)) path = alt.string();
    }
    return path;
}

static pll_partition_t* build_pll_partition(
    pll_utree_t* tree,
    const std::vector<std::string>& tip_names,
    const std::vector<std::string>& tip_seqs,
    size_t sites,
    int states,
    int rate_cats,
    const std::vector<double>& pi,
    const std::vector<double>& rate_weights,
    const std::vector<double>& rate_rates,
    const std::vector<double>& Q,
    unsigned int attrib)
{
    if (tip_names.size() != tip_seqs.size())
        throw std::runtime_error("MSA names/rows size mismatch.");
    if ((unsigned int)tip_names.size() != tree->tip_count) {
        throw std::runtime_error("Tree tip count does not match MSA names.");
    }
    if (pi.size() != (size_t)states)
        throw std::runtime_error("Frequency vector size mismatch.");
    if (rate_weights.size() != (size_t)rate_cats)
        throw std::runtime_error("Rate weights size mismatch.");
    const size_t states_sz = static_cast<size_t>(states);
    if (Q.size() != states_sz * states_sz)
        throw std::runtime_error("Q matrix size mismatch.");
    if (states != 4)
        throw std::runtime_error("PLL verification currently supports only 4-state models.");

    unsigned int tip_count = tree->tip_count;
    unsigned int inner_count = tree->inner_count;
    unsigned int branch_count = tip_count + inner_count - 1;

    pll_partition_t* partition = pll_partition_create(
        tip_count,
        inner_count,
        (unsigned int)states,
        (unsigned int)sites,
        1,
        branch_count,
        (unsigned int)rate_cats,
        inner_count,
        attrib);
    if (!partition)
        throw std::runtime_error("pll_partition_create failed.");

    std::vector<float> freqs(states);
    for (int i = 0; i < states; ++i) freqs[i] = static_cast<float>(pi[i]);
    pll_set_frequencies(partition, 0, freqs.data());

    std::array<unsigned int, 6> idx_a = {0, 0, 0, 1, 1, 2};
    std::array<unsigned int, 6> idx_b = {1, 2, 3, 2, 3, 3};
    float subst_params[6];
    for (size_t i = 0; i < 6; ++i) {
        unsigned int a = idx_a[i];
        unsigned int b = idx_b[i];
        double denom = pi[b];
        if (denom <= 1e-12) denom = 1e-12;
        subst_params[i] = static_cast<float>(Q[a * states_sz + b] / denom);
    }
    pll_set_subst_params(partition, 0, subst_params);

    std::vector<float> cat_rates(rate_cats);
    std::vector<float> cat_weights(rate_cats);
    for (int rc = 0; rc < rate_cats; ++rc) {
        cat_rates[rc] = static_cast<float>(rate_rates[rc]);
        cat_weights[rc] = static_cast<float>(rate_weights[rc]);
    }
    pll_set_category_rates(partition, cat_rates.data());
    pll_set_category_weights(partition, cat_weights.data());

    std::unordered_map<std::string, unsigned int> label_to_index;
    label_to_index.reserve(tip_count * 2);
    for (unsigned int i = 0; i < tip_count; ++i) {
        pll_unode_t* node = tree->nodes[i];
        if (!node || !node->label) {
            throw std::runtime_error("Tip node lacks label or pointer.");
        }
        label_to_index[node->label] = node->clv_index;
    }

    for (size_t i = 0; i < tip_names.size(); ++i) {
        if (tip_seqs[i].size() != sites)
            throw std::runtime_error("MSA row length mismatch.");
        auto it = label_to_index.find(tip_names[i]);
        if (it == label_to_index.end())
            throw std::runtime_error("MSA tip not found in tree: " + tip_names[i]);
        pll_set_tip_states(partition, it->second, pll_map_nt, tip_seqs[i].c_str());
    }

    return partition;
}

static void verify_tree_loglik_pll_vs_gpu(
    pll_utree_t* rtree,
    pll_partition_t* part,
    unsigned int attrib,
    const DeviceTree& dev_tree)
{
    (void)attrib;
    if (!rtree || !part)
        throw std::runtime_error("PLL structures not provided.");
    if (!g_gpu_tree_for_verify || !g_gpu_host_for_verify ||
        !g_gpu_pi_for_verify || !g_gpu_rate_w_for_verify) {
        throw std::runtime_error("GPU structures not available for verification.");
    }

    const unsigned int tip_count = rtree->tip_count;
    const unsigned int inner_count = rtree->inner_count;
    const unsigned int nodes_count = tip_count + inner_count;
    const unsigned int branch_count = nodes_count - 1;
    const unsigned int rate_cats = dev_tree.rate_cats;
    if (rate_cats == 0)
        throw std::runtime_error("Device tree has zero rate categories.");

    std::vector<pll_unode_t*> travbuffer(nodes_count);
    unsigned int traversal_size = 0;
    pll_unode_t* root = rtree->vroot;
    if (!root)
        throw std::runtime_error("PLL tree has no root.");

    if (!pll_utree_traverse(root,
                            PLL_TREE_TRAVERSE_POSTORDER,
                            cb_full_traversal,
                            travbuffer.data(),
                            &traversal_size)) {
        throw std::runtime_error("PLL traversal failed.");
    }

    std::vector<float> branch_lengths(branch_count);
    std::vector<unsigned int> matrix_indices(branch_count);
    std::vector<pll_operation_t> operations(inner_count);
    unsigned int matrix_count = 0;
    unsigned int ops_count = 0;

    pll_utree_create_operations(
        travbuffer.data(),
        traversal_size,
        branch_lengths.data(),
        matrix_indices.data(),
        operations.data(),
        &matrix_count,
        &ops_count);

    if (ops_count == 0)
        throw std::runtime_error("No operations generated for PLL partial update.");

    std::vector<unsigned int> params_indices(rate_cats, 0u);
    pll_update_prob_matrices(part,
                             params_indices.data(),
                             matrix_indices.data(),
                             branch_lengths.data(),
                             matrix_count);
                             
    const auto start_pll = std::chrono::steady_clock::now();
    pll_update_partials(part, operations.data(), ops_count);
    float pll_logl = pll_compute_root_loglikelihood(
        part,
        root->clv_index,
        root->scaler_index,
        params_indices.data(),
        nullptr);
    const auto end_pll = std::chrono::steady_clock::now();
    const double pll_ms = std::chrono::duration<double, std::milli>(end_pll - start_pll).count();
    cudaEvent_t start, stop;
    float gpu_ms_kernel = 0.0f;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    const auto start_gpu = std::chrono::steady_clock::now();

    cudaEventRecord(start);
    double gpu_logl = EvaluateTreeLogLikelihood(
        dev_tree,
        *g_gpu_tree_for_verify,
        *g_gpu_host_for_verify,
        *g_gpu_pi_for_verify,
        *g_gpu_rate_w_for_verify,
        0);
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);

    // Get elapsed time (milliseconds)
    cudaEventElapsedTime(&gpu_ms_kernel, start, stop);
    const auto end_gpu = std::chrono::steady_clock::now();
    const double gpu_ms = std::chrono::duration<double, std::milli>(end_gpu - start_gpu).count();
    printf("GPU kernel time = %.3f ms\n", gpu_ms);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    

    
    const double pll_logl_d = static_cast<double>(pll_logl);
    const double abs_diff = std::fabs(pll_logl_d - gpu_logl);
    const double denom = std::max(std::fabs(pll_logl_d), 1e-12);
    const double rel_diff = abs_diff / denom;

    std::cout << "[PLL ] logL = " << pll_logl_d << ", time = " << pll_ms << " ms\n";
    std::cout << "[SELF] logL = " << gpu_logl << ", time = " << gpu_ms << " ms\n";
    std::cout << "[Speed Up GPU / CPU] = " << std::fixed << std::setprecision(3) << (pll_ms / gpu_ms) << "\n";
    std::cout << "[DIFF] abs = " << abs_diff << ", rel = " << std::setprecision(5) <<  rel_diff << "\n";
}


void print_tree_rec(const TreeBuildResult& T, int node_id, int depth)
{
    if (node_id < 0) return;
    const TreeNode &nd = T.nodes[node_id];

    // 縮排
    for (int i = 0; i < depth; ++i) std::cout << "  ";

    // 節點基本資訊
    std::cout << "[" << nd.id << "]";
    if (nd.is_tip) {
        std::cout << " (tip: " << nd.name << ")";
    } else {
        std::cout << " (inner)";
    }

    if (nd.parent >= 0) {
        std::cout << "  len=" << nd.branch_length_to_parent
                  << "  parent=" << nd.parent;
    } else {
        std::cout << "  <ROOT>";
    }
    std::cout << "\n";

    // 如果是 internal node，就印左右小孩
    if (!nd.is_tip) {
        if (nd.left  >= 0) print_tree_rec(T, nd.left,  depth + 1);
        if (nd.right >= 0) print_tree_rec(T, nd.right, depth + 1);
    }
}

void print_tree_structure(const TreeBuildResult& T)
{
    std::cout << "==== Tree structure (indented) ====\n";
    if (T.root_id < 0) {
        std::cout << "No root_id set!\n";
        return;
    }
    print_tree_rec(T, T.root_id, 0);
    std::cout << "===================================\n";
}

static void normalize_vector(std::vector<double>& vec) {
    double sum = 0.0;
    for (double v : vec) sum += v;
    if (sum <= 0.0) return;
    for (double& v : vec) v /= sum;
}

static std::vector<double> build_mixture_weights(const parse::ModelConfig& model, int rate_cats) {
    std::vector<double> weights;
    if (rate_cats <= 0) return weights;
    if ((int)model.rate_weights.size() == rate_cats) {
        weights = model.rate_weights;
    } else {
        weights.assign(rate_cats, 1.0 / rate_cats);
    }

    double sum = 0.0;
    for (double v : weights) sum += v;
    if (sum <= 0.0) {
        weights.assign(rate_cats, 1.0 / rate_cats);
        return weights;
    }
    for (double& v : weights) v /= sum;
    return weights;
}

static std::vector<double> build_gamma_rate_categories(double alpha, int rate_cats) {
    std::vector<double> rates(rate_cats, 1.0);
    if (rate_cats <= 1 || alpha <= 0.0) return rates;
    std::vector<float> gamma_tmp(rate_cats);
    int status = pll_compute_gamma_cats(static_cast<float>(alpha), rate_cats, gamma_tmp.data(), PLL_GAMMA_RATES_MEAN);
    if (status != PLL_SUCCESS) {
        throw std::runtime_error("pll_compute_gamma_cats failed.");
    }
    for (int i = 0; i < rate_cats; ++i) {
        rates[i] = static_cast<double>(gamma_tmp[i]);
    }
    return rates;
}

static std::vector<double> build_gtr_q_matrix(
    int states,
    const parse::ModelConfig& model,
    const std::vector<double>& pi)
{
    std::vector<double> Q(states * states, 0.0);
    if (states != 4 || model.rates.size() < 6 || pi.size() != 4)
        return Q;

    auto set_pair = [&](int i, int j, double rate) {
        Q[i * states + j] = rate * pi[j];
        Q[j * states + i] = rate * pi[i];
    };

    set_pair(0, 1, model.rates[0]); // A-C
    set_pair(0, 2, model.rates[1]); // A-G
    set_pair(0, 3, model.rates[2]); // A-T
    set_pair(1, 2, model.rates[3]); // C-G
    set_pair(1, 3, model.rates[4]); // C-T
    set_pair(2, 3, model.rates[5]); // G-T

    for (int i = 0; i < states; ++i) {
        double row_sum = 0.0;
        for (int j = 0; j < states; ++j) {
            if (i != j) row_sum += Q[i * states + j];
        }
        Q[i * states + i] = -row_sum;
    }

    // 1) 算 μ = - Σ_i π_i * Q_ii
    double mu = 0.0;
    for (int i = 0; i < states; ++i) {
        mu -= pi[i] * Q[i * states + i];  // 注意：Q_ii 是負的
    }

    // 2) 整個 Q 都除以同一個 μ
    for (int idx = 0; idx < states * states; ++idx) {
        Q[idx] /= mu;
    }

    return Q;
}

int main(int argc, char** argv) {
    std::string config_path = read_config_path_from_args(argc, argv);
    bool verify_mode = has_flag(argc, argv, "--verify-pll");

    try {
        auto inputs = parse::load_inputs_from_config(config_path);
        const auto& alignment = inputs.alignment;
        if (alignment.names.empty())
            throw std::runtime_error("Alignment from config contains no sequences.");
        if (alignment.sites == 0)
            throw std::runtime_error("Alignment from config contains zero sites.");

        const auto& msa_names = alignment.names;
        const auto& rows = alignment.sequences;
        size_t sites = alignment.sites;
        std::string newick = normalize_newick(inputs.tree);
        const auto& model = inputs.config.model;
        int states = model.states;
        if (states <= 0)
            throw std::runtime_error("Model states must be positive.");
        int rate_cats = model.ncat;
        if (rate_cats <= 0)
            throw std::runtime_error("Model ncat must be positive.");
        bool per_rate_scaling = model.per_rate_scaling;

        std::vector<double> pi = model.freqs;
        if ((int)pi.size() != states) pi.assign(states, 1.0 / states);
        normalize_vector(pi);
        std::vector<double> rate_weights = build_mixture_weights(model, rate_cats);
        std::vector<double> rate_rates = build_gamma_rate_categories(model.alpha, rate_cats);
        std::vector<double> Q = build_gtr_q_matrix(states, model, pi);
        size_t free_before, total_before;
        cudaMemGetInfo(&free_before, &total_before);

        const auto start_gpu = std::chrono::steady_clock::now();
        auto res = BuildAllToGPU(msa_names, rows, newick, Q, pi, rate_rates, rate_weights, sites, states, rate_cats, per_rate_scaling);
        
        g_gpu_tree_for_verify = &res.tree;
        g_gpu_host_for_verify = &res.hostPack;
        g_gpu_pi_for_verify = &pi;
        g_gpu_rate_w_for_verify = &rate_weights;
        
        std::cout << "Uploaded. N=" << res.dev.N << ", tips=" << res.dev.tips
                  << ", per_node_elems=" << res.dev.per_node_elems() << "\n";
        print_tree_structure(res.tree);
        if (verify_mode) {
            pll_utree_t* rtree = pll_utree_parse_newick_string_rooted(newick.c_str());
            if (!rtree) throw std::runtime_error("PLL tree parsing failed.");
            set_missing_branch_length(rtree, 0.000001f);
            unsigned int pll_attrib = PLL_ATTRIB_ARCH_CPU;
            pll_partition_t* part = build_pll_partition(
                rtree,
                msa_names,
                rows,
                sites,
                states,
                rate_cats,
                pi,
                rate_weights,
                rate_rates,
                Q,
                pll_attrib);
            verify_tree_loglik_pll_vs_gpu(rtree, part, pll_attrib, res.dev);
            size_t free_after, total_after;
            cudaMemGetInfo(&free_after, &total_after);
            size_t used = free_before - free_after;
            printf("Peak memory delta = %zu GB\n", used / (1024 * 1024 * 1024));
            pll_partition_destroy(part);
            pll_utree_destroy(rtree, nullptr);
            free_device_tree(res.dev);
            return 0;
        }

        

        double logL = EvaluateTreeLogLikelihood(
            res.dev,
            res.tree,
            res.hostPack,
            pi,
            rate_weights,
            0
        );
        const auto end_gpu = std::chrono::steady_clock::now();
        const double gpu_ms = std::chrono::duration<double, std::milli>(end_gpu - start_gpu).count();
        std::cout << "time = " << gpu_ms << " ms\n";
        std::cout << "Tree log-likelihood = " << logL << "\n";
        free_device_tree(res.dev);

    } catch(const std::exception& e){
        std::cerr << "ERR: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
