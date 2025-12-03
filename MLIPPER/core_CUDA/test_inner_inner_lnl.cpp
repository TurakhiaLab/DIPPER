// test_inner_inner_lnl.cpp — Run PLL (float) and your GPU (double) with identical inputs
#include <cstdio>
#include <vector>
#include <random>
#include <chrono>
#include <cmath>
#include <cctype>
#include <cstring>
#include <string>
#include <cuda_runtime.h>
#include <libpll/pll.h>

#include "core_likelihood.cuh"  // Param / Likelihood_Inner_Inner

// ===== Usage =====
static void print_usage(const char* prog) {
    std::fprintf(stderr,
        "Usage: %s [--states 4|5|20] [--sites N] [--rate_cats R]\n"
        "Defaults: --states 5  --sites 8192  --rate_cats 1\n"
        "Examples:\n"
        "  %s --states 4 --sites 10000 --rate_cats 4\n"
        "  %s --states 5\n"
        "  %s --states 20 --sites 4096 --rate_cats 8\n",
        prog, prog, prog, prog);
}

int main(int argc, char** argv) {
    using namespace core_likelihood;

    // ===== flags =====
    unsigned states    = 5;      // 4:DNA, 5:DNA+gap, 20:protein
    unsigned rate_cats = 1;
    unsigned sites     = 8192;

    for (int i=1; i<argc; ++i) {
        std::string a = argv[i];
        auto need = [&](int &idx)->const char*{
            if (idx+1 >= argc) { print_usage(argv[0]); std::exit(1); }
            return argv[++idx];
        };
        if      (a == "--states")     states    = (unsigned)std::stoul(need(i));
        else if (a == "--sites")      sites     = (unsigned)std::stoul(need(i));
        else if (a == "--rate_cats")  rate_cats = (unsigned)std::stoul(need(i));
        else if (a == "-h" || a == "--help") { print_usage(argv[0]); return 0; }
        else {
            std::fprintf(stderr, "Unknown flag: %s\n", a.c_str());
            print_usage(argv[0]); return 1;
        }
    }
    if (!(states==4 || states==5 || states==20)) {
        std::fprintf(stderr, "[ERROR] --states must be 4, 5, or 20 (got %u)\n", states);
        return 1;
    }

    // ===== 1) 準備輸入 (float 版本供 PLL；double 版本供 GPU) =====
    // 1.1 隨機且每 (site, rate) 正規化的 CLV
    std::vector<float>  parent_clv_f(sites * rate_cats * states);
    std::vector<float>  child_clv_f (sites * rate_cats * states);
    {
        std::mt19937_64 rng(42);
        std::uniform_real_distribution<float> U(0.1f, 1.0f);
        auto fill_norm = [&](std::vector<float>& clv){
            for (unsigned i = 0; i < sites; ++i) {
                for (unsigned r = 0; r < rate_cats; ++r) {
                    float sum = 0.0f;
                    float* row = &clv[(i*rate_cats + r)*states];
                    for (unsigned s = 0; s < states; ++s) { row[s] = U(rng); sum += row[s]; }
                    const float inv = 1.0f / sum;
                    for (unsigned s = 0; s < states; ++s) row[s] *= inv;
                }
            }
        };
        fill_norm(parent_clv_f);
        fill_norm(child_clv_f);
    }
    std::vector<double> parent_clv_d(parent_clv_f.begin(), parent_clv_f.end());
    std::vector<double> child_clv_d (child_clv_f.begin(),  child_clv_f.end());

    // 1.2 P 矩陣（每 rate 單位矩陣，方便 smoke 對比）
    std::vector<float>  pmatrix_f(rate_cats * states * states, 0.0f);
    for (unsigned r = 0; r < rate_cats; ++r)
        for (unsigned j = 0; j < states; ++j)
            pmatrix_f[r*states*states + j*states + j] = 1.0f;
    std::vector<double> pmatrix_d(pmatrix_f.begin(), pmatrix_f.end());

    // 1.3 頻率（均勻）
    std::vector<float>  freqs_f(states, 1.0f / float(states));
    std::vector<double> freqs_d(freqs_f.begin(), freqs_f.end());

    // 1.4 rate 權重（全 1）
    std::vector<float>  rate_w_f(rate_cats, 1.0f);
    std::vector<double> rate_w_d(rate_w_f.begin(), rate_w_f.end());

    // 1.5 pattern weights
    std::vector<unsigned int> patt_w(sites, 1u);

    // 1.6 mixtures（只有一個 mixture → 全 0）
    std::vector<unsigned int> freqs_indices(rate_cats, 0u);

    // 1.7 invariants：先關閉
    float  inv_prop_f = 0.0f;
    double inv_prop_d = 0.0;
    const float*  invar_prop_ptr_f = &inv_prop_f;
    const double* invar_prop_ptr_d = &inv_prop_d;
    std::vector<int> invar_indices(sites, -1);

    // 1.8 scaler：先給 nullptr/零 buffer（與 tip-inner smoke 一致）
    const unsigned int* parent_scaler_pll = nullptr;
    const unsigned int* child_scaler_pll  = nullptr;
    std::vector<unsigned int> parent_scaler_gpu(sites * rate_cats, 0u);
    std::vector<unsigned int> child_scaler_gpu (sites * rate_cats, 0u);

    // ===== 2) libpll (float) =====
    std::vector<float*> freqs_ptrs_f(rate_cats, freqs_f.data());
    std::vector<float> per_site_ll(sites, 0.0f);
    const unsigned int attrib = 0u; // 不開其他旗標

    auto t0 = std::chrono::steady_clock::now();
    float pll_ll = pll_core_edge_loglikelihood_ii(
        states,
        sites,
        rate_cats,
        parent_clv_f.data(),
        parent_scaler_pll,
        child_clv_f.data(),
        child_scaler_pll,
        pmatrix_f.data(),
        freqs_ptrs_f.data(),
        rate_w_f.data(),
        patt_w.data(),
        invar_prop_ptr_f,
        invar_indices.data(),
        freqs_indices.data(),
        per_site_ll.data(),
        attrib);
    auto t1 = std::chrono::steady_clock::now();
    const double ms_pll = std::chrono::duration<double, std::milli>(t1 - t0).count();

    if (!std::isfinite(pll_ll)) {
        std::fprintf(stderr, "[ERROR] libpll total logL not finite: %.9g\n", (double)pll_ll);
        return 2;
    }

    // ===== 3) GPU (double) =====
    Param p{};
    p.states            = states;
    p.rate_cats         = rate_cats;
    p.sites             = sites;
    p.per_rate_scaling  = true;
    p.invar_proportion  = 0.0;

    Likelihood_Inner_Inner inner_inner;
    inner_inner.Initialize(p);

    inner_inner.ConstructionOnGpu(
        parent_clv_d.data(),
        parent_scaler_gpu.data(),
        child_clv_d.data(),
        child_scaler_gpu.data(),
        pmatrix_d.data(),
        freqs_d.data(),
        rate_w_d.data(),
        patt_w.data(),
        invar_indices.data(),
        p
    );

    cudaEvent_t ev_start, ev_stop;
    cudaEventCreate(&ev_start);
    cudaEventCreate(&ev_stop);
    cudaDeviceSynchronize();
    cudaEventRecord(ev_start, 0);

    inner_inner.ComputeLikelihood(p);

    cudaEventRecord(ev_stop, 0);
    cudaEventSynchronize(ev_stop);
    float ms_gpu = 0.0f;
    cudaEventElapsedTime(&ms_gpu, ev_start, ev_stop);
    cudaEventDestroy(ev_start);
    cudaEventDestroy(ev_stop);
    cudaDeviceSynchronize();

    const double gpu_ll = inner_inner.total_loglik;

    // ===== 4) Summary =====
    std::printf("[CFG ] states=%u, sites=%u, rate_cats=%u\n", states, sites, rate_cats);
    std::printf("[GPU(double)] logL = %.12f, time = %.3f ms\n", gpu_ll, ms_gpu);
    std::printf("[PLL(float)] logL = %.12f, time = %.3f ms\n", (double)pll_ll, ms_pll);
    std::printf("[DIFF] abs = %.6e\n", std::fabs(gpu_ll - (double)pll_ll));
    if (ms_pll > 0.0)
        std::printf("[Speed Up] GPU / PLL = %.3f\n", ms_pll / ms_gpu);
    else
        std::printf("[TIME] GPU / PLL = inf (PLL ms ~ 0)\n");

    inner_inner.CleanUp();
    return 0;
}
