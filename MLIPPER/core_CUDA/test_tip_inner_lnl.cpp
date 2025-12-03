// test_tip_inner_lnl.cpp — Run PLL (float) and your GPU (double) with identical inputs
#include <cstdio>
#include <vector>
#include <random>
#include <chrono>
#include <cassert>
#include <cmath>
#include <cctype>
#include <cstring>
#include <string>
#include <cuda_runtime.h>
#include <libpll/pll.h>

#include "core_likelihood.cuh"  // Param / Likelihood_Tip_Inner

// ========== Usage ==========
static void print_usage(const char* prog) {
    std::fprintf(stderr,
        "Usage: %s [--states 4|5|20] [--sites N] [--rate_cats R]\n"
        "Defaults: --states 5  --sites 8192  --rate_cats 1\n"
        "Examples:\n"
        "  %s --states 4 --sites 10000 --rate_cats 4   (DNA)\n"
        "  %s --states 5                                 (DNA + gap)\n"
        "  %s --states 20 --sites 4096 --rate_cats 8     (Protein)\n",
        prog, prog, prog, prog);
}

// ---------- tipmap helpers ----------
// DNA(4): A,C,G,T  (no dedicated GAP state)
static void make_tipmap_DNA4(std::vector<pll_state_t>& tipmap_pll,
                             std::vector<unsigned>&    tipmap_gpu)
{
    tipmap_pll.assign(256, (pll_state_t)0);
    tipmap_gpu.assign(256, 0u);

    auto set_mask = [&](unsigned char ch, unsigned m){
        tipmap_pll[ch] = (pll_state_t)m;
        tipmap_gpu[ch] = m;
    };

    // bit layout: A(0),C(1),G(2),T(3)
    const unsigned A = 1u << 0;
    const unsigned C = 1u << 1;
    const unsigned G = 1u << 2;
    const unsigned T = 1u << 3;
    const unsigned N = A | C | G | T; // ambiguous

    set_mask('A', A); set_mask('C', C); set_mask('G', G); set_mask('T', T);
    set_mask('a', A); set_mask('c', C); set_mask('g', G); set_mask('t', T);
    set_mask('N', N); set_mask('n', N);

    // treat '-' as unknown in 4-state (avoid -inf); maps to N
    set_mask('-', N);

    // IUPAC subset (optional)
    set_mask('R', A|G); set_mask('r', A|G);
    set_mask('Y', C|T); set_mask('y', C|T);
    set_mask('S', G|C); set_mask('s', G|C);
    set_mask('W', A|T); set_mask('w', A|T);
    set_mask('K', G|T); set_mask('k', G|T);
    set_mask('M', A|C); set_mask('m', A|C);
}

// DNA+GAP(5): A,C,G,T,GAP('-' as state 4)
static void make_tipmap_DNA5(std::vector<pll_state_t>& tipmap_pll,
                             std::vector<unsigned>&    tipmap_gpu)
{
    tipmap_pll.assign(256, (pll_state_t)0);
    tipmap_gpu.assign(256, 0u);
    auto set_mask = [&](unsigned char ch, unsigned m){
        tipmap_pll[ch] = (pll_state_t)m;
        tipmap_gpu[ch] = m;
    };

    const unsigned A = 1u << 0;
    const unsigned C = 1u << 1;
    const unsigned G = 1u << 2;
    const unsigned T = 1u << 3;
    const unsigned GAP = 1u << 4;
    const unsigned N = A | C | G | T; // N 不包含 GAP

    set_mask('A', A); set_mask('C', C); set_mask('G', G); set_mask('T', T);
    set_mask('a', A); set_mask('c', C); set_mask('g', G); set_mask('t', T);
    set_mask('N', N); set_mask('n', N);

    // 專用 gap 狀態
    set_mask('-', GAP);

    // IUPAC (對 4 個核酸位；不把 GAP 混進來)
    set_mask('R', A|G); set_mask('r', A|G);
    set_mask('Y', C|T); set_mask('y', C|T);
    set_mask('S', G|C); set_mask('s', G|C);
    set_mask('W', A|T); set_mask('w', A|T);
    set_mask('K', G|T); set_mask('k', G|T);
    set_mask('M', A|C); set_mask('m', A|C);
}

// Protein(20): 20 AA + ambiguous
static void make_tipmap_AA20(std::vector<pll_state_t>& tipmap_pll,
                             std::vector<unsigned>&    tipmap_gpu)
{
    tipmap_pll.assign(256, (pll_state_t)0);
    tipmap_gpu.assign(256, 0u);

    auto set_mask = [&](unsigned char ch, unsigned m){
        tipmap_pll[ch] = (pll_state_t)m;
        tipmap_gpu[ch] = m;
    };
    auto bit = [](int i){ return 1u << i; }; // 0..19

    // Order: A C D E F G H I K L M N P Q R S T V W Y
    set_mask('A', bit(0));  set_mask('a', bit(0));
    set_mask('C', bit(1));  set_mask('c', bit(1));
    set_mask('D', bit(2));  set_mask('d', bit(2));
    set_mask('E', bit(3));  set_mask('e', bit(3));
    set_mask('F', bit(4));  set_mask('f', bit(4));
    set_mask('G', bit(5));  set_mask('g', bit(5));
    set_mask('H', bit(6));  set_mask('h', bit(6));
    set_mask('I', bit(7));  set_mask('i', bit(7));
    set_mask('K', bit(8));  set_mask('k', bit(8));
    set_mask('L', bit(9));  set_mask('l', bit(9));
    set_mask('M', bit(10)); set_mask('m', bit(10));
    set_mask('N', bit(11)); set_mask('n', bit(11));
    set_mask('P', bit(12)); set_mask('p', bit(12));
    set_mask('Q', bit(13)); set_mask('q', bit(13));
    set_mask('R', bit(14)); set_mask('r', bit(14));
    set_mask('S', bit(15)); set_mask('s', bit(15));
    set_mask('T', bit(16)); set_mask('t', bit(16));
    set_mask('V', bit(17)); set_mask('v', bit(17));
    set_mask('W', bit(18)); set_mask('w', bit(18));
    set_mask('Y', bit(19)); set_mask('y', bit(19));

    const unsigned ALL20 = (1u << 20) - 1u;
    set_mask('X', ALL20); set_mask('x', ALL20);    // unknown → any AA
    set_mask('-', ALL20);                          // gap 當未知（避免 -inf）
    set_mask('B', bit(2) | bit(11)); set_mask('b', bit(2) | bit(11)); // D/N
    set_mask('Z', bit(3) | bit(13)); set_mask('z', bit(3) | bit(13)); // E/Q
    set_mask('J', bit(7) | bit(9));  set_mask('j', bit(7) | bit(9));  // I/L
}

// ---------- tipchars makers ----------
static void make_tipchars_DNA(std::vector<unsigned char>& tipchars, bool with_gap=false)
{
    if (!with_gap) {
        static const unsigned char bases[5] = {'A','C','G','T','N'};
        for (size_t i=0;i<tipchars.size();++i) tipchars[i] = bases[i % 5];
    } else {
        static const unsigned char bases[6] = {'A','C','G','T','N','-'};
        for (size_t i=0;i<tipchars.size();++i) tipchars[i] = bases[i % 6];
    }
}

static void make_tipchars_AA(std::vector<unsigned char>& tipchars)
{
    static const unsigned char aas[21] =
        {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','X'};
    for (size_t i=0;i<tipchars.size();++i) tipchars[i] = aas[i % 21];
}

int main(int argc, char** argv) {
    using namespace core_likelihood;

    // ===== flags =====
    unsigned states    = 5;      // default: DNA+gap
    unsigned rate_cats = 1;
    unsigned sites     = 8192;

    for (int i=1; i<argc; ++i) {
        std::string a = argv[i];
        auto need_val = [&](int &idx)->const char*{
            if (idx+1 >= argc) { print_usage(argv[0]); std::exit(1); }
            return argv[++idx];
        };
        if (a == "--states") {
            states = static_cast<unsigned>(std::stoul(need_val(i)));
        } else if (a == "--sites") {
            sites = static_cast<unsigned>(std::stoul(need_val(i)));
        } else if (a == "--rate_cats") {
            rate_cats = static_cast<unsigned>(std::stoul(need_val(i)));
        } else if (a == "-h" || a == "--help") {
            print_usage(argv[0]); return 0;
        } else {
            std::fprintf(stderr, "Unknown flag: %s\n", a.c_str());
            print_usage(argv[0]); return 1;
        }
    }
    if (!(states==4 || states==5 || states==20)) {
        std::fprintf(stderr, "[ERROR] --states must be 4, 5, or 20 (got %u)\n", states);
        return 1;
    }

    // ========== 1) Inputs ==========
    // 1.1 tipmap
    std::vector<pll_state_t> tipmap_pll;
    std::vector<unsigned>    tipmap_gpu;
    if (states == 4)      make_tipmap_DNA4(tipmap_pll, tipmap_gpu);
    else if (states == 5) make_tipmap_DNA5(tipmap_pll, tipmap_gpu);
    else                  make_tipmap_AA20(tipmap_pll, tipmap_gpu);
    const unsigned tipmap_size = 256;

    // 1.2 tipchars
    std::vector<unsigned char> tipchars(sites);
    if (states == 20) make_tipchars_AA(tipchars);
    else              make_tipchars_DNA(tipchars, /*with_gap=*/(states==5));

    // 1.3 parent_clv (random, normalized per site-rate)
    std::vector<double> parent_clv_d(sites * rate_cats * states);
    {
        std::mt19937_64 rng(42);
        std::uniform_real_distribution<double> U(0.1, 1.0);
        for (unsigned i = 0; i < sites; ++i) {
            for (unsigned r = 0; r < rate_cats; ++r) {
                double sum = 0.0;
                double* row = &parent_clv_d[(i*rate_cats + r)*states];
                for (unsigned s = 0; s < states; ++s) { row[s] = U(rng); sum += row[s]; }
                const double inv = 1.0 / sum;
                for (unsigned s = 0; s < states; ++s) row[s] *= inv;
            }
        }
    }

    // 1.4 P matrices = identity per rate
    std::vector<double> pmatrix_d(rate_cats * states * states, 0.0);
    for (unsigned r = 0; r < rate_cats; ++r)
        for (unsigned j = 0; j < states; ++j)
            pmatrix_d[r*states*states + j*states + j] = 1.0;

    // 1.5 frequencies = uniform
    std::vector<double> freqs_d(states, 1.0 / static_cast<double>(states));
    double* freqs_ptrs_d[1] = { freqs_d.data() }; // 如果 GPU API 需要 **

    // 1.6 rate_weights = all 1
    std::vector<double> rate_w_d(rate_cats, 1.0);

    // 1.7 pattern_weights = all 1
    std::vector<unsigned> patt_w(sites, 1u);

    // 1.8 invariants OFF
    std::vector<int> invar_idx(sites, -1);
    const double inv_prop = 0.0;

    // 1.9 mixtures index (single mixture -> 0)
    std::vector<unsigned> freqs_indices(rate_cats, 0u);

    // 1.A scaler
    const unsigned* parent_scaler_pll = nullptr;
    std::vector<unsigned> parent_scaler_gpu(sites * rate_cats, 0u);

    // ========== 2) libpll (float) ==========
    auto to_float = [](const std::vector<double>& v){
        std::vector<float> out(v.size());
        for (size_t i=0;i<v.size();++i) out[i] = static_cast<float>(v[i]);
        return out;
    };
    std::vector<float> parent_clv_f = to_float(parent_clv_d);
    std::vector<float> pmatrix_f    = to_float(pmatrix_d);
    std::vector<float> freqs_f      = to_float(freqs_d);
    std::vector<float> rate_w_f     = to_float(rate_w_d);
    float* freqs_ptrs_f[1] = { freqs_f.data() };

    float inv_prop_arr[1]          = { 0.0f };
    const float* invar_proportion_f = inv_prop_arr;
    const int*   invar_indices_f    = nullptr;
    const unsigned attrib_pll        = 0u;

    std::vector<float> per_site_ll(sites, 0.0f);

    auto t0 = std::chrono::steady_clock::now();
    float pll_ll = pll_core_edge_loglikelihood_ti(
        states,
        sites,
        rate_cats,
        parent_clv_f.data(),
        parent_scaler_pll,
        tipchars.data(),
        tipmap_pll.data(),
        tipmap_size,
        pmatrix_f.data(),
        freqs_ptrs_f,
        rate_w_f.data(),
        patt_w.data(),
        invar_proportion_f,
        invar_indices_f,
        freqs_indices.data(),
        per_site_ll.data(),
        attrib_pll);
    auto t1 = std::chrono::steady_clock::now();
    const double ms_pll = std::chrono::duration<double, std::milli>(t1 - t0).count();

    if (!std::isfinite(pll_ll)) {
        std::fprintf(stderr,
                     "[ERROR] libpll total logL not finite: %.9g\n"
                     "Check: tipmap/tipchars and states/layout.\n",
                     pll_ll);
        return 2;
    }

    // ========== 3) GPU (double) ==========
    Param p{};
    p.states            = states;
    p.rate_cats         = rate_cats;
    p.sites             = sites;
    p.per_rate_scaling  = false;
    p.invar_proportion  = 0.0;

    Likelihood_Tip_Inner tip_inner;
    tip_inner.Initialize(p);

    tip_inner.ConstructionOnGpu(
        parent_clv_d.data(),
        parent_scaler_gpu.data(),            // zero buffer
        tipchars.data(),                     // ASCII
        tipmap_gpu.data(),                   // 256 entries, unsigned*
        pmatrix_d.data(),                    // identity
        freqs_d.data(),                      // uniform
        rate_w_d.data(),                     // all 1
        patt_w.data(),                       // all 1
        invar_idx.data(),                    // all -1
        p);

    cudaEvent_t ev_start, ev_stop;
    cudaEventCreate(&ev_start);
    cudaEventCreate(&ev_stop);
    cudaDeviceSynchronize();
    cudaEventRecord(ev_start, 0);

    tip_inner.ComputeLikelihood(p);

    cudaEventRecord(ev_stop, 0);
    cudaEventSynchronize(ev_stop);
    float ms_gpu = 0.0f;
    cudaEventElapsedTime(&ms_gpu, ev_start, ev_stop);
    cudaEventDestroy(ev_start);
    cudaEventDestroy(ev_stop);
    cudaDeviceSynchronize();

    const double gpu_ll = tip_inner.total_loglik;

    // ========== 4) Summary ==========
    std::printf("[CFG ] states=%u, sites=%u, rate_cats=%u\n", states, sites, rate_cats);
    std::printf("[GPU ] logL = %.12f, time = %.3f ms\n", gpu_ll, ms_gpu);
    std::printf("[PLLf] logL = %.12f, time = %.3f ms\n", static_cast<double>(pll_ll), ms_pll);
    std::printf("[DIFF] abs = %.6e\n", std::fabs(gpu_ll - static_cast<double>(pll_ll)));
    if (ms_pll > 0.0)
        std::printf("[Speed Up] GPU / PLLf = %.3f\n",  ms_pll / ms_gpu);
    else
        std::printf("[TIME] GPU / PLLf = inf (PLLf ms ~ 0)\n");

    tip_inner.CleanUp();
    return 0;
}
