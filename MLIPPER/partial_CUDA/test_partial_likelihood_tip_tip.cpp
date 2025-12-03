#include <cstdio>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cuda_runtime.h>
#include "partial_likelihood.cuh"
#include <libpll/pll.h>

#ifndef CUDA_CHECK
#define CUDA_CHECK(expr) do { \
  cudaError_t _err = (expr);  \
  if (_err != cudaSuccess) {  \
    fprintf(stderr, "[CUDA] %s failed: %s (%s:%d)\n", #expr, cudaGetErrorString(_err), __FILE__, __LINE__); \
    std::exit(1); \
  } \
} while(0)
#endif

static double max_abs_diff_mixed(const std::vector<float>& a, const std::vector<double>& b) { // ★修正處：混合型別比較
    double m = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        m = std::max(m, std::fabs((double)a[i] - b[i]));
    }
    return m;
}

static double cpu_loglik_from_parent_clv_double(const double* parent,
                                                size_t sites, int states, int rate_cats,
                                                const std::vector<double>& freqs) {
    const size_t span = (size_t)states * (size_t)rate_cats;
    const double inv_rc = 1.0 / (double)rate_cats;
    double total = 0.0;
    for (size_t n = 0; n < sites; ++n) {
        const double* site_vec = parent + n * span;
        double site_like = 0.0;
        for (int r = 0; r < rate_cats; ++r) {
            const double* v = site_vec + r * states;
            double dot = 0.0;
            for (int x = 0; x < states; ++x) dot += freqs[x] * v[x];
            site_like += inv_rc * dot;
        }
        if (site_like <= 0) site_like = 1e-300;
        total += std::log(site_like);
    }
    return total;
}

// ★修正處：給 float* 的 logL（轉成 double 再算）
static double cpu_loglik_from_parent_clv_float(const float* parent,
                                               size_t sites, int states, int rate_cats,
                                               const std::vector<double>& freqs) {
    const size_t span = (size_t)states * (size_t)rate_cats;
    const double inv_rc = 1.0 / (double)rate_cats;
    double total = 0.0;
    for (size_t n = 0; n < sites; ++n) {
        const float* site_vec = parent + n * span;
        double site_like = 0.0;
        for (int r = 0; r < rate_cats; ++r) {
            const float* v = site_vec + r * states;
            double dot = 0.0;
            for (int x = 0; x < states; ++x) dot += freqs[x] * (double)v[x];
            site_like += inv_rc * dot;
        }
        if (site_like <= 0) site_like = 1e-300;
        total += std::log(site_like);
    }
    return total;
}

int main() {
    partial_likelihood::Param P;
    P.sites = 4;
    P.states = 4;             // DNA
    P.rate_cats = 4;
    P.per_rate_scaling = false;

    const unsigned int tipmap_size = 4;   // DNA
    const unsigned int log2_stride = 2;   // stride = 4
    const size_t span = (size_t)P.states * (size_t)P.rate_cats;
    const size_t parent_count = P.sites * span;

    std::vector<double> freqs(P.states, 1.0 / P.states);

    std::vector<unsigned char> h_left(P.sites), h_right(P.sites);
    for (size_t i = 0; i < P.sites; ++i) {
        h_left[i]  = (unsigned char)(i % tipmap_size);
        h_right[i] = (unsigned char)((i * 3) % tipmap_size);
    }

    // GPU/測試用 double 版 lookup（保留）
    const size_t nblocks = ((size_t)1u << log2_stride) * tipmap_size;
    const size_t lookup_elems = nblocks * span;
    std::vector<double> h_lookup_d(lookup_elems, 0.0);
    for (unsigned int j = 0; j < tipmap_size; ++j) {
        for (unsigned int k = 0; k < tipmap_size; ++k) {
            const size_t base = (((size_t)j << log2_stride) + k) * span;
            for (size_t t = 0; t < span; ++t) {
                h_lookup_d[base + t] = 1000.0 * j + 10.0 * k + (double)t;
            }
        }
    }

    // ★修正處：libpll 是 float 版 → 建 float 版 lookup 給 libpll
    std::vector<float> h_lookup_f(lookup_elems, 0.0f);
    for (size_t i = 0; i < lookup_elems; ++i) h_lookup_f[i] = (float)h_lookup_d[i];

    // ==== GPU ====
    cudaStream_t stream = 0;
    partial_likelihood::Partial_Likelihood_Tip_Tip tt;
    tt.ConstructionOnGpu(
        P,
        h_left.data(),
        h_right.data(),
        h_lookup_d.data(),   // GPU 仍用 double 版 lookup
        lookup_elems,
        stream);
    tt.UpdatePartialLikelihood(P, 0);

    std::vector<double> gpu_parent(parent_count);
    CUDA_CHECK(cudaMemcpy(gpu_parent.data(), tt.d_parent_clv,
                          parent_count * sizeof(double), cudaMemcpyDeviceToHost));

    // ==== CPU（libpll, float 版）====
    std::vector<float> cpu_parent_f(parent_count, 0.0f);              // ★修正處：float
    std::vector<unsigned int> cpu_scaler(P.sites, 0);
    const pll_state_t tipmap[4] = {1u, 2u, 4u, 8u};

    pll_core_update_partial_tt(
        (unsigned int)P.states,
        (unsigned int)P.sites,
        (unsigned int)P.rate_cats,
        cpu_parent_f.data(),                 // ★ float*
        cpu_scaler.data(),
        h_left.data(),
        h_right.data(),
        tipmap,                              // ★ const unsigned int*
        tipmap_size,
        h_lookup_f.data(),                   // ★ float*
        false                     // ★ 關閉縮放
    );

    // ==== 比對 CLV（float vs double）====
    const double diff = max_abs_diff_mixed(cpu_parent_f, gpu_parent); // ★修正處
    printf("[TT] max_abs_diff(parent CLV, GPU vs libpll CPU) = %.3e\n", diff);

    // ==== logL（各自算）====
    const double logL_gpu = cpu_loglik_from_parent_clv_double(gpu_parent.data(),
                            P.sites, P.states, P.rate_cats, freqs);
    const double logL_cpu = cpu_loglik_from_parent_clv_float(cpu_parent_f.data(),
                            P.sites, P.states, P.rate_cats, freqs);    // ★修正處
    printf("[TT] logL(GPU) = %.12f\n", logL_gpu);
    printf("[TT] logL(CPU) = %.12f\n", logL_cpu);
    printf("[TT] |Δ logL|   = %.3e\n", std::fabs(logL_gpu - logL_cpu));

    tt.CleanUp();
    return 0;
}
