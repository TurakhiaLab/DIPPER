#include <cstdio>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cuda_runtime.h>
#include "partial_likelihood.cuh"
#include <libpll/pll.h>

// ----------------- CUDA error helper -----------------
#ifndef CUDA_CHECK
#define CUDA_CHECK(expr) do { \
  cudaError_t _err = (expr);  \
  if (_err != cudaSuccess) {  \
    fprintf(stderr, "[CUDA] %s failed: %s (%s:%d)\n", #expr, cudaGetErrorString(_err), __FILE__, __LINE__); \
    std::exit(1); \
  } \
} while(0)
#endif

// ----------------- diff & logL helpers -----------------
static double max_abs_diff_mixed(const std::vector<float>& a, const std::vector<double>& b) {
    double m = 0.0;
    const size_t n = std::min(a.size(), b.size());
    for (size_t i = 0; i < n; ++i) m = std::max(m, std::fabs((double)a[i] - b[i]));
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

// ----------------- main test -----------------
int main() {
    // 基本參數
    partial_likelihood::Param P;
    P.sites = 8196;
    P.states = 4;            // DNA
    P.rate_cats = 4;
    P.per_rate_scaling = false; // 本測試先不做縮放

    const unsigned int tipmap_size = 4;   // DNA: A,C,G,T -> {1,2,4,8}
    const size_t span = (size_t)P.states * (size_t)P.rate_cats;
    const size_t parent_count = (size_t)P.sites * span;

    // root/頻率：均勻
    std::vector<double> freqs(P.states, 1.0 / P.states);

    // 左 tip chars
    std::vector<unsigned char> h_left(P.sites);
    for (size_t i = 0; i < P.sites; ++i) h_left[i] = (unsigned char)((i * 2) % tipmap_size);
    std::vector<unsigned char> h_left_bits(P.sites);
    for (size_t i = 0; i < P.sites; ++i) {
        unsigned idx = (unsigned)h_left[i];   // 0..3
        h_left_bits[i] = (unsigned char)(1u << idx);
    }

    // tipmap（bitmask）
    const unsigned tipmap[4] = {1u, 2u, 4u, 8u};

    // 右 child CLV（double for GPU）
    // 佈局: [site][rate][state]
    std::vector<double> h_right_clv_d(parent_count);
    for (size_t n = 0; n < P.sites; ++n) {
        for (int r = 0; r < P.rate_cats; ++r) {
            for (int x = 0; x < P.states; ++x) {
                // 可預測的值：100*n + 10*r + x + 1
                h_right_clv_d[n*span + r*P.states + x] = 100.0*n + 10.0*r + (double)x + 1.0;
            }
        }
    }

    // 左/右 P 矩陣（double for GPU）：[rate][from][to]，row-major（from*k + to）
    const size_t pmat_elems = (size_t)P.rate_cats * P.states * P.states;
    std::vector<double> h_left_mat_d(pmat_elems), h_right_mat_d(pmat_elems);
    for (int r = 0; r < P.rate_cats; ++r) {
        for (int from = 0; from < P.states; ++from) {
            for (int to = 0; to < P.states; ++to) {
                size_t idx = (size_t)r*P.states*P.states + (size_t)from*P.states + to;
                // 左邊：依靠 tip mask → 只會累加允許的 from；給個可預測數值
                h_left_mat_d[idx]  = 0.01 * (double)(r+1) + 0.1 * from + 0.001 * to;
                // 右邊：乘上 right_clv；給另一種可預測數值
                h_right_mat_d[idx] = 0.02 * (double)(r+1) + 0.05 * to  + 0.002 * from;
            }
        }
    }

    // ===== GPU（double）=====
    cudaStream_t stream = 0;
    partial_likelihood::Partial_Likelihood_Tip_Inner ti;
    ti.ConstructionOnGpu(
        P,
        h_left.data(),
        tipmap,                 // tipmap on host
        h_right_clv_d.data(),
        h_left_mat_d.data(),
        h_right_mat_d.data(),
        stream);
    ti.UpdatePartialLikelihood(P, stream);

    std::vector<double> gpu_parent(parent_count, 0.0);
    CUDA_CHECK(cudaMemcpy(gpu_parent.data(), ti.d_parent_clv,
                          parent_count * sizeof(double),
                          cudaMemcpyDeviceToHost));

    // ===== CPU baseline（libpll，float 版）=====
    // 建 float 緩衝
    
    std::vector<float> right_clv_f(parent_count);
    std::vector<float> left_mat_f(pmat_elems), right_mat_f(pmat_elems);
    for (size_t i = 0; i < parent_count; ++i) right_clv_f[i] = (float)h_right_clv_d[i];
    for (size_t i = 0; i < pmat_elems; ++i) {
        left_mat_f[i]  = (float)h_left_mat_d[i];
        right_mat_f[i] = (float)h_right_mat_d[i];
    }

    std::vector<float> cpu_parent_f(parent_count, 0.0f);
    std::vector<unsigned int> parent_scaler(P.sites, 0); // 不縮放，仍給零陣（與 API 對齊）
    // right_scaler 可為 nullptr（右子不縮放）
    const unsigned int* right_scaler = nullptr;
    const pll_state_t tipmap_CPU[4] = {1u, 2u, 4u, 8u};
    // 呼叫 libpll：pll_core_update_partial_ti(...)
    // 簽名（float 版）大致為：
    // pll_core_update_partial_ti(states, sites, rate_cats,
    //   parent_clv, parent_scaler,
    //   left_tipchars,
    //   right_clv,
    //   left_matrix, right_matrix,
    //   right_scaler,
    //   tipmap, tipmap_size,
    //   attrib)
    pll_core_update_partial_ti(
        (unsigned int)P.states,
        (unsigned int)P.sites,
        (unsigned int)P.rate_cats,
        cpu_parent_f.data(),
        parent_scaler.data(),
        h_left_bits.data(),
        right_clv_f.data(),
        left_mat_f.data(),
        right_mat_f.data(),
        right_scaler,
        tipmap_CPU,
        tipmap_size,
        false);

    // ===== 比對 CLV =====
    double diff = max_abs_diff_mixed(cpu_parent_f, gpu_parent);
    printf("[TI] max_abs_diff(parent CLV, GPU vs libpll CPU) = %.3e\n", diff);

    // ===== log-likelihood（雙方各算）=====
    double logL_gpu = cpu_loglik_from_parent_clv_double(
        gpu_parent.data(), P.sites, P.states, P.rate_cats, freqs);

    double logL_cpu = cpu_loglik_from_parent_clv_float(
        cpu_parent_f.data(), P.sites, P.states, P.rate_cats, freqs);

    printf("[TI] logL(GPU) = %.12f\n", logL_gpu);
    printf("[TI] logL(CPU) = %.12f\n", logL_cpu);
    printf("[TI] |Δ logL|   = %.3e\n", std::fabs(logL_gpu - logL_cpu));

    ti.CleanUp();
    return 0;
}
