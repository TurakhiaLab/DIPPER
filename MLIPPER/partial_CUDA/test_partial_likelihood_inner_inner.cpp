#include <cstdio>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cuda_runtime.h>
#include <libpll/pll.h>
#include "partial_likelihood.cuh"

// ---------- 小工具 ----------
static double max_abs_diff_mixed(const std::vector<float>& a, const std::vector<double>& b) {
    const size_t n = std::min(a.size(), b.size());
    double m = 0.0;
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

// ---------- 主程式：II 對拍 ----------
int main() {
    // 參數（可改）
    partial_likelihood::Param P;
    P.sites = 16;
    P.states = 4;        // DNA
    P.rate_cats = 4;
    P.per_rate_scaling = false; // 本測試不做縮放

    const size_t span = (size_t)P.states * (size_t)P.rate_cats;
    const size_t clv_elems = (size_t)P.sites * span;
    const size_t pmat_elems = (size_t)P.rate_cats * (size_t)P.states * (size_t)P.states;

    // 均勻頻率
    std::vector<double> freqs(P.states, 1.0 / P.states);

    // ===== 準備 GPU 輸入（double）=====
    // 佈局: CLV = [site][rate][state] -> n*span + r*states + k
    std::vector<double> h_left_clv_d(clv_elems), h_right_clv_d(clv_elems);
    for (size_t n = 0; n < P.sites; ++n) {
        for (int r = 0; r < P.rate_cats; ++r) {
            for (int k = 0; k < P.states; ++k) {
                const size_t idx = n*span + r*P.states + k;
                h_left_clv_d[idx]  =  0.1 * (double)(n+1) + 0.01 * (double)(r+1) + 0.001 * (double)(k+1);
                h_right_clv_d[idx] =  0.2 * (double)(n+1) + 0.02 * (double)(r+1) + 0.002 * (double)(k+1);
            }
        }
    }

    // P 矩陣（double），佈局: [rate][row=j][col=k] -> r*states*states + j*states + k（row-major）
    std::vector<double> h_left_mat_d(pmat_elems), h_right_mat_d(pmat_elems);
    for (int r = 0; r < P.rate_cats; ++r) {
        for (int j = 0; j < P.states; ++j) {
            for (int k = 0; k < P.states; ++k) {
                const size_t idx = (size_t)r*P.states*P.states + (size_t)j*P.states + k;
                // 給可預測值（非隨機）
                h_left_mat_d[idx]  = 0.05 * (r+1) + 0.1 * j + 0.001 * k;
                h_right_mat_d[idx] = 0.04 * (r+1) + 0.07 * j + 0.002 * k;
            }
        }
    }

    // ===== GPU 執行 =====
    cudaStream_t stream = 0;
    partial_likelihood::Partial_Likelihood_Inner_Inner ii;
    ii.ConstructionOnGpu(
        P,
        h_left_clv_d.data(),
        h_right_clv_d.data(),
        h_left_mat_d.data(),
        h_right_mat_d.data(),
        stream);
    ii.UpdatePartialLikelihood(P, stream);

    std::vector<double> gpu_parent(clv_elems, 0.0);
    {
        cudaError_t err = cudaMemcpy(gpu_parent.data(), ii.d_parent_clv,
                                     clv_elems * sizeof(double),
                                     cudaMemcpyDeviceToHost);
        if (err != cudaSuccess) {
            fprintf(stderr, "[CUDA] memcpy D2H parent failed: %s\n", cudaGetErrorString(err));
            return 1;
        }
    }

    // ===== CPU baseline（libpll，float 版）=====
    // libpll API 使用 float，所以轉成 float 測 baseline
    std::vector<float> left_clv_f(clv_elems), right_clv_f(clv_elems);
    std::vector<float> left_mat_f(pmat_elems), right_mat_f(pmat_elems);
    for (size_t i = 0; i < clv_elems; ++i) {
        left_clv_f[i]  = (float)h_left_clv_d[i];
        right_clv_f[i] = (float)h_right_clv_d[i];
    }
    for (size_t i = 0; i < pmat_elems; ++i) {
        left_mat_f[i]  = (float)h_left_mat_d[i];
        right_mat_f[i] = (float)h_right_mat_d[i];
    }

    std::vector<float> cpu_parent_f(clv_elems, 0.0f);
    // 不做縮放 → scaler 皆為 nullptr，attrib 用 NONE
    const unsigned int* left_scaler  = nullptr;
    const unsigned int* right_scaler = nullptr;

    // 注意：II 的 libpll 介面不需要 tipmap/left_tipchars
    // 原型通常為：
    // pll_core_update_partial_ii(states, sites, rate_cats,
    //     parent_clv, parent_scaler,  // parent_scaler 可給 nullptr 或零陣（這裡給 nullptr）
    //     left_clv, right_clv,
    //     left_matrix, right_matrix,
    //     left_scaler, right_scaler,
    //     attrib);
    pll_core_update_partial_ii(
        (unsigned int)P.states,
        (unsigned int)P.sites,
        (unsigned int)P.rate_cats,
        cpu_parent_f.data(),
        /* parent_scaler */ nullptr,
        left_clv_f.data(),
        right_clv_f.data(),
        left_mat_f.data(),
        right_mat_f.data(),
        left_scaler,
        right_scaler,
        false);

    // ===== 比對 CLV =====
    const double diff = max_abs_diff_mixed(cpu_parent_f, gpu_parent);
    printf("[II] max_abs_diff(parent CLV, GPU vs libpll CPU) = %.3e\n", diff);

    // ===== 以均勻頻率計 logL（雙方一致的聚合方式）=====
    const double logL_gpu = cpu_loglik_from_parent_clv_double(gpu_parent.data(),
                             P.sites, P.states, P.rate_cats, freqs);
    const double logL_cpu = cpu_loglik_from_parent_clv_float(cpu_parent_f.data(),
                             P.sites, P.states, P.rate_cats, freqs);
    printf("[II] logL(GPU) = %.12f\n", logL_gpu);
    printf("[II] logL(CPU) = %.12f\n", logL_cpu);
    printf("[II] |Δ logL|   = %.3e\n", std::fabs(logL_gpu - logL_cpu));

    ii.CleanUp();
    return 0;
}
