#include <cstdio>
#include <vector>
#include <cmath>
#include <chrono>
#include <cuda_runtime.h>
#include "core_likelihood.cuh"
#include <libpll/pll.h>

#ifndef CUDA_CHECK
#define CUDA_CHECK(call) do { \
  cudaError_t _e = (call); \
  if (_e != cudaSuccess) { \
    fprintf(stderr, "[CUDA] %s failed: %s (%d) @ %s:%d\n", \
            #call, cudaGetErrorString(_e), int(_e), __FILE__, __LINE__); \
    std::exit(1); \
  } \
} while(0)
#endif

int main() {
  using namespace core_likelihood;

  const unsigned sites     = 8192;
  const int      states    = 20;
  const int      rate_cats = 4;

  // double inputs (GPU)
  std::vector<double> root_clv_d(size_t(sites) * size_t(rate_cats) * size_t(states));
  std::vector<unsigned> site_scaler(sites);
  std::vector<double> frequencies_d(states, 0.25);
  std::vector<double> rate_weights_d(rate_cats, 1.0 / double(rate_cats));
  std::vector<unsigned> pattern_weights(sites, 1u);
  std::vector<int> invar_indices(sites, -1);

  for (unsigned i = 0; i < sites; ++i) {
    double s = 0.0;
    for(int j = 0; j < rate_cats; ++j) {
      for (int k = 0; k < states; ++k) {
          double v = 0.5 + 0.5 * std::sin((i * states * rate_cats + j * states + k) * 0.01);
          root_clv_d[(size_t)i * rate_cats * states + j * states + k] = v;
          s += v;
      }
      for (int k = 0; k < states; ++k) {
          root_clv_d[(size_t)i * rate_cats * states + j * states + k] /= s;
      }
    }
    for (int k = 0; k < states; ++k) root_clv_d[i*states + k] /= s;
    site_scaler[i] = i % 3;
  }

  // GPU path
  Param p;
  p.sites = sites;
  p.states = states;
  p.rate_cats = rate_cats;
  p.invar_proportion = 0.0;
  p.per_rate_scaling = false;

  Likelihood_Root lr;
  lr.Initialize(p);


  lr.ConstructionOnGpu(
      root_clv_d.data(),
      frequencies_d.data(),
      rate_weights_d.data(),
      pattern_weights.data(),
      site_scaler.data(),
      invar_indices.data(),
      p);
      
  cudaEvent_t start, stop;
  CUDA_CHECK(cudaEventCreate(&start));
  CUDA_CHECK(cudaEventCreate(&stop));
  CUDA_CHECK(cudaEventRecord(start));

  lr.ComputeLikelihood(p);

  CUDA_CHECK(cudaEventRecord(stop));
  CUDA_CHECK(cudaEventSynchronize(stop));
  float ms_gpu = 0.f;
  CUDA_CHECK(cudaEventElapsedTime(&ms_gpu, start, stop));
  CUDA_CHECK(cudaEventDestroy(start));
  CUDA_CHECK(cudaEventDestroy(stop));

  double gpu_ll = lr.total_loglik[0];

  // Convert to float for PLL float build
  std::vector<float> root_clv_f(root_clv_d.begin(), root_clv_d.end());
  std::vector<float> frequencies_f(frequencies_d.begin(), frequencies_d.end());
  std::vector<float> rate_weights_f(rate_weights_d.begin(), rate_weights_d.end());

  // Build float* const* frequencies (1 mixture)
  float* mix0_ptr = frequencies_f.data();
  float* const* frequencies_pp = &mix0_ptr;
  std::vector<unsigned> freqs_indices(rate_cats, 0u);

  // PLL float
  auto t1 = std::chrono::high_resolution_clock::now();

  double pll_ll = pll_core_root_loglikelihood(
      (unsigned)states,
      (unsigned)sites,
      (unsigned)rate_cats,
      (const float*)root_clv_f.data(),
      (const unsigned*)site_scaler.data(),
      (float* const*)frequencies_pp,
      (const float*)rate_weights_f.data(),
      (const unsigned*)pattern_weights.data(),
      (const float*)nullptr,
      (const int*)invar_indices.data(),
      (const unsigned*)freqs_indices.data(),
      (float*)nullptr,
      0u);

  auto t2 = std::chrono::high_resolution_clock::now();
  double ms_pll = std::chrono::duration<double, std::milli>(t2 - t1).count();

  printf("[GPU ] logL = %.12f, time = %.3f ms\n", gpu_ll, ms_gpu);
  printf("[PLLf] logL = %.12f, time = %.3f ms\n", pll_ll, ms_pll);
  printf("[DIFF] abs = %.6e\n", std::fabs(gpu_ll - pll_ll));
  printf("[TIME] GPU / PLLf = %.3f\n", ms_gpu / ms_pll);
  lr.CleanUp();
  return 0;
}
