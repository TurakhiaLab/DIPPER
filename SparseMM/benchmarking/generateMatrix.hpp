#include <cstdlib>            // std::rand
#include <cuda_fp16.h>

const int N = 3200; // dimension of square matrix
const int runs = 100;

void generate(__half *hA, __half *hB) {
    srand(55);
    // generate matrix with forced 2:4 sparsity
    for (size_t i = 0; i < N * N; i++)                                            
        if (i%2) {
            hA[i] = static_cast<__half>(static_cast<float>(std::rand() % 10));
        } else {
            hA[i] = 0;
        }
    // generate dense matrix
    for (size_t i = 0; i < N * N; i++)
        hB[i] = static_cast<__half>(static_cast<float>(std::rand() % 10));
}

void generate(float *hA, float *hB) {
    srand(55);
    // generate matrix with forced 2:4 sparsity
    for (size_t i = 0; i < N * N; i++)                                            
        if (i%2) {
            hA[i] = static_cast<float>(std::rand() % 10);
        } else {
            hA[i] = 0;
        }
    // generate dense matrix
    for (size_t i = 0; i < N * N; i++)                                            
        hB[i] = static_cast<float>(std::rand() % 10);
}
