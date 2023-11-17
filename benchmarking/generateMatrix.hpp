#include <cstdlib>            // std::rand
#include <cuda_fp16.h>

const int N = 320; // dimension of square matrix
void generate(__half *hA, __half *hB) {
    srand(55);
    for (int i = 0; i < N * N; i++)                                            // generate matrix with forced 2:4 sparsity
        if (i%2) {
            hA[i] = static_cast<__half>(static_cast<float>(std::rand() % 10));
        } else {
            hA[i] = 0;
        }
    for (int i = 0; i < N * N; i++)                                            // generate dense matrix
        hB[i] = static_cast<__half>(static_cast<float>(std::rand() % 10));
}
