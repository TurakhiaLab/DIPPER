#include <cstdlib>            // std::rand
#include <cuda_fp16.h>

const int N = 32000; // dimension of square matrix
void generate(__half *hA, __half *hB) {
    srand(55);
    // generate matrix with forced 2:4 sparsity
    for (int i = 0; i < N * N; i++)                                            
        if (i%2) {
            hA[i] = static_cast<__half>(static_cast<float>(std::rand() % 10));
        } else {
            hA[i] = 0;
        }
    // generate dense matrix
    for (int i = 0; i < N * N; i++)
        hB[i] = static_cast<__half>(static_cast<float>(std::rand() % 10));
}
