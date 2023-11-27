#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "generateMatrix.hpp"

__global__ void matrixMultiply(const float* A, const float* B, float* C, int n) {
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (row < n && col < n) {
        float sum = 0.0f;
        for (int i = 0; i < n; ++i) {
            float a = A[row * n + i];
            float b = B[i * n + col];
            sum += a*b;
        }
        C[row * n + col] = sum;
    }
}

int main() {
    int matrix_size = N * N;

    float *h_A = new float[matrix_size];
    float *h_B = new float[matrix_size]; 
    float *result = new float[matrix_size];

    generate(h_A, h_B);

    for (int iterations = 0; iterations < runs; iterations++) {
        float *d_A, *d_B, *d_C;
        cudaMalloc((void**)&d_A, matrix_size * sizeof(float));
        cudaMalloc((void**)&d_B, matrix_size * sizeof(float));
        cudaMalloc((void**)&d_C, matrix_size * sizeof(float));

        cudaMemcpy(d_A, h_A, matrix_size * sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_B, h_B, matrix_size * sizeof(float), cudaMemcpyHostToDevice);
        
        dim3 blockSize(16, 16);
        dim3 gridSize((N + blockSize.x - 1) / blockSize.x, (N + blockSize.y - 1) / blockSize.y);

        //for (int iterations = 0; iterations < runs; iterations++) {
            matrixMultiply<<<gridSize, blockSize>>>(d_A, d_B, d_C, N);
        //}
        
        cudaDeviceSynchronize();

        cudaMemcpy(result, d_C, matrix_size * sizeof(float), cudaMemcpyDeviceToHost);

        cudaFree(d_A);
        cudaFree(d_B);
        cudaFree(d_C);
    }
    return 0;
}
