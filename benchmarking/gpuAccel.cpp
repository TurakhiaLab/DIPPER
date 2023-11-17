#include <cuda_runtime_api.h> // cudaMalloc, cudaMemcpy, etc.
#include "generateMatrix.hpp"


__global__ 
void matrixMultiply(const __half* A, const __half* B, __half* C, int n) {
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (row < n && col < n) {
        __half sum = __float2half(0.0f);
        for (int i = 0; i < n; ++i) {
            __half a = A[row * n + i];
            __half b = B[i * n + col];
            sum = __hadd2(sum, __hmul2(a, b));
        }
        C[row * n + col] = sum;
    }
}

int main() {
    int matrix_size = N * N;

    __half* h_A[matrix_size];
    __half* h_B[matrix_size]; 
    __half* result = new __half[matrix_size]();

    generateMatrix(h_A, h_B);
    
    __half *d_A, *d_B, *d_C;
    cudaMalloc((void**)&d_A, matrix_size * sizeof(__half));
    cudaMalloc((void**)&d_B, matrix_size * sizeof(__half));
    cudaMalloc((void**)&d_C, matrix_size * sizeof(__half));

    cudaMemcpy(d_A, h_A, matrix_size * sizeof(__half), cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, h_B, matrix_size * sizeof(__half), cudaMemcpyHostToDevice);
    
    dim3 blockSize(16, 16);
    dim3 gridSize((N + blockSize.x - 1) / blockSize.x, (N + blockSize.y - 1) / blockSize.y);

    matrixMultiply<<<gridSize, blockSize>>>(d_A, d_B, d_C, N);

    cudaMemcpy(result, d_C, matrix_size * sizeof(__half), cudaMemcpyDeviceToHost);

    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);

    return 0;
}