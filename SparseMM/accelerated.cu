/*
Matrix multiplication algorithm using GPU acceleration (naive approach)

Compute e^R using Maclaurin series approximation for 100 iterations. 
*/

#include <iostream>
#include <vector>
#include <cstdlib> // for random number generation
#include <ctime>   // for seeding random number generator

using namespace std;

const int N = 10000; // matrix size N*N

// generate a random sparse 2D matrix
void generateSparseMatrix(int *matrix, double sparsity) {

    for (int i = 0; i < N*N; i++) {
        double randomValue = (double)rand() / RAND_MAX; // Generate a random value between 0 and 1

        if (randomValue > sparsity) {
            // Set a non-zero value in the matrix
            matrix[i] = rand() % 100; // You can adjust the range of non-zero values as needed
        }
        
    }
}

// compute the factorial of a number
long factorial(const int n) {
    long f = 1;
    for (int i=1; i<=n; i++)
        f *= i;
    return f;
}

// CUDA kernel for matrix multiplication
__global__ void matrixMultiply(int* A, int* B, int* C, int n, int block_size) {
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int row = blockIdx.y * blockDim.y + threadIdx.y;

    int sum = 0;
    if (row < n && col < n) {
        for (int i = 0; i < n; i++) {
            sum += A[row * n + i] * B[i * n + col];
        }
        C[row * n + col] = sum;
    }
}

// Function to perform matrix multiplication with grid-striding
void performMatrixMultiplication(int* d_A, int* d_B, int* d_C, int N) {
    dim3 blockDim(16, 16); // Number of threads per block
    dim3 gridDim((N + blockDim.x - 1) / blockDim.x, (N + blockDim.y - 1) / blockDim.y); // Number of blocks

    // Launch the CUDA kernel for matrix multiplication in a grid-striding loop
    for (int i = 0; i < N; i += blockDim.x) {
        for (int j = 0; j < N; j += blockDim.y) {
            matrixMultiply<<<gridDim, blockDim>>>(d_A + i * N, d_B + j, d_C + i * N + j, N, blockDim.x);
            cudaDeviceSynchronize(); // Ensure all previous kernel launches complete before the next one
        }
    }
}

// add two square matrices together
__global__
void mAdd(int *A, int *B, int *result) {
    int indexWithinTheGrid = threadIdx.x + blockIdx.x * blockDim.x;
    int gridStride = gridDim.x * blockDim.x;

    for (int i = indexWithinTheGrid; i < N*N; i += gridStride) {
            result[i] = A[i] + B[i];
    }
}

// divide a vector by a scalar.
__global__
void sDivide(float f, int *matrix) {
    int indexWithinTheGrid = threadIdx.x + blockIdx.x * blockDim.x;
    int gridStride = gridDim.x * blockDim.x;

    for (int i = indexWithinTheGrid; i < N; i += gridStride) {
        matrix[i] = matrix[i]/f;
    }
}

// copy one vector to another
__global__
void mCopy(int *to, int *from) {
    int indexWithinTheGrid = threadIdx.x + blockIdx.x * blockDim.x;
    int gridStride = gridDim.x * blockDim.x;

    for (int i = indexWithinTheGrid; i < N*N; i += gridStride) {
        to[i] = from[i];
    }
}


int main() {
    srand(55); // Seed the random number generator with the current time

    double sparsity = 0.7; // Set the sparsity factor (higher value = more sparse)

    int *h_matrix; // Host vector
    int *d_matrix_original; // Device vector
    int *d_matrix_product; // Device vector
    int *temp_product;

    int *output; // output of computation

    // Allocate memory for host vector
    h_matrix = new int[N*N]();
    generateSparseMatrix(h_matrix, sparsity);

    // Allocate memory for device vectors
    cudaMalloc((void**)&d_matrix_original, N * N * sizeof(int));
    cudaMalloc((void**)&d_matrix_product, N * N * sizeof(int));
    cudaMalloc((void**)&temp_product, N * N * sizeof(int));
    cudaMalloc((void**)&output, N * N * sizeof(int));

    // Copy data from host to device
    cudaMemcpy(d_matrix_original, h_matrix, N * N * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_matrix_product, h_matrix, N * N * sizeof(int), cudaMemcpyHostToDevice);

    // Define thread block and grid dimensions
    dim3 blockDim(256); // Number of threads per block
    dim3 gridDim((N*N + blockDim.x - 1) / blockDim.x); // Number of blocks

    ///////////////////////////////////////////////////////////////////////////
    // compute series
    int *ones = new int[N*N];
    for (int i = 0; i < N*N; i++) {
        ones[i] = 1;
    }

    mAdd<<<gridDim, blockDim>>>(d_matrix_original, ones, output);
    mCopy<<<gridDim, blockDim>>>(d_matrix_product, d_matrix_original);

    for (int i = 2; i <= 100; i++) { // compute up to 100 iterations of the series

        performMatrixMultiplication(temp_product, d_matrix_original, d_matrix_product, N);
        mCopy<<<gridDim, blockDim>>>(d_matrix_product, temp_product);

        sDivide<<<gridDim, blockDim>>>(factorial(i), temp_product);
        mAdd<<<gridDim, blockDim>>>(temp_product, output, output);
    }
    // cudaDeviceSynchronize();
    ///////////////////////////////////////////////////////////////////////////

    // Copy the result from device to host

    cudaMemcpy(h_matrix, output, N *N  * sizeof(int), cudaMemcpyDeviceToHost);
    cudaFree(d_matrix_original);
    cudaFree(d_matrix_product);
    cudaFree(temp_product);
    cudaFree(output);
    delete[] ones;
    cudaDeviceSynchronize();
    // Print the result vector
    for (int i = 0; i < N; i++) {
        std::cout << h_matrix[i] << " ";
    }
    std::cout << std::endl;

    // Free device and host memory
    
    delete[] h_matrix;
    
   

    return 0;
}
