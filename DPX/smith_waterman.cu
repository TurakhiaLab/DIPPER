#include <iostream>
#include <vector>
#include <cuda_runtime.h>

const int match_score = 3;
const int mismatch_penalty = -3;
const int gap_penalty = -2;


// CUDA kernel for Smith-Waterman
__global__ 
void smithWatermanKernel(const char* seq1, const char* seq2, int* score_matrix, 
                        int* traceback_matrix, int m, int n) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    int iter = __vimax3_s32(i, j, 0); // max iterations needed for cell

    if (i < m && j < n) {
        for (int k = 0; k < n; k++) {
            if (i == 0 || j == 0) {
                score_matrix[i * n + j] = 0;
                traceback_matrix[i * n + j] = 0;
            } else {
                if (k > iter) {
                    break;
                } 
                else {
                    int match = (seq1[i - 1] == seq2[j - 1]) ? match_score : mismatch_penalty;

                    int diagonal = score_matrix[(i - 1) * n + (j - 1)] + match;
                    int up = score_matrix[(i - 1) * n + j] + gap_penalty;
                    int left = score_matrix[i * n + (j - 1)] + gap_penalty;

                    int max_score = __vimax3_s32(diagonal, up, left);

                    score_matrix[i * n + j] = max_score;

                    if (max_score == diagonal) {
                        traceback_matrix[i * n + j] = 1;  // diagonal
                    } else if (max_score == up) {
                        traceback_matrix[i * n + j] = 2;  // up
                    } else {
                        traceback_matrix[i * n + j] = 3;  // left
                    }
                } 
            } 
        }
    }
}

// Host function to run Smith-Waterman on GPU
void smithWatermanGPU(const std::vector<char>& seq1, const std::vector<char>& seq2,
                      int* score_matrix, int* traceback_matrix, int m, int n) {
    char *d_seq1, *d_seq2;
    int *d_score_matrix, *d_traceback_matrix;

    // Allocate device memory
    cudaMalloc((void**)&d_seq1, m * sizeof(char));
    cudaMalloc((void**)&d_seq2, n * sizeof(char));
    cudaMalloc((void**)&d_score_matrix, m * n * sizeof(int));
    cudaMalloc((void**)&d_traceback_matrix, m * n * sizeof(int));

    // Copy data from host to device
    cudaMemcpy(d_seq1, seq1.data(), m * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(d_seq2, seq2.data(), n * sizeof(char), cudaMemcpyHostToDevice);

    // Set grid and block dimensions
    dim3 blockSize(16, 16);
    dim3 gridSize((m + blockSize.x - 1) / blockSize.x, (n + blockSize.y - 1) / blockSize.y);

    // Launch the kernel
    smithWatermanKernel<<<gridSize, blockSize>>>(d_seq1, d_seq2, d_score_matrix, d_traceback_matrix, m, n);

    // Copy results from device to host
    cudaMemcpy(score_matrix, d_score_matrix, m * n * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(traceback_matrix, d_traceback_matrix, m * n * sizeof(int), cudaMemcpyDeviceToHost);

    // Free device memory
    cudaFree(d_seq1);
    cudaFree(d_seq2);
    cudaFree(d_score_matrix);
    cudaFree(d_traceback_matrix);
}

int main() {
    // Example sequences
    std::vector<char> seq1 = {'A', 'C', 'G', 'T', 'A', 'C', 'C', 'A', 'G', 'T'};
    std::vector<char> seq2 = {'A', 'C', 'C', 'A', 'G', 'T', 'A', 'C', 'G', 'T'};

    // Get the size of the sequences
    int m = seq1.size() + 1;
    int n = seq2.size() + 1;

    // Allocate memory for the score and traceback matrices
    std::vector<int> score_matrix(m * n, 0);
    std::vector<int> traceback_matrix(m * n, 0);

    // Run Smith-Waterman on the GPU
    smithWatermanGPU(seq1, seq2, score_matrix.data(), traceback_matrix.data(), m, n);

    // Print the resulting score matrix and traceback matrix

    std::cout << " \t\t";

    for (char ch : seq2) {
        std::cout << ch << '\t';
    }

    std::cout << std::endl;

    for (int i = 0; i < m; i++) {
        if (i > 0) 
            std::cout << seq1.at(i-1) << '\t';
        else 
            std::cout << '\t';
        for (int j = 0; j < n; j++) {
            std::cout << score_matrix[i * n + j] << '\t';
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;

    std::cout << " \t\t";

    for (char ch : seq2) {
        std::cout << ch << '\t';
    }

    std::cout << std::endl;

    for (int i = 0; i < m; i++) {
        if (i > 0) 
            std::cout << seq1.at(i-1) << '\t';
        else 
            std::cout << '\t';
        for (int j = 0; j < n; j++) {
            std::cout << traceback_matrix[i * n + j] << '\t';
        }
        std::cout << std::endl;
    }

    return 0;
}