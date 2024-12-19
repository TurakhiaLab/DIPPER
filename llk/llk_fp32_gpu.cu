#ifndef CORE_LIKELIHOOD_FP32
#include "core_likelihood_fp32.hpp"
#endif

#include "core_likelihood_fp32.cpp"
#include "tree_fp32.cpp"
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <cublas_v2.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unordered_map>
#include <cmath> // For exp function on CPU
#include <device_launch_parameters.h>

#include <mma.h>
#include <cassert>

using namespace nvcuda;

// error checking
#define CHECK_CUDA(call)                                                                                    \
    {                                                                                                       \
        cudaError_t err = call;                                                                             \
        if (err != cudaSuccess)                                                                             \
        {                                                                                                   \
            fprintf(stderr, "CUDA Error: %s (code %d), line %d\n", cudaGetErrorString(err), err, __LINE__); \
            exit(err);                                                                                      \
        }                                                                                                   \
    }

#define CHECK_CUSOLVER(call)                                                 \
    {                                                                        \
        cusolverStatus_t err = call;                                         \
        if (err != CUSOLVER_STATUS_SUCCESS)                                  \
        {                                                                    \
            fprintf(stderr, "cuSOLVER Error: %d, line %d\n", err, __LINE__); \
            exit(err);                                                       \
        }                                                                    \
    }

#define CHECK_CUBLAS(call)                                                 \
    {                                                                      \
        cublasStatus_t err = call;                                         \
        if (err != CUBLAS_STATUS_SUCCESS)                                  \
        {                                                                  \
            fprintf(stderr, "cuBLAS Error: %d, line %d\n", err, __LINE__); \
            exit(err);                                                     \
        }                                                                  \
    }


typedef float exp_f;
typedef float llk_f;
typedef double scaling_f;

constexpr int STATE_COUNT = 5;
constexpr int MAX_CHILD = 3; // 2
// constexpr int IDENTIFIER_LENGTH = 100;
constexpr int NODE_COUNT = 72;        // 72
constexpr int SEQUENCE_LENGTH = 5794; // 5794
constexpr int siteNum = 1;

// GPU NODE STRUCTURE
struct GPUNode
{
    int parentIndex;
    int childrenIndices[MAX_CHILD];
    int numChildren;
    // char identifier[IDENTIFIER_LENGTH];
    llk_f branchLength;
    llk_f bottom[STATE_COUNT][SEQUENCE_LENGTH];
    scaling_f scaleVector[SEQUENCE_LENGTH];
    llk_f probabilityMatrix[STATE_COUNT][STATE_COUNT];
    char sequence[SEQUENCE_LENGTH];
};

// GPU TREE STRUCTURE
struct GPUTree
{
    GPUNode nodes[NODE_COUNT];
    exp_f rateMatrix[STATE_COUNT * STATE_COUNT];
    int order[NODE_COUNT];
};

__global__ void float_to_double(exp_f *dest, const llk_f *src, int count)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < count)
    {
        dest[idx] = static_cast<exp_f>(src[idx]);
    }
}

__global__ void double_to_float(llk_f *dest, const exp_f *src, int count)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < count)
    {
        dest[idx] = static_cast<llk_f>(src[idx]);
    }
}

__global__ void exp_diagonal_elements(exp_f *d_diag, int n)
{
    int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx < n)
    {
        d_diag[idx] = exp(d_diag[idx]);
    }
}

__global__ void scale_matrix(exp_f *d_matrix, exp_f *d_dest_matrix, llk_f scalar, int n)
{
    int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx < n * n)
    {
        d_dest_matrix[idx] = d_matrix[idx] * scalar;
    }
}


void copyNode(utility::Node *srcNode, GPUTree &destTree, int parentIndex, int &nodeIndex, std::unordered_map<std::string, int> &indexMap)
{
    if (!srcNode) {
        return;
    }

    // Assign current node index and update global node index
    int currentIndex = nodeIndex++;
    assert(currentIndex < NODE_COUNT);

    // Create a GPUNode and copy data
    GPUNode &gpuNode = destTree.nodes[currentIndex];
    gpuNode.parentIndex = parentIndex;

    gpuNode.numChildren = srcNode->children.size();
    gpuNode.branchLength = srcNode->branchLength;

    std::string seq = utility::seqs[srcNode->identifier].second;
    for (int i = 0; i < seq.size(); ++i)
    {
        gpuNode.sequence[i] = seq[i];
    }

    indexMap[srcNode->identifier] = currentIndex;

    // std::cout << "Node " << srcNode->identifier << " mapped to index " << currentIndex << std::endl;

    // Copy children nodes
    for (int i = 0; i < srcNode->children.size(); ++i)
    {

        gpuNode.childrenIndices[i] = nodeIndex;

        copyNode(srcNode->children[i], destTree, currentIndex, nodeIndex, indexMap);
    }
}

void copyTree(utility::Tree &sourceTree, GPUTree &destTree)
{
    // Map to store the node identifier to index mapping for the destination tree
    for (int i = 0; i < STATE_COUNT; ++i)
    {
        for (int j = 0; j < STATE_COUNT; ++j)
        {
            destTree.rateMatrix[i * STATE_COUNT + j] = static_cast<exp_f>(utility::rate_matrix[i][j]);
        }
    }

    std::unordered_map<std::string, int> indexMap;

    // Start node index from 0
    int nodeIndex = 0;

    // Start the recursive copy from the root
    copyNode(sourceTree.root, destTree, -1, nodeIndex, indexMap);
}

__device__ void initialize_leaf(GPUNode &node)
{
    for (int i = 0; i < SEQUENCE_LENGTH; i++)
    {
        for (int j = 0; j < STATE_COUNT; j++)
        {
            node.bottom[j][i] = 0;
        }

        char nucleotide = node.sequence[i];

        if (nucleotide == 'A')
        {
            node.bottom[0][i] = 1;
        }
        else if (nucleotide == 'C')
        {
            node.bottom[1][i] = 1;
        }
        else if (nucleotide == 'G')
        {
            node.bottom[2][i] = 1;
        }
        else if (nucleotide == 'T')
        {
            node.bottom[3][i] = 1;
        }
        else
        {
            if (STATE_COUNT == 5)
            {
                node.bottom[4][i] = 1;
            }
        }
        node.scaleVector[i] = 1.0;
    }
}

__device__ void initialize_leaf_per_site(GPUNode &node, int seqPosStart)
{
    for (int j = 0; j < STATE_COUNT; j++)
    {
        node.bottom[j][seqPosStart] = 0;
    }

    char nucleotide = node.sequence[seqPosStart];

    if (nucleotide == 'A')
    {
        node.bottom[0][seqPosStart] = 1;
    }
    else if (nucleotide == 'C')
    {
        node.bottom[1][seqPosStart] = 1;
    }
    else if (nucleotide == 'G')
    {
        node.bottom[2][seqPosStart] = 1;
    }
    else if (nucleotide == 'T')
    {
        node.bottom[3][seqPosStart] = 1;
    }
    else
    {
        if (STATE_COUNT == 5)
        {
            node.bottom[4][seqPosStart] = 1;
        }
    }
    node.scaleVector[seqPosStart] = 1.0;
}

__global__ void non_recursive_dfs(GPUNode *nodes, int sequenceLength)
{
    int seqPosStart = (blockIdx.x * blockDim.x + threadIdx.x) * siteNum;
    if (seqPosStart >= sequenceLength)
        return;

    int stack[NODE_COUNT];
    float results[STATE_COUNT][siteNum] = {0};                 // Initialize with zero
    float childResults[MAX_CHILD][STATE_COUNT][siteNum] = {0}; // Initialize with zero
    bool ready[NODE_COUNT] = {0};                              // Initialize with false

    int stackSize = 0;
    stack[stackSize++] = 0; // root in stack

    while (stackSize > 0)
    {
        int nodeIdx = stack[stackSize - 1]; // peek
        GPUNode &node = nodes[nodeIdx];

        if (!ready[nodeIdx])
        {
            if (node.numChildren == 0)
            { // leaf
                initialize_leaf(node);
                ready[nodeIdx] = true;
            }
            else
            {
                bool childrenReady = true;
                for (int i = 0; i < node.numChildren; ++i)
                {
                    int childIdx = node.childrenIndices[i];
                    if (!ready[childIdx])
                    {
                        stack[stackSize++] = childIdx; // child in stack
                        childrenReady = false;
                    }
                }
                if (childrenReady)
                {
                    ready[nodeIdx] = true;
                }
                else
                {
                    continue;
                }
            }
        }

        if (ready[nodeIdx])
        {
            if (node.numChildren > 0)
            {
                for (int s = 0; s < siteNum; ++s)
                {
                    for (int i = 0; i < STATE_COUNT; ++i)
                    {
                        results[i][s] = 1.0f;
                    }
                }

                for (int i = 0; i < node.numChildren; ++i)
                {
                    int childIdx = node.childrenIndices[i];

                    for (int s = 0; s < siteNum; ++s)
                    {
                        int seqPos = seqPosStart + s;
                        if (seqPos >= sequenceLength)
                            break;

                        for (int j = 0; j < STATE_COUNT; ++j)
                        {
                            childResults[i][j][s] = 0.0f; // Ensure initialization
                            for (int k = 0; k < STATE_COUNT; ++k)
                            {
                                childResults[i][j][s] += nodes[childIdx].probabilityMatrix[j][k] * nodes[childIdx].bottom[k][seqPos];
                            }
                        }
                    }
                }
                for (int s = 0; s < siteNum; ++s)
                {
                    int seqPos = seqPosStart + s;
                    if (seqPos >= sequenceLength)
                    {
                        break;
                    }

                    float sumLikelihoods = 0.0f;

                    for (int i = 0; i < STATE_COUNT; ++i)
                    {
                        for (int j = 0; j < node.numChildren; ++j)
                        {
                            if (childResults[j][i][s] == 0)
                            {
                                continue;
                            }
                            results[i][s] *= childResults[j][i][s];
                        }
                        sumLikelihoods += results[i][s];
                    }

                    node.scaleVector[seqPos] = sumLikelihoods; // Store scaling factor
                    for (int i = 0; i < STATE_COUNT; ++i)
                    {
                        node.bottom[i][seqPos] = results[i][s] / sumLikelihoods; // Normalize result
                    }
                }
            }
            stackSize--;
        }
    }
}
__global__ void felsenstein_pruning_kernel(GPUNode *nodes, int *order, int sequenceLength, int nodeCount)
{
    int seqPosStart = (blockIdx.x * blockDim.x + threadIdx.x) * siteNum;
    if (seqPosStart >= sequenceLength)
        return;

    const int MAX_STACK_SIZE = 15; // Maximum depth of the tree
    const int MAX_SITE_NUM = siteNum; // Number of sites processed per thread
    float bottomValuesStack[MAX_STACK_SIZE][STATE_COUNT][MAX_SITE_NUM];
    float results[STATE_COUNT][MAX_SITE_NUM]; // Stores intermediate results
    int sp = 0; // Stack pointer

    for (int idx = 0; idx < nodeCount; idx++)
    {
        int nodeIdx = order[idx];
        GPUNode &node = nodes[nodeIdx];

        if (node.numChildren == 0)
        { // Leaf node
            for (int s = 0; s < siteNum; ++s)
            {
                int seqPos = seqPosStart + s;
                if (seqPos >= sequenceLength)
                    break;

                // Initialize bottomValues for leaf nodes
                for (int j = 0; j < STATE_COUNT; j++)
                {
                    bottomValuesStack[sp][j][s] = 0.0f;
                }

                char nucleotide = node.sequence[seqPos];

                if (nucleotide == 'A')
                {
                    bottomValuesStack[sp][0][s] = 1.0f;
                }
                else if (nucleotide == 'C')
                {
                    bottomValuesStack[sp][1][s] = 1.0f;
                }
                else if (nucleotide == 'G')
                {
                    bottomValuesStack[sp][2][s] = 1.0f;
                }
                else if (nucleotide == 'T')
                {
                    bottomValuesStack[sp][3][s] = 1.0f;
                }
                else
                {
                    if (STATE_COUNT == 5)
                    {
                        bottomValuesStack[sp][4][s] = 1.0f;
                    }
                }

                // Set scaling factor to 1.0 in global memory
                node.scaleVector[seqPos] = 1.0;
            }
            // Push bottomValues onto the stack
            sp++;
        }
        else
        { // Internal node
            // Pop children's bottomValues from the stack
            sp -= node.numChildren;

            for (int s = 0; s < siteNum; ++s)
            {
                int seqPos = seqPosStart + s;
                if (seqPos >= sequenceLength)
                    break;

                // Initialize results to 1.0
                for (int j = 0; j < STATE_COUNT; ++j)
                {
                    results[j][s] = 1.0f;
                }
                
                // Compute the product of children's likelihoods
                for (int i = 0; i < node.numChildren; ++i)
                {
                    int childIdx = node.childrenIndices[i];

                    for (int j = 0; j < STATE_COUNT; ++j)
                    {
                        float sum = 0.0f;
                        for (int k = 0; k < STATE_COUNT; ++k)
                        {
                            sum += nodes[childIdx].probabilityMatrix[j][k] * bottomValuesStack[sp + i][k][s];
                        }
                        results[j][s] *= sum;
                    }
                }

                float sumLikelihoods = 0.0f;
                
                int exponents[STATE_COUNT];
                int minExp = INT_MAX;
                int maxExp = INT_MIN;

                for (int j = 0; j < STATE_COUNT; ++j)
                {
                    sumLikelihoods += results[j][s];
                    frexpf(results[j][s], &exponents[j]);
                    maxExp = max(exponents[j], maxExp);
                }
                scaling_f scaleFactor = pow(2,maxExp); // 2^{maxExp}

                // Store scaling factor in global memory

                // node.scaleVector[seqPos] = scaleFactor;
                node.scaleVector[seqPos] *= sumLikelihoods;

                // Normalize and store bottomValues in the stack
                for (int j = 0; j < STATE_COUNT; ++j)
                {
                    bottomValuesStack[sp][j][s] = results[j][s] / sumLikelihoods;
                    // bottomValuesStack[sp][j][s] = results[j][s] / scaleFactor;

                }
                

            }
            // Push the current node's bottomValues onto the stack
            sp++;
        }
    }

    // After processing all nodes, write the root's bottomValues to global memory
    int rootIdx = order[nodeCount - 1];
    GPUNode &rootNode = nodes[rootIdx];

    for (int s = 0; s < siteNum; ++s)
    {
        int seqPos = seqPosStart + s;
        if (seqPos >= sequenceLength)
            break;

        // 'scaleVector' is already stored; copy 'bottom' values to global memory
        for (int i = 0; i < STATE_COUNT; ++i)
        {
            rootNode.bottom[i][seqPos] = bottomValuesStack[sp - 1][i][s];
        }
    }
}



__global__ void felsenstein_pruning_kernel1(GPUNode *nodes, int *order, int sequenceLength, int nodeCount)
{
    int seqPosStart = (blockIdx.x * blockDim.x + threadIdx.x) * siteNum;
    if (seqPosStart >= sequenceLength)
        return;

    // each thread processes a site
    for (int idx = 0; idx < nodeCount; idx++)
    {

        int nodeIdx = order[idx];
        GPUNode &node = nodes[nodeIdx];

        if (node.numChildren == 0)
        { // Leaf node
            for (int s = 0; s < siteNum; ++s)
            {
                int seqPos = seqPosStart + s;
                if (seqPos >= sequenceLength)
                    break;
                initialize_leaf_per_site(node, seqPos);
            }
        }
        else
        {                                                              // Internal node
            float results[STATE_COUNT][siteNum] = {0};                 // Initialize with zero
            float childResults[MAX_CHILD][STATE_COUNT][siteNum] = {0}; // Initialize with zero

            for (int s = 0; s < siteNum; ++s)
            {
                for (int i = 0; i < STATE_COUNT; ++i)
                {
                    results[i][s] = 1.0f;
                }
            }

            for (int i = 0; i < node.numChildren; ++i)
            {
                int childIdx = node.childrenIndices[i];

                for (int s = 0; s < siteNum; ++s)
                {
                    int seqPos = seqPosStart + s;
                    if (seqPos >= sequenceLength)
                        break;

                    for (int j = 0; j < STATE_COUNT; ++j)
                    {
                        childResults[i][j][s] = 0.0f; // Ensure initialization
                        for (int k = 0; k < STATE_COUNT; ++k)
                        {
                            childResults[i][j][s] += nodes[childIdx].probabilityMatrix[j][k] * nodes[childIdx].bottom[k][seqPos];
                        }
                    }
                }
            }

//             for (int s = 0; s < siteNum; ++s)
//             {
//                 int seqPos = seqPosStart + s;
//                 if (seqPos >= sequenceLength)
//                     break;

//                 float sumLikelihoods = 0.0f;

//                 for (int i = 0; i < STATE_COUNT; ++i)
//                 {
//                     for (int j = 0; j < node.numChildren; ++j)
//                     {
//                         if (childResults[j][i][s] == 0)
//                         {
//                             continue;
//                         }
//                         results[i][s] *= childResults[j][i][s];
//                     }
//                     sumLikelihoods += results[i][s];
//                 }

//                 node.scaleVector[seqPos] = sumLikelihoods; // Store scaling factor
                
//                 for (int i = 0; i < STATE_COUNT; ++i)
//                 {
//                     node.bottom[i][seqPos] = results[i][s] / sumLikelihoods; // Normalize result
//                 }
//             }
// //

            for (int s = 0; s < siteNum; ++s)
            {
                int seqPos = seqPosStart + s;
                if (seqPos >= sequenceLength)
                    break;

                int exponents[STATE_COUNT];
                int minExp = INT_MAX;
                int maxExp = INT_MIN;

                for (int i = 0; i < STATE_COUNT; ++i)
                {

                    for (int j = 0; j < node.numChildren; ++j)
                    {
                        if (childResults[j][i][s] == 0.0f)
                        {
                            continue;
                        }
                        results[i][s] *= childResults[j][i][s];
                    }
                    frexpf(results[i][s], &exponents[i]);
                    maxExp = max(exponents[i], maxExp);
                }
                scaling_f scaleFactor = pow(2,maxExp); // 2^{maxExp}
                // scaling_f scaleFactor = ldexpf(1.0f, -maxExp); 
                // node.scaleVector[seqPos] = scaleFactor;
                node.scaleVector[seqPos] = scaleFactor;
                if (threadIdx.x == 0 && blockIdx.x == 0)
                {
                    printf("seqPos = %d\n", seqPos);
                    printf("maxExp = %d, scaleFactor = %d\n", maxExp, scaleFactor);
                    for (int i = 0; i < STATE_COUNT; ++i)
                    {
                        printf("Before scaling - results[%d][%d] = %e, exponent = %d\n", i, s, results[i][s], exponents[i]);
                    }
                }

                // 对结果进行缩放并存储到 bottom 中
                for (int i = 0; i < STATE_COUNT; ++i)
                {
                    results[i][s] /= scaleFactor;
                    
                    node.bottom[i][seqPos] = results[i][s];

                    // 添加调试输出（仅在第一个线程中）
                    if (threadIdx.x == 0 && blockIdx.x == 0)
                    {
                        printf("After scaling - results[%d][%d] = %e\n", i, s, results[i][s]);
                    }
                }
            }



        }
    }
}

exp_f *d_rateMatrix;
exp_f *d_matrixV;
exp_f *d_eigenvalues;

// CUDA kernel to scale and exponentiate eigenvalues
__global__ void scale_and_exp_eigenvalues(exp_f *d_eigenvalues, exp_f *d_matrixExpD, exp_f scalar, int n)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n)
    {
        exp_f exp_val = exp(scalar * d_eigenvalues[idx]);
        d_matrixExpD[idx * n + idx] = exp_val;
    }
}

__host__ void precompute_eigen_decomposition(int n, cusolverDnHandle_t solverHandle)
{

    exp_f *d_matrixA;

    CHECK_CUDA(cudaMalloc((void **)&d_matrixA, n * n * sizeof(exp_f)));

    // Copy the rate matrix to d_matrixA
    CHECK_CUDA(cudaMemcpy(d_matrixA, d_rateMatrix, n * n * sizeof(exp_f), cudaMemcpyDeviceToDevice));

    // Compute eigenvalues and eigenvectors
    int lwork = 0, *devInfo;

    CHECK_CUDA(cudaMalloc((void **)&devInfo, sizeof(int)));
    CHECK_CUSOLVER(cusolverDnSsyevd_bufferSize(solverHandle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, n, d_matrixA, n, d_eigenvalues, &lwork));

    exp_f *d_work;
    CHECK_CUDA(cudaMalloc((void **)&d_work, lwork * sizeof(exp_f)));
    CHECK_CUSOLVER(cusolverDnSsyevd(solverHandle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, n, d_matrixA, n, d_eigenvalues, d_work, lwork, devInfo));
    
    CHECK_CUDA(cudaDeviceSynchronize());

    // Copy eigenvectors to d_matrixV
    CHECK_CUDA(cudaMemcpy(d_matrixV, d_matrixA, n * n * sizeof(exp_f), cudaMemcpyDeviceToDevice));

    // Free resources
    CHECK_CUDA(cudaFree(d_matrixA));
    CHECK_CUDA(cudaFree(d_work));
    CHECK_CUDA(cudaFree(devInfo));
}

__global__ void computeMatrixExponential(exp_f *d_eigenvalues, GPUNode *d_nodes, exp_f *d_matrixV)
{
    int matrixIndex = blockIdx.x;
    int threadIndex = threadIdx.x;
    int i = threadIndex / STATE_COUNT;
    int j = threadIndex % STATE_COUNT;

    if (matrixIndex < NODE_COUNT)
    {
        __shared__ exp_f matrixV[STATE_COUNT][STATE_COUNT];
        __shared__ exp_f matrixExpD[STATE_COUNT][STATE_COUNT];
        __shared__ exp_f matrixVExpDVt[STATE_COUNT][STATE_COUNT];

        if (threadIndex < STATE_COUNT * STATE_COUNT)
        {
            matrixExpD[i][j] = 0.0;
            matrixV[i][j] = 0.0;
            matrixVExpDVt[i][j] = 0.0;
            // finalResult[i][j] = 0.0;
        }
        __syncthreads();

        // Load V matrix and eigenvalues into shared memory
        if (threadIndex < STATE_COUNT * STATE_COUNT)
        {
            matrixV[i][j] = d_matrixV[i * STATE_COUNT + j];
        }
        __syncthreads();

        if (threadIndex < STATE_COUNT)
        {
            matrixExpD[threadIndex][threadIndex] = exp(d_nodes[matrixIndex + 1].branchLength * d_eigenvalues[threadIndex]);
        }
        __syncthreads();

        // Perform matrix multiplication: V * exp(sD)
        exp_f sum = 0.0;
        for (int k = 0; k < STATE_COUNT; ++k)
        {
            sum += matrixExpD[i][k] * matrixV[k][j];
            // sum += matrixV[i][k] * (k == j ? matrixExpD[k] : 0.0);
        }

        if (threadIndex < STATE_COUNT * STATE_COUNT)
        {
            matrixVExpDVt[i][j] = sum;
        }

        __syncthreads();

        // Perform matrix multiplication: (V * exp(sD)) * V^T
        sum = 0.0;
        for (int k = 0; k < STATE_COUNT; ++k)
        {
            sum += matrixV[k][i] * matrixVExpDVt[k][j]; // Transpose multiplication
        }
        if (threadIndex < STATE_COUNT * STATE_COUNT)
        {
            d_nodes[matrixIndex + 1].probabilityMatrix[i][j] = sum;
        }

        __syncthreads();

        // // Write the result back to global memory
        // if (threadIndex < STATE_COUNT * STATE_COUNT) {
        //     d_nodes[matrixIndex + 1].probabilityMatrix[i][j] = finalResult[i][j];
        // }
        // __syncthreads();

        // if (threadIndex == 0 && matrixIndex == 0) {
        //     printf("matrixV for block %d:\n", matrixIndex);
        //     for (int i = 0; i < STATE_COUNT; ++i) {
        //         for (int j = 0; j < STATE_COUNT; ++j) {
        //             printf("%f ", matrixV[i][j]);
        //         }
        //         printf("\n");
        //     }
        //     printf("matrixExpD result for block %d:\n", matrixIndex);
        //     for (int i = 0; i < STATE_COUNT; ++i) {
        //         for (int j = 0; j < STATE_COUNT; ++j) {
        //             printf("%f ", matrixExpD[i][j]);
        //         }
        //         printf("\n");
        //     }
        //     printf("matrixVExpDVt result for block %d:\n", matrixIndex);
        //     for (int i = 0; i < STATE_COUNT; ++i) {
        //         for (int j = 0; j < STATE_COUNT; ++j) {
        //             printf("%f ", matrixVExpDVt[i][j]);
        //         }
        //         printf("\n");
        //     }
        //     printf("Transpose V result for block %d:\n", matrixIndex);
        //     for (int i = 0; i < STATE_COUNT; ++i) {
        //         for (int j = 0; j < STATE_COUNT; ++j) {
        //             printf("%f ", matrixVT[i][j]);
        //         }
        //         printf("\n");
        //     }
        //     printf("FINAL result for block %d:\n", matrixIndex);
        //     for (int i = 0; i < STATE_COUNT; ++i) {
        //         for (int j = 0; j < STATE_COUNT; ++j) {
        //             printf("%f ", finalResult[i][j]);
        //         }
        //         printf("\n");
        //     }

        // }
        // __syncthreads();
    }
}

__host__ void compute_exponential(GPUNode *d_nodes)
{

    int threads_per_block = STATE_COUNT * STATE_COUNT;
    int blocks_per_grid = NODE_COUNT;

    computeMatrixExponential<<<blocks_per_grid, threads_per_block>>>(d_eigenvalues, d_nodes, d_matrixV);
    CHECK_CUDA(cudaPeekAtLastError());
    CHECK_CUDA(cudaDeviceSynchronize());
}


// Modify the function to return the depth of the subtree rooted at nodeIndex
int postorderTraversal(const GPUTree &tree, int nodeIndex, std::vector<int> &order)
{
    const GPUNode &node = tree.nodes[nodeIndex];

    int maxChildDepth = 0; // To keep track of the maximum depth among child subtrees

    // Recursively traverse child nodes
    for (int i = 0; i < node.numChildren; ++i)
    {
        int childIndex = node.childrenIndices[i];
        int childDepth = postorderTraversal(tree, childIndex, order);
        maxChildDepth = std::max(maxChildDepth, childDepth);
    }

    order.push_back(nodeIndex); // Add current node to post-order traversal

    return maxChildDepth + 1; // Depth of current node is max child depth + 1
}

void computePostorder(GPUTree &tree, int rootIndex)
{
    std::vector<int> order;
    int maxTreeDepth = postorderTraversal(tree, rootIndex, order);

    // Copy the traversal order into the tree's order array
    for (size_t i = 0; i < order.size(); ++i)
    {
        tree.order[i] = order[i];
    }

    // Print out the maximum tree depth
    std::cout << "Maximum tree depth: " << maxTreeDepth << std::endl;
}

// void postorderTraversal(const GPUTree &tree, int nodeIndex, std::vector<int> &order)
// {
//     GPUNode node = tree.nodes[nodeIndex];

//     for (int i = 0; i < node.numChildren; ++i)
//     {
//         int childIndex = node.childrenIndices[i];
//         postorderTraversal(tree, childIndex, order);
//     }

//     order.push_back(nodeIndex);
// }

// void computePostorder(GPUTree &tree, int rootIndex)
// {
//     std::vector<int> order;
//     postorderTraversal(tree, rootIndex, order);

//     // add to order
//     for (int i = 0; i < order.size(); ++i)
//     {
//         tree.order[i] = order[i];
//     }
// }



#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

void matrix_exp(float bl, const std::vector<std::vector<float>>& rate_matrix, std::vector<std::vector<float>>& mat_out)
{
    Eigen::MatrixXf A(rate_matrix.size(), rate_matrix.size());
    for (size_t i = 0; i < rate_matrix.size(); i++)
    {
        for (size_t j = 0; j < rate_matrix[i].size(); j++)
        {
            A(i, j) = rate_matrix[i][j];
        }
    }

    Eigen::MatrixXf expA = (A * bl).exp();

    mat_out.resize(rate_matrix.size());
    for (size_t i = 0; i < rate_matrix.size(); i++)
    {
        mat_out[i].resize(rate_matrix.size());
        for (size_t j = 0; j < rate_matrix.size(); j++)
        {
            mat_out[i][j] = expA(i, j);
        }
    }
}


void computeProbabilityMatrices(GPUTree &h_gpuTree)
{
    std::vector<std::vector<float>> rate_matrix(STATE_COUNT, std::vector<float>(STATE_COUNT));
    for (int i = 0; i < STATE_COUNT; i++)
    {
        for (int j = 0; j < STATE_COUNT; j++)
        {
            rate_matrix[i][j] = h_gpuTree.rateMatrix[i * STATE_COUNT + j];
        }
    }

    for (int i = 0; i < NODE_COUNT; i++)
    {
        GPUNode &node = h_gpuTree.nodes[i];
        for (int i = 0; i < SEQUENCE_LENGTH; i++)
        {
            node.scaleVector[i] = 1.0;
        }

        std::vector<std::vector<float>> prob_matrix;

        matrix_exp(node.branchLength, rate_matrix, prob_matrix);

        for (int m = 0; m < STATE_COUNT; m++)
        {
            for (int n = 0; n < STATE_COUNT; n++)
            {
                node.probabilityMatrix[m][n] = prob_matrix[m][n];
            }
        }
    }
}
void top_down(GPUTree &h_gpuTree)
{
    computeProbabilityMatrices(h_gpuTree);

    GPUNode *d_nodes;
    CHECK_CUDA(cudaMalloc((void **)&d_rateMatrix, STATE_COUNT * STATE_COUNT * sizeof(exp_f)));
    CHECK_CUDA(cudaMemcpy(d_rateMatrix, h_gpuTree.rateMatrix, STATE_COUNT * STATE_COUNT * sizeof(exp_f), cudaMemcpyHostToDevice));

    CHECK_CUDA(cudaMalloc(&d_nodes, sizeof(GPUNode) * NODE_COUNT));
    CHECK_CUDA(cudaMemcpy(d_nodes, h_gpuTree.nodes, sizeof(GPUNode) * NODE_COUNT, cudaMemcpyHostToDevice));

    int *order;
    CHECK_CUDA(cudaMalloc((void **)&order, sizeof(int) * NODE_COUNT));
    CHECK_CUDA(cudaMemcpy(order, h_gpuTree.order, sizeof(int) * NODE_COUNT, cudaMemcpyHostToDevice));

    CHECK_CUDA(cudaMalloc((void **)&d_matrixV, STATE_COUNT * STATE_COUNT * sizeof(exp_f)));
    CHECK_CUDA(cudaMalloc((void **)&d_eigenvalues, STATE_COUNT * sizeof(exp_f)));

    // Initialize cuSOLVER and cuBLAS handles
    cusolverDnHandle_t solverHandle;
    cublasHandle_t blasHandle;
    CHECK_CUSOLVER(cusolverDnCreate(&solverHandle));
    CHECK_CUBLAS(cublasCreate(&blasHandle));

    precompute_eigen_decomposition(STATE_COUNT, solverHandle);

    auto start = std::chrono::high_resolution_clock::now();

    // compute_exponential(d_nodes);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Matrix exponential execution took " << duration.count() << " microseconds." << std::endl;

    CHECK_CUSOLVER(cusolverDnDestroy(solverHandle));
    CHECK_CUBLAS(cublasDestroy(blasHandle));

    CHECK_CUDA(cudaFree(d_rateMatrix));
    CHECK_CUDA(cudaFree(d_matrixV));
    CHECK_CUDA(cudaFree(d_eigenvalues));

    start = std::chrono::high_resolution_clock::now();

    int threadsPerBlock = 32;                                                                           
    int blocksPerGrid = (SEQUENCE_LENGTH + threadsPerBlock * siteNum - 1) / (threadsPerBlock * siteNum);
    // int blocksPerGrid = (SEQUENCE_LENGTH + threadsPerBlock - 1) / threadsPerBlock;

    std::cout << "Threads per block: " << threadsPerBlock << std::endl;
    std::cout << "Blocks per grid: " << blocksPerGrid << std::endl;

    // less optimized
    // non_recursive_dfs<<<blocksPerGrid, threadsPerBlock>>>(d_nodes, SEQUENCE_LENGTH);

    // optimized
    felsenstein_pruning_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_nodes, order, SEQUENCE_LENGTH, NODE_COUNT);

    cudaError_t error = cudaPeekAtLastError();
    if (error != cudaSuccess)
    {
        printf("CUDA error: %s\n", cudaGetErrorString(error));
    }
    cudaDeviceSynchronize();

    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Tree Traversal GPU execution took " << duration.count() << " microseconds." << std::endl;

    CHECK_CUDA(cudaMemcpy(h_gpuTree.nodes, d_nodes, sizeof(GPUNode) * NODE_COUNT, cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaFree(d_nodes));
    CHECK_CUDA(cudaFree(order));
}

double postorder_traversal_scaling(size_t node_index, GPUTree &h_gpuTree)
{
    double lk = 0.0;
    for (int i = 0; i < h_gpuTree.nodes[node_index].numChildren; i++)
    {
        int child_index = h_gpuTree.nodes[node_index].childrenIndices[i];
        lk += postorder_traversal_scaling(child_index, h_gpuTree);
    }

    for (int i = 0; i < SEQUENCE_LENGTH; i++)
    {
        lk += log(h_gpuTree.nodes[node_index].scaleVector[i]);
    }

    return lk;
}

std::vector<float> pi;

void felsenstein_pruning(GPUTree &h_gpuTree)
{

    std::cout << "Starting tree traversal on GPU" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    top_down(h_gpuTree);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "GPU execution took " << duration.count() << " microseconds." << std::endl;

    // // display the first leaf node bottom
    // for (int i = 0; i < NODE_COUNT; i ++) {
    //     if (h_gpuTree.nodes[i].numChildren != 0) {
    //         std::cout << "Leaf Node " << i << std::endl;
    //         for (int j = 0; j < SEQUENCE_LENGTH; j++) {
    //             for (int k = 0; k < STATE_COUNT; k++) {
    //                 // check if it is nan
    //                 if (std::isnan(h_gpuTree.nodes[i].bottom[k][j])) {
    //                     std::cout << "NAN on site " << j << " state " << k << std::endl;
    //                     break;
    //                 }
    //             }
    //         }
    //     }
    // }

    for (int i = 0; i < STATE_COUNT; ++i)
    {
        std::cout << pi[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "Starting scaling traversal" << std::endl;

    // for (size_t i = 0; i < 20; i++)
    // {

    //     for (size_t j = 0; j < STATE_COUNT; j++)
    //     {
    //         printf("%.40f ", h_gpuTree.nodes[0].bottom[j][i]);
    //         // std::cout << i << ":" << h_gpuTree.nodes[0].bottom[j][i] << "\t";
    //     }
    //     std::cout << std::endl;
    // }

    std::vector<std::vector<float>> bottom(STATE_COUNT, std::vector<float>(SEQUENCE_LENGTH));

    for (size_t i = 0; i < SEQUENCE_LENGTH; i++)
    {
        for (size_t j = 0; j < STATE_COUNT; j++)
        {
            bottom[j][i] = h_gpuTree.nodes[0].bottom[j][i];
        }
    }

    double lk = 0;

    double root_lk = 0.0;

    double sum = 0.0;
    for (size_t i = 0; i < SEQUENCE_LENGTH; i++)
    {
        float lk_node = 0;
        for (size_t j = 0; j < STATE_COUNT; j++)
        {
            // std::cout << bottom[i][j] << "\t" << pi[j] << "\t";
            lk_node += bottom[j][i] * pi[j];
            sum += bottom[j][i];
        }

        lk += (double)log(lk_node);
    }
    printf("sum: %.50f\n", sum);

    printf("likelihood of Root llkROot: %.20f\n", root_lk);
    printf("likelihood of Root GPU: %lf\n", lk);

    // Adding scaling factor
    lk += postorder_traversal_scaling(0, h_gpuTree);

    printf("likelihood of GPU: %lf\n", lk);
}

int main(int argc, char **argv)
{
    if (argc != 3)
        std::cerr << "Usage: " << argv[0] << " [alignment file] [newick tree file]" << std::endl;

    utility::msa_seq(argv[1]); // Store MSA data into data-structure

    utility::subs_param.resize(10); // Set Subs parameter

    for (size_t i = 0; i < 10; i++)
        utility::subs_param[i] = 1;
    utility::rate_matrix_calc(); // Find the rate matrix

    if (1) // Print matrix_exp Matrix
    {
        for (size_t i = 0; i < utility::rate_matrix.size(); i++)
        {
            for (size_t j = 0; j < utility::rate_matrix[0].size(); j++) std::cout << utility::rate_matrix[i][j] << "\t";
            std::cout << "\n";
        }
        for (auto &p: utility::pi) std::cout << p << "\n";
    }

    std::ifstream newickTree(argv[2]);
    std::string newickString;
    newickTree >> newickString;
    utility::Tree tree(newickString);
    std::cerr << "after tree" << std::endl;
    fflush(stderr);

    int node_cnt = tree.allNodes.size();
    int seq_len = utility::seqs.begin()->second.second.length();
    int max_child = 0;

    for (auto &n : tree.allNodes)
        max_child = std::max(max_child, (int)n.second->children.size());

    std::cerr << "NODE_COUNT: " << node_cnt << std::endl;
    std::cerr << "SEQUENCE_LENGTH: " << seq_len << std::endl;
    std::cerr << "MAX_CHILD: " << max_child << std::endl;

    pi = utility::pi;

    GPUTree *tree_gpu1 = new GPUTree();

    // initializeTree(tree_gpu);
    copyTree(tree, *tree_gpu1);
    computePostorder(*tree_gpu1, 0);

    // for (int i = 0; i < NODE_COUNT; ++i) {
    //     std::cout << "Node " << tree_gpu1->order[i] << " is at position " << i << " in postorder traversal.\n";
    // }

    auto start = std::chrono::high_resolution_clock::now();

    felsenstein_pruning(*tree_gpu1);

    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "Felsensteins pruning GPU execution took " << duration.count() << " microseconds." << std::endl;

    // Display results
    // // print root node
    // printf("Root Node\n");
    // printf("branchLength: %f\n", tree_gpu->nodes[0].branchLength);

    // printf("Bottom Matrix:\n");
    // for (int i = 0; i < 10; ++i) {
    //     for (int j = 0; j < STATE_COUNT; ++j) {
    //         printf("%f ", tree_gpu->nodes[0].bottom[j][i]);
    //     }
    //     printf("\n");
    // }

    // // print rate matrix
    // for (int i = 0; i < STATE_COUNT; ++i) {
    //     for (int j = 0; j < STATE_COUNT; ++j) {
    //         printf("%f ", tree_gpu1->rateMatrix[i*STATE_COUNT+j]);
    //     }
    //     printf("\n");
    // }

    // for (int i = 0; i < 2; ++i) {
    //     printf("Node %d\n", i);
    //     printf("branchLength: %f\n", tree_gpu1->nodes[i].branchLength);
    //     printf("Probability Matrix:\n");
    //     for (int j = 0; j < STATE_COUNT; ++j) {
    //         for (int k = 0; k < STATE_COUNT; ++k) {
    //             printf("%f ", tree_gpu1->nodes[i].probabilityMatrix[j][k]);
    //         }
    //         printf("\n");
    //     }
    // }
    //     // printf("Bottom Matrix:\n");
    //     // for (int j = 0; j < STATE_COUNT; ++j) {
    //     //     for (int k = 0; k < SEQUENCE_LENGTH; ++k) {
    //     //         printf("%f ", *tree_gpu.nodes[i].bottom[j][k]);
    //     //     }
    //     //     printf("\n");
    //     // }
    //     // printf("\n");
    // }
    // utility::bottom_up(tree);
    // std::cout << "\n";
    // utility::top_down(tree);
    // std::cout << "\n";
    // utility::marginal(tree);
    
    delete tree_gpu1;

    // utility::fitch(tree);

    // if (PRINT_INFERENCE)
    // {
    //     std::string tip="Mon";
    //     auto search = tree.allNodes.find(tip);
    //     if (search != tree.allNodes.end()) utility::printLikelihoodInference(tree, search->second);
    //     else {std::cerr<<"Tip not found";}
    // }

    return 0;
}

// __host__ void symmetric_matrix_exponential(exp_f scalar, int n) {

//     cusolverDnHandle_t solverHandle;
//     cublasHandle_t blasHandle;
//     CHECK_CUSOLVER(cusolverDnCreate(&solverHandle));
//     CHECK_CUBLAS(cublasCreate(&blasHandle));

//     exp_f* d_matrixA, *d_matrixV, *d_eigenvalues, *d_matrixExpD, *d_matrixVExpDVt;
//     CHECK_CUDA(cudaMalloc((void**)&d_matrixA, n * n * sizeof(exp_f)));
//     CHECK_CUDA(cudaMalloc((void**)&d_matrixV, n * n * sizeof(exp_f)));
//     CHECK_CUDA(cudaMalloc((void**)&d_eigenvalues, n * sizeof(exp_f)));
//     CHECK_CUDA(cudaMalloc((void**)&d_matrixExpD, n * n * sizeof(exp_f)));
//     CHECK_CUDA(cudaMalloc((void**)&d_matrixVExpDVt, n * n * sizeof(exp_f)));

//     // Scale the matrix by the scalar on device
//     int threads_per_block = 256;
//     int blocks_per_grid = (n * n + threads_per_block - 1) / threads_per_block;
//     scale_matrix<<<blocks_per_grid, threads_per_block>>>(d_rateMatrix, d_matrixA, scalar, n);
//     CHECK_CUDA(cudaPeekAtLastError());
//     CHECK_CUDA(cudaDeviceSynchronize());

//     // Compute eigenvalues and eigenvectors
//     int lwork = 0, *devInfo;
//     CHECK_CUDA(cudaMalloc((void**)&devInfo, sizeof(int)));
//     CHECK_CUSOLVER(cusolverDnSsyevd_bufferSize(solverHandle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, n, d_matrixA, n, d_eigenvalues, &lwork));
//     exp_f* d_work;
//     CHECK_CUDA(cudaMalloc((void**)&d_work, lwork * sizeof(exp_f)));
//     CHECK_CUSOLVER(cusolverDnSsyevd(solverHandle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, n, d_matrixA, n, d_eigenvalues, d_work, lwork, devInfo));
//     CHECK_CUDA(cudaDeviceSynchronize());

//     auto start = std::chrono::high_resolution_clock::now();

//     // Copy eigenvectors to d_matrixV
//     CHECK_CUDA(cudaMemcpy(d_matrixV, d_matrixA, n * n * sizeof(exp_f), cudaMemcpyDeviceToDevice));

//     // Exponentiate eigenvalues on device using CUDA kernel

//     exp_diagonal_elements<<<blocks_per_grid, threads_per_block>>>(d_eigenvalues, n);
//     CHECK_CUDA(cudaPeekAtLastError());
//     CHECK_CUDA(cudaDeviceSynchronize());

//     // Construct exp(D)
//     CHECK_CUDA(cudaMemset(d_matrixExpD, 0, n * n * sizeof(exp_f)));
//     for (int i = 0; i < n; i++) {
//         CHECK_CUDA(cudaMemcpy(&d_matrixExpD[i * n + i], &d_eigenvalues[i], sizeof(exp_f), cudaMemcpyDeviceToDevice));
//     }

//     // Calculate V * exp(D)
//     const exp_f alpha = 1.0;
//     const exp_f beta = 0.0;
//     CHECK_CUBLAS(cublasSgemm(blasHandle, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, &alpha, d_matrixV, n, d_matrixExpD, n, &beta, d_matrixVExpDVt, n));

//     // Calculate (V * exp(D)) * V^T
//     CHECK_CUBLAS(cublasSgemm(blasHandle, CUBLAS_OP_N, CUBLAS_OP_T, n, n, n, &alpha, d_matrixVExpDVt, n, d_matrixV, n, &beta, d_matrixA, n));

//     // Copy the result back to host
//     CHECK_CUDA(cudaMemcpy(d_tempProbMatrix, d_matrixA, n * n * sizeof(exp_f), cudaMemcpyDeviceToDevice));

//     auto end = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//     std::cout << "scaling took " << duration.count() << " microseconds." << std::endl;

//     // Free resources
//     CHECK_CUDA(cudaFree(d_matrixA));
//     CHECK_CUDA(cudaFree(d_matrixV));
//     CHECK_CUDA(cudaFree(d_eigenvalues));
//     CHECK_CUDA(cudaFree(d_matrixExpD));
//     CHECK_CUDA(cudaFree(d_matrixVExpDVt));
//     CHECK_CUDA(cudaFree(d_work));
//     CHECK_CUDA(cudaFree(devInfo));
//     CHECK_CUSOLVER(cusolverDnDestroy(solverHandle));
//     CHECK_CUBLAS(cublasDestroy(blasHandle));

// }
