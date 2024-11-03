// clustering.cu

// #include <cuda_runtime.h>
// #include <curand.h>
// #include <curand_kernel.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <ctime>
#include <stdexcept>
#include <unordered_set>

#define THREADS_PER_BLOCK 256

// Function to check for CUDA errors
// #define cudaCheckError(ans) { gpuAssert((ans), __FILE__, __LINE__); }
// inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true){
//     if (code != cudaSuccess){
//         fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
//         if (abort) exit(code);
//     }
// }

// // Distance function (Euclidean distance in this case)
// __device__ float distance(float x1, float x2){
//     return fabsf(x1 - x2);
// }

// // Kernel to compute distances and assign clusters
// __global__ void assignClusters(float *E, int size, float a, float b, int *clusterAssignments){
//     int idx = blockIdx.x * blockDim.x + threadIdx.x;
//     if(idx >= size) return;

//     float distToA = distance(E[idx], a);
//     float distToB = distance(E[idx], b);

//     if(distToA < distToB)
//         clusterAssignments[idx] = 0; // Assign to cluster Ea
//     else
//         clusterAssignments[idx] = 1; // Assign to cluster Eb
// }

// // Host function to perform clustering recursively
// void cluster(float *E, int size, int k, std::vector<std::vector<float>> &clusters){
//     if(size <= k){
//         // Base case: add the cluster to the list
//         clusters.push_back(std::vector<float>(E, E + size));
//         return;
//     }

//     // Randomly select two pivot points a and b
//     srand(time(NULL));
//     int idxA = rand() % size;
//     int idxB = rand() % size;
//     while(idxB == idxA){
//         idxB = rand() % size;
//     }
//     float a = E[idxA];
//     float b = E[idxB];

//     // Allocate memory on device
//     float *d_E;
//     int *d_clusterAssignments;
//     cudaCheckError(cudaMalloc((void**)&d_E, size * sizeof(float)));
//     cudaCheckError(cudaMalloc((void**)&d_clusterAssignments, size * sizeof(int)));

//     // Copy data to device
//     cudaCheckError(cudaMemcpy(d_E, E, size * sizeof(float), cudaMemcpyHostToDevice));

//     // Launch kernel to assign clusters
//     int blocks = (size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
//     assignClusters<<<blocks, THREADS_PER_BLOCK>>>(d_E, size, a, b, d_clusterAssignments);
//     cudaDeviceSynchronize();

//     // Copy cluster assignments back to host
//     int *clusterAssignments = new int[size];
//     cudaCheckError(cudaMemcpy(clusterAssignments, d_clusterAssignments, size * sizeof(int), cudaMemcpyDeviceToHost));

//     // Separate elements into two clusters
//     std::vector<float> Ea, Eb;
//     for(int i = 0; i < size; i++){
//         if(i == idxA || i == idxB) continue; // Exclude a and b
//         if(clusterAssignments[i] == 0)
//             Ea.push_back(E[i]);
//         else
//             Eb.push_back(E[i]);
//     }

//     // Free device memory
//     cudaCheckError(cudaFree(d_E));
//     cudaCheckError(cudaFree(d_clusterAssignments));
//     delete[] clusterAssignments;

//     // Recursively cluster Ea and Eb
//     if(!Ea.empty()){
//         cluster(Ea.data(), Ea.size(), k, clusters);
//     }
//     if(!Eb.empty()){
//         cluster(Eb.data(), Eb.size(), k, clusters);
//     }
// }
std::vector<int> distCal(const std::vector<int> &cluster, int index)
{
    std::vector<int> distCalDiffMatrix;
    int ref = cluster[index]; // Store the reference value from cluster

    for (size_t iter = 0; iter < cluster.size(); ++iter)
    {
        distCalDiffMatrix.push_back(std::abs(ref - cluster[iter])); // Calculate the difference
    }

    return distCalDiffMatrix;
}

std::vector<int> clusterAlgo(const std::vector<int> &distMatrix1, const std::vector<int> &distMatrix2, const std::vector<int> &cluster, int *&boundary)
{
    std::vector<int> cluster1, cluster2;

    // Ensure both matrices have the same size
    if (distMatrix1.size() != distMatrix2.size())
    {
        throw std::invalid_argument("Distance matrices must have the same size");
    }

    for (size_t i = 0; i < distMatrix1.size(); ++i)
    {
        if (std::abs(distMatrix1[i]) <= std::abs(distMatrix2[i]))
        {
            cluster1.push_back(cluster[i]);
        }
        else
        {
            cluster2.push_back(cluster[i]);
        }
    }
    *boundary = cluster1.size(); // Combine clusters into a single result vector
    cluster1.insert(cluster1.end(), cluster2.begin(), cluster2.end());

    return cluster1;
}

struct treeNode
{
    int nodeNum;
    int child1;
    int child2;
};

// void cluster(int *cInstr, int numCluster, int totalSize, int *clusterMap, int *dataset)
// {
//     for (int idx = 0; idx < totalSize; idx++)
//     {
//         if (clusterMap[idx] >= 0)
//         {
//             for (int clusterIdx = 0; clusterIdx < 3 * (numCluster); clusterIdx += 3)
//                 if (cInstr[clusterIdx] == clusterMap[idx])
//                 {
//                     int distance1 = std::abs(cInstr[clusterIdx + 1] - dataset[idx]);
//                     int distance2 = std::abs(cInstr[clusterIdx + 2] - dataset[idx]);
//                     if (distance1 < distance2)
//                         clusterMap[idx] = cInstr[clusterIdx] * 2 + 1;
//                     else
//                         clusterMap[idx] = cInstr[clusterIdx] * 2 + 2;
//                 }
//         }
//     }
// }

void getTwoRandomIndices(int *clusterMap, int clusterSize, int searchIndex, treeNode *node)
{
    static int result[2] = {-1, -1}; // Static array to store the results
    std::unordered_set<int> uniqueIndices;
    const int MAX_ITERATIONS = 10; // Maximum number of iterations to prevent infinite loops
    int iterations = 0;

    // Seed the random number generator
    std::srand(std::time(nullptr));

    while (uniqueIndices.size() < 2 && iterations < MAX_ITERATIONS)
    {
        iterations++;
        // Randomly get a number within the size of the array
        int randomIndex = std::rand() % clusterSize;

        // Start iterating from the random index
        for (int i = 0; i < clusterSize; i++)
        {
            int indx = (randomIndex + i) % clusterSize; // Wrap around if needed
            if (clusterMap[indx] == searchIndex)
            {
                uniqueIndices.insert(indx);
                break; // Break the inner loop and start over for the next index
            }
        }
    }

    // Error handling for insufficient unique indices
    if (uniqueIndices.size() == 2)
    {
        auto it = uniqueIndices.begin();
        node->child1 = *it++;
        node->child2 = *it;
        node->nodeNum = searchIndex;
    }
}
void invalidateExtraOccurrences(int *arr, int size)
{
    // Find the maximum number in the array
    int maxNum = *std::max_element(arr, arr + size);

    // Allocate memory for the counts array
    int resultSize = maxNum + 1;
    int *counts = new int[resultSize](); // Initialize with zeros

    // Count occurrences in a single pass
    for (int i = 0; i < size; ++i)
    {
        counts[arr[i]]++;
    }

    // Invalidate numbers with less than 4 occurrences in a single pass
    for (int i = 0; i < size; ++i)
    {
        if (counts[arr[i]] < 4)
        {
            arr[i] = -arr[i];
        }
    }

    // Free dynamically allocated memory
    delete[] counts;
}

// int main()
// {
//     int *boundary = new int;
//     srand(time(NULL));

// #define clusterSize 21

//     int clusters[clusterSize] = {1, 2, 3, 98, 5, 6, 67, 89, 91, 100, 102, 10, 32, 45, 6, 7, 432, 54, 23, 66, 43};
//     int clusterMap[clusterSize] = {0};
//     treeNode *nodes[10];
//     int n = -1;
//     for (int i = 0; i < 10; i++)
//     {
//         nodes[i] = new treeNode();
//     }
//     getTwoRandomIndices(clusterMap, clusterSize, 0, nodes[0]);
//     int cInstr[3] = {0, nodes[0]->child1, nodes[0]->child2};
//     cluster(cInstr, 1, clusterSize, clusterMap, clusters);
//     invalidateExtraOccurrences(clusterMap, clusterSize);
//     while (n++ < clusterSize - 1)
//     {
//         printf("index %d value %d \n", n, clusterMap[n]);
//     }

//     getTwoRandomIndices(clusterMap, clusterSize, nodes[0]->nodeNum * 2 + 1, nodes[1]);
//     getTwoRandomIndices(clusterMap, clusterSize, nodes[0]->nodeNum * 2 + 2, nodes[2]);
//     int cInstr2[6] = {nodes[0]->nodeNum * 2 + 1, nodes[1]->child1, nodes[1]->child2,
//                       nodes[0]->nodeNum * 2 + 2, nodes[2]->child1, nodes[2]->child2};
//     cluster(cInstr2, 2, clusterSize, clusterMap, clusters);
//     invalidateExtraOccurrences(clusterMap, clusterSize);
//     n = -1;
//     while (n++ < clusterSize - 1)
//     {
//         printf("index %d value %d \n", n, clusterMap[n]);
//     }

//     getTwoRandomIndices(clusterMap, clusterSize, nodes[1]->nodeNum * 2 + 1, nodes[3]);
//     getTwoRandomIndices(clusterMap, clusterSize, nodes[1]->nodeNum * 2 + 2, nodes[4]);
//     getTwoRandomIndices(clusterMap, clusterSize, nodes[2]->nodeNum * 2 + 1, nodes[5]);
//     getTwoRandomIndices(clusterMap, clusterSize, nodes[2]->nodeNum * 2 + 2, nodes[6]);
//     int cInstr3[12] = {nodes[1]->nodeNum * 2 + 1, nodes[3]->child1, nodes[3]->child2,
//                        nodes[1]->nodeNum * 2 + 2, nodes[4]->child1, nodes[4]->child2,
//                        nodes[2]->nodeNum * 2 + 2, nodes[5]->child1, nodes[5]->child2,
//                        nodes[2]->nodeNum * 2 + 2, nodes[6]->child1, nodes[6]->child2};
//     cluster(cInstr3, 4, clusterSize, clusterMap, clusters);
//     invalidateExtraOccurrences(clusterMap, clusterSize);

//     n = -1;
//     while (n++ < clusterSize - 1)
//     {
//         printf("index %d value %d \n", n, clusterMap[n]);
//     }
//     return 0;
// }

// CUDA kernel
__global__ void clusterKernel(int *cInstr, int numCluster, int totalSize, int *clusterMap, int *dataset)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < totalSize)
    {
        if (clusterMap[idx] >= 0)
        {
            for (int clusterIdx = 0; clusterIdx < 3 * (numCluster); clusterIdx += 3)
            {
                if (cInstr[clusterIdx] == clusterMap[idx])
                {
                    int distance1 = abs(cInstr[clusterIdx + 1] - dataset[idx]);
                    int distance2 = abs(cInstr[clusterIdx + 2] - dataset[idx]);
                    if (distance1 < distance2)
                        clusterMap[idx] = cInstr[clusterIdx] * 2 + 1;
                    else
                        clusterMap[idx] = cInstr[clusterIdx] * 2 + 2;
                }
            }
        }
    }
}

// Host function to call the CUDA kernel
void cluster(int *cInstr, int numCluster, int totalSize, int *clusterMap, int *dataset)
{
    int *d_cInstr, *d_clusterMap, *d_dataset;
    
    // Allocate device memory
    cudaMalloc((void**)&d_cInstr, 3 * numCluster * sizeof(int));
    cudaMalloc((void**)&d_clusterMap, totalSize * sizeof(int));
    cudaMalloc((void**)&d_dataset, totalSize * sizeof(int));

    // Copy data from host to device
    cudaMemcpy(d_cInstr, cInstr, 3 * numCluster * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_clusterMap, clusterMap, totalSize * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_dataset, dataset, totalSize * sizeof(int), cudaMemcpyHostToDevice);

    // Launch kernel
    int threadsPerBlock = 256;
    int blocksPerGrid = (totalSize + threadsPerBlock - 1) / threadsPerBlock;
    clusterKernel<<<blocksPerGrid, threadsPerBlock>>>(d_cInstr, numCluster, totalSize, d_clusterMap, d_dataset);

    // Copy result back to host
    cudaMemcpy(clusterMap, d_clusterMap, totalSize * sizeof(int), cudaMemcpyDeviceToHost);

    // Free device memory
    cudaFree(d_cInstr);
    cudaFree(d_clusterMap);
    cudaFree(d_dataset);
}

// Function to initialize GPU memory for dataset (call this once at the start of your program)
int* initializeGPUDataset(int* dataset, int size)
{
    int* d_dataset;
    cudaMalloc((void**)&d_dataset, size * sizeof(int));
    cudaMemcpy(d_dataset, dataset, size * sizeof(int), cudaMemcpyHostToDevice);
    return d_dataset;
}

// Function to free GPU memory (call this at the end of your program)
void freeGPUMemory(int* d_dataset)
{
    cudaFree(d_dataset);
}

int main()
{
    srand(time(NULL));

    #define clusterSize 21
    
        int clusters[clusterSize] = {1, 2, 3, 98, 5, 6, 67, 89, 91, 100, 102, 10, 32, 45, 6, 7, 432, 54, 23, 66, 43};
        int clusterMap[clusterSize] = {0};
        treeNode *nodes[10];
        int n = -1;
        for (int i = 0; i < 10; i++)
        {
            nodes[i] = new treeNode();
        }

    // Initialize GPU memory for dataset
    int* d_clusters = initializeGPUDataset(clusters, clusterSize);

    getTwoRandomIndices(clusterMap, clusterSize, 0, nodes[0]);
    int cInstr[3] = {0, nodes[0]->child1, nodes[0]->child2};
    cluster(cInstr, 1, clusterSize, clusterMap, d_clusters);
    invalidateExtraOccurrences(clusterMap, clusterSize);

        while (n++ < clusterSize - 1)
    {
        printf("index %d value %d \n", n, clusterMap[n]);
    }

    getTwoRandomIndices(clusterMap, clusterSize, nodes[0]->nodeNum * 2 + 1, nodes[1]);
    getTwoRandomIndices(clusterMap, clusterSize, nodes[0]->nodeNum * 2 + 2, nodes[2]);
    int cInstr2[6] = {nodes[0]->nodeNum * 2 + 1, nodes[1]->child1, nodes[1]->child2,
                      nodes[0]->nodeNum * 2 + 2, nodes[2]->child1, nodes[2]->child2};
    cluster(cInstr2, 2, clusterSize, clusterMap, d_clusters);
    invalidateExtraOccurrences(clusterMap, clusterSize);

    n = -1;
        while (n++ < clusterSize - 1)
    {
        printf("index %d value %d \n", n, clusterMap[n]);
    }

    getTwoRandomIndices(clusterMap, clusterSize, nodes[1]->nodeNum * 2 + 1, nodes[3]);
    getTwoRandomIndices(clusterMap, clusterSize, nodes[1]->nodeNum * 2 + 2, nodes[4]);
    getTwoRandomIndices(clusterMap, clusterSize, nodes[2]->nodeNum * 2 + 1, nodes[5]);
    getTwoRandomIndices(clusterMap, clusterSize, nodes[2]->nodeNum * 2 + 2, nodes[6]);
    int cInstr3[12] = {nodes[1]->nodeNum * 2 + 1, nodes[3]->child1, nodes[3]->child2,
                       nodes[1]->nodeNum * 2 + 2, nodes[4]->child1, nodes[4]->child2,
                       nodes[2]->nodeNum * 2 + 2, nodes[5]->child1, nodes[5]->child2,
                       nodes[2]->nodeNum * 2 + 2, nodes[6]->child1, nodes[6]->child2};
    cluster(cInstr3, 4, clusterSize, clusterMap, d_clusters);
    invalidateExtraOccurrences(clusterMap, clusterSize);
   
    n = -1;
        while (n++ < clusterSize - 1)
    {
        printf("index %d value %d \n", n, clusterMap[n]);
    }

    // Free GPU memory
    freeGPUMemory(d_clusters);

    return 0;
}