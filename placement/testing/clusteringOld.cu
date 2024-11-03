// clustering.cu

#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>

#define THREADS_PER_BLOCK 256

// Function to check for CUDA errors
#define cudaCheckError(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true){
    if (code != cudaSuccess){
        fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}

// Distance function (Euclidean distance in this case)
__device__ float distance(float x1, float x2){
    return fabsf(x1 - x2);
}

// Kernel to compute distances and assign clusters
__global__ void assignClusters(float *E, int size, float a, float b, int *clusterAssignments){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= size) return;
    
    float distToA = distance(E[idx], a);
    float distToB = distance(E[idx], b);
    
    if(distToA < distToB)
        clusterAssignments[idx] = 0; // Assign to cluster Ea
    else
        clusterAssignments[idx] = 1; // Assign to cluster Eb
}

// Host function to perform clustering recursively
void cluster(float *E, int size, int k, std::vector<std::vector<float>> &clusters){
    if(size <= k){
        // Base case: add the cluster to the list
        clusters.push_back(std::vector<float>(E, E + size));
        return;
    }
    
    // Randomly select two pivot points a and b
    srand(time(NULL));
    int idxA = rand() % size;
    int idxB = rand() % size;
    while(idxB == idxA){
        idxB = rand() % size;
    }
    float a = E[idxA];
    float b = E[idxB];
    
    // Allocate memory on device
    float *d_E;
    int *d_clusterAssignments;
    cudaCheckError(cudaMalloc((void**)&d_E, size * sizeof(float)));
    cudaCheckError(cudaMalloc((void**)&d_clusterAssignments, size * sizeof(int)));
    
    // Copy data to device
    cudaCheckError(cudaMemcpy(d_E, E, size * sizeof(float), cudaMemcpyHostToDevice));
    
    // Launch kernel to assign clusters
    int blocks = (size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    assignClusters<<<blocks, THREADS_PER_BLOCK>>>(d_E, size, a, b, d_clusterAssignments);
    cudaDeviceSynchronize();
    
    // Copy cluster assignments back to host
    int *clusterAssignments = new int[size];
    cudaCheckError(cudaMemcpy(clusterAssignments, d_clusterAssignments, size * sizeof(int), cudaMemcpyDeviceToHost));
    
    // Separate elements into two clusters
    std::vector<float> Ea, Eb;
    for(int i = 0; i < size; i++){
        if(i == idxA || i == idxB) continue; // Exclude a and b
        if(clusterAssignments[i] == 0)
            Ea.push_back(E[i]);
        else
            Eb.push_back(E[i]);
    }
    
    // Free device memory
    cudaCheckError(cudaFree(d_E));
    cudaCheckError(cudaFree(d_clusterAssignments));
    delete[] clusterAssignments;
    
    // Recursively cluster Ea and Eb
    if(!Ea.empty()){
        cluster(Ea.data(), Ea.size(), k, clusters);
    }
    if(!Eb.empty()){
        cluster(Eb.data(), Eb.size(), k, clusters);
    }
}



int main(){
    // Example data array E
    int N = 10000; // Size of E
    int k = 1000;  // Maximum cluster size
    float *E = new float[N];
    
    // Initialize E with random data
    srand(time(NULL));
    for(int i = 0; i < N; i++){
        E[i] = static_cast<float>(rand()) / RAND_MAX;
    }
    
    // Vector to store clusters
    std::vector<std::vector<float>> clusters = {1,2,3,4,5,6,67,89,91,100,102,103};

    std::vector<std::vector<float>> distMatrix1=  distanceCal(clusters,8);
    std::vector<std::vector<float>> distMatrix2=  distanceCal(clusters,1);
    
    // // Perform clustering
    // cluster(E, N, k, clusters);
    
    // // Output results
    // std::cout << "Number of clusters formed: " << clusters.size() << std::endl;
    // for(size_t i = 0; i < clusters.size(); i++){
    //     std::cout << "Cluster " << i+1 << " size: " << clusters[i].size() << std::endl;
    // }
    
    // // Clean up
    // delete[] E;
    return 0;
}