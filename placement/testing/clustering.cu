#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <ctime>
#include <stdexcept>
#include <unordered_set>
#include <vector>
#include <random>
#include <algorithm>


#define THREADS_PER_BLOCK 256

const int clusterSize = 10000000;
const int threshold = 1000000;

struct treeNode {
    int nodeNum;
    int nodechild1;
    int nodechild2;
};

// CUDA kernel for the cluster function
__global__ void clusterKernel(int *cInstr, int numCluster, int clusterSize, int *clusterMap, int *dataset, int MAX_LEVELS, int * stopFlag, int* sharedCount ) {
    //Randomization to select 2 random 
   

    // __syncthreads(); 
    //Clustering phase
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < clusterSize) {
        if (clusterMap[idx] >= 0) {
            bool clusterFound = false;
            for (int clusterIdx = 0; clusterIdx < 3 * (numCluster); clusterIdx += 3) {
                if (cInstr[clusterIdx] == clusterMap[idx]) {
                    int distance1 = abs(dataset[cInstr[clusterIdx + 1]] - dataset[idx]);
                    int distance2 = abs(dataset[cInstr[clusterIdx + 2]] - dataset[idx]);
                    clusterMap[idx] = cInstr[clusterIdx] * 2 + (distance1 < distance2 ? 1 : 2);
                    clusterFound = true;
                    break;
                }
            }
            if (!clusterFound) {
                printf("Warning: No matching cluster found for clusterMap index %d with value %d\n", idx, clusterMap[idx]);
            }
        }
    }
    
    __syncthreads();

    // Invalidation step
    if(idx==0)
        stopFlag[0]=0;

    __syncthreads();
    

    if(idx < (1<<(MAX_LEVELS+1)))
        sharedCount[clusterMap[idx]] =0;
   
    __syncthreads();
    
    if (idx < clusterSize) 
    {
        if(clusterMap[idx]>0){
            atomicAdd(&sharedCount[clusterMap[idx]], 1);  
                if(sharedCount[clusterMap[idx]] <=threshold)
                    clusterMap[idx] = -clusterMap[idx] ;
        }
        else
            atomicAdd(&stopFlag[0], 1); 
    }


}

// Function to handle CUDA errors
void checkCudaError(cudaError_t error, const char *file, int line) {
    if (error != cudaSuccess) {
        printf("CUDA error at %s:%d: %s\n", file, line, cudaGetErrorString(error));
        exit(EXIT_FAILURE);
    }
}

#define CHECK_CUDA_ERROR(error) checkCudaError(error, __FILE__, __LINE__)

void getTwoRandomIndices(int *clusterMap, int clusterSize, int searchIndex, treeNode *node)
{

    if (!clusterMap || !node || clusterSize <= 0)
    {
        throw std::invalid_argument("Invalid arguments passed to getTwoRandomIndices");
    }

    std::unordered_set<int> uniqueIndices;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, clusterSize - 1); // Range from 0 to clusterSize - 1
    // int n = -1;
    // while (n++ < clusterSize - 1)
    // {
    //     printf("index %d value %d \n", n, clusterMap[n]);
    // }
    for (int i = 0; i < clusterSize; i++)
    {
        int indx = (dis(gen) + i) % clusterSize;
        if (clusterMap[indx] == searchIndex && uniqueIndices.find(indx) == uniqueIndices.end())
        {
            uniqueIndices.insert(indx);
            if (uniqueIndices.size() == 2)
                break;
        }
    }
    // if )
    // {
    //     throw std::runtime_error("Not enough unique indices found for the search index");
    // }
    if (uniqueIndices.size() >= 2)
    {

        auto it = uniqueIndices.begin();
        node->nodechild1 = *it++;
        node->nodechild2 = *it;
        node->nodeNum = searchIndex;
    }
}

void processClusterLevels(int *clusterMap, int clusterSize, treeNode *nodes[], int *d_dataset, int MAX_LEVELS) {
    if (!clusterMap || !nodes || !d_dataset || clusterSize <= 0) {
        throw std::invalid_argument("Invalid arguments passed to processClusterLevels");
    }

    int nodeIndex = 0;
    int *stopFlag = new int;
    int size = (1<<(MAX_LEVELS+1));
    int *sharedCount = new int[size];
    for(int i=0;i<size;i++)
        sharedCount[i]=0;

    int *d_cInstr, *d_clusterMap, *d_stopFlag,*d_sharedCount;
    for (int level = 0; level < MAX_LEVELS; level++) {
        int nodesInThisLevel = 1 << level;
        int totalInstructions = nodesInThisLevel * 3;

        int *cInstr = new int[totalInstructions];
        int instrIndex = 0;

        for (int i = 0; i < nodesInThisLevel; i++) {
            int parentIndex = (nodeIndex - 1) / 2;
            int baseClusterIndex = (level == 0) ? 0 : nodes[parentIndex]->nodeNum * 2 + i % 2 + 1;

            try {
                getTwoRandomIndices(clusterMap, clusterSize, baseClusterIndex, nodes[nodeIndex]);
            } catch (const std::exception &e) {
                printf("Error in getTwoRandomIndices: %s\n", e.what());
                delete[] cInstr;
                return;
            }

            cInstr[instrIndex++] = baseClusterIndex;
            cInstr[instrIndex++] = nodes[nodeIndex]->nodechild1;
            cInstr[instrIndex++] = nodes[nodeIndex]->nodechild2;
            nodeIndex++;
        }

        try {

            CHECK_CUDA_ERROR(cudaMalloc(&d_stopFlag, sizeof(int)));
            CHECK_CUDA_ERROR(cudaMemcpy(d_stopFlag, stopFlag, sizeof(int), cudaMemcpyHostToDevice));

            CHECK_CUDA_ERROR(cudaMalloc(&d_sharedCount, sizeof(int)));
            CHECK_CUDA_ERROR(cudaMemcpy(d_clusterMapCount, sharedCount, sizeof(int), cudaMemcpyHostToDevice));
            
            CHECK_CUDA_ERROR(cudaMalloc(&d_cInstr, 3 * nodesInThisLevel * sizeof(int)));
            CHECK_CUDA_ERROR(cudaMalloc(&d_clusterMap, clusterSize * sizeof(int)));
        
            // Copy data to device
            CHECK_CUDA_ERROR(cudaMemcpy(d_cInstr, cInstr, 3 * nodesInThisLevel * sizeof(int), cudaMemcpyHostToDevice));
            CHECK_CUDA_ERROR(cudaMemcpy(d_clusterMap, clusterMap, clusterSize * sizeof(int), cudaMemcpyHostToDevice));
        
            // Launch kernel
            int blocksPerGrid = (clusterSize + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
            clusterKernel<<<blocksPerGrid, THREADS_PER_BLOCK>>>(d_cInstr, nodesInThisLevel, clusterSize, d_clusterMap, d_dataset,MAX_LEVELS, d_stopFlag, d_sharedCount);
        
            CHECK_CUDA_ERROR(cudaMemcpy(clusterMap, d_clusterMap, clusterSize* sizeof(int), cudaMemcpyDeviceToHost));
            CHECK_CUDA_ERROR(cudaMemcpy(stopFlag, d_stopFlag, sizeof(int), cudaMemcpyDeviceToHost));
          
            // clusterGPU(cInstr, nodesInThisLevel, clusterSize, clusterMap, d_dataset);

            // int n = -1;
            // while (n++ < clusterSize - 1)
            // {
            //     printf("cluster value index %d value %d \n", n, clusterMap[n]);
            // }

            if(stopFlag[0] != 0)
                break;

            // Check for kernel launch errors
            CHECK_CUDA_ERROR(cudaGetLastError());
        
            // Copy result back to host
          




            // clusterGPU(cInstr, nodesInThisLevel, clusterSize, clusterMap, d_dataset);
            // if (invalidateExtraOccurrences(clusterMap, clusterSize)) {
            //     delete[] cInstr;
            //     return;
            // }
        } catch (const std::exception &e) {
            printf("Error in cluster or invalidateExtraOccurrences: %s\n", e.what());
            delete[] cInstr;
            return;
        }

        delete[] cInstr;
        CHECK_CUDA_ERROR(cudaFree(d_sharedCount));
        CHECK_CUDA_ERROR(cudaFree(d_cInstr));
        CHECK_CUDA_ERROR(cudaFree(d_stopFlag));
    }
}



int main() {
    int *boundary = new int;
    int clusterSizeVar = (int)clusterSize;
    int MAX_LEVELS = 0;
    while (clusterSizeVar >>= 1) ++MAX_LEVELS;

    std::vector<int> clustersVec(clusterSize);
    std::generate(clustersVec.begin(), clustersVec.end(), []() { return rand() % 10000 + 1; });

    int *clusters = new int[clusterSize];
    std::copy(clustersVec.begin(), clustersVec.end(), clusters);

    int *clusterMap = new int[clusterSize]();

    treeNode **nodes = new treeNode*[1 << MAX_LEVELS];
    for (int i = 0; i < 1 << MAX_LEVELS; i++) {
        nodes[i] = new treeNode();
    }

    // Allocate and copy dataset to GPU
    int *d_dataset;
    CHECK_CUDA_ERROR(cudaMalloc(&d_dataset, clusterSize * sizeof(int)));
    CHECK_CUDA_ERROR(cudaMemcpy(d_dataset, clusters, clusterSize * sizeof(int), cudaMemcpyHostToDevice));

    processClusterLevels(clusterMap, clusterSize, nodes, d_dataset, MAX_LEVELS);


    // int n = -1;
    // while (n++ < clusterSize - 1)
    // {
    //     printf("cluster value %d index %d value %d \n", clusters[n], n, clusterMap[n]);
    // }

    // Clean up
    for (int i = 0; i < 1 << MAX_LEVELS; i++) {
        delete nodes[i];
    }
    delete[] nodes;
    delete[] clusters;
    delete[] clusterMap;
    delete boundary;

    // Free GPU memory
    CHECK_CUDA_ERROR(cudaFree(d_dataset));

    return 0;
}