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
#include <vector>
#include <random>
#include <algorithm>

#define THREADS_PER_BLOCK 256

const int clusterSize = 1000000;
// Adjust this based on how many levels you want
// const int MAX_NODES = 21;
const int threshold = 10000;

struct treeNode
{
    int nodeNum;
    int nodechild1;
    int nodechild2;
    int valuechild1;
    int valuechild2;
};

void cluster(int *cInstr, int numCluster, int totalSize, int *clusterMap, int *dataset)
{
    if (!cInstr || !clusterMap || !dataset)
    {
        throw std::invalid_argument("Null pointer passed to cluster function");
    }
    if (numCluster <= 0 || totalSize <= 0)
    {
        throw std::invalid_argument("Invalid cluster or size parameters");
    }

    for (int idx = 0; idx < totalSize; idx++)
    {
        if (clusterMap[idx] >= 0)
        {
            bool clusterFound = false;
            for (int clusterIdx = 0; clusterIdx < 3 * (numCluster); clusterIdx += 3)
            {
                if (cInstr[clusterIdx] == clusterMap[idx])
                {
                    int distance1 = std::abs(dataset[cInstr[clusterIdx + 1]] - dataset[idx]);
                    int distance2 = std::abs(dataset[cInstr[clusterIdx + 2]] - dataset[idx]);
                    clusterMap[idx] = cInstr[clusterIdx] * 2 + (distance1 < distance2 ? 1 : 2);
                    clusterFound = true;
                    break;
                }
            }
            if (!clusterFound)
            {
                printf("Warning: No matching cluster found for clusterMap index %d with value %d\n", idx, clusterMap[idx]);
            }
        }
    }
}

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
    int n = -1;
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

int invalidateExtraOccurrences(int *arr, int size)
{
    if (!arr || size <= 0)
    {
        throw std::invalid_argument("Invalid arguments passed to invalidateExtraOccurrences");
    }

    int maxNum = *std::max_element(arr, arr + size);
    if (maxNum < 0)
    {
        printf("Built clusters with give max sub-tree size\n");
        // break;
        return 1;

        // return;
    }

    int resultSize = maxNum + 1;
    int *counts = new int[resultSize]();

    for (int i = 0; i < size; ++i)
    {
        if (arr[i] >= 0 && arr[i] < resultSize)
        {
            counts[arr[i]]++;
        }
    }

    for (int i = 0; i < size; ++i)
    {
        if (arr[i] >= 0 && arr[i] < resultSize && counts[arr[i]] < threshold)
        {
            arr[i] = -arr[i];
        }
    }

    delete[] counts;
    return 0;
}

void processClusterLevels(int *clusterMap, int clusterSize, treeNode *nodes[], int *clusters, int MAX_LEVELS)
{
    if (!clusterMap || !nodes || !clusters || clusterSize <= 0)
    {
        throw std::invalid_argument("Invalid arguments passed to processClusterLevels");
    }

    int nodeIndex = 0;

    for (int level = 0; level < MAX_LEVELS; level++)
    {
        int nodesInThisLevel = 1 << level;
        int totalInstructions = nodesInThisLevel * 3;

        // if (nodeIndex + nodesInThisLevel > MAX_NODES)
        // {
        //     printf("Warning: Not enough nodes available for level %d\n", level);
        // //break;
        // }

        int *cInstr = new int[totalInstructions];
        int instrIndex = 0;

        for (int i = 0; i < nodesInThisLevel; i++)
        {
            int parentIndex = (nodeIndex - 1) / 2;
            int baseClusterIndex = (level == 0) ? 0 : nodes[parentIndex]->nodeNum * 2 + i % 2 + 1;

            try
            {
                getTwoRandomIndices(clusterMap, clusterSize, baseClusterIndex, nodes[nodeIndex]);
                // printf("nodes are %d, index=%d and %d, index=%d\n", clusters[nodes[nodeIndex]->nodechild1], nodes[nodeIndex]->nodechild1, clusters[nodes[nodeIndex]->nodechild2], nodes[nodeIndex]->nodechild2);
            }
            catch (const std::exception &e)
            {
                printf("Error in getTwoRandomIndices: %s\n", e.what());
                delete[] cInstr;
                return;
            }

            cInstr[instrIndex++] = baseClusterIndex;
            cInstr[instrIndex++] = nodes[nodeIndex]->nodechild1;
            cInstr[instrIndex++] = nodes[nodeIndex]->nodechild2;

            nodeIndex++;
        }

        try
        {
            cluster(cInstr, nodesInThisLevel, clusterSize, clusterMap, clusters);
            if (invalidateExtraOccurrences(clusterMap, clusterSize))
            {
                delete[] cInstr;
                return;
            }
        }
        catch (const std::exception &e)
        {
            printf("Error in cluster or invalidateExtraOccurrences: %s\n", e.what());
            delete[] cInstr;
            return;
        }

        // for (int n = 0; n < clusterSize; n++)
        // {
        //     printf("index %d value %d \n", n, clusterMap[n]);
        // }

        delete[] cInstr;
    }
}

int *convertVectorToIntArray(const std::vector<int> &vec)
{
    int *arr = new int[vec.size()];
    std::copy(vec.begin(), vec.end(), arr);
    return arr;
}

int main()
{
    int *boundary = new int;
    int clusterSizeVar = (int)clusterSize;
    int MAX_LEVELS;
    while (clusterSizeVar >>= 1)
        ++MAX_LEVELS;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(1, 10000); // Range from 1 to 10000

    // Create and fill the vector with random numbers
    std::vector<int> clustersVec(clusterSize);
    std::generate(clustersVec.begin(), clustersVec.end(), [&]()
                  { return dis(gen); });

    int *clusters = convertVectorToIntArray(clustersVec);

    // int clusters[clusterSize] = {1, 2, 3, 98, 5, 6, 67, 89, 91, 100, 102, 10, 32, 45, 6, 7, 432, 54, 23, 66, 43};
    // int clusters[clusterSize] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
    int clusterMap[clusterSize] = {0};
    treeNode *nodes[1 << MAX_LEVELS];
    for (int i = 0; i < 1 << MAX_LEVELS; i++)
    {
        nodes[i] = new treeNode();
    }

    processClusterLevels(clusterMap, clusterSize, nodes, clusters, MAX_LEVELS);

    // getTwoRandomIndices(clusterMap, clusterSize, 0, nodes[0]);
    // int cInstr[3] = {0, nodes[0]->child1, nodes[0]->child2};
    // cluster(cInstr, 1, clusterSize, clusterMap, clusters);
    // invalidateExtraOccurrences(clusterMap, clusterSize);
    // int n = -1;
    // while (n++ < clusterSize - 1)
    // {
    //     printf("index %d value %d \n", n, clusterMap[n]);
    // }

    // getTwoRandomIndices(clusterMap, clusterSize, nodes[0]->nodeNum * 2 + 1, nodes[1]);
    // getTwoRandomIndices(clusterMap, clusterSize, nodes[0]->nodeNum * 2 + 2, nodes[2]);
    // int cInstr2[6] = {nodes[0]->nodeNum * 2 + 1, nodes[1]->child1, nodes[1]->child2,
    //                   nodes[0]->nodeNum * 2 + 2, nodes[2]->child1, nodes[2]->child2};
    // cluster(cInstr2, 2, clusterSize, clusterMap, clusters);
    // invalidateExtraOccurrences(clusterMap, clusterSize);
    // n = -1;
    // while (n++ < clusterSize - 1)
    // {
    //     printf("index %d value %d \n", n, clusterMap[n]);
    // }

    // getTwoRandomIndices(clusterMap, clusterSize, nodes[1]->nodeNum * 2 + 1, nodes[3]);
    // getTwoRandomIndices(clusterMap, clusterSize, nodes[1]->nodeNum * 2 + 2, nodes[4]);
    // getTwoRandomIndices(clusterMap, clusterSize, nodes[2]->nodeNum * 2 + 1, nodes[5]);
    // getTwoRandomIndices(clusterMap, clusterSize, nodes[2]->nodeNum * 2 + 2, nodes[6]);
    // int cInstr3[12] = {nodes[1]->nodeNum * 2 + 1, nodes[3]->child1, nodes[3]->child2,
    //                    nodes[1]->nodeNum * 2 + 2, nodes[4]->child1, nodes[4]->child2,
    //                    nodes[2]->nodeNum * 2 + 2, nodes[5]->child1, nodes[5]->child2,
    //                    nodes[2]->nodeNum * 2 + 2, nodes[6]->child1, nodes[6]->child2};
    // cluster(cInstr3, 4, clusterSize, clusterMap, clusters);
    // invalidateExtraOccurrences(clusterMap, clusterSize);
    int n = -1;
    // while (n++ < clusterSize - 1)
    // {
    //     printf("cluster value %d index %d value %d \n", clusters[n], n, clusterMap[n]);
    // }
    delete[] clusters;
    return 0;
}