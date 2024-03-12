#ifndef MASH_CUH
#include "mash.cuh"
#endif

#include <stdio.h>
#include <queue>
#include <thrust/sort.h>
#include <thrust/scan.h>
#include <thrust/binary_search.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <chrono>
#include <iostream>

/* Note: d_aggseqLengths is the aggregated compressed length of original string (h_seqLengths)
Ex: ["dog", "mouse", "cat"] 
h_seqLengths -> [3, 5, 3] 
d_aggseqLengths -> [1, 2, 3] */
void GpuSketch::DeviceArrays::allocateDeviceArrays(uint32_t ** h_compressedSeqs, uint32_t * h_seqLengths, size_t numSequences, Param& params)
{
    cudaError_t err;

    d_numSequences = numSequences;

    uint32_t kmerSize = params.kmerSize;
    size_t hashListLength = 0;   

    // Allocate memory
    err = cudaMalloc(&d_aggseqLengths, numSequences*sizeof(uint32_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_seqLengths, numSequences*sizeof(uint32_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    /* Flatten data */
    uint32_t * h_aggseqLengths = new uint32_t[numSequences];
    uint32_t flatStringLength=0;
    for (size_t i =0; i<numSequences; i++) flatStringLength+= (h_seqLengths[i]+15)/16;
    uint32_t * h_flattenCompressSeqs = new uint32_t[flatStringLength];
    flatStringLength=0;
    for (size_t i =0; i<numSequences; i++) 
    {
        uint32_t flatStringLengthLocal = (h_seqLengths[i]+15)/16;
        hashListLength += h_seqLengths[i] - kmerSize + 1;
        flatStringLength+=flatStringLengthLocal;
        for (size_t j=0; j<flatStringLengthLocal;j++)  
        {
            h_flattenCompressSeqs[j] = h_compressedSeqs[i][j];
            // if (i==9) printf("%u\n",h_flattenCompressSeqs[j]); 
        }
        h_flattenCompressSeqs += flatStringLengthLocal;
        h_aggseqLengths[i] = flatStringLength;
    }

    h_flattenCompressSeqs -= flatStringLength;
    //printf("%d", flatStringLength);

    err = cudaMalloc(&d_hashList, hashListLength*sizeof(uint32_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }


    err = cudaMalloc(&d_compressedSeqs, flatStringLength*sizeof(uint32_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    err = cudaMalloc(&d_mashDist, (numSequences*(numSequences -1)/2)*sizeof(float));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
        exit(1);
    }

    // Transfer data
    err = cudaMemcpy(d_aggseqLengths, h_aggseqLengths, numSequences*sizeof(uint32_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) 
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
        exit(1);
    }

    err = cudaMemcpy(d_seqLengths, h_seqLengths, numSequences*sizeof(uint32_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) 
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
        exit(1);
    }

    err = cudaMemcpy(d_compressedSeqs, h_flattenCompressSeqs, flatStringLength*sizeof(uint32_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) 
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
        exit(1);
    }

    cudaDeviceSynchronize();
}

void GpuSketch::DeviceArrays::deallocateDeviceArrays(){
    cudaFree(d_compressedSeqs);
    cudaFree(d_aggseqLengths);
    cudaFree(d_seqLengths);
    cudaFree(d_hashList);
    // cudaFree(d_mashDist);
}

__device__ uint32_t fmix32 ( uint32_t h )
{
  h ^= h >> 16;
  h *= 0x85ebca6b;
  h ^= h >> 13;
  h *= 0xc2b2ae35;
  h ^= h >> 16;

  return h;
}

__device__ uint32_t MurmurHash3_x86_32 ( uint32_t key, int len, uint32_t seed)
{
    uint32_t h1 = seed;

    const uint32_t c1 = 0xcc9e2d51;
    const uint32_t c2 = 0x1b873593;

    uint32_t k1 = key;

    k1 *= c1;
    k1 = (k1 << 15) | (k1 >> (32 - 15));
    k1 *= c2;
    h1 ^= k1;
    h1 = (h1 << 13) | (h1 >> (32 - 13));
    h1 = h1*5+0xe6546b64;

    h1 ^= len;
    h1 = fmix32(h1);

    return h1;
} 

__global__ void sketchConstructionSerial
(
    uint32_t * d_compressedSeqs,
    uint32_t * d_aggseqLengths,
    uint32_t * d_seqLengths,
    size_t d_numSequences,
    uint32_t * d_hashList,
    uint32_t kmerSize
){
    int tx = threadIdx.x;
    int bx = blockIdx.x;


    uint32_t kmer = 0;
    uint32_t mask = (1<<2*kmerSize) - 1;

    uint32_t * hashList = d_hashList;
    uint32_t * compressedSeqs = d_compressedSeqs;

    //printf("hashList pointer in device%p\n", d_hashList);

    if (tx==0 && bx==0)
    {
        for (size_t i=0; i<d_numSequences; i++)
        {
            uint32_t seqLength = d_seqLengths[i];
            
            //if (i==9)printf("%ld:\t", i);

            for (size_t j=0; j<=seqLength-kmerSize; j++)
            {
                uint32_t index = j/16;
                uint32_t shift1 = 2*(j%16);
                if (shift1>0)
                {
                    uint32_t shift2 = 32-shift1;
                    kmer = ((compressedSeqs[index] >> shift1) | (compressedSeqs[index+1] << shift2)) & mask;
                }
                else
                {   
                    kmer = compressedSeqs[index] & mask;
                }
                uint32_t hash = MurmurHash3_x86_32  (kmer, 30, 53);
                hashList[j] = hash;
                //if (i==9) printf("(%u, %u)\t",kmer, hash);
            }
            //if (i==9) printf("\n");
            hashList += seqLength-kmerSize+1;
            compressedSeqs += (seqLength+15)/16;

        }
    }
    //printf("hashList pointer in device%p\n", d_hashList);

}


__global__ void sketchConstruction
(
    uint32_t * d_compressedSeqs,
    uint32_t * d_aggseqLengths,
    uint32_t * d_seqLengths,
    uint32_t * d_prefixHashlist,
    uint32_t * d_prefixCompressed,
    size_t d_numSequences,
    uint32_t * d_hashList,
    uint32_t kmerSize
){
    size_t tx = threadIdx.x;
    size_t bx = blockIdx.x;
    size_t threads_per_block = blockDim.x;
    size_t blocks_per_grid = gridDim.x;

    uint32_t kmer = 0;
    uint32_t mask = (1<<2*kmerSize) - 1;

    //printf("hashList pointer in device%p\n", d_hashList);

    if (bx >= d_numSequences) return;

    for (size_t j = bx; j < d_numSequences; j+=blocks_per_grid) 
    {
        uint32_t seqLength = d_seqLengths[j];
        uint32_t * hashList = d_hashList + d_prefixHashlist[j];
        uint32_t * compressedSeqs = d_compressedSeqs + d_prefixCompressed[j];
        
        for (size_t i = tx; i <= seqLength - kmerSize; i += threads_per_block) 
        {
            
            uint32_t index = i/16;
            uint32_t shift1 = 2*(i%16);
            if (shift1>0)
            {
                uint32_t shift2 = 32-shift1;
                kmer = ((compressedSeqs[index] >> shift1) | (compressedSeqs[index+1] << shift2)) & mask;
            }
            else
            {   
                kmer = compressedSeqs[index] & mask;
            }
            uint32_t hash = MurmurHash3_x86_32  (kmer, 30, 53);
            //if (j==0) printf("(%u, %u)\n",kmer, hash);
            hashList[i] = hash;
        }

    }
}


__device__ void swap(uint32_t &a, uint32_t &b) {
    uint32_t temp = a;
    a = b;
    b = temp;
}


__device__ int partition(uint32_t *arr, int low, int high) {
    uint32_t pivot = arr[high];
    int i = (low - 1);
    
    for (int j = low; j <= high- 1; j++) {
        if (arr[j] < pivot) {
            i++;
            swap(arr[i], arr[j]);
        }
    }
    swap(arr[i + 1], arr[high]);
    return (i + 1);
}



// Device function for insertion sort
__device__ void insertionSort(uint32_t *arr, int left, int right) {
    for (int i = left + 1; i <= right; i++) {
        int key = arr[i];
        int j = i - 1;
        
        while (j >= left && arr[j] > key) {
            arr[j + 1] = arr[j];
            j = j - 1;
        }
        arr[j + 1] = key;
    }
}



__device__ int findMedian(uint32_t *arr, int n) {
    insertionSort(arr, 0, n - 1);
    return arr[n / 2];
}



__device__ int quickSelect(uint32_t *arr, int left, int right, size_t k) {
    if (left == right) return arr[left];
    
    int pi = partition(arr, left, right);

    int kth = pi - left + 1;
    if (k == kth) return arr[pi];
    else if (k < kth) return quickSelect(arr, left, pi - 1, k);
    else return quickSelect(arr, pi + 1, right, k - kth);
}



__global__ void introSelectKernel(uint32_t *arr, int n, int k) {
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        *arr = quickSelect(arr, 0, n - 1, k);
    }
}


__global__ void pruneHashes
(
    uint32_t * d_prefixHashlist,
    size_t d_numSequences,
    uint32_t * d_hashList,
    uint32_t * d_hashListPruned
){
    size_t tx = threadIdx.x;
    size_t bx = blockIdx.x;
    size_t threads_per_block = blockDim.x;
    size_t blocks_per_grid = gridDim.x;

    if (bx >= d_numSequences) return;

    for (size_t j = bx; j < d_numSequences; j+=blocks_per_grid) 
    {
        uint32_t * hashList = d_hashList + d_prefixHashlist[j];
        
        for (size_t i = tx; i < 1000; i += threads_per_block) 
        {
            d_hashListPruned[j * 1000 + i] = hashList[i];
        }

    }
}



void GpuSketch::selectNthSmallestOnGpu
(
    uint32_t* d_hashList, 
    uint32_t * h_seqLengths, 
    size_t d_numSequences, 
    Param& params
) {
    std::vector<cudaStream_t> streams(d_numSequences); // tune the number of streams needed, can reuse

    uint32_t * hashList = d_hashList;
    for (size_t i = 0; i < d_numSequences; i++)
    {
        cudaStreamCreate(&streams[i]);
        uint32_t numKmers = (h_seqLengths[i] - params.kmerSize + 1);
        introSelectKernel<<<params.numBlocks, params.blockSize, 0, streams[i]>>>(hashList, numKmers, 1000);
        hashList += numKmers;
    }

    // Wait for all streams to complete
    for (size_t i = 0; i < d_numSequences; ++i) {
        cudaStreamSynchronize(streams[i]);
        cudaStreamDestroy(streams[i]);
    }
}




void GpuSketch::sketchConstructionOnGpu
(
    uint32_t * d_compressedSeqs,
    uint32_t * d_aggseqLengths,
    uint32_t * d_seqLengths,
    size_t d_numSequences,
    uint32_t * d_hashList,
    uint32_t * h_seqLengths,
    Param& params
){
    auto timerStart = std::chrono::high_resolution_clock::now();

    cudaError_t err;

    // prefix-sum of d_seqLengths using thrust
    uint32_t * d_prefixHashlist;
    uint32_t * d_prefixCompressed;

    int bytes = d_numSequences * sizeof(uint32_t);

    cudaMalloc(&d_prefixHashlist, bytes);
    cudaMalloc(&d_prefixCompressed, bytes);
    

    thrust::device_ptr<uint32_t> dev_seqLengths(d_seqLengths);
    thrust::device_ptr<uint32_t> dev_prefixHash(d_prefixHashlist);
    thrust::device_ptr<uint32_t> dev_prefixComp(d_prefixCompressed);

    thrust::transform(thrust::device,
        dev_seqLengths, dev_seqLengths + d_numSequences, dev_prefixComp, 
        [] __device__ (const uint32_t& x) -> uint32_t { 
            return (x + 15) / 16;
        }
    );

    const uint32_t kmerSize = params.kmerSize; // Extract kmerSize

    thrust::transform(
        thrust::device,
        dev_seqLengths, dev_seqLengths + d_numSequences, dev_prefixHash,
        [kmerSize] __device__ (const uint32_t& x) -> uint32_t {
            return x - kmerSize + 1;
        }
    );

    thrust::exclusive_scan(dev_prefixHash, dev_prefixHash + d_numSequences, dev_prefixHash);
    thrust::exclusive_scan(dev_prefixComp, dev_prefixComp + d_numSequences, dev_prefixComp);
    auto timerEnd = std::chrono::high_resolution_clock::now();

    std::chrono::nanoseconds time = timerEnd - timerStart;
    std::cout << "Time to create prefix array: " << time.count() << "ns\n";

    

    timerStart = std::chrono::high_resolution_clock::now();
    // Serial kernel call
    //sketchConstructionSerial<<<params.numBlocks, params.blockSize>>>(d_compressedSeqs, d_aggseqLengths, d_seqLengths, d_numSequences, d_hashList, kmerSize);

    // New kernel call
    sketchConstruction<<<params.numBlocks, params.blockSize>>>(d_compressedSeqs, d_aggseqLengths, d_seqLengths, d_prefixHashlist, d_prefixCompressed, d_numSequences, d_hashList, kmerSize);
    
    cudaDeviceSynchronize();

    timerEnd = std::chrono::high_resolution_clock::now();
    time = timerEnd - timerStart;
    std::cout << "Time to generate hashes: " << time.count() << "ns\n";

    cudaFree(d_prefixCompressed);
    

    timerStart = std::chrono::high_resolution_clock::now();

    // // Old implementation -- for correctness
    // uint32_t * hashList = d_hashList;
    // for (size_t i = 0; i < d_numSequences; i++)
    // {
    //     thrust::device_ptr<uint32_t> hashPtr(hashList);
    //     uint32_t numKmers = (h_seqLengths[i] - kmerSize + 1);   
    //     thrust::sort(hashPtr, hashPtr + numKmers);
    //     hashList += numKmers;
    // }

    selectNthSmallestOnGpu(d_hashList, h_seqLengths, d_numSequences, params);

    timerEnd = std::chrono::high_resolution_clock::now();
    time = timerEnd - timerStart;
    std::cout << "Time to sort: " << time.count() << "ns\n";

    uint32_t * d_hashListPruned;
    bytes = d_numSequences * 1000 * sizeof(uint32_t);
    cudaMalloc(&d_hashListPruned, bytes);
    uint32_t * h_pruned = (uint32_t*)malloc(bytes);
    pruneHashes<<<params.numBlocks, params.blockSize>>>(d_prefixHashlist, d_numSequences, d_hashList, d_hashListPruned);

    cudaMemcpy(h_pruned, d_hashListPruned, bytes, cudaMemcpyDeviceToHost);
    for (int i = 0; i < d_numSequences; i++) {
        printf("i       hashList[i] (%d)\n", i);
        for (int j = 0; j < 1000; j++) {
            printf("%d       %d\n", j, h_pruned[i * 1000 + j]);

        }
    }

    // prevent mem leak for now
    cudaFree(d_hashListPruned);
    free(h_pruned);

    cudaFree(d_prefixHashlist);
    

    

}

__device__ float mashDistance
(
    uint32_t * A,
    uint32_t * B,
    uint32_t kmerSize,
    uint32_t sketchSize
){
    uint32_t unionPtr = 0, APtr = 0, BPtr = 0;
    float inter=0, uni=0;

    while (true)
    {
        if ((APtr >= sketchSize) || (BPtr >= sketchSize) || (unionPtr >= sketchSize)) break;

        if (A[APtr]==B[BPtr]) 
        {   
            inter++; uni++;
            APtr++; BPtr++; unionPtr++;
        } 
        else if (A[APtr]>B[BPtr])
        {
            uni++;
            BPtr++; unionPtr++;
        }
        else
        {
            uni++;
            APtr++; unionPtr++;
        }
    }

    while (unionPtr<sketchSize-1 && (APtr<sketchSize-1 || BPtr<sketchSize-1))
    {
        if(APtr<sketchSize-1) {unionPtr++; APtr++; uni++;}
        if(BPtr<sketchSize-1) {unionPtr++; BPtr++; uni++;}
    }

    if (unionPtr<sketchSize-1)
    {
        printf("Error: Not enough hashes to build %u size union sketch\n", sketchSize);
    }
    float jaccardEstimate = (inter/uni);

    float mashDist = (log(2.0*jaccardEstimate/(1.0+jaccardEstimate)))/kmerSize;

    return mashDist;

}

__global__ void mashDistConstruction
(
    uint32_t * d_hashList,
    uint32_t * d_seqLengths,
    size_t d_numSequences,
    float * d_mashDist,
    uint32_t kmerSize,
    uint32_t sketchSize
){
    int tx = threadIdx.x;
    int bx = blockIdx.x;

    uint32_t * hashList = d_hashList;


    if (tx==0 && bx==0)
    {
        uint32_t mashDistCount=0;
        for (size_t i=0; i<d_numSequences; i++)
        {
            uint32_t * hashListStartIndex = hashList + d_seqLengths[i] - kmerSize + 1;
            for (size_t j=i+1; j<d_numSequences; j++)
            {
                float mashDist = mashDistance(hashList, hashListStartIndex, kmerSize, sketchSize);
                d_mashDist[mashDistCount++] = mashDist;
                hashListStartIndex += d_seqLengths[j] - kmerSize + 1;
            }
            hashList += d_seqLengths[i] - kmerSize + 1;
        }
    }
}


void GpuSketch::mashDistConstructionOnGpu
(   
    uint32_t * d_hashList,
    uint32_t * d_seqLengths,
    size_t d_numSequences,
    float * d_mashDist,
    uint32_t * h_seqLengths,
    Param& params
){

    mashDistConstruction<<<params.numBlocks, params.blockSize>>>(d_hashList, d_seqLengths, d_numSequences, d_mashDist, params.kmerSize, params.sketchSize);

    cudaDeviceSynchronize();

}

void GpuSketch::DeviceArrays::printSketchValues(int numValues, uint32_t * h_seqLengths) 
{
    uint32_t * h_hashList = new uint32_t[numValues];

    uint32_t * hashList = d_hashList;

    cudaError_t err;

    for (size_t j = 0; j < d_numSequences; j++)
    {

        err = cudaMemcpy(h_hashList, hashList, numValues*sizeof(uint32_t), cudaMemcpyDeviceToHost);
        if (err != cudaSuccess) {
            fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
            exit(1);
        }

        printf("i\thashList[i] (%zu)\n",j);
        for (int i=0; i<numValues; i++) {
            printf("%i\t%u\n", i, h_hashList[i]);
        }
        hashList += h_seqLengths[j] - 15 + 1;
       
    }

}

void GpuSketch::DeviceArrays::printMashDist(uint32_t h_numSequences) 
{
    
    float * h_mashDist = new float[h_numSequences*(h_numSequences - 1)/2];

    float * mashDist = d_mashDist;

    cudaError_t err;


    err = cudaMemcpy(h_mashDist, mashDist, (h_numSequences*(h_numSequences - 1)/2)*sizeof(uint32_t), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
        exit(1);
    }

    for (int i=0; i<h_numSequences; i++) 
    {
        printf("1.0\t");
        for (int j=i+1; j<h_numSequences; j++) 
        {
            printf("%f\t", *h_mashDist);
            h_mashDist++;
        }
        printf("\n");        
    }
       

}