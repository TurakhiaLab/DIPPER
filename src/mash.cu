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
#include <cub/cub.cuh>

/* Note: d_aggseqLengths is the aggregated compressed length of original string (h_seqLengths)
Ex: ["dog", "mouse", "cat"] 
h_seqLengths -> [3, 5, 3] 
d_aggseqLengths -> [1, 2, 3] */
void GpuSketch::DeviceArrays::allocateDeviceArrays(uint64_t ** h_compressedSeqs, uint64_t * h_seqLengths, size_t numSequences, Param& params)
{
    cudaError_t err;

    d_numSequences = numSequences;

    uint64_t kmerSize = params.kmerSize;
    size_t hashListLength = 0;   

    // Allocate memory
    err = cudaMalloc(&d_aggseqLengths, numSequences*sizeof(uint64_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_seqLengths, numSequences*sizeof(uint64_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    /* Flatten data */
    uint64_t * h_aggseqLengths = new uint64_t[numSequences];
    uint64_t flatStringLength=0;
    for (size_t i =0; i<numSequences; i++) flatStringLength+= (h_seqLengths[i]+15)/16;
    uint64_t * h_flattenCompressSeqs = new uint64_t[flatStringLength];
    flatStringLength=0;
    for (size_t i =0; i<numSequences; i++) 
    {
        uint64_t flatStringLengthLocal = (h_seqLengths[i]+15)/16;
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

    // err = cudaMalloc(&d_hashList, hashListLength*sizeof(uint64_t));
    err = cudaMalloc(&d_hashList, 1000*numSequences*sizeof(uint64_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    printf("%p", d_hashList);


    err = cudaMalloc(&d_compressedSeqs, flatStringLength*sizeof(uint64_t));
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
    err = cudaMemcpy(d_aggseqLengths, h_aggseqLengths, numSequences*sizeof(uint64_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) 
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
        exit(1);
    }

    err = cudaMemcpy(d_seqLengths, h_seqLengths, numSequences*sizeof(uint64_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) 
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
        exit(1);
    }

    err = cudaMemcpy(d_compressedSeqs, h_flattenCompressSeqs, flatStringLength*sizeof(uint64_t), cudaMemcpyHostToDevice);
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

__device__ uint64_t MurmurHash3_x64_128_updated(const uint64_t key, const int len, const uint64_t seed /* , void* out */) {
    uint64_t h1 = seed;
    uint64_t h2 = seed;

    const uint64_t c1 = 0x87c37b91114253d5;
    const uint64_t c2 = 0x4cf5ad432745937f;

    // Masks for processing parts of the key based on len
    const uint64_t masks[9] = {
        0x0000000000000000, // 0 (not used)
        0x00000000000000FF, // 1
        0x000000000000FFFF, // 2
        0x0000000000FFFFFF, // 3
        0x00000000FFFFFFFF, // 4
        0x000000FFFFFFFFFF, // 5
        0x0000FFFFFFFFFFFF, // 6
        0x00FFFFFFFFFFFFFF, // 7
        0xFFFFFFFFFFFFFFFF  // 8
    };

    uint64_t k1 = key & (len < 8 ? masks[len] : 0xFFFFFFFFFFFFFFFF);
    uint64_t k2 = len > 8 ? (key >> (64 - (len - 8) * 8)) & masks[len - 8] : 0;

    // Process k1
    k1 *= c1; k1 = (k1 << 31) | (k1 >> 33); k1 *= c2; h1 ^= k1;
    h1 = (h1 << 27) | (h1 >> 37); h1 += h2; h1 = h1 * 5 + 0x52dce729;

    // Process k2
    k2 *= c2; k2 = (k2 << 33) | (k2 >> 31); k2 *= c1; h2 ^= k2;
    h2 = (h2 << 31) | (h2 >> 33); h2 += h1; h2 = h2 * 5 + 0x38495ab5;

    // Finalization
    h1 ^= len; h2 ^= len;

    h1 += h2;
    h2 += h1;

    h1 ^= h1 >> 33;
    h1 *= 0xff51afd7ed558ccd;
    h1 ^= h1 >> 33;
    h1 *= 0xc4ceb9fe1a85ec53;
    h1 ^= h1 >> 33;

    h2 ^= h2 >> 33;
    h2 *= 0xff51afd7ed558ccd;
    h2 ^= h2 >> 33;
    h2 *= 0xc4ceb9fe1a85ec53;
    h2 ^= h2 >> 33;

    h1 += h2;
    h2 += h1;

    // Total 128-bit hash
    // ((uint64_t*)out)[0] = h1;
    // ((uint64_t*)out)[1] = h2;

    // pick one of the 64-bit hashes
    return h2; // which one is picked shouldnt affect in theory
}

/*
Original serial reference
__global__ void sketchConstructionSerial
(
    uint64_t * d_compressedSeqs,
    uint64_t * d_aggseqLengths,
    uint64_t * d_seqLengths,
    size_t d_numSequences,
    uint64_t * d_hashList,
    uint64_t kmerSize
){
    int tx = threadIdx.x;
    int bx = blockIdx.x;


    uint64_t kmer = 0;
    uint64_t mask = (1<<2*kmerSize) - 1;

    uint64_t * hashList = d_hashList;
    uint64_t * compressedSeqs = d_compressedSeqs;

    //printf("hashList pointer in device%p\n", d_hashList);

    if (tx==0 && bx==0)
    {
        for (size_t i=0; i<d_numSequences; i++)
        {
            uint64_t seqLength = d_seqLengths[i];
            
            //if (i==9)printf("%ld:\t", i);

            for (size_t j=0; j<=seqLength-kmerSize; j++)
            {
                uint64_t index = j/16;
                uint64_t shift1 = 2*(j%16);
                if (shift1>0)
                {
                    uint64_t shift2 = 32-shift1;
                    kmer = ((compressedSeqs[index] >> shift1) | (compressedSeqs[index+1] << shift2)) & mask;
                }
                else
                {   
                    kmer = compressedSeqs[index] & mask;
                }
                uint64_t hash = MurmurHash3_x86_32  (kmer, 30, 53);
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
*/


__global__ void sketchConstruction(
    uint64_t * d_compressedSeqs,
    uint64_t * d_seqLengths,
    uint64_t * d_prefixCompressed,
    size_t d_numSequences,
    uint64_t * d_hashList,
    uint64_t kmerSize
) {
    extern __shared__ uint64_t stored[];

    typedef cub::BlockRadixSort<uint64_t, 500, 3> BlockRadixSort;

    __shared__ typename BlockRadixSort::TempStorage temp_storage;

    int tx = threadIdx.x;
    int bx = blockIdx.x;

    for (size_t i = bx; i < d_numSequences; i += gridDim.x) {
        stored[tx] = 0xFFFFFFFFFFFFFFFF;
        stored[tx + 500] = 0xFFFFFFFFFFFFFFFF;

        uint64_t seqLength = d_seqLengths[i];
        uint64_t * compressedSeqs = d_compressedSeqs + d_prefixCompressed[i];

        uint64_t kmer = 0;
        uint64_t mask = (1 << (2 * kmerSize)) - 1;

        for (size_t j = tx; j <= seqLength - kmerSize; j += blockDim.x) {
            uint64_t index = j/16;
            uint64_t shift1 = 2*(j%16);
            if (shift1>0) {
                uint64_t shift2 = 32-shift1;
                kmer = ((compressedSeqs[index] >> shift1) | (compressedSeqs[index+1] << shift2)) & mask;
            }
            else {   
                kmer = compressedSeqs[index] & mask;
            }
            uint64_t hash = MurmurHash3_x64_128_updated(kmer, 30, 53);
            // Combine stored and computed to sort and rank
            uint64_t keys[3];

            keys[0] = (tx < 500) ? stored[tx] : 0xFFFFFFFFFFFFFFFF;
            keys[1] = (tx < 500) ? stored[tx + 500] : 0xFFFFFFFFFFFFFFFF;
            keys[2] = hash;

            __syncthreads();
            BlockRadixSort(temp_storage).Sort(keys);

            // if (i < 3) { 
            //     if (tx == 0) {
            //         printf("Sequence %lu, Length: %lu, First K-mer: %lu, First Hash: %lu\n", i, seqLength, kmer, keys[2]);
            //     }
            // }

            // Move top 1000 hashes back to stored
            __syncthreads();
            if (tx < 333) {
                stored[3*tx] = keys[0];
                stored[3*tx + 1] = keys[1];
                stored[3*tx + 2] = keys[2];
            } else if (tx == 333) {
                stored[999] = keys[0];
            }
            // if (i == 1) {
            //     if (tx < 10) {
            //         printf("Thread %lu, Value %lu\n", tx, stored[tx]);
            //     }
            // }
            __syncthreads();

        }

        // if (tx == 0) {
        //     printf("i       hashList[i] (%d)\n", i);
        //     for (int j = 0; j < 10; j++) {
        //         printf("%lu       %lu\n", j, stored[j]);
        //     }
        // }

        // Result writing back to global memory.
        if (tx < 500) {
            d_hashList[i*1000 + tx] = stored[tx];
            d_hashList[i*1000 + tx + 500] = stored[tx + 500];
        }

         
        if (tx == 0) {
            printf("i       hashList[i] (%d), %ld\n", i, d_numSequences);
            for (int j = 0; j < 20; j++) {
                printf("%lu       %lu\n", j, d_hashList[i * 1000 + j]);
            }
        }
       
    }

}


void GpuSketch::sketchConstructionOnGpu
(
    uint64_t * d_compressedSeqs,
    uint64_t * d_aggseqLengths,
    uint64_t * d_seqLengths,
    size_t d_numSequences,
    uint64_t * d_hashList,
    uint64_t * h_seqLengths,
    Param& params
){

    // prefix-sum of d_seqLengths using thrust
    uint64_t * d_prefixCompressed;

    int bytes = d_numSequences * sizeof(uint64_t);
    cudaMalloc(&d_prefixCompressed, bytes);
    
    thrust::device_ptr<uint64_t> dev_seqLengths(d_seqLengths);
    thrust::device_ptr<uint64_t> dev_prefixComp(d_prefixCompressed);

    auto timerStart = std::chrono::high_resolution_clock::now();

    thrust::transform(thrust::device,
        dev_seqLengths, dev_seqLengths + d_numSequences, dev_prefixComp, 
        [] __device__ (const uint64_t& x) -> uint64_t { 
            return (x + 15) / 16;
        }
    );

    const uint64_t kmerSize = params.kmerSize; // Extract kmerSize
    thrust::exclusive_scan(dev_prefixComp, dev_prefixComp + d_numSequences, dev_prefixComp);
    auto timerEnd = std::chrono::high_resolution_clock::now();

    std::chrono::nanoseconds time = timerEnd - timerStart;
    std::cout << "Time to create prefix array: " << time.count() << "ns\n";

    timerStart = std::chrono::high_resolution_clock::now();

    // New kernel call
    int memSize = 25000;
    int threads = 500; // hard code 500 for structure
    sketchConstruction<<<1 /*params.numBlocks*/, threads, memSize>>>(d_compressedSeqs, d_seqLengths, d_prefixCompressed, d_numSequences, d_hashList, kmerSize);
    cudaDeviceSynchronize();

    timerEnd = std::chrono::high_resolution_clock::now();
    time = timerEnd - timerStart;
    std::cout << "Time to generate hashes: " << time.count() << "ns\n";

    // cudaFree(d_prefixCompressed);

    // bytes = d_numSequences * 1000 * sizeof(uint64_t);
    // // uint64_t * h_pruned = (uint64_t*)malloc(bytes);
    // uint64_t* h_pruned = new uint64_t[d_numSequences * 1000];
    // uint64_t * hashlist= d_hashList;


    // cudaError_t err = cudaMemcpy(h_pruned, hashlist, d_numSequences * 1000 * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    // if(err!= cudaSuccess)
    // {
    //     std::cout<<"Error; " << err << std::endl;
    // }
    // for (int i = 0; i < d_numSequences; i++) {
    //     printf("i       hashList[i] (%d)\n", i);
    //     for (int j = 0; j < 10; j++) {
    //         printf("%lu       %lu\n", j, h_pruned[i * 1000 + j]);

    //     }
    // }
   

    // prevent mem leak for now
    // free(h_pruned);

    
}

__device__ float mashDistance
(
    uint64_t * A,
    uint64_t * B,
    uint64_t kmerSize,
    uint64_t sketchSize
){
    uint64_t unionPtr = 0, APtr = 0, BPtr = 0;
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
    uint64_t * d_hashList,
    uint64_t * d_seqLengths,
    size_t d_numSequences,
    float * d_mashDist,
    uint64_t kmerSize,
    uint64_t sketchSize
){
    int tx = threadIdx.x;
    int bx = blockIdx.x;

    uint64_t * hashList = d_hashList;


    if (tx==0 && bx==0)
    {
        uint64_t mashDistCount=0;
        for (size_t i=0; i<d_numSequences; i++)
        {
            uint64_t * hashListStartIndex = hashList + d_seqLengths[i] - kmerSize + 1;
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
    uint64_t * d_hashList,
    uint64_t * d_seqLengths,
    size_t d_numSequences,
    float * d_mashDist,
    uint64_t * h_seqLengths,
    Param& params
){

    mashDistConstruction<<<params.numBlocks, params.blockSize>>>(d_hashList, d_seqLengths, d_numSequences, d_mashDist, params.kmerSize, params.sketchSize);

    cudaDeviceSynchronize();

}

void GpuSketch::DeviceArrays::printSketchValues(int numValues) 
{
    uint64_t * h_hashList = new uint64_t[1000*30];

    printf("%p", d_hashList);


    uint64_t * hashList = d_hashList;

    cudaError_t err;

    printf("%d", d_numSequences*1000);

    err = cudaMemcpy(h_hashList, hashList, 30*100*sizeof(uint64_t), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
        exit(1);
    }

    // printf("i\thashList[i] (%zu)\n");
    for (int i=0; i<numValues; i++) {
        printf("%i\t%ld\n", i, h_hashList[i]);
    }

}




void GpuSketch::DeviceArrays::printMashDist(uint64_t h_numSequences) 
{
    
    float * h_mashDist = new float[h_numSequences*(h_numSequences - 1)/2];

    float * mashDist = d_mashDist;

    cudaError_t err;


    err = cudaMemcpy(h_mashDist, mashDist, (h_numSequences*(h_numSequences - 1)/2)*sizeof(uint64_t), cudaMemcpyDeviceToHost);
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