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
    for (size_t i =0; i<numSequences; i++) flatStringLength+= (h_seqLengths[i]+31)/32;
    uint64_t * h_flattenCompressSeqs = new uint64_t[flatStringLength];
    flatStringLength=0;
    for (size_t i =0; i<numSequences; i++) 
    {
        uint64_t flatStringLengthLocal = (h_seqLengths[i]+31)/32;
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

    err = cudaMalloc(&d_prefixCompressed, d_numSequences*sizeof(uint64_t));
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

    // generate prefix array    
    thrust::device_ptr<uint64_t> dev_seqLengths(d_seqLengths);
    thrust::device_ptr<uint64_t> dev_prefixComp(d_prefixCompressed);

    thrust::transform(thrust::device,
        dev_seqLengths, dev_seqLengths + d_numSequences, dev_prefixComp, 
        [] __device__ (const uint64_t& x) -> uint64_t { 
            return (x + 31) / 32;
        }
    );

    thrust::exclusive_scan(dev_prefixComp, dev_prefixComp + d_numSequences, dev_prefixComp);

    cudaDeviceSynchronize();
}

void GpuSketch::DeviceArrays::deallocateDeviceArrays(){
    cudaFree(d_compressedSeqs);
    cudaFree(d_aggseqLengths);
    cudaFree(d_seqLengths);
    cudaFree(d_hashList);
    cudaFree(d_prefixCompressed);
    // cudaFree(d_mashDist);
}

#define BIG_CONSTANT(x) (x##LLU)

__device__ inline uint64_t rotl64 ( uint64_t x, int8_t r )
{
    return (x << r) | (x >> (64 - r));
}

#define ROTL64(x,y)    rotl64(x,y)

__device__ uint64_t getblock64 ( const uint64_t * p, int i )
{
    return p[i];
}

__device__ uint64_t fmix64 ( uint64_t k )
{
    k ^= k >> 33;
    k *= BIG_CONSTANT(0xff51afd7ed558ccd);
    k ^= k >> 33;
    k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
    k ^= k >> 33;

    return k;
}

// First hashing function using raw sequence
__device__ void MurmurHash3_x64_128_MASH ( void * key, const int len,
                           const uint32_t seed, void * out)
{
    const uint8_t * data = (const uint8_t*)key;
    const int nblocks = len / 16;

    uint64_t h1 = seed;
    uint64_t h2 = seed;

    const uint64_t c1 = BIG_CONSTANT(0x87c37b91114253d5);
    const uint64_t c2 = BIG_CONSTANT(0x4cf5ad432745937f);

    //----------
    // body

    const uint64_t * blocks = (const uint64_t *)(data);

    for(int i = 0; i < nblocks; i++)
    {
        uint64_t k1 = getblock64(blocks,i*2+0);
        uint64_t k2 = getblock64(blocks,i*2+1);

        k1 *= c1; k1  = ROTL64(k1,31); k1 *= c2; h1 ^= k1;

        h1 = ROTL64(h1,27); h1 += h2; h1 = h1*5+0x52dce729;

        k2 *= c2; k2  = ROTL64(k2,33); k2 *= c1; h2 ^= k2;

        h2 = ROTL64(h2,31); h2 += h1; h2 = h2*5+0x38495ab5;
    }

    //----------
    // tail

    const uint8_t * tail = (const uint8_t*)(data + nblocks*16);

    uint64_t k1 = 0;
    uint64_t k2 = 0;

    switch(len & 15)
    {
        case 15: k2 ^= ((uint64_t)tail[14]) << 48;
        case 14: k2 ^= ((uint64_t)tail[13]) << 40;
        case 13: k2 ^= ((uint64_t)tail[12]) << 32;
        case 12: k2 ^= ((uint64_t)tail[11]) << 24;
        case 11: k2 ^= ((uint64_t)tail[10]) << 16;
        case 10: k2 ^= ((uint64_t)tail[ 9]) << 8;
        case  9: k2 ^= ((uint64_t)tail[ 8]) << 0;
               k2 *= c2; k2  = ROTL64(k2,33); k2 *= c1; h2 ^= k2;

        case  8: k1 ^= ((uint64_t)tail[ 7]) << 56;
        case  7: k1 ^= ((uint64_t)tail[ 6]) << 48;
        case  6: k1 ^= ((uint64_t)tail[ 5]) << 40;
        case  5: k1 ^= ((uint64_t)tail[ 4]) << 32;
        case  4: k1 ^= ((uint64_t)tail[ 3]) << 24;
        case  3: k1 ^= ((uint64_t)tail[ 2]) << 16;
        case  2: k1 ^= ((uint64_t)tail[ 1]) << 8;
        case  1: k1 ^= ((uint64_t)tail[ 0]) << 0;
               k1 *= c1; k1  = ROTL64(k1,31); k1 *= c2; h1 ^= k1;
    };

    //----------
    // finalization

    h1 ^= len; h2 ^= len;

    h1 += h2;
    h2 += h1;

    h1 = fmix64(h1);
    h2 = fmix64(h2);

    h1 += h2;
    h2 += h1;

    ((uint64_t*)out)[0] = h1;
    ((uint64_t*)out)[1] = h2;
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

__device__ void decompress(uint64_t compressedSeq, char * decompressedSeq) {
    static const char lookupTable[4] = {'A', 'C', 'G', 'T'};
    for (int i = 31; i >= 0; i--) {
        uint64_t twoBitVal = (compressedSeq >> (2 * i)) & 0x3;
        decompressedSeq[i] = lookupTable[twoBitVal];
    }
}

__global__ void sketchConstruction(
    uint64_t * d_compressedSeqs,
    uint64_t * d_seqLengths,
    uint64_t * d_prefixCompressed,
    size_t d_numSequences,
    uint64_t * d_hashList,
    uint64_t kmerSize
) {
    extern __shared__ uint64_t stored[];

    typedef cub::BlockRadixSort<uint64_t, 512, 3> BlockRadixSort;

    __shared__ typename BlockRadixSort::TempStorage temp_storage;

    int tx = threadIdx.x;
    int bx = blockIdx.x;
    int stride = gridDim.x;

    uint64_t kmer = 0;
    uint64_t mask = (1 << (2 * kmerSize)) - 1;
    uint8_t out[16];
    char decompressedSeq[32];

    for (size_t i = bx; i < d_numSequences; i += stride) {
        //if (tx == 0) printf("Started block %d\n", i);
        stored[tx] = 0xFFFFFFFFFFFFFFFF; // reset block's stored values
        stored[tx + 500] = 0xFFFFFFFFFFFFFFFF;

        uint64_t  * hashList = d_hashList + 1000*i;

        uint64_t seqLength = d_seqLengths[i];
        uint64_t * compressedSeqs = d_compressedSeqs + d_prefixCompressed[i];

        size_t j = tx;
        size_t iteration = 0;
        while (iteration <= seqLength - kmerSize) {
            uint64_t keys[3];
            if (j <= seqLength - kmerSize) {
                uint64_t index = j/32;
                uint64_t shift1 = 2*(j%32);
                if (shift1>0) {
                    uint64_t shift2 = 64-shift1;
                    kmer = ((compressedSeqs[index] >> shift1) | (compressedSeqs[index+1] << shift2)); //& mask;
                }
                else {   
                    kmer = compressedSeqs[index];// & mask;
                }

                decompress(kmer, decompressedSeq);

                // convert to char representation and call w/ original
                MurmurHash3_x64_128_MASH(decompressedSeq, kmerSize, 42, out);
                uint64_t hash = *((uint64_t *)out);

                // Combine stored and computed to sort and rank
                keys[0] = (tx < 500) ? stored[tx] : 0xFFFFFFFFFFFFFFFF;
                keys[1] = (tx < 500) ? stored[tx + 500] : 0xFFFFFFFFFFFFFFFF;
                keys[2] = hash;
            } else {
                keys[0] = (tx < 500) ? stored[tx] : 0xFFFFFFFFFFFFFFFF;
                keys[1] = (tx < 500) ? stored[tx + 500] : 0xFFFFFFFFFFFFFFFF;
                keys[2] = 0xFFFFFFFFFFFFFFFF;
            }

            BlockRadixSort(temp_storage).Sort(keys);

            // // Move top 1000 hashes back to stored
            if (tx < 333) {
                stored[3*tx] = keys[0];
                stored[3*tx + 1] = keys[1];
                stored[3*tx + 2] = keys[2];
            } else if (tx == 333) {
                stored[999] = keys[0];
            }

            __syncthreads();

            iteration += blockDim.x;
            j += blockDim.x;

        }

        // if (tx == 0) {
        //     printf("i       hashList[i] (%d)\n", i);
        //     for (int j = 0; j < 10; j++) {
        //         printf("%lu       %lu\n", j, stored[j]);
        //     }
        // }
    
        // Result writing back to global memory.
        if (tx < 500) {
            // d_hashList[i*1000 + tx] = stored[tx];
            // d_hashList[i*1000 + tx + 500] = stored[tx + 500];
            hashList[tx] = stored[tx];
            hashList[tx + 500] = stored[tx + 500];
        }

        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            printf("CUDA Error: %s\n", cudaGetErrorString(err));
        }

        // if (tx == 0) {
        //     printf("i       hashList[i] (%d), %ld, %p\n", i, d_numSequences, hashList);
        //     for (int j = 0; j < 20; j++) {
        //         printf("%lu       %lu\n", j, d_hashList[i * 1000 + j]);
        //     }
        // }
       
    }

}


void GpuSketch::sketchConstructionOnGpu
(
    uint64_t * d_compressedSeqs,
    uint64_t * d_prefixCompressed,
    uint64_t * d_aggseqLengths,
    uint64_t * d_seqLengths,
    size_t d_numSequences,
    uint64_t * d_hashList,
    uint64_t * h_seqLengths,
    Param& params
){

    const uint64_t kmerSize = params.kmerSize; // Extract kmerSize

    auto timerStart = std::chrono::high_resolution_clock::now();

    int threadsPerBlock = 512;
    int blocksPerGrid = (d_numSequences + threadsPerBlock - 1) / threadsPerBlock;
    size_t sharedMemorySize = sizeof(uint64_t) * (2000);
    sketchConstruction<<<1, threadsPerBlock, sharedMemorySize>>>(
        d_compressedSeqs, d_seqLengths, d_prefixCompressed, d_numSequences, d_hashList, kmerSize
    );

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("CUDA Error: %s\n", cudaGetErrorString(err));
    }
    cudaDeviceSynchronize();

    auto timerEnd = std::chrono::high_resolution_clock::now();
    auto time = timerEnd - timerStart;
    std::cout << "Time to generate hashes: " << time.count() << "ns\n";

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
    uint64_t * h_hashList = new uint64_t[1000*d_numSequences];


    uint64_t * hashList = d_hashList;

    cudaError_t err;

    //printf("Total Hashes: %d", d_numSequences*1000);

    err = cudaMemcpy(h_hashList, hashList, d_numSequences*1000*sizeof(uint64_t), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
        exit(1);
    }

    // printf("i\thashList[i] (%zu)\n");
    for (int j = 0; j < d_numSequences; j++) {
        printf("Sequence (%d)\n", j);
        for (int i=0; i<numValues; i++) {
            printf("%i\t%lu\n", i, h_hashList[j*1000 + i]);
        }
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