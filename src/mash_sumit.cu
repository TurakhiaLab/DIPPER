#ifndef MASH_CUH
#include "mash_sumit.cuh"
#endif

#include <stdio.h>
#include <queue>
#include <thrust/sort.h>
#include <thrust/scan.h>
#include <thrust/binary_search.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <cub/cub.cuh>

/* Note: d_aggseqLengths is the aggregated compressed length of original string (h_seqLengths)
Ex: ["dog", "mouse", "cat"] 
h_seqLengths -> [3, 5, 3] 
d_aggseqLengths -> [1, 2, 3] */
void GpuSketch::DeviceArrays::allocateDeviceArrays(char ** h_seqsIn, uint64_t ** h_seqsLenIn, size_t numSequences, Param& params)
{

    int device;
    cudaGetDevice(&device);

    struct cudaDeviceProp props;
    cudaGetDeviceProperties(&props, device);
    // printf("Shared memory per block (Kbytes) %.1f\n",(float)(props.sharedMemPerBlock)/1024.0);

    cudaError_t err;
    uint64_t kmerSize = params.kmerSize;
    size_t hashListLength = 0;   
    char * h_seqs = *h_seqsIn;
    uint64_t * h_seqsLen = *h_seqsLenIn;

    d_numSequences = numSequences;
    // Allocate memory
    err = cudaMalloc(&d_aggrSeqsLen, numSequences*sizeof(uint64_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed d_aggseqLengths!\n");
        exit(1);
    }
    err = cudaMalloc(&d_seqsLen, numSequences*sizeof(uint64_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed d_seqLengths!\n");
        exit(1);
    }

    uint64_t * aggrSeqsLen = new uint64_t[numSequences];
    uint64_t totalLen = 0;
    aggrSeqsLen[0] = 0;
    for (int i = 1; i < numSequences; i++)
    {
        aggrSeqsLen[i]=aggrSeqsLen[i-1] + h_seqsLen[i-1];
    }
    totalLen = aggrSeqsLen[numSequences-1] + h_seqsLen[numSequences-1];

    std::cout << "Length of sequence transferred to GPU: " << totalLen << std::endl;
    
    err = cudaMalloc(&d_seqs, totalLen*sizeof(char));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed d_compressedSeqs!\n");
        exit(1);
    }


    // Transfer data
    err = cudaMemcpy(d_aggrSeqsLen, aggrSeqsLen, numSequences*sizeof(uint64_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) 
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed d_aggseqLengths!\n");
        exit(1);
    }

    err = cudaMemcpy(d_seqsLen, h_seqsLen, numSequences*sizeof(uint64_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) 
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed d_seqLengths!\n");
        exit(1);
    }

    err = cudaMemcpy(d_seqs, h_seqs, totalLen*sizeof(char), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) 
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed d_compressedSeqs!\n");
        exit(1);
    }

    cudaDeviceSynchronize();
}

void GpuSketch::DeviceArrays::deallocateDeviceArrays(){
    cudaFree(d_seqs);
    cudaFree(d_aggrSeqsLen);
    cudaFree(d_seqsLen);
    cudaFree(d_hashList);
    cudaFree(d_hashListPruned);
    // cudaFree(d_mashDist);
}

#define BIG_CONSTANT(x) (x)
__device__ uint64_t getblock64 ( const uint64_t * p, int i )
{
  return p[i];
}

__device__ uint64_t rotl64 ( uint64_t x, int8_t r )
{
  return (x << r) | (x >> (64 - r));
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

__device__ void MurmurHash3_x64_128 ( const void * key, const int len,
                           const uint32_t seed, void * out )
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

    k1 *= c1; k1  = rotl64(k1,31); k1 *= c2; h1 ^= k1;

    h1 = rotl64(h1,27); h1 += h2; h1 = h1*5+0x52dce729;

    k2 *= c2; k2  = rotl64(k2,33); k2 *= c1; h2 ^= k2;

    h2 = rotl64(h2,31); h2 += h1; h2 = h2*5+0x38495ab5;
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
            k2 *= c2; k2  = rotl64(k2,33); k2 *= c1; h2 ^= k2;

    case  8: k1 ^= ((uint64_t)tail[ 7]) << 56;
    case  7: k1 ^= ((uint64_t)tail[ 6]) << 48;
    case  6: k1 ^= ((uint64_t)tail[ 5]) << 40;
    case  5: k1 ^= ((uint64_t)tail[ 4]) << 32;
    case  4: k1 ^= ((uint64_t)tail[ 3]) << 24;
    case  3: k1 ^= ((uint64_t)tail[ 2]) << 16;
    case  2: k1 ^= ((uint64_t)tail[ 1]) << 8;
    case  1: k1 ^= ((uint64_t)tail[ 0]) << 0;
            k1 *= c1; k1  = rotl64(k1,31); k1 *= c2; h1 ^= k1;
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

template <int BLOCK_THREADS, int ITEMS_PER_THREAD>
__global__ void sketchConstruction
(
    char ** d_seqsIn,
    uint64_t ** d_aggrSeqsLenIn,
    uint64_t ** d_seqsLenIn,
    size_t d_numSequences,
    uint64_t * d_hashList,
    uint64_t kmerSize,
    uint64_t sketchSize,
    uint64_t seed
){

    char * d_seqs = *d_seqsIn;
    uint64_t * d_aggrSeqsLen = *d_aggrSeqsLenIn;
    uint64_t * d_seqsLen = *d_seqsLenIn;

    size_t tx = threadIdx.x;
    size_t bx = blockIdx.x;
    size_t gx = gridDim.x;

    extern __shared__ uint64_t sketches[];
    uint64_t hash;

    typedef cub::BlockRadixSort<uint64_t, BLOCK_THREADS, ITEMS_PER_THREAD> BlockRadixSort;

    __shared__ typename BlockRadixSort::TempStorage temp_storage;
    __shared__ typename BlockRadixSort::TempStorage temp_storage2;

    for (size_t j = bx; j < d_numSequences; j += gx) 
    {
        uint64_t seqLen = d_seqsLen[j];
        char * currentSeq = d_seqs + d_aggrSeqsLen[j];
        
        for (size_t s = tx; s <= seqLen - kmerSize; s += sketchSize)
        {
            uint64_t hashPerThread [ITEMS_PER_THREAD]; 
            int index = 0;  
            for (size_t i = s; i <= s + sketchSize; i += blockDim.x) 
            {
                if (i <= seqLen - kmerSize)
                {
                    char * kmer = currentSeq + i;
                    char hashChar[16];
                    MurmurHash3_x64_128  (kmer, kmerSize, seed, hashChar);
                    // sketches[sketchSize + (i%sketchSize)] = *((uint64_t *)hashChar);
                    hashPerThread[index++] = *((uint64_t *)hashChar);
                }
            }
            __syncthreads();

            BlockRadixSort(temp_storage).Sort(hashPerThread);

            for (size_t i = s; i <= s + sketchSize; i += blockDim.x) 
            {
                if (i <= seqLen - kmerSize)
                {
                    // sketches[sketchSize + (i%sketchSize)] = *((uint64_t *)hashChar);
                    sketches[sketchSize + (i%sketchSize)] = hashPerThread[--index];
                }
            }

            __syncthreads();

            

        } 
    }
    

}


void GpuSketch::sketchConstructionOnGpu
(
    char * d_seqs,
    uint64_t * d_aggrSeqsLen,
    uint64_t * d_seqsLen,
    size_t d_numSequences,
    uint64_t * d_hashList,
    uint64_t * d_hashListPruned,
    Param& params,
    uint64_t seed
){

    cudaError_t err;

    //Constructing Kmer and Hashes for each sequence
    sketchConstruction<128, (1280/128)><<<params.numBlocks, params.blockSize>>>(&d_seqs, &d_aggrSeqsLen, &d_seqsLen, d_numSequences, d_hashList, params.kmerSize, params.sketchSize, seed);
    cudaDeviceSynchronize();

}

