#include "mash_placement.cuh"

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

void MashPlacement::MashDeviceArrays::allocateDeviceArrays(uint64_t ** h_compressedSeqs, uint64_t * h_seqLengths, size_t num, Param& params)
{
    cudaError_t err;
    this->totalNumSequences = num;
    this->backboneSize = params.backboneSize;

    /* Allocate device arrays */
    size_t maxLength = 0;
    for (size_t i = 0; i < num; i++) {
        if (h_seqLengths[i] > maxLength)
            maxLength = h_seqLengths[i];
    }
    size_t maxLengthCompressed = (maxLength + 31) / 32;

    err = cudaMalloc(&d_seqLengths, params.batchSize*sizeof(uint64_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: d_seqLengths cudaMalloc failed!\n");
        exit(1);
    }
    // printf("Allocated for d_seqLengths: %zu %p\n", params.batchSize, d_seqLengths);
    
    err = cudaMalloc(&d_hashList, params.sketchSize*params.batchSize*sizeof(uint64_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: d_hashList cudaMalloc failed!\n");
        exit(1);
    }
    // printf("Allocated for d_hashList: %zu %p\n", params.sketchSize*params.batchSize, d_hashList);

    err = cudaMalloc(&d_compressedSeqs, params.batchSize*maxLengthCompressed*sizeof(uint64_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: d_compressedSeqs cudaMalloc failed!\n");
        exit(1);
    }
    // printf("Allocated for d_compressedSeqs: %zu %p\n", params.batchSize*maxLengthCompressed, d_compressedSeqs);

    err = cudaMalloc(&d_hashListConst, params.sketchSize*params.backboneSize*sizeof(uint64_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: d_hashListConst cudaMalloc failed!\n");
        exit(1);
    }
    std::cerr << "Allocated for d_hashListConst: "<< d_hashListConst << std::endl;

    err = cudaMalloc(&d_hashListBackbone, params.sketchSize*params.backboneSize*sizeof(uint64_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: d_hashListBackbone cudaMalloc failed!\n");
        exit(1);
    }

    err = cudaMalloc(&d_leafID_map, totalNumSequences*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    cudaDeviceSynchronize();
}

void MashPlacement::MashDeviceArrays::deallocateDeviceArrays(){
    cudaFree(d_compressedSeqs);
    cudaFree(d_seqLengths);
    cudaFree(d_hashList);
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


__device__ void decompress(uint64_t compressedSeq, uint64_t kmerSize, char * decompressedSeq_fwd, char * decompressedSeq_rev) {
    static const char lookupTable[4] = {'A', 'C', 'G', 'T'};
    for (int i = kmerSize - 1; i >= 0; i--) {
        uint64_t twoBitVal = (compressedSeq >> (2 * i)) & 0x3;
        decompressedSeq_fwd[i] = lookupTable[twoBitVal];
        decompressedSeq_rev[kmerSize - 1 - i] = lookupTable[3 - twoBitVal];
    }
}

__device__ int memcmp_device(const char* kmer_fwd, const char* kmer_rev, int kmerSize) {
    for (int i = 0; i < kmerSize; i++) {
        if (kmer_fwd[i] < kmer_rev[i]) {
            return -1;
        }
        if (kmer_fwd[i] > kmer_rev[i]) {
            return 1;
        }
    }
    return 0;
}

__global__ void sketchConstruction(
    uint64_t * d_compressedSeqs,
    uint64_t * d_seqLengths,
    size_t maxLengthCompressed,
    size_t numSequences,
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
    char kmer_fwd[32];
    char kmer_rev[32];
    uint8_t out[16];
    
    for (size_t i = bx; i < numSequences; i += stride) {
        // if (tx == 0) printf("Started block %d\n", i);
        stored[tx] = 0xFFFFFFFFFFFFFFFF; // reset block's stored values
        stored[tx + 500] = 0xFFFFFFFFFFFFFFFF;

        // printf("hashList Pointer: %p\n", d_hashList);
        // printf("Sequence length: %lu\n", d_seqLengths[i]);
        // printf("compressedSeqs Pointer: %p\n", d_compressedSeqs);
        // printf("maxLengthCompressed: %lu\n", maxLengthCompressed);
        
        uint64_t  * hashList = d_hashList + 1000*i;
        uint64_t seqLength = d_seqLengths[i];
        uint64_t * compressedSeqs = d_compressedSeqs + (uint64_t)maxLengthCompressed*i;
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

                decompress(kmer, kmerSize, kmer_fwd, kmer_rev);

                // if ((i == 0) && (tx == 0)) {
                //     for (char c : kmer_fwd) printf("%c", c);   
                
                //     printf("\n");
                // }

                // convert to char representation and call w/ original
                MurmurHash3_x64_128_MASH( (memcmp_device(kmer_fwd, kmer_rev, kmerSize) <= 0) 
                    ? kmer_fwd : kmer_rev, kmerSize, 42, out);
                
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
  
            __syncthreads();
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
    
        // Result writing back to global memory.
        if (tx < 500) {
            // d_hashList[i*1000 + tx] = stored[tx];
            // d_hashList[i*1000 + tx + 500] = stored[tx + 500];
            hashList[tx] = stored[tx];
            hashList[tx + 500] = stored[tx + 500];
        }

        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            // printf("hashList Pointer: %p\n", d_hashList);
            // printf("Sequence length: %lu\n", d_seqLengths[i]);
            // printf("compressedSeqs Pointer: %p\n", d_compressedSeqs);
            // printf("maxLengthCompressed: %lu\n", maxLengthCompressed);
            printf("CUDA Error: %s\n", cudaGetErrorString(err));
        }
       
    }

}

__global__ void rearrangeHashList(
    int numSequences,
    int sketchSize,
    uint64_t * original,
    uint64_t * target
){
    int tx = threadIdx.x, bx = blockIdx.x;
    int bs = blockDim.x;
    int idx = tx+bs*bx;
    if(idx>=numSequences) return;
    for(int i=0;i<sketchSize;i++){
        target[i*numSequences+idx] = original[idx*sketchSize + i];
    }
}

void MashPlacement::MashDeviceArrays::sketchConstructionOnGpu(Param& params, uint64_t** h_compressedSeqs, uint64_t * seqLengths, uint64_t numSequences){
    
    cudaError_t err;

    h_hashList = new uint64_t[params.sketchSize*numSequences];

    size_t maxLength = 0;
    for (size_t i = 0; i < numSequences; i++) {
        if (seqLengths[i] > maxLength)
            maxLength = seqLengths[i];
    }
    size_t maxLengthCompressed = (maxLength + 31) / 32;
    const uint64_t kmerSize = params.kmerSize; // Extract kmerSize
    auto timerStart = std::chrono::high_resolution_clock::now();

    uint64_t localBatchSize = params.batchSize;
    for (int i=0; i < numSequences; i+= localBatchSize) {
        std::cerr << "Processing batch " << i <<std::endl;
        uint64_t * h_seqLengths = new uint64_t[localBatchSize];
        uint64_t * h_flattenCompressSeqs = new uint64_t[localBatchSize*maxLengthCompressed];
        if (i+localBatchSize > numSequences) {
            localBatchSize = numSequences - i;
        }
        for (auto j=i; j<i+localBatchSize && j<numSequences; j++) {
            for (size_t k=0; k<(seqLengths[j]+31)/32;k++)  
            {
                h_flattenCompressSeqs[(j-i)*maxLengthCompressed+k] = h_compressedSeqs[j][k];
                
            }
            h_seqLengths[j-i] = seqLengths[j];  
        }

        err = cudaMemcpy(d_compressedSeqs, h_flattenCompressSeqs, localBatchSize*maxLengthCompressed*sizeof(uint64_t), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) 
        {
            fprintf(stderr, "Gpu_ERROR: h_flattenCompressSeqs cudaMemCpy failed!\n");
            exit(1);
        }
        // std::cerr << "Copied to d_compressedSeqs: " << localBatchSize*maxLengthCompressed << " " << d_compressedSeqs << std::endl;
        // printf("Copied to d_compressedSeqs: %zu %p\n", localBatchSize*maxLengthCompressed, d_compressedSeqs);

        err = cudaMemcpy(d_seqLengths, h_seqLengths, localBatchSize*sizeof(uint64_t), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) 
        {
            fprintf(stderr, "Gpu_ERROR: h_seqLengths cudaMemCpy failed!\n");
            exit(1);
        }
        // std::cerr << "Copied to d_seqLengths: " << localBatchSize << " " << d_seqLengths << std::endl;
        // printf("Copied to d_seqLengths: %zu %p\n", localBatchSize, d_seqLengths);
    
        int threadsPerBlock = 512;
        int blocksPerGrid = 1024;
        size_t sharedMemorySize = sizeof(uint64_t) * (2000);
        sketchConstruction<<<blocksPerGrid, threadsPerBlock, sharedMemorySize>>>(
            d_compressedSeqs, d_seqLengths, maxLengthCompressed, localBatchSize, d_hashList, kmerSize
        );

        cudaDeviceSynchronize();
        cudaMemcpy(h_hashList,d_hashList,params.sketchSize*localBatchSize*sizeof(uint64_t),cudaMemcpyDeviceToHost);
        if (err != cudaSuccess){
            fprintf(stderr, "Gpu_ERROR: h_hashList cudaMalloc failed!\n");
            exit(1);
        }
        h_hashList += params.sketchSize*localBatchSize;
    }

    h_hashList -= params.sketchSize*numSequences;

    /* Rearrange only for backbone tree */
    uint64_t * temp_hashList;
    err = cudaMalloc(&temp_hashList, params.sketchSize*params.backboneSize*sizeof(uint64_t));
    if (err != cudaSuccess){
        fprintf(stderr, "Gpu_ERROR: temp_hashList cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMemcpy(d_hashListBackbone, h_hashList, params.sketchSize*params.backboneSize*sizeof(uint64_t),cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: d_hashListBackbone cudaMemcpy failed!\n");
        exit(1);
    }
    int threadsPerBlock = 512;
    int blocksPerGrid = 1024;
    rearrangeHashList <<<blocksPerGrid, threadsPerBlock >>>(
        params.backboneSize,
        int(params.sketchSize),
        d_hashListBackbone,
        temp_hashList
    );
    std::swap(d_hashListBackbone, temp_hashList);
    cudaFree(temp_hashList);
    
    err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("CUDA Error: %s\n", cudaGetErrorString(err));
    }
    cudaDeviceSynchronize();

    auto timerEnd = std::chrono::high_resolution_clock::now();
    auto time = timerEnd - timerStart;

    // // printf("i\thashList[i] (%zu)\n");
    // for (int j = 0; j < numSequences; j++) {
    //     fprintf(stderr, "Sequence (%d)\n", j);
    //     for (int i=950; i<960; i++) {
    //         fprintf(stderr, "%i\t%lu\n", i, h_hashList[j*params.sketchSize+i]);
    //     }
    // }
}

__global__ void mashDistConstruction(
    int rowId,
    uint64_t * d_hashList,
    double * d_mashDist,
    uint64_t kmerSize,
    uint64_t sketchSize,
    int numSequences
) {
    int tx = threadIdx.x, bx = blockIdx.x, bs = blockDim.x;
    int idx = tx+bx*bs;
    if(idx>=rowId) return;
    int uni = 0, bPos = rowId, inter = 0;
    uint64_t aval, bval;
    for(int i=idx; uni < sketchSize; i+=numSequences, uni++){
        aval = d_hashList[i];
        while(uni < sketchSize && bPos < numSequences * sketchSize){
            bval = d_hashList[bPos];
            // printf("%ull %ull\n",aval,bval);
            if(bval > aval) break;
            if(bval < aval) uni++;
            else inter++;
            bPos += numSequences;
        }
        if(uni >= sketchSize) break;
    }
    double jaccardEstimate = max(double(inter),1.0)/uni;
    d_mashDist[idx] = min(1.0, abs(log(2.0*jaccardEstimate/(1.0+jaccardEstimate))/kmerSize));
}

__global__ void mashDistConstructionRangeForClustering(
    int rowId,
    uint64_t * d_hashListBackbone,
    uint64_t * d_hashListConst,
    double * d_mashDist,
    uint64_t kmerSize,
    uint64_t sketchSize,
    int numSequences,
    int st,
    int ed
) {
    int tx = threadIdx.x, bx = blockIdx.x, bs = blockDim.x;
    int idx = tx+bx*bs;
    if(idx>ed-st) return;
    idx += st;
    int uni = 0, bPos = rowId*sketchSize, inter = 0;
    uint64_t aval, bval;
    for(int i=idx; uni < sketchSize; i+=numSequences, uni++){
        aval = d_hashListBackbone[i];
        while(uni < sketchSize && bPos < rowId*sketchSize + sketchSize){
            bval = d_hashListConst[bPos];
            if(bval > aval) break;
            if(bval < aval) uni++;
            else inter++;
            bPos += 1;
        }
        if(uni >= sketchSize) break;
    }
    double jaccardEstimate = max(double(inter),1.0)/uni;
    d_mashDist[idx] = min(1.0, abs(log(2.0*jaccardEstimate/(1.0+jaccardEstimate))/kmerSize));
}


__global__ void mashDistConstructionRange(
    int rowId,
    uint64_t * d_hashList,
    double * d_mashDist,
    uint64_t kmerSize,
    uint64_t sketchSize,
    int numSequences,
    int st,
    int ed
) {
    int tx = threadIdx.x, bx = blockIdx.x, bs = blockDim.x;
    int idx = tx+bx*bs;
    if(idx>ed-st) return;
    idx += st;
    int uni = 0, bPos = rowId, inter = 0;
    uint64_t aval, bval;
    for(int i=idx; uni < sketchSize; i+=numSequences, uni++){
        aval = d_hashList[i];
        while(uni < sketchSize && bPos < numSequences * sketchSize){
            bval = d_hashList[bPos];
            // printf("%ull %ull\n",aval,bval);
            if(bval > aval) break;
            if(bval < aval) uni++;
            else inter++;
            bPos += numSequences;
        }
        if(uni >= sketchSize) break;
    }
    double jaccardEstimate = max(double(inter),1.0)/uni;
    d_mashDist[idx] = min(1.0, abs(log(2.0*jaccardEstimate/(1.0+jaccardEstimate))/kmerSize));
}


__global__ void mashDistConstructionSpecialID(
    int rowId,
    uint64_t * d_hashListBackbone,
    uint64_t * d_hashListConst,
    double * d_mashDist,
    uint64_t kmerSize,
    uint64_t sketchSize,
    int backboneSize,
    int numToConstruct,
    int * d_id,
    int * d_leafMap
) {
    int tx = threadIdx.x, bx = blockIdx.x, bs = blockDim.x;
    int idx = tx+bx*bs;
    if(idx>=numToConstruct) return;
    
    // if (idx == 0) {
    //     printf("d_leafMask now: ");
    //     for (int z=0;z<numToConstruct;z++){
    //         printf("%d\t", d_id[z]);
    //     }
    //     printf("\n");
    //     printf("leafMap now: ");
    //     for (int z=0;z<numToConstruct;z++){
    //         printf("%d\t", d_leafMap[z]);
    //     }
    //     printf("\n");
    // }

    int mapIdx = d_leafMap[idx];
    idx = d_id[idx];
    if(idx==-1) return;
    int uni = 0, bPos = rowId, inter = 0;
    uint64_t aval, bval;
    int idx_const = idx;
    if (idx > backboneSize) idx_const = mapIdx;
    // if (idx == 618) {
    //     for (int i=idx_const; i<sketchSize*backboneSize; i+=backboneSize) {
    //         printf("%d %llu\n", i, d_hashListConst[i]);
    //     }
    //     for (int i=rowId; i<sketchSize*backboneSize; i+=backboneSize) {
    //         printf("%d %llu\n", i, d_hashListConst[i]);
    //     }
    // }
    for(int i=idx_const; uni < sketchSize; i+=backboneSize, uni++){
        if (idx > backboneSize) aval = d_hashListConst[i];
        else aval = d_hashListBackbone[i];
        while(uni < sketchSize && bPos < backboneSize * sketchSize){
            bval = d_hashListConst[bPos];
            if(bval > aval) break;
            if(bval < aval) uni++;
            else inter++;
            bPos += backboneSize;
        }
        if(uni >= sketchSize) break;
    }
    double jaccardEstimate = max(double(inter),1.0)/uni;
    d_mashDist[idx] = min(1.0, abs(log(2.0*jaccardEstimate/(1.0+jaccardEstimate))/kmerSize));
    // printf("idx: %d new idx: %d mapIdx: %d, rowID %d d_mashDist %lf backbone %d idx_const %d\n", idx, d_id[idx], mapIdx, rowId, d_mashDist[idx], backboneSize, idx_const);
}

__global__ void mashDistConstructionSpecialIDClustering(
    int rowId,
    uint64_t * d_hashListBackbone,
    uint64_t * d_hashListConst,
    double * d_mashDist,
    uint64_t kmerSize,
    uint64_t sketchSize,
    int backboneSize,
    int numToConstruct,
    int * d_id,
    int * d_leafID_map
) {
    int tx = threadIdx.x, bx = blockIdx.x, bs = blockDim.x;
    int idx = tx+bx*bs;
    if(idx>=numToConstruct) return;
    idx = d_id[idx];
    if(idx==-1) return;
    int uni = 0, bPos = 0, inter = 0;
    uint64_t aval, bval;
    int idx_const = idx;

    for (int i=idx_const*sketchSize; i<idx_const*sketchSize+sketchSize; i++) {
        if (idx > backboneSize) aval = d_hashListConst[i];
        else                    aval = d_hashListBackbone[i];
        
        while (bPos < sketchSize && uni < sketchSize) {
            bval = d_hashListConst[d_leafID_map[rowId]*sketchSize + bPos];
            if (bval > aval) break;
            if (bval < aval) uni++;
            else inter++;
            bPos += 1;
        }
        if(uni >= sketchSize) break;
    }
    
    double jaccardEstimate = max(double(inter),1.0)/uni;
    d_mashDist[idx] = min(1.0, abs(log(2.0*jaccardEstimate/(1.0+jaccardEstimate))/kmerSize));
}



void MashPlacement::MashDeviceArrays::distConstructionOnGpuForBackbone(Param& params, int rowId, double* d_mashDist) const{
    int threadNum = 256, blockNum = (this->backboneSize+threadNum-1)/threadNum;
    mashDistConstruction <<<blockNum, threadNum>>> (
        rowId, 
        this->d_hashListBackbone, 
        d_mashDist, 
        params.kmerSize, 
        params.sketchSize, 
        this->backboneSize
    );
}

void MashPlacement::MashDeviceArrays::distConstructionOnGpu(Param& params, int rowId, double* d_mashDist) const{
    int threadNum = 256, blockNum = (this->backboneSize+threadNum-1)/threadNum;
    mashDistConstruction <<<blockNum, threadNum>>> (
        rowId, 
        d_hashList, 
        d_mashDist, 
        params.kmerSize, 
        params.sketchSize, 
        this->backboneSize
    );
}


void MashPlacement::MashDeviceArrays::distRangeConstructionOnGpu(Param& params, int rowId, double* d_mashDist, int l, int r, bool clustering) const{
    
    int threadNum = 1024, blockNum = (r-l+1+threadNum-1)/threadNum;
    if (!clustering) { 
        mashDistConstructionRange <<<blockNum, threadNum>>> (
            rowId, 
            d_hashListBackbone, 
            d_mashDist, 
            params.kmerSize, 
            params.sketchSize, 
            this->backboneSize,
            l,
            r
        );
    } else {
        mashDistConstructionRangeForClustering <<<blockNum, threadNum>>> (
            rowId, 
            d_hashListBackbone, 
            d_hashListConst,
            d_mashDist, 
            params.kmerSize, 
            params.sketchSize, 
            this->backboneSize,
            l,
            r
        );
    }
}



void 
MashPlacement::MashDeviceArrays::distSpecialIDConstructionOnGpu(
        Param& params, 
        int rowId, 
        double* d_mashDist, 
        int numToConstruct, 
        int * d_id,
        int * d_leafMap
    ) const {
    
    int threadNum = 256, blockNum = (numToConstruct+threadNum-1)/threadNum;
    mashDistConstructionSpecialID <<<blockNum, threadNum>>> (
        rowId, 
        d_hashListBackbone,
        d_hashListConst, 
        d_mashDist, 
        params.kmerSize, 
        params.sketchSize, 
        backboneSize,
        numToConstruct,
        d_id,
        d_leafMap
    );
    // cudaDeviceSynchronize();
}


void MashPlacement::MashDeviceArrays::printSketchValues(int numValues) 
{
    uint64_t * h_hashList = new uint64_t[1000*totalNumSequences];


    uint64_t * hashList = d_hashList;

    cudaError_t err;

    //printf("Total Hashes: %d", numSequences*1000);

    err = cudaMemcpy(h_hashList, hashList, totalNumSequences*1000*sizeof(uint64_t), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
        exit(1);
    }

    // printf("i\thashList[i] (%zu)\n");
    for (int j = 0; j < totalNumSequences; j++) {
        fprintf(stderr, "Sequence (%d)\n", j);
        for (int i=0; i<numValues; i++) {
            fprintf(stderr, "%i\t%lu\n", i, h_hashList[i*totalNumSequences+j]);
        }
    }
    

}