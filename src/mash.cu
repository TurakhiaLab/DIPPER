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

__global__ void sketchConstruction
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

    // printf("hashList pointer in device%p\n", d_hashList);

    if (tx==0 && bx==0)
    {
        for (size_t i=0; i<d_numSequences; i++)
        {
            uint32_t seqLength = d_seqLengths[i];
            
            // if (i==9)printf("%ld:\t", i);

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
                // if (i==9) printf("(%u, %u)\t",kmer, hash);
            }
            // if (i==9) printf("\n");
            hashList += seqLength-kmerSize+1;
            compressedSeqs += (seqLength+15)/16;

        }
    }
    // printf("hashList pointer in device%p\n", d_hashList);

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

    sketchConstruction<<<params.numBlocks, params.blockSize>>>(d_compressedSeqs, d_aggseqLengths, d_seqLengths, d_numSequences, d_hashList, params.kmerSize);

    uint32_t * hashList = d_hashList;
    for (size_t i = 0; i < d_numSequences; i++)
    {
        thrust::device_ptr<uint32_t> hashPtr(hashList);
        uint32_t numKmers = (h_seqLengths[i] - 15 + 1);   
        thrust::sort(hashPtr, hashPtr + numKmers);
        hashList += numKmers;
    }

    cudaDeviceSynchronize();
}

__device__ float mashDistance //local to each thread
(
    uint32_t * A,
    uint32_t * B,
    uint32_t kmerSize,
    uint32_t sketchSize,
    uint32_t * seqA,
    uint32_t * seqB
){

    //outer loop 
    //union, 0, <sketchSize, +=32, outer loop needs to keep updating the arrays
    //inner loop
    //subUnion, 0, condition min(32,sketchSize-union), ++

    uint32_t offset = 32*blockDim.x, APtrLocal = 32*threadIdx.x, BPtrLocal = offset + 32*threadIdx.x;
    uint32_t APtrGlobal = 0, BPtrGlobal = 0;
    float inter=0, uni=0;
    // for(uint32_t i = 0; i < sketchSize; i++){
    //     if (threadIdx.x == 1 && blockIdx.x == 0) printf("tx: %d, bx: %d, A[%d] is %u, B[%d] is %u\n",threadIdx.x, blockIdx.x, i, A[i], i, B[i]);
    // }
    //printf("tx: %d, bx: %d, A[0] is %u, B[0] is %u\n",threadIdx.x, blockIdx.x, *A, *B);

    //printf("start intersection value is %f, union value is %f, sketchSize %u\n",inter, uni, sketchSize);
    for(uint32_t u = 0; u < sketchSize; u += 32){ //union always increments by 32 in inner loop
        
        //fetch next 32 OR leftover elements from global memory
        for(uint32_t i = 0; i < min(32, sketchSize - u); i++){
            seqA[APtrLocal + i] = A[APtrGlobal + i];
            seqB[BPtrLocal + i] = B[BPtrGlobal + i];
            //if (threadIdx.x == 1 && blockIdx.x == 0)
            //    printf("union: %d, in mashDistance shared mem init: seqA[%d] is %u, seqB[%d] is %u, A[%d] is %u, B[%d] is %u\n", u, APtrLocal+i, seqA[APtrLocal + i], BPtrLocal+i, seqB[BPtrLocal + i], APtrGlobal+i, A[APtrGlobal+i], BPtrGlobal+i, B[BPtrGlobal+i]);
        }


        for(uint32_t su = 0; su < min(32,sketchSize - u); su++){ // subUnion (union of 32 hashes in seq A and B)

            //printf("sequence A: %u and sequence B: %u\n", seqA[APtrLocal], seqB[BPtrGlobal]);
            if (seqA[APtrLocal] == seqB[BPtrLocal]) 
            {   
                inter++; uni++;
                // if ((threadIdx.x == 1 && blockIdx.x == 0) && (APtrLocal == 1 && BPtrLocal == 1))
                //     printf("tx %d bx %d, sequence A: %u and sequence B: %u, inter: %f\n", threadIdx.x, blockIdx.x, seqA[APtrLocal], seqB[BPtrLocal], inter);
                APtrLocal++; BPtrLocal++;
                APtrGlobal++; BPtrGlobal++;
            } 
            else if (seqA[APtrLocal] > seqB[BPtrLocal])
            {
                uni++;
                BPtrLocal++;
                BPtrGlobal++;
            }
            else
            {
                uni++;
                APtrLocal++;
                APtrGlobal++;
            }
        }

        APtrLocal = 32*threadIdx.x, BPtrLocal = offset + 32*threadIdx.x; //reset APtrLocal and BPtrLocal

    }
  
    //printf("intersection value is %f, union value is %f\n",inter, uni);
    float jaccardEstimate = (inter/uni);
    float mashDist = (-1)*(log(2.0*jaccardEstimate/(1.0+jaccardEstimate)))/kmerSize;
    return mashDist;

}


__device__ int* matrixRowComputation
(
    int* dim2info,
    int index1D,
    int d_numSequences
){
    int left = 0;
    int right = d_numSequences - 1;
    int n = d_numSequences;
    int row_index = 0;
    int col_index = 0;
    //int numElementsPrior = 0;

    while(left<=right){
        int mid = left + (right-left)/2;
        int numElementsPrior = mid*n - (mid*(mid+1))/2;
        int numElementsCurrent = (mid+1)*n - ((mid+1)*(mid+2))/2;
        if(numElementsPrior < (index1D + 1)){
            if (numElementsCurrent < (index1D + 1)){
                left = mid + 1;
            }
            else{
                row_index = mid;
                col_index = (index1D - numElementsPrior) + row_index + 1;
                dim2info[0] = row_index;
                dim2info[1] = col_index;
                return dim2info;
            }

        }

        else if (numElementsPrior >= (index1D + 1)){
            right = mid - 1;
        }

    }
}


__global__ void mashDistConstruction //shared memory allocated here is common to a thread block
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

    int bs = blockDim.x;
    int gs = gridDim.x;

    int tid = bs*bx + tx; //global thread ID
    int threads = bs*gs;

    extern __shared__ uint32_t seqA[];
    extern __shared__ uint32_t seqB[];

    uint32_t * hashList = d_hashList;
    uint32_t mashDistCount=0;
    uint32_t totalElements = d_numSequences * (d_numSequences - 1)/2;
    //printf("number of sequences is %d\n",d_numSequences);
    uint32_t* firstSeq = d_hashList;
    uint32_t* secondSeq = d_hashList; 
    int* dim2info = new int[2];
    
        //full distance matrix
        /*for(int i = tid; i < totalElements; i += threads){
            int row_index = i/d_numSequences;
            int col_index = i%d_numSequences;
            firstSeq = d_hashList + sketchSize*row_index;
            secondSeq = d_hashList + sketchSize*col_index;
            // if (threadIdx.x == 1 && blockIdx.x == 0)
            //     printf("totalElements: %d, iteration: %d, first sequence: %u and second sequence is %u, threadIdx is %d, blockIdx is %d, global thread ID is %d\n", totalElements, i, *firstSeq, *secondSeq, tx, bx, tid);
            //if (col_index >= row_index){
                float mashDist = mashDistance(firstSeq, secondSeq ,kmerSize, sketchSize, seqA, seqB); //2D
                d_mashDist[i] = mashDist; //1D
            //}
        }*/

        //upper triangle 
        for(int i = tid; i < totalElements; i += threads){
            //printf("going to compute matrix row and column for i=%d\n",i);
            dim2info = matrixRowComputation(dim2info, i, d_numSequences);
            //printf("computed matrix row and column for i=%d\n",i);
            int row_index = dim2info[0];
            int col_index = dim2info[1];
            //printf("i is %d, row index is %d, col index is %d\n", i, row_index, col_index);
            firstSeq = d_hashList + sketchSize*row_index;
            secondSeq = d_hashList + sketchSize*col_index;
            // if (threadIdx.x == 1 && blockIdx.x == 0)
            //     printf("totalElements: %d, iteration: %d, first sequence: %u and second sequence is %u, threadIdx is %d, blockIdx is %d, global thread ID is %d\n", totalElements, i, *firstSeq, *secondSeq, tx, bx, tid);
            //if (col_index >= row_index){
            float mashDist = mashDistance(firstSeq, secondSeq ,kmerSize, sketchSize, seqA, seqB); //2D
            //printf("distance for iteration i=%d is %f\n",i,mashDist);
            d_mashDist[i] = mashDist; //1D
            //}
        }
}


__global__ void pruneHashList
(
    uint32_t * d_hashList,
    uint32_t * d_seqLengths,
    size_t d_numSequences,
    float * d_mashDist,
    uint32_t kmerSize,
    uint32_t sketchSize
){
    uint32_t * hashList = d_hashList;
    uint32_t hashListIdx1 = 0;

    uint32_t * hashListPruned = d_hashList;
    //uint32_t * hashListPruned = new uint32_t[sketchSize*d_numSequences*sizeof(uint32_t)];

    //uint32_t hashListIdx2 = hashList + d_seqLengths[i] - kmerSize + 1;
    uint32_t hashListIdx2 = 0;
    for(uint32_t i=0;i<d_numSequences;i++)
    {
        for(uint32_t j=0;j<sketchSize;j++) //would work if sketchSize is less than seqLength
        {
            hashListPruned[i*sketchSize+j] = hashList[hashListIdx1];
            //hashList[i*params.sketchSize+j] = hashList[hashListIdx1];
            hashListIdx1++;
        }

        hashListIdx2 += d_seqLengths[i] - kmerSize + 1;
        hashListIdx1 = hashListIdx2;

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

    pruneHashList<<<1, 1>>>(d_hashList, d_seqLengths, d_numSequences, d_mashDist, params.kmerSize, params.sketchSize);
    //mashDistConstruction<<<params.numBlocks, params.blockSize>>>(d_hashList, d_seqLengths, d_numSequences, d_mashDist, params.kmerSize, params.sketchSize);
    mashDistConstruction<<<params.numBlocks, params.blockSize, 2*32*params.blockSize*sizeof(uint32_t)>>>(d_hashList, d_seqLengths, d_numSequences, d_mashDist, params.kmerSize, params.sketchSize);

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
    
    float * h_mashDist = new float[h_numSequences*(h_numSequences - 1)/2]; //NC2

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