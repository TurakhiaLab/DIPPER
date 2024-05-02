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
void GpuSketch::DeviceArrays::allocateDeviceArrays(uint64_t ** h_compressedSeqs, uint64_t * h_seqLengths, size_t numSequences, Param& params)
{

    int device;
    cudaGetDevice(&device);

    struct cudaDeviceProp props;
    cudaGetDeviceProperties(&props, device);
    // printf("Shared memory per block (Kbytes) %.1f\n",(float)(props.sharedMemPerBlock)/1024.0);

    cudaError_t err;

    d_numSequences = numSequences;

    uint64_t kmerSize = params.kmerSize;
    size_t hashListLength = 0;   

    // Allocate memory
    err = cudaMalloc(&d_aggseqLengths, numSequences*sizeof(uint64_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed d_aggseqLengths!\n");
        exit(1);
    }
    err = cudaMalloc(&d_seqLengths, numSequences*sizeof(uint64_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed d_seqLengths!\n");
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
    
    err = cudaMalloc(&d_hashList, hashListLength*sizeof(uint64_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed d_hashList!\n");
        exit(1);
    }

    err = cudaMalloc(&d_hashListPruned, hashListLength*sizeof(uint64_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed d_hashListPruned!\n");
        exit(1);
    }

    err = cudaMalloc(&d_compressedSeqs, flatStringLength*sizeof(uint64_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed d_compressedSeqs!\n");
        exit(1);
    }

    err = cudaMalloc(&d_mashDist, (numSequences*(numSequences -1)/2)*sizeof(float));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed d_mashDist!\n");
        exit(1);
    }

    // Transfer data
    err = cudaMemcpy(d_aggseqLengths, h_aggseqLengths, numSequences*sizeof(uint64_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) 
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed d_aggseqLengths!\n");
        exit(1);
    }

    err = cudaMemcpy(d_seqLengths, h_seqLengths, numSequences*sizeof(uint64_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) 
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed d_seqLengths!\n");
        exit(1);
    }

    err = cudaMemcpy(d_compressedSeqs, h_flattenCompressSeqs, flatStringLength*sizeof(uint64_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) 
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed d_compressedSeqs!\n");
        exit(1);
    }

    cudaDeviceSynchronize();
}

void GpuSketch::DeviceArrays::deallocateDeviceArrays(){
    cudaFree(d_compressedSeqs);
    cudaFree(d_aggseqLengths);
    cudaFree(d_seqLengths);
    cudaFree(d_hashList);
    cudaFree(d_hashListPruned);
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

__device__ void MurmurHash3_x64_128_updated(const uint64_t key, const int len, const uint32_t seed, void* out) {
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

    ((uint64_t*)out)[0] = h1;
    ((uint64_t*)out)[1] = h2;
}

__global__ void sketchConstruction
(
    uint64_t * d_compressedSeqs,
    uint64_t * d_aggseqLengths,
    uint64_t * d_seqLengths,
    uint64_t * d_prefixHashlist,
    uint64_t * d_prefixCompressed,
    size_t d_numSequences,
    uint64_t * d_hashList,
    uint64_t kmerSize,
    uint64_t seed
){
    size_t tx = threadIdx.x;
    size_t bx = blockIdx.x;

    uint64_t kmer = 0;
    uint64_t mask = (1<<2*kmerSize) - 1;

    //printf("hashList pointer in device%p\n", d_hashList);

    if (bx >= d_numSequences) return;

    char hashChar[16];
    uint64_t hash;
    for (size_t j = bx; j < d_numSequences; j += gridDim.x) 
    {
        uint64_t seqLength = d_seqLengths[j];
        uint64_t * hashList = d_hashList + d_prefixHashlist[j];
        uint64_t * compressedSeqs = d_compressedSeqs + d_prefixCompressed[j];
        
        for (size_t i = tx; i <= seqLength - kmerSize; i += blockDim.x) 
        {
            
            uint64_t index = i/32;
            uint64_t shift1 = 2*(i%32);
            if (shift1>0)
            {
                uint64_t shift2 = 64-shift1;
                kmer = ((compressedSeqs[index] >> shift1) | (compressedSeqs[index+1] << shift2)) & mask;
            }
            else
            {   
                kmer = compressedSeqs[index] & mask;
            }
            MurmurHash3_x64_128_updated  (kmer, (kmerSize+3/4), seed, hashChar);
            hash = *((uint64_t *)hashChar);
            // hash = MurmurHash3_x86_32  (kmer, 30, seed);
            //if (j==0) printf("(%u, %u)\n",kmer, hash);
            hashList[i] = hash;
        }

    }
    

}

__global__ void pruneHashList
(
    uint64_t * d_hashList,
    uint64_t * d_hashListPruned,
    uint64_t * d_prefixHashlist,
    uint64_t * d_seqLengths,
    size_t d_numSequences,
    uint64_t kmerSize,
    uint32_t sketchSize
){

    size_t tx = threadIdx.x;
    size_t bx = blockIdx.x;

    uint64_t * hashList = d_hashList;


    for(uint64_t i=bx;i<d_numSequences;i+=gridDim.x)
    {
        hashList = d_hashList + d_prefixHashlist[i];
        //__syncthreads();
        if(d_seqLengths[i] < sketchSize){
            printf("Error: Not enough hashes to build %u size union sketch\n", sketchSize);
        }
        for(uint64_t j=tx;j<sketchSize;j+=blockDim.x) //would work if sketchSize is less than seqLength
        {
            d_hashListPruned[i*sketchSize+j] = hashList[j];
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
    uint64_t * d_hashListPruned,
    uint64_t * h_seqLengths,
    Param& params,
    uint64_t seed
){

    cudaError_t err;

    // prefix-sum of d_seqLengths using thrust
    uint64_t * d_prefixHashlist;
    uint64_t * d_prefixCompressed;

    int bytes = d_numSequences * sizeof(uint64_t);

    err = cudaMalloc(&d_prefixHashlist, bytes);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    err = cudaMalloc(&d_prefixCompressed, bytes);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    thrust::device_ptr<uint64_t> dev_seqLengths(d_seqLengths);
    thrust::device_ptr<uint64_t> dev_prefixHash(d_prefixHashlist);
    thrust::device_ptr<uint64_t> dev_prefixComp(d_prefixCompressed);

    thrust::transform(thrust::device,
        dev_seqLengths, dev_seqLengths + d_numSequences, dev_prefixComp, 
        [] __device__ (const uint64_t& x) -> uint64_t { 
            return (x + 31) / 32;
        }
    );

    const uint64_t kmerSize = params.kmerSize; // Extract kmerSize

    thrust::transform(
        thrust::device,
        dev_seqLengths, dev_seqLengths + d_numSequences, dev_prefixHash,
        [kmerSize] __device__ (const uint64_t& x) -> uint64_t {
            return x - kmerSize + 1;
        }
    );

    //Find out start index of sequences
    thrust::exclusive_scan(dev_prefixHash, dev_prefixHash + d_numSequences, dev_prefixHash);
    thrust::exclusive_scan(dev_prefixComp, dev_prefixComp + d_numSequences, dev_prefixComp);

    //Constructing Kmer and Hashes for each sequence
    sketchConstruction<<<params.numBlocks, params.blockSize>>>(d_compressedSeqs, d_aggseqLengths, d_seqLengths, d_prefixHashlist, d_prefixCompressed, d_numSequences, d_hashList, kmerSize, seed);
    cudaDeviceSynchronize();

    uint64_t * hashList = d_hashList;

    //Sorting the hashList in lexicographical manner
    for (size_t i = 0; i < d_numSequences; i++) 
    {
        thrust::device_ptr<uint64_t> hashPtr(hashList);
        uint64_t numKmers = (h_seqLengths[i] - kmerSize + 1);   
        thrust::sort(hashPtr, hashPtr + numKmers);
        hashList += numKmers;
    }

    //Prune the hashList to create sketch of sketchSize
    pruneHashList<<<params.numBlocks, params.blockSize>>>(d_hashList, d_hashListPruned, d_prefixHashlist, d_seqLengths, d_numSequences, params.kmerSize, params.sketchSize);
    cudaDeviceSynchronize();

    cudaFree(d_prefixHashlist);
    cudaFree(d_prefixCompressed);

}

__device__ float mashDistance //local to each thread
(
    uint64_t * A,
    uint64_t * B,
    uint64_t kmerSize,
    uint32_t sketchSize,
    uint64_t * seqA,
    uint64_t * seqB
){



    //Shared array used by each thread in a block is 33 elements wide to reduce bank conflicts
    uint64_t offset = 33*blockDim.x, APtrLocal = 33*threadIdx.x, BPtrLocal = offset + 33*threadIdx.x;
    uint64_t APtrGlobal = 0, BPtrGlobal = 0;
    float inter=0, uni=0;
   

    printf("start intersection value is %f, union value is %f, sketchSize %u\n",inter, uni, sketchSize);

    //outer loop on union
    for(uint64_t u = 0; u < sketchSize; u += 32){ //union always increments by 32 for every 32 elements
        
        //fetch next 32 OR leftover elements from global memory
        for(uint64_t i = 0; i < min((uint64_t)32, sketchSize - u); i++){
            seqA[APtrLocal + i] = A[APtrGlobal + i];
            seqB[BPtrLocal + i] = B[BPtrGlobal + i];
        }

        //inner loop
        for(uint64_t su = 0; su < min((uint64_t)32,sketchSize - u); su++){ // subUnion (union of 32 hashes in seq A and B)

            //printf("sequence A: %u and sequence B: %u\n", seqA[APtrLocal], seqB[BPtrGlobal]);
            if (seqA[APtrLocal] == seqB[BPtrLocal]) 
            {   
                inter++; uni++;
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

        APtrLocal = 33*threadIdx.x, BPtrLocal = offset + 33*threadIdx.x; //reset APtrLocal and BPtrLocal

    }
  
    float jaccardEstimate = (inter/uni);

    float mashDist = abs((log(2.0*jaccardEstimate/(1.0+jaccardEstimate)))/kmerSize);

    printf("Printing values/////:::: %f\t%f\t%f\n", mashDist, inter, uni);

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
    uint64_t * d_hashList,
    uint64_t * d_seqLengths,
    size_t d_numSequences,
    float * d_mashDist,
    uint64_t kmerSize,
    uint32_t sketchSize
){
    int tx = threadIdx.x;
    int bx = blockIdx.x;

    int bs = blockDim.x;
    int gs = gridDim.x;

    int tid = bs*bx + tx; //global thread ID
    int threads = bs*gs;

    extern __shared__ uint64_t seqA[]; //shared memory for storing first sequence for each thread
    extern __shared__ uint64_t seqB[]; //shared memory for storing second sequence for each thread

    uint64_t totalElements = (d_numSequences * (d_numSequences - 1))/2; //total elements in the distance matrix are NC2
    uint64_t* firstSeq = d_hashList;
    uint64_t* secondSeq = d_hashList; 
    int* dim2info = new int[2];

        //computing elements in the upper triangle
        for(int i = tid; i < totalElements; i += threads){
            dim2info = matrixRowComputation(dim2info, i, d_numSequences); //compute the row and column index in the 2D matrix for the sequences for which we need to calculate the distance as per their placement in the flattened distance matrix
            int row_index = dim2info[0];
            int col_index = dim2info[1];
            firstSeq = d_hashList + sketchSize*row_index; //index of first sequence in the hashList
            secondSeq = d_hashList + sketchSize*col_index; //index of second sequence in the hashList
            float mashDist = mashDistance(firstSeq, secondSeq ,kmerSize, sketchSize, seqA, seqB); //Calculating mash distance for the pair of sequences
            d_mashDist[i] = mashDist; //Filling the 1D distance matrix
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

    mashDistConstruction<<<params.numBlocks, params.blockSize, 2*33*params.blockSize*sizeof(uint32_t)>>>(d_hashList, d_seqLengths, d_numSequences, d_mashDist, params.kmerSize, params.sketchSize);

    cudaDeviceSynchronize();

}

void GpuSketch::DeviceArrays::printSketchValues(uint64_t numValues) 
{
    uint64_t * h_hashList = new uint64_t[numValues];

    uint64_t * hashList = d_hashList;

    cudaError_t err;

    err = cudaMemcpy(h_hashList, hashList, numValues*d_numSequences*sizeof(uint64_t), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
        exit(1);
    }

    for (size_t j = 0; j < d_numSequences; j++)
    {

        printf("i\thashList[i] (%zu)\n",j);
        for (int i=0; i<numValues; i++) {
            printf("%i\t%ld\n", i, h_hashList[i]);
        }
        h_hashList += numValues;
       
    }

}

void GpuSketch::DeviceArrays::printMashDist(uint64_t h_numSequences, std::vector<std::string> seqs) 
{
    
    float * h_mashDist = new float[h_numSequences*(h_numSequences - 1)/2]; //NC2

    float * mashDist = d_mashDist;

    cudaError_t err;


    err = cudaMemcpy(h_mashDist, mashDist, (h_numSequences*(h_numSequences - 1)/2)*sizeof(float), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
        exit(1);
    }

    printf("%ld\n", h_numSequences);
    for (int i=0; i<h_numSequences; i++) 
    {
        std::cout << seqs[i] << "\t";
        for (int j=0; j<=i; j++) 
        {
            printf("1.0\t");
        }        
        for (int j=i+1; j<h_numSequences; j++) 
        {
            printf("%f\t", *h_mashDist);
            h_mashDist++;
        }
        printf("\n");        
    }
       

}