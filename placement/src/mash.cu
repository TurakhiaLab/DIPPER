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

#define THREADS_PER_BLOCK 256

// Function to handle CUDA errors
void checkCudaError(cudaError_t error, const char *file, int line) {
    if (error != cudaSuccess) {
        printf("CUDA error at %s:%d: %s\n", file, line, cudaGetErrorString(error));
        exit(EXIT_FAILURE);
    }
}



#define CHECK_CUDA_ERROR(error) checkCudaError(error, __FILE__, __LINE__)


__device__ void mashDistConstructionRowWise(
    int rowId,
    uint64_t * d_hashList,
    double * d_mashDist,
    uint64_t kmerSize,
    uint64_t sketchSize,
    int numSequences
) {
    for(int bx =0 ;bx <(numSequences+THREADS_PER_BLOCK-1)/THREADS_PER_BLOCK;bx++){
        for(int tx=0; tx < THREADS_PER_BLOCK ;tx++){
        // int tx = threadIdx.x, bx = blockIdx.x, bs = blockDim.x;
            int indx = tx+bx*(numSequences+THREADS_PER_BLOCK-1)/THREADS_PER_BLOCK;
        // if(idx>=rowId) return;
            int uni = 0, bPos = rowId, inter = 0;
            uint64_t aval, bval;
            for(int i=indx; uni < sketchSize; i+=numSequences, uni++){
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
            double jaccardEstimate = double(inter)/uni;
            d_mashDist[indx] = min(1.0, abs((log(2.0*jaccardEstimate/(1.0+jaccardEstimate)))/kmerSize));
        }
    }
    // printf("\n");
}


__device__ double jukesCantor(const uint64_t* d_compressedSeqs, const uint64_t* d_seqLengths, const uint64_t* d_prefixCompressed, int i, int j, int numSequences) {
    // Check for valid input indices
    if (i < 0 || i >= numSequences || j < 0 || j >= numSequences) {
        printf("Error: Invalid sequence indices i=%d, j=%d (numSequences=%d)\n", i, j, numSequences);
        return -1.0; // Return an error value
    }

    int mismatches = 0;
    uint64_t seqLengthI = d_seqLengths[i];
    uint64_t seqLengthJ = d_seqLengths[j];
    
    // Check for valid sequence lengths
    if (seqLengthI == 0 || seqLengthJ == 0) {
        printf("Error: Zero length sequence detected for i=%d or j=%d\n", i, j);
        return -1.0; // Return an error value
    }
    
    // Ensure we're comparing the shorter sequence length
    int totalBases = min(seqLengthI, seqLengthJ);
    
    // Calculate starting positions for each sequence in the flattened array
    uint64_t startI = d_prefixCompressed[i];
    uint64_t startJ = d_prefixCompressed[j];
    
    // Calculate how many uint64_t are needed to store each sequence
    int numUint64I = (seqLengthI + 31) / 32;
    int numUint64J = (seqLengthJ + 31) / 32;
    
    for (int k = 0; k < min(numUint64I, numUint64J); ++k) {
        // Check for array bounds
        // if (startI + k >= numSequences || startJ + k >= numSequences) {
        //     printf("Error: Array index out of bounds in jukesCantor\n");
        //     return -1.0; // Return an error value
        // }
        
        uint64_t seqI = d_compressedSeqs[startI + k];
        uint64_t seqJ = d_compressedSeqs[startJ + k];
        uint64_t xor_result = seqI ^ seqJ;
        
        // Count mismatches in this uint64_t
        mismatches += __popcll(xor_result & 0x5555555555555555ULL);
    }
    
    // Adjust mismatches if we processed more bases than totalBases
    int processedBases = min(numUint64I, numUint64J) * 32;
    if (processedBases > totalBases) {
        int excessBases = processedBases - totalBases;
        uint64_t mask = (1ULL << (2 * excessBases)) - 1;
        
        // // Check for array bounds
        // if (startI + min(numUint64I, numUint64J) - 1 >= numSequences || 
        //     startJ + min(numUint64I, numUint64J) - 1 >= numSequences) {
        //     printf("Error: Array index out of bounds when adjusting mismatches\n");
        //     return -1.0; // Return an error value
        // }
        
        uint64_t lastSeqI = d_compressedSeqs[startI + min(numUint64I, numUint64J) - 1] & ~mask;
        uint64_t lastSeqJ = d_compressedSeqs[startJ + min(numUint64I, numUint64J) - 1] & ~mask;
        uint64_t lastXor = lastSeqI ^ lastSeqJ;
        int lastMismatches = __popcll(lastXor & 0x5555555555555555ULL);
        mismatches -= lastMismatches;
    }
    
    // Check for negative mismatches (shouldn't happen, but just in case)
    if (mismatches < 0) {
        printf("Error: Negative mismatch count in jukesCantor\n");
        return -1.0; // Return an error value
    }
    
    // Calculate p (proportion of sites that differ)
    double p = static_cast<double>(mismatches) / (2 * totalBases);
    
    // Check if p is too large for the Jukes-Cantor model
    if (p >= 0.75) {
        return 1e10; // A very large number instead of infinity
    }
    
    // Calculate and return the Jukes-Cantor distance
    double distance = -0.75 * log1p(-(4.0 / 3.0) * p);
    
    // Check for NaN or infinity
    if (isnan(distance) || isinf(distance)) {
        printf("Error: Invalid Jukes-Cantor distance calculated\n");
        return -1.0; // Return an error value
    }
    
    return distance;
}

__global__ void distanceMatrixTester(uint64_t kmerSize, uint64_t sketchSize,uint64_t  * d_hashList,  int numSequences, double *d_mashDist){

    for(int i =0; i<numSequences;i++){

    
        mashDistConstructionRowWise(i,   d_hashList, 
            d_mashDist, 
            kmerSize, 
            sketchSize, 
            numSequences);

        // if(i==0){
        //     for (int j=0; j<numSequences; j++) 
        //     {
        //         std::cout<<seqs[j]<<'\t';
        //     }
        // }
        printf("\n");  

        for (int j=0; j<numSequences; j++) 
        {
            printf("%f\t", d_mashDist[j]);
        }
    }
}
__global__ void clusterKernelGPU( uint64_t kmerSize, uint64_t sketchSize, int *d_clusterMap, int numSequences, uint64_t *d_compressedSeqs, uint64_t *d_seqLengths, uint64_t *d_prefixCompressed, int MAX_LEVELS, int *d_stopFlag, int *d_sharedCount, int *d_cInstr, int nodesInThisLevel, int *stopFlag, uint64_t  * d_hashList, double *d_mashDist1, double* d_mashDist2)
{
    int size = (1 << (MAX_LEVELS + 1));
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < numSequences) {
        if (d_clusterMap[idx] >= 0) {
            bool clusterFound = false;          
            for (int clusterIdx = 0; clusterIdx < 3 * nodesInThisLevel; clusterIdx += 3) {
                if (d_cInstr[clusterIdx] == d_clusterMap[idx]) {


                    mashDistConstructionRowWise(d_cInstr[clusterIdx+1],
                    d_hashList, 
                    d_mashDist1, 
                    kmerSize, 
                    sketchSize, 
                    numSequences);
                    mashDistConstructionRowWise(d_cInstr[clusterIdx+2],
                    d_hashList, 
                    d_mashDist2, 
                    kmerSize, 
                    sketchSize, 
                    numSequences);

                    // double distance1 =1, distance2=2;
                    
                    // Error checking for jukesCantor function calls
                    if (d_cInstr[clusterIdx+1] >= numSequences || d_cInstr[clusterIdx+2] >= numSequences) {
                        printf("Error: Invalid index in d_cInstr at clusterIdx %d\n", clusterIdx);
                        return;
                    }
                    double distance1 = d_mashDist1[idx];
                    double distance2 = d_mashDist2[idx];
                    // printf("%f\t %f\n" ,distance1, distance2);
                    
                    // distance1 = jukesCantor(d_compressedSeqs, d_seqLengths, d_prefixCompressed, d_cInstr[clusterIdx+1], idx, numSequences);
                    // distance2 = jukesCantor(d_compressedSeqs, d_seqLengths, d_prefixCompressed, d_cInstr[clusterIdx+2], idx, numSequences);
                    
                    if (isinf(distance1) || isinf(distance2) || isnan(distance1) || isnan(distance2)) {
                        printf("Error: Invalid distance calculated for idx %d\n", idx);
                        return;
                    }

                    int newClusterIdx = d_cInstr[clusterIdx] * 2 + (distance1 <= distance2 ? 1 : 2);
                    
                    // Check if the new cluster index is valid
                    if (newClusterIdx >= size) {
                        printf("Error: Invalid cluster index %d for idx %d\n", newClusterIdx, idx);
                        return;
                    }
                    
                    d_clusterMap[idx] = newClusterIdx;
                    clusterFound = true;
                    break;
                }
            }
            if (!clusterFound) {
                printf("Warning: No matching cluster found for clusterMap index %d with value %d\n", idx, d_clusterMap[idx]);
            }
        }
    }


    if (idx == 0) {
        *d_stopFlag = 0;
    }

    if (idx < size) {
        d_sharedCount[idx] = 0;
    }

    __syncthreads();

    if (idx < numSequences) {
        if (d_clusterMap[idx] > 0) {
            if (d_clusterMap[idx] >= size) {
                printf("Error: Invalid d_clusterMap value %d for idx %d\n", d_clusterMap[idx], idx);
                return;
            }
            atomicAdd(&d_sharedCount[d_clusterMap[idx]], 1);
            if (d_sharedCount[d_clusterMap[idx]] <= 20) {
                d_clusterMap[idx] = -d_clusterMap[idx];
            }
        }
    }
}


void MashPlacement::MashDeviceArrays::processClusterLevels(int *clusterMap, int numSequences, treeNode *nodes[], int MAX_LEVELS, MashPlacement::Param& params) {
    int nodeIndex = 0;
    int *stopFlag = new int;
    int sharedCount = 0;
    int *d_cInstr, *d_clusterMap, *d_stopFlag, *d_sharedCount;

    // std::cout << "Starting processClusterLevels with numSequences: " << numSequences << ", MAX_LEVELS: " << MAX_LEVELS << std::endl;

    double *d_mashDist;
    double * d_mashDist1,*d_mashDist2;
    cudaMalloc(&d_mashDist, numSequences*sizeof(double));
    // distanceMatrixTester<<< 1, 1>>>(params.kmerSize, params.sketchSize,d_hashList,  numSequences,  d_mashDist);
    
    CHECK_CUDA_ERROR(cudaMalloc(&d_stopFlag, sizeof(int)));

    for (int level = 0; level < MAX_LEVELS; level++) {
        // std::cout << "Processing level " << level << std::endl;

        int nodesInThisLevel = 1 << level;
        int totalInstructions = nodesInThisLevel * 3;

        // std::cout << "Allocating cInstr array with size: " << 3 * nodesInThisLevel << std::endl;
        int *cInstr = new int[3 * nodesInThisLevel];
        if (!cInstr) {
            std::cerr << "Error: Failed to allocate cInstr array" << std::endl;
            return;
        }

        int instrIndex = 0;
        for (int i = 0; i < nodesInThisLevel; i++) {
            int parentIndex = (nodeIndex - 1) / 2;
            int baseClusterIndex = (level == 0) ? 0 : nodes[parentIndex]->nodeNum * 2 + i % 2 + 1;
            
            try {
                getTwoRandomIndices(clusterMap, numSequences, baseClusterIndex, nodes[nodeIndex]);
            } catch (const std::exception &e) {
                std::cerr << "Error in getTwoRandomIndices: " << e.what() << std::endl;
                delete[] cInstr;
                return;
            }

            cInstr[instrIndex++] = baseClusterIndex;
            cInstr[instrIndex++] = nodes[nodeIndex]->nodechild1;
            cInstr[instrIndex++] = nodes[nodeIndex]->nodechild2;
            nodeIndex++;
        }

        try {
            // std::cout << "Copying data to device and launching kernel" << std::endl;

            CHECK_CUDA_ERROR(cudaMemcpy(d_stopFlag, stopFlag, sizeof(int), cudaMemcpyHostToDevice));
            CHECK_CUDA_ERROR(cudaMalloc(&d_sharedCount, sizeof(int)));
            CHECK_CUDA_ERROR(cudaMemcpy(d_sharedCount, &sharedCount, sizeof(int), cudaMemcpyHostToDevice));
            CHECK_CUDA_ERROR(cudaMalloc(&d_cInstr, 3 * nodesInThisLevel * sizeof(int)));
            CHECK_CUDA_ERROR(cudaMalloc(&d_clusterMap, numSequences * sizeof(int)));

            CHECK_CUDA_ERROR(cudaMemcpy(d_cInstr, cInstr, 3 * nodesInThisLevel * sizeof(int), cudaMemcpyHostToDevice));
            CHECK_CUDA_ERROR(cudaMemcpy(d_clusterMap, clusterMap, numSequences * sizeof(int), cudaMemcpyHostToDevice));

            int blocksPerGrid = (numSequences + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
            // std::cout << "Launching kernel with " << blocksPerGrid << " blocks" << std::endl;
            
            cudaMalloc(&d_mashDist1, numSequences*sizeof(double));
            cudaMalloc(&d_mashDist2, numSequences*sizeof(double));

            clusterKernelGPU<<< blocksPerGrid, THREADS_PER_BLOCK>>>(
                params.kmerSize, params.sketchSize,d_clusterMap, numSequences, d_compressedSeqs, d_seqLengths, d_prefixCompressed,
                MAX_LEVELS, d_stopFlag, d_sharedCount, d_cInstr, nodesInThisLevel, stopFlag, d_hashList, d_mashDist1, d_mashDist2
            );

            CHECK_CUDA_ERROR(cudaGetLastError());
            CHECK_CUDA_ERROR(cudaDeviceSynchronize());

            // std::cout << "Kernel execution completed" << std::endl;

            CHECK_CUDA_ERROR(cudaMemcpy(clusterMap, d_clusterMap, numSequences * sizeof(int), cudaMemcpyDeviceToHost));
            CHECK_CUDA_ERROR(cudaMemcpy(stopFlag, d_stopFlag, sizeof(int), cudaMemcpyDeviceToHost));

            // std::cout << "Data copied back to host" << std::endl;

            // for (int i = 0; i < numSequences; i++) {
                // std::cout << "index " << i << " clusterMapValue " << clusterMap[i] << std::endl;
            // }

            if (*stopFlag != 0) {
                // std::cout << "Stop flag set, breaking out of loop" << std::endl;
                break;
            }

        } catch (const std::exception &e) {
            std::cerr << "Error in cluster processing: " << e.what() << std::endl;
            delete[] cInstr;
            return;
        }

        delete[] cInstr;
        CHECK_CUDA_ERROR(cudaFree(d_cInstr));
        CHECK_CUDA_ERROR(cudaFree(d_mashDist1));
        CHECK_CUDA_ERROR(cudaFree(d_mashDist2));
        CHECK_CUDA_ERROR(cudaFree(d_sharedCount));
    }

    CHECK_CUDA_ERROR(cudaFree(d_stopFlag));
    // DEBUG_PRINT("processClusterLevels completed successfully");
}

void MashPlacement::MashDeviceArrays::allocateDeviceArrays(uint64_t ** h_compressedSeqs, uint64_t * h_seqLengths, size_t num, Param& params)
{
    cudaError_t err;

    numSequences = int(num);

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

    // printf("%p", d_hashList);


    err = cudaMalloc(&d_compressedSeqs, flatStringLength*sizeof(uint64_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    err = cudaMalloc(&d_prefixCompressed, numSequences*sizeof(uint64_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
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
        dev_seqLengths, dev_seqLengths + numSequences, dev_prefixComp, 
        [] __device__ (const uint64_t& x) -> uint64_t { 
            return (x + 31) / 32;
        }
    );

    thrust::exclusive_scan(dev_prefixComp, dev_prefixComp + numSequences, dev_prefixComp);

    free(h_flattenCompressSeqs);
    free(h_aggseqLengths);
    cudaDeviceSynchronize();
}

void MashPlacement::MashDeviceArrays::deallocateDeviceArrays(){
    cudaFree(d_compressedSeqs);
    cudaFree(d_aggseqLengths);
    cudaFree(d_seqLengths);
    cudaFree(d_prefixCompressed);
    cudaFree(d_hashList);
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
    uint64_t * d_prefixCompressed,
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

void MashPlacement::MashDeviceArrays::sketchConstructionOnGpu(Param& params){
    const uint64_t kmerSize = params.kmerSize; // Extract kmerSize
    auto timerStart = std::chrono::high_resolution_clock::now();

    int threadsPerBlock = 512;
    int blocksPerGrid = (numSequences + threadsPerBlock - 1) / threadsPerBlock;
    size_t sharedMemorySize = sizeof(uint64_t) * (2000);
    sketchConstruction<<<1, threadsPerBlock, sharedMemorySize>>>(
        d_compressedSeqs, d_seqLengths, d_prefixCompressed, numSequences, d_hashList, kmerSize
    );

    uint64_t * temp_hashList;
    auto err = cudaMalloc(&temp_hashList, numSequences*int(params.sketchSize)*sizeof(uint64_t));
    if (err != cudaSuccess){
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    rearrangeHashList <<<blocksPerGrid, threadsPerBlock >>>(
        numSequences,
        int(params.sketchSize),
        d_hashList,
        temp_hashList
    );
    std::swap(d_hashList, temp_hashList);
    cudaFree(temp_hashList);

    
    err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("CUDA Error: %s\n", cudaGetErrorString(err));
    }
    cudaDeviceSynchronize();
    h_hashList = new uint64_t[1000*numSequences];
    cudaMemcpy(h_hashList,d_hashList,1000*numSequences*sizeof(uint64_t),cudaMemcpyDeviceToHost);
    auto timerEnd = std::chrono::high_resolution_clock::now();
    auto time = timerEnd - timerStart;
    // std::cout << "Time to generate hashes: " << time.count() << "ns\n";

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
    double jaccardEstimate = double(inter)/uni;
    d_mashDist[idx] = min(1.0, abs((log(2.0*jaccardEstimate/(1.0+jaccardEstimate)))/kmerSize));
}

void MashPlacement::MashDeviceArrays::distConstructionOnGpu(Param& params, int rowId, double* d_mashDist) const{
    // if(rowId%100==0) cudaMemcpy(d_hashList,h_hashList,1000*numSequences*sizeof(uint64_t),cudaMemcpyHostToDevice);
    int threadNum = 256, blockNum = (numSequences+threadNum-1)/threadNum;
    mashDistConstruction <<<threadNum, blockNum>>> (
        rowId, 
        d_hashList, 
        d_mashDist, 
        params.kmerSize, 
        params.sketchSize, 
        numSequences
    );
}


void MashPlacement::MashDeviceArrays::printSketchValues(int numValues) 
{
    uint64_t * h_hashList = new uint64_t[1000*numSequences];


    uint64_t * hashList = d_hashList;

    cudaError_t err;

    //printf("Total Hashes: %d", numSequences*1000);

    err = cudaMemcpy(h_hashList, hashList, numSequences*1000*sizeof(uint64_t), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
        exit(1);
    }

    // printf("i\thashList[i] (%zu)\n");
    for (int j = 0; j < numSequences; j++) {
        fprintf(stderr, "Sequence (%d)\n", j);
        for (int i=0; i<numValues; i++) {
            fprintf(stderr, "%i\t%lu\n", i, h_hashList[i*numSequences+j]);
        }
    }
    

}


__global__ void mashDistConstructionRowPrint(
    int rowId,
    uint64_t * d_hashList,
    double * d_mashDist,
    uint64_t kmerSize,
    uint64_t sketchSize,
    int numSequences
) {
    int tx = threadIdx.x, bx = blockIdx.x, bs = blockDim.x;
    int idx = tx+bx*bs;
    // if(idx>=rowId) return;
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
    double jaccardEstimate = double(inter)/uni;
    d_mashDist[idx] = min(1.0, abs((log(2.0*jaccardEstimate/(1.0+jaccardEstimate)))/kmerSize));
}

void MashPlacement::MashDeviceArrays::distConstructionOnGpuRowPrint(Param& params, int rowId, double* d_mashDist) const{
    // if(rowId%100==0) cudaMemcpy(d_hashList,h_hashList,1000*numSequences*sizeof(uint64_t),cudaMemcpyHostToDevice);
    int threadNum = 256, blockNum = (numSequences+threadNum-1)/threadNum;
    mashDistConstructionRowPrint <<<threadNum, blockNum>>> (
        rowId, 
        d_hashList, 
        d_mashDist, 
        params.kmerSize, 
        params.sketchSize, 
        numSequences
    );
}


void MashPlacement::MashDeviceArrays::printMashDist(Param& params,uint64_t h_numSequences, std::vector<std::string> seqs) 
{
    
    cudaError_t err;
    
    double *d_mashDist;
    err = cudaMalloc(&d_mashDist, numSequences*sizeof(double));
 
    for (int i=0; i<h_numSequences; i++) 
    { 
        double * h_mashDist = new double[h_numSequences];
        err = cudaMemcpy(h_mashDist, d_mashDist, (h_numSequences)*sizeof(uint64_t), cudaMemcpyDeviceToHost);
        distConstructionOnGpuRowPrint(params, i,d_mashDist);
        if (err != cudaSuccess) {
            fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
            exit(1);
        }

        // if(i==0){
        //     for (int j=0; j<h_numSequences; j++) 
        //     {
        //         std::cout<<seqs[j]<<'\t';
        //     }
        // }
        printf("\n");  

        for (int j=0; j<h_numSequences; j++) 
        {
            printf("%f\t", h_mashDist[j]);
        }
    
        delete []    h_mashDist;
    }
    
    printf("\n");
}