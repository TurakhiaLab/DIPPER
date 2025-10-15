#include "mash_placement.cuh"

#include <stdio.h>
#include <queue>
#include <tbb/parallel_for.h>
#include <chrono>
#include <cassert>
#include <iostream>

#define BIG_CONSTANT(x) (x##LLU)

// CPU Implementation start here


inline void exclusive_scan(uint64_t* data, size_t n) {
    uint64_t runningSum = 0;
    for (size_t i = 0; i < n; ++i) {
        uint64_t temp = data[i];
        data[i] = runningSum;
        runningSum += temp;
    }
}

void MashPlacement::MashDeviceArrays::allocateDeviceArrays(uint64_t ** h_compressedSeqs, uint64_t * h_seqLengths, size_t num, Param& params)
{
    numSequences = int(num);

    uint64_t kmerSize = params.kmerSize;
    size_t hashListLength = 0;   

    // Allocate memory
    d_aggseqLengths = new uint64_t[numSequences];  
    d_seqLengths = new uint64_t[numSequences]; 

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
    params.flatStringLength = flatStringLength;

    d_hashList = new uint64_t[1000*numSequences];
    d_compressedSeqs = new uint64_t[flatStringLength];
    d_prefixCompressed = new uint64_t[numSequences];

    // Transfer data
    for (size_t i = 0; i < numSequences; i++) {
        d_aggseqLengths[i] = h_aggseqLengths[i];
        d_seqLengths[i] = h_seqLengths[i];
    }
    for (size_t i = 0; i < flatStringLength; i++) {
        d_compressedSeqs[i] = h_flattenCompressSeqs[i];
    }

    // generate prefix array    
    for (int i = 0; i < numSequences; i++) {
        d_prefixCompressed[i] = (d_seqLengths[i] + 31) / 32;
    }
    exclusive_scan(d_prefixCompressed, numSequences);

}

void MashPlacement::MashDeviceArrays::deallocateDeviceArrays(){
    delete [] d_compressedSeqs;
    delete [] d_aggseqLengths;
    delete [] d_seqLengths;
    delete [] d_prefixCompressed;
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

inline uint64_t rotl64_cpu(uint64_t x, int8_t r) {
    return (x << r) | (x >> (64 - r));
}

inline uint64_t getblock64_cpu(const uint64_t *p, int i) {
    return p[i];
}

inline uint64_t fmix64_cpu(uint64_t k) {
    k ^= k >> 33;
    k *= BIG_CONSTANT(0xff51afd7ed558ccd);
    k ^= k >> 33;
    k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
    k ^= k >> 33;
    return k;
}


void MurmurHash3_x64_128_MASH_CPU( void * key, const int len, const uint32_t seed, void * out)
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
        uint64_t k1 = getblock64_cpu(blocks,i*2+0);
        uint64_t k2 = getblock64_cpu(blocks,i*2+1);

        k1 *= c1; k1  = rotl64_cpu(k1,31); k1 *= c2; h1 ^= k1;

        h1 = rotl64_cpu(h1,27); h1 += h2; h1 = h1*5+0x52dce729;

        k2 *= c2; k2  = rotl64_cpu(k2,33); k2 *= c1; h2 ^= k2;

        h2 = rotl64_cpu(h2,31); h2 += h1; h2 = h2*5+0x38495ab5;
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
               k2 *= c2; k2  = rotl64_cpu(k2,33); k2 *= c1; h2 ^= k2;

        case  8: k1 ^= ((uint64_t)tail[ 7]) << 56;
        case  7: k1 ^= ((uint64_t)tail[ 6]) << 48;
        case  6: k1 ^= ((uint64_t)tail[ 5]) << 40;
        case  5: k1 ^= ((uint64_t)tail[ 4]) << 32;
        case  4: k1 ^= ((uint64_t)tail[ 3]) << 24;
        case  3: k1 ^= ((uint64_t)tail[ 2]) << 16;
        case  2: k1 ^= ((uint64_t)tail[ 1]) << 8;
        case  1: k1 ^= ((uint64_t)tail[ 0]) << 0;
               k1 *= c1; k1  = rotl64_cpu(k1,31); k1 *= c2; h1 ^= k1;
    };

    //----------
    // finalization

    h1 ^= len; h2 ^= len;

    h1 += h2;
    h2 += h1;

    h1 = fmix64_cpu(h1);
    h2 = fmix64_cpu(h2);

    h1 += h2;
    h2 += h1;

    ((uint64_t*)out)[0] = h1;
    ((uint64_t*)out)[1] = h2;
}

// CPU version of decompress
void decompressCPU(
    uint64_t compressedSeq, 
    uint64_t kmerSize,
    char* decompressedSeq_fwd, 
    char* decompressedSeq_rev) 
    {
    static const char lookupTable[4] = {'A', 'C', 'G', 'T'};
    for (int i = kmerSize - 1; i >= 0; i--) {
        uint64_t twoBitVal = (compressedSeq >> (2 * i)) & 0x3;
        decompressedSeq_fwd[i] = lookupTable[twoBitVal];
        decompressedSeq_rev[kmerSize - 1 - i] = lookupTable[3 - twoBitVal];
    }
}

// CPU version of memcmp_device
int memcmpCPU(const char* kmer_fwd, const char* kmer_rev, int kmerSize) {
    for (int i = 0; i < kmerSize; i++) {
        if (kmer_fwd[i] < kmer_rev[i]) return -1;
        if (kmer_fwd[i] > kmer_rev[i]) return 1;
    }
    return 0;
}


void sketchConstruction(
    uint64_t* d_compressedSeqs,
    uint64_t* d_seqLengths,       // bases (not words)
    uint64_t* d_prefixCompressed, // offsets in words into d_compressedSeqs
    size_t    numSequences,
    uint64_t* d_hashList,         // output: numSequences * 1000
    uint64_t  kmerSize            // in bases
) {

    tbb::parallel_for(tbb::blocked_range<int>(0, numSequences), [&](tbb::blocked_range<int> range){ 
    for (int i = range.begin(); i < range.end(); ++i) {
    
        uint64_t stored[2000];
        const int threadsPerBlock = 512;
        int blocksPerGrid = (numSequences + threadsPerBlock - 1) / threadsPerBlock;
        std::vector<uint64_t>keys (threadsPerBlock * 3);
    // for (size_t i = 0; i < numSequences; ++i) {
        int bx = i;
        //if (tx == 0) printf("Started block %d\n", i);
        for (int tx = 0; tx < threadsPerBlock; ++tx) {
            stored[tx]       = 0xFFFFFFFFFFFFFFFF; // reset block's stored values
            stored[tx + 500] = 0xFFFFFFFFFFFFFFFF;
        }

        uint64_t seqLength = d_seqLengths[i];
        uint64_t * compressedSeqs = d_compressedSeqs + d_prefixCompressed[i];

        size_t iteration = 0;
        while (iteration <= seqLength - kmerSize) {
            for (int j = iteration; j < iteration+threadsPerBlock; ++j) {
                int tx = j - iteration;
                if (j <= seqLength - kmerSize) {
                    uint64_t kmer = 0;
                    char kmer_fwd[32] = {0};
                    char kmer_rev[32] = {0};
                    uint8_t out[16];
                    uint64_t index = j/32;
                    uint64_t shift1 = 2*(j%32);
                    if (shift1>0) {
                        uint64_t shift2 = 64-shift1;
                        kmer = ((compressedSeqs[index] >> shift1) | (compressedSeqs[index+1] << shift2)); //& mask;
                    }
                    else {   
                        kmer = compressedSeqs[index];// & mask;
                    }
                    decompressCPU(kmer, kmerSize, kmer_fwd, kmer_rev);
                    MurmurHash3_x64_128_MASH_CPU(
                        (memcmpCPU(kmer_fwd, kmer_rev, static_cast<int>(kmerSize)) <= 0) ? kmer_fwd : kmer_rev,
                        static_cast<int>(kmerSize),
                        42,
                        out
                    );
                    uint64_t hash = *((uint64_t *)out);
                    // Combine stored and computed to sort and rank
                    keys[3*tx+0] = (tx < 500) ? stored[tx] : 0xFFFFFFFFFFFFFFFF;
                    keys[3*tx+1] = (tx < 500) ? stored[tx + 500] : 0xFFFFFFFFFFFFFFFF;
                    keys[3*tx+2] = hash;
                }
                else {
                    keys[3*tx+0] = (tx < 500) ? stored[tx] : 0xFFFFFFFFFFFFFFFF;
                    keys[3*tx+1] = (tx < 500) ? stored[tx + 500] : 0xFFFFFFFFFFFFFFFF;
                    keys[3*tx+2] = 0xFFFFFFFFFFFFFFFF;
                }
            }

            // if (i == 0) {
            //     printf("%d, Before Sort: ", iteration);
            //     for (int kk = 0; kk < 3; ++kk) printf("%22llu ", keys[kk]);
            //         printf("\n");
            // }
            std::sort(keys.begin(), keys.end());
            // if (i == 0) {
            //     printf("%d, After Sort:  ", iteration);
            //     for (int kk = 0; kk < 3; ++kk) printf("%22llu ", keys[kk]);
            //         printf("\n");
            // }
            // Move top 1000 hashes back to stored
            for (int j = iteration; j <= iteration+threadsPerBlock; ++j) {
                int tx = j - iteration;
                if (tx < 333) {
                    stored[3*tx] =     keys[3*tx+0];
                    stored[3*tx + 1] = keys[3*tx+1];
                    stored[3*tx + 2] = keys[3*tx+2];
                } else if (tx == 333) {
                    stored[999] = keys[3*tx+0];
                }
            }
            iteration += threadsPerBlock;
        }
    
        // Result writing back to global memory.
        for (int tx = 0; tx < threadsPerBlock; ++tx) {
            if (tx < 500) {
                d_hashList[1000*i+tx] = stored[tx];
                d_hashList[1000*i+tx + 500] = stored[tx + 500];
            }
        }
    // }
        keys.clear();
    }
    });
    return;

}


void rearrangeHashList(
    int numSequences,
    int sketchSize,
    uint64_t* original,
    uint64_t* target
) {
    for (int seq = 0; seq < numSequences; seq++) {
        for (int i = 0; i < sketchSize; i++) {
            target[i * numSequences + seq] =
                original[seq * sketchSize + i];
        }
    }
    return;
}

void MashPlacement::MashDeviceArrays::sketchConstructionOnGpu(
    Param &params) {

    const uint64_t kmerSize = params.kmerSize;
    const int sketchSize = params.sketchSize;
    auto timerStart = std::chrono::high_resolution_clock::now();

    // Temporary CPU buffers for device data
    // std::vector<uint64_t> compressedSeqsHost(params.flatStringLength); // adjust size if needed
    // std::vector<uint64_t> seqLengthsHost(numSequences);
    // std::vector<uint64_t> prefixCompressedHost(numSequences);
    // std::vector<uint64_t> hashListHost(numSequences * 1000);
    // cudaMemcpy(compressedSeqsHost.data(), d_compressedSeqs, sizeof(uint64_t) * compressedSeqsHost.size(), cudaMemcpyDeviceToHost);
    // cudaMemcpy(seqLengthsHost.data(), d_seqLengths, sizeof(uint64_t) * seqLengthsHost.size(), cudaMemcpyDeviceToHost);
    // cudaMemcpy(prefixCompressedHost.data(), d_prefixCompressed, sizeof(uint64_t) * prefixCompressedHost.size(), cudaMemcpyDeviceToHost);
    // CPU version of sketchConstruction

    /*
    sketchConstruction(
        compressedSeqsHost.data(),
        seqLengthsHost.data(),
        prefixCompressedHost.data(),
        numSequences,
        hashListHost.data(),
        kmerSize
    );
    */

    sketchConstruction(
        d_compressedSeqs,
        d_seqLengths,
        d_prefixCompressed,
        numSequences,
        d_hashList,
        kmerSize
    );

    // Rearrange hashes on CPU
    uint64_t * temp_hashList = new uint64_t[numSequences * sketchSize];
    rearrangeHashList(numSequences, sketchSize, d_hashList, temp_hashList);
    std::swap(d_hashList, temp_hashList);
    delete [] temp_hashList;

    // for (int i = 0; i < numSequences; i++) {
    //     printf("Sequence %d: \n", i);
    //     for (int j = 0; j < 20; j++) {
    //         printf("%llu,", d_hashList[j * numSequences + i]);
    //         if (j % 5 == 4) printf("\n");
    //     }
    // }

    // Copy back to device memory if needed
    // q_ct1.memcpy(d_hashList, temp_hashList.data(), sizeof(uint64_t) * temp_hashList.size()).wait();
	// std::memcpy(d_hashList, temp_hashList.data(), sizeof(uint64_t) * temp_hashList.size());
    // cudaMemcpy(d_hashList, temp_hashList.data(), sizeof(uint64_t) * temp_hashList.size(), cudaMemcpyHostToDevice);

    auto timerEnd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = timerEnd - timerStart;
    std::cout << "Time to generate hashes (CPU sequential): " << elapsed.count() << " s\n";
}

void mashDistConstruction(
    int rowId,
    uint64_t * d_hashList,
    double * d_mashDist,
    uint64_t kmerSize,
    uint64_t sketchSize,
    int numSequences
) {
    tbb::parallel_for(tbb::blocked_range<int>(0, rowId), [&](tbb::blocked_range<int> range){ 
    for (int idx = range.begin(); idx < range.end(); ++idx) {
    // for (int idx = 0; idx < rowId; ++idx) {
        int uni = 0, bPos = rowId, inter = 0;
        uint64_t aval, bval;
        for(int i=idx; uni < sketchSize; i+=numSequences, uni++){
            aval = d_hashList[i];
            while(uni < sketchSize && bPos < numSequences * sketchSize){
                bval = d_hashList[bPos];
                if(bval > aval) break;
                if(bval < aval) uni++;
                else inter++;
                bPos += numSequences;
            }
            if(uni >= sketchSize) break;
        }
        assert(uni==1000);
        double jaccardEstimate = max(double(inter),1.0)/uni;
        d_mashDist[idx] = min(1.0, abs(log(2.0*jaccardEstimate/(1.0+jaccardEstimate))/kmerSize));
    }
    });
}

void MashPlacement::MashDeviceArrays::distConstructionOnGpu(Param& params, int rowId, double* d_mashDist) const{
    
    mashDistConstruction (
        rowId, 
        d_hashList, 
        d_mashDist, 
        params.kmerSize, 
        params.sketchSize, 
        numSequences
    );
            
}
