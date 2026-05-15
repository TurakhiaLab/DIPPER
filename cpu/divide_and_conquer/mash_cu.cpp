/**
 * cpu/divide_and_conquer/mash_cu.cpp
 *
 * CPU DC MASH: allocate DC arrays (batch/backbone hash lists, leaf map),
 * MurmurHash3, sketch construction for DC, distConstruction for backbone/range/special-ID.
 * Mirrors GPU DC mash but runs on host with TBB.
 */

#include "../mash_placement.hpp"

#include <stdio.h>
#include <queue>
#include <chrono>
#include <iostream>
#include <tbb/parallel_for.h>
#include <cmath>

/** Allocate DC MASH host arrays (batch/backbone hash lists, leaf map). */
void MashPlacement::MashDeviceArraysDC::allocateDeviceArraysDC(uint64_t **h_compressedSeqs, uint64_t *h_seqLengths, size_t num, Param &params)
{
    this->totalNumSequences = num;
    this->backboneSize = params.backboneSize;
    /* Allocate device arrays */
    size_t maxLength = 0;
    for (size_t i = 0; i < num; i++)
    {
        if (h_seqLengths[i] > maxLength)
            maxLength = h_seqLengths[i];
    }
    size_t maxLengthCompressed = params.isProtein ? (maxLength + 7) / 8 : (maxLength + 31) / 32;

    d_seqLengths = new uint64_t[params.batchSize];
    d_hashList = new uint64_t[params.sketchSize * params.batchSize];
    d_compressedSeqs = new uint64_t[params.batchSize * maxLengthCompressed];
    d_hashListConst = new uint64_t[params.sketchSize * params.backboneSize];
    d_hashListBackbone = new uint64_t[params.sketchSize * params.backboneSize];
    d_leafID_map = new int[totalNumSequences];
    return;
}

void MashPlacement::MashDeviceArraysDC::deallocateDeviceArraysDC()
{
    // cudaFree(d_compressedSeqs);
    // cudaFree(d_seqLengths);
    // cudaFree(d_hashList);
    delete[] d_compressedSeqs;
    delete[] d_seqLengths;
    delete[] d_hashList;
}

#define BIG_CONSTANTDC(x) (x##LLU)

inline uint64_t rotl64DC(uint64_t x, int8_t r)
{
    return (x << r) | (x >> (64 - r));
}

#define ROTL64DC(x, y) rotl64DC(x, y)

uint64_t getblock64DC(const uint64_t *p, int i)
{
    return p[i];
}

uint64_t fmix64DC(uint64_t k)
{
    k ^= k >> 33;
    k *= BIG_CONSTANTDC(0xff51afd7ed558ccd);
    k ^= k >> 33;
    k *= BIG_CONSTANTDC(0xc4ceb9fe1a85ec53);
    k ^= k >> 33;

    return k;
}

// First hashing function using raw sequence
void MurmurHash3_x64_128_MASHDC(void *key, const int len, const uint32_t seed, void *out)
{
    const uint8_t *data = (const uint8_t *)key;
    const int nblocks = len / 16;

    uint64_t h1 = seed;
    uint64_t h2 = seed;

    const uint64_t c1 = BIG_CONSTANTDC(0x87c37b91114253d5);
    const uint64_t c2 = BIG_CONSTANTDC(0x4cf5ad432745937f);

    //----------
    // body

    const uint64_t *blocks = (const uint64_t *)(data);

    for (int i = 0; i < nblocks; i++)
    {
        uint64_t k1 = getblock64DC(blocks, i * 2 + 0);
        uint64_t k2 = getblock64DC(blocks, i * 2 + 1);

        k1 *= c1;
        k1 = ROTL64DC(k1, 31);
        k1 *= c2;
        h1 ^= k1;

        h1 = ROTL64DC(h1, 27);
        h1 += h2;
        h1 = h1 * 5 + 0x52dce729;

        k2 *= c2;
        k2 = ROTL64DC(k2, 33);
        k2 *= c1;
        h2 ^= k2;

        h2 = ROTL64DC(h2, 31);
        h2 += h1;
        h2 = h2 * 5 + 0x38495ab5;
    }

    //----------
    // tail

    const uint8_t *tail = (const uint8_t *)(data + nblocks * 16);

    uint64_t k1 = 0;
    uint64_t k2 = 0;

    switch (len & 15)
    {
    case 15:
        k2 ^= ((uint64_t)tail[14]) << 48;
    case 14:
        k2 ^= ((uint64_t)tail[13]) << 40;
    case 13:
        k2 ^= ((uint64_t)tail[12]) << 32;
    case 12:
        k2 ^= ((uint64_t)tail[11]) << 24;
    case 11:
        k2 ^= ((uint64_t)tail[10]) << 16;
    case 10:
        k2 ^= ((uint64_t)tail[9]) << 8;
    case 9:
        k2 ^= ((uint64_t)tail[8]) << 0;
        k2 *= c2;
        k2 = ROTL64DC(k2, 33);
        k2 *= c1;
        h2 ^= k2;

    case 8:
        k1 ^= ((uint64_t)tail[7]) << 56;
    case 7:
        k1 ^= ((uint64_t)tail[6]) << 48;
    case 6:
        k1 ^= ((uint64_t)tail[5]) << 40;
    case 5:
        k1 ^= ((uint64_t)tail[4]) << 32;
    case 4:
        k1 ^= ((uint64_t)tail[3]) << 24;
    case 3:
        k1 ^= ((uint64_t)tail[2]) << 16;
    case 2:
        k1 ^= ((uint64_t)tail[1]) << 8;
    case 1:
        k1 ^= ((uint64_t)tail[0]) << 0;
        k1 *= c1;
        k1 = ROTL64DC(k1, 31);
        k1 *= c2;
        h1 ^= k1;
    };

    //----------
    // finalization

    h1 ^= len;
    h2 ^= len;

    h1 += h2;
    h2 += h1;

    h1 = fmix64DC(h1);
    h2 = fmix64DC(h2);

    h1 += h2;
    h2 += h1;

    ((uint64_t *)out)[0] = h1;
    ((uint64_t *)out)[1] = h2;
}

void decompressDC(uint64_t compressedSeq, uint64_t kmerSize, char *decompressedSeq_fwd, char *decompressedSeq_rev)
{
    static const char lookupTable[4] = {'A', 'C', 'G', 'T'};
    for (int i = kmerSize - 1; i >= 0; i--)
    {
        uint64_t twoBitVal = (compressedSeq >> (2 * i)) & 0x3;
        decompressedSeq_fwd[i] = lookupTable[twoBitVal];
        decompressedSeq_rev[kmerSize - 1 - i] = lookupTable[3 - twoBitVal];
    }
}

int memcmp_deviceDC(const char *kmer_fwd, const char *kmer_rev, int kmerSize)
{
    for (int i = 0; i < kmerSize; i++)
    {
        if (kmer_fwd[i] < kmer_rev[i])
        {
            return -1;
        }
        if (kmer_fwd[i] > kmer_rev[i])
        {
            return 1;
        }
    }
    return 0;
}

static void decompress_protein_kmer_cpuDC(
    uint64_t *compressedSeqs, uint64_t j, uint64_t kmerSize, char *kmer_fwd)
{
    static const char dec[21] = {
        'A','C','G','T','D','E','F','H','I','K','L','M','N','P','Q','R','S','V','W','Y','X'
    };
    for (uint64_t p = 0; p < kmerSize; p++)
    {
        uint64_t pos = j + p;
        uint64_t wi = pos / 8;
        int sh = (int)((pos % 8) * 8);
        uint64_t val = (compressedSeqs[wi] >> sh) & 0xFFULL;
        if (val > 20) val = 20;
        kmer_fwd[p] = dec[val];
    }
}

void sketchConstructionDC(
    uint64_t *d_compressedSeqs,
    uint64_t *d_seqLengths,
    size_t maxLengthCompressed,
    size_t numSequences,
    uint64_t *d_hashList,
    uint64_t kmerSize,
    bool isProtein)
{

    tbb::parallel_for(tbb::blocked_range<int>(0, numSequences), [&](tbb::blocked_range<int> range)
                      { 
    for (int i = range.begin(); i < range.end(); ++i) {
        uint64_t stored[2000];
        const int threadsPerBlock = 512;
        int blocksPerGrid = (numSequences + threadsPerBlock - 1) / threadsPerBlock;
        std::vector<uint64_t>keys (threadsPerBlock * 3);
        for (int tx = 0; tx < threadsPerBlock; ++tx) {
            stored[tx]       = 0xFFFFFFFFFFFFFFFF; // reset block's stored values
            stored[tx + 500] = 0xFFFFFFFFFFFFFFFF;
        }
        uint64_t seqLength = d_seqLengths[i];
        uint64_t * compressedSeqs = d_compressedSeqs + (uint64_t)maxLengthCompressed*i;
        // size_t j = tx;
        size_t iteration = 0;
        while (iteration <= seqLength - kmerSize) {
            for (int j = iteration; j < iteration+threadsPerBlock; ++j) {
                int tx = j - iteration;
                if (j <= seqLength - kmerSize) {
                    uint64_t kmer = 0;
                    char kmer_fwd[32] = {0};
                    char kmer_rev[32] = {0};
                    uint8_t out[16];
                    if (isProtein) {
                        decompress_protein_kmer_cpuDC(compressedSeqs, j, kmerSize, kmer_fwd);
                        MurmurHash3_x64_128_MASHDC(
                            kmer_fwd,
                            static_cast<int>(kmerSize),
                            42,
                            out);
                    } else {
                        uint64_t index = j/32;
                        uint64_t shift1 = 2*(j%32);
                        if (shift1>0) {
                            uint64_t shift2 = 64-shift1;
                            kmer = ((compressedSeqs[index] >> shift1) | (compressedSeqs[index+1] << shift2));
                        }
                        else {
                            kmer = compressedSeqs[index];
                        }
                        decompressDC(kmer, kmerSize, kmer_fwd, kmer_rev);
                        MurmurHash3_x64_128_MASHDC(
                            (memcmp_deviceDC(kmer_fwd, kmer_rev, static_cast<int>(kmerSize)) <= 0) ? kmer_fwd : kmer_rev,
                            static_cast<int>(kmerSize),
                            42,
                            out
                        );
                    }
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
            std::sort(keys.begin(), keys.end());
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
    } });
}

void rearrangeHashListDC(
    int numSequences,
    int sketchSize,
    uint64_t *original,
    uint64_t *target)
{
    tbb::parallel_for(tbb::blocked_range<int>(0, numSequences), [&](tbb::blocked_range<int> range) {
    for (int idx = range.begin(); idx < range.end(); ++idx) {
        for(int i=0;i<sketchSize;i++){
            target[i*numSequences+idx] = original[idx*sketchSize + i];
        }
    } 
    });
    // if(idx>=numSequences) return;
    // for(int i=0;i<sketchSize;i++){
    //     target[i*numSequences+idx] = original[idx*sketchSize + i];
    // }
}

void MashPlacement::MashDeviceArraysDC::sketchConstructionOnGpuDC(Param &params, uint64_t **h_compressedSeqs, uint64_t *seqLengths, uint64_t numSequences)
{

    h_hashList = new uint64_t[params.sketchSize * numSequences];

    size_t maxLength = 0;
    for (size_t i = 0; i < numSequences; i++)
    {
        if (seqLengths[i] > maxLength)
            maxLength = seqLengths[i];
    }
    size_t maxLengthCompressed = params.isProtein ? (maxLength + 7) / 8 : (maxLength + 31) / 32;
    const uint64_t kmerSize = params.kmerSize; // Extract kmerSize
    auto timerStart = std::chrono::high_resolution_clock::now();

    uint64_t localBatchSize = params.batchSize;
    for (int i = 0; i < numSequences; i += localBatchSize)
    {
        std::cerr << "Processing batch " << i << std::endl;
        uint64_t *h_seqLengths = new uint64_t[localBatchSize];
        uint64_t *h_flattenCompressSeqs = new uint64_t[localBatchSize * maxLengthCompressed];
        if (i + localBatchSize > numSequences)
        {
            localBatchSize = numSequences - i;
        }
        for (auto j = i; j < i + localBatchSize && j < numSequences; j++)
        {
            size_t nWords = params.isProtein ? (seqLengths[j] + 7) / 8 : (seqLengths[j] + 31) / 32;
            for (size_t k = 0; k < nWords; k++)
            {
                h_flattenCompressSeqs[(j - i) * maxLengthCompressed + k] = h_compressedSeqs[j][k];
            }
            h_seqLengths[j - i] = seqLengths[j];
        }
        for (int copy = 0; copy < localBatchSize * maxLengthCompressed; ++copy)
            d_compressedSeqs[copy] = h_flattenCompressSeqs[copy];
        for (int copy = 0; copy < localBatchSize; ++copy)
            d_seqLengths[copy] = h_seqLengths[copy];

        int threadsPerBlock = 512;
        int blocksPerGrid = 1024;
        // size_t sharedMemorySize = sizeof(uint64_t) * (2000);
        sketchConstructionDC(
            d_compressedSeqs, d_seqLengths, maxLengthCompressed, localBatchSize, d_hashList, kmerSize,
            params.isProtein);

        for (int copy = 0; copy < params.sketchSize * localBatchSize; ++copy)
            h_hashList[copy] = d_hashList[copy];

        h_hashList += params.sketchSize * localBatchSize;
    }

    h_hashList -= params.sketchSize * numSequences;

    /* Rearrange only for backbone tree */
    uint64_t *temp_hashList;
    temp_hashList = new uint64_t[params.sketchSize * params.backboneSize];
    for (int copy = 0; copy < params.sketchSize * params.backboneSize; ++copy)
        d_hashListBackbone[copy] = h_hashList[copy];

    int threadsPerBlock = 512;
    int blocksPerGrid = 1024;
    rearrangeHashListDC(
        params.backboneSize,
        int(params.sketchSize),
        d_hashListBackbone,
        temp_hashList);
    std::swap(d_hashListBackbone, temp_hashList);
    delete[] temp_hashList;
    // cudaFree(temp_hashList);

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

void mashDistConstructionDC(
    int rowId,
    uint64_t *d_hashList,
    double *d_mashDist,
    uint64_t kmerSize,
    uint64_t sketchSize,
    int numSequences)
{
    // int tx = threadIdx.x, bx = blockIdx.x, bs = blockDim.x;
    // int idx = tx+bx*bs;
    // if(idx>=rowId) return;
    tbb::parallel_for(tbb::blocked_range<int>(0, rowId), [&](tbb::blocked_range<int> range)
                      {
    for (int idx_ = range.begin(); idx_ < range.end(); ++idx_) {
        int uni = 0, bPos = rowId, inter = 0;
        uint64_t aval, bval;
        for(int i=idx_; uni < sketchSize; i+=numSequences, uni++){
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
        double jaccardEstimate = std::max(double(inter),1.0)/uni;
        d_mashDist[idx_] = std::min(1.0, fabs(log(2.0*jaccardEstimate/(1.0+jaccardEstimate))/kmerSize));
    } });
}

void mashDistConstructionRangeForClusteringDC(
    int rowId,
    uint64_t *d_hashListBackbone,
    uint64_t *d_hashListConst,
    double *d_mashDist,
    uint64_t kmerSize,
    uint64_t sketchSize,
    int numSequences,
    int st,
    int ed)
{
    // int tx = threadIdx.x, bx = blockIdx.x, bs = blockDim.x, gs= gridDim.x;
    // int idx = tx+bx*bs;
    tbb::parallel_for(tbb::blocked_range<int>(0, ed - st), [&](tbb::blocked_range<int> range)
                      {
    for (int idx_ = range.begin(); idx_ < range.end(); ++idx_) {
    // for (int idx_=idx; idx_<=ed-st; idx_+=gs*bs){
        if (idx_>ed-st) return;
        idx_ += st;
        int uni = 0, bPos = rowId*sketchSize, inter = 0;
        uint64_t aval, bval;
        for(int i=idx_; uni < sketchSize; i+=numSequences, uni++){
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
        double jaccardEstimate = std::max(double(inter),1.0)/uni;
        d_mashDist[idx_] = std::min(1.0, fabs(log(2.0*jaccardEstimate/(1.0+jaccardEstimate))/kmerSize));
    } });
}

void mashDistConstructionRangeDC(
    int rowId,
    uint64_t *d_hashList,
    double *d_mashDist,
    uint64_t kmerSize,
    uint64_t sketchSize,
    int numSequences,
    int st,
    int ed)
{
    // int tx = threadIdx.x, bx = blockIdx.x, bs = blockDim.x, gs = gridDim.x;
    // int idx = tx+bx*bs;
    tbb::parallel_for(tbb::blocked_range<int>(0, ed - st), [&](tbb::blocked_range<int> range)
                      {
    for (int idx_ = range.begin(); idx_ < range.end(); ++idx_) {
    // for (int idx_=idx; idx_<=ed-st; idx_+=gs*bs){
        if (idx_>ed-st) return;
        idx_+=st;
        int uni = 0, bPos = rowId, inter = 0;
        uint64_t aval, bval;
        for(int i=idx_; uni < sketchSize; i+=numSequences, uni++){
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
        double jaccardEstimate = std::max(double(inter),1.0)/uni;
        d_mashDist[idx_] = std::min(1.0, fabs(log(2.0*jaccardEstimate/(1.0+jaccardEstimate))/kmerSize));
    } });
}

void mashDistConstructionSpecialIDDC(
    int rowId,
    uint64_t *d_hashListBackbone,
    uint64_t *d_hashListConst,
    double *d_mashDist,
    uint64_t kmerSize,
    uint64_t sketchSize,
    int backboneSize,
    int numToConstruct,
    int *d_id,
    int *d_leafMap)
{
    // int tx = threadIdx.x, bx = blockIdx.x, bs = blockDim.x, gs = gridDim.x;
    // int idx = tx+bx*bs;
    // if(idx>=numToConstruct) return;

    tbb::parallel_for(tbb::blocked_range<int>(0, numToConstruct), [&](tbb::blocked_range<int> range)
                      {
    for (int idx_ = range.begin(); idx_ < range.end(); ++idx_) {
    // for (int idx_=idx; idx_<numToConstruct; idx_+=bs*gs){
        if (idx_>=numToConstruct) return;
        
        int mapIdx = d_leafMap[idx_];
        idx_ = d_id[idx_];
        if(idx_==-1) return;
        int uni = 0, bPos = rowId, inter = 0;
        uint64_t aval, bval;
        int idx_const = idx_;
        if (idx_ > backboneSize) idx_const = mapIdx;
    
        for(int i=idx_const; uni < sketchSize; i+=backboneSize, uni++){
            if (idx_ > backboneSize) aval = d_hashListConst[i];
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
        double jaccardEstimate = std::max(double(inter),1.0)/uni;
        d_mashDist[idx_] = std::min(1.0, fabs(log(2.0*jaccardEstimate/(1.0+jaccardEstimate))/kmerSize));
    } });
}

void mashDistConstructionSpecialIDClusteringDC(
    int rowId,
    uint64_t *d_hashListBackbone,
    uint64_t *d_hashListConst,
    double *d_mashDist,
    uint64_t kmerSize,
    uint64_t sketchSize,
    int backboneSize,
    int numToConstruct,
    int *d_id,
    int *d_leafID_map)
{

    tbb::parallel_for(tbb::blocked_range<int>(0, numToConstruct), [&](tbb::blocked_range<int> range)
                      {
    for (int idx_ = range.begin(); idx_ < range.end(); ++idx_) {
        int idx_1 = d_id[idx_];
        if (idx_1 != -1) {
            int uni = 0, bPos = 0, inter = 0;
            uint64_t aval, bval;
            int idx_const = idx_1;
            for (int i=idx_const*sketchSize; i<idx_const*sketchSize+sketchSize; i++) {
                if (idx_1 > backboneSize) aval = d_hashListConst[i];
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
            double jaccardEstimate = std::max(double(inter),1.0)/uni;
            d_mashDist[idx_1] = std::min(1.0, fabs(log(2.0*jaccardEstimate/(1.0+jaccardEstimate))/kmerSize));
        }
    } });
    /*
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
    */
}

void MashPlacement::MashDeviceArraysDC::distConstructionOnGpuForBackboneDC(Param &params, int rowId, double *d_mashDist) const
{
    // int threadNum = 256, blockNum = (this->backboneSize+threadNum-1)/threadNum;
    int threadNum = 1024, blockNum = 1024;
    mashDistConstructionDC(
        rowId,
        this->d_hashListBackbone,
        d_mashDist,
        params.kmerSize,
        params.sketchSize,
        this->backboneSize);
}

void MashPlacement::MashDeviceArraysDC::distConstructionOnGpuDC(Param &params, int rowId, double *d_mashDist) const
{
    int threadNum = 256, blockNum = (this->backboneSize + threadNum - 1) / threadNum;
    mashDistConstructionDC(
        rowId,
        d_hashList,
        d_mashDist,
        params.kmerSize,
        params.sketchSize,
        this->backboneSize);
}

void MashPlacement::MashDeviceArraysDC::distRangeConstructionOnGpuDC(Param &params, int rowId, double *d_mashDist, int l, int r, bool clustering) const
{

    // int threadNum = 1024, blockNum = (r-l+1+threadNum-1)/threadNum;
    int threadNum = 1024, blockNum = 1024;
    if (!clustering)
    {
        mashDistConstructionRangeDC (
            rowId,
            d_hashListBackbone,
            d_mashDist,
            params.kmerSize,
            params.sketchSize,
            this->backboneSize,
            l,
            r);
    }
    else
    {
        mashDistConstructionRangeForClusteringDC (
            rowId,
            d_hashListBackbone,
            d_hashListConst,
            d_mashDist,
            params.kmerSize,
            params.sketchSize,
            this->backboneSize,
            l,
            r);
    }
}

void MashPlacement::MashDeviceArraysDC::distSpecialIDConstructionOnGpuDC(
    Param &params,
    int rowId,
    double *d_mashDist,
    int numToConstruct,
    int *d_id,
    int *d_leafMap) const
{

    // int threadNum = 256, blockNum = (numToConstruct+threadNum-1)/threadNum;
    int threadNum = 1024, blockNum = 1024;
    mashDistConstructionSpecialIDDC (
        rowId,
        d_hashListBackbone,
        d_hashListConst,
        d_mashDist,
        params.kmerSize,
        params.sketchSize,
        backboneSize,
        numToConstruct,
        d_id,
        d_leafMap);
    // cudaDeviceSynchronize();
}

void MashPlacement::MashDeviceArraysDC::printSketchValuesDC(int numValues)
{
    uint64_t *h_hashList = new uint64_t[1000 * totalNumSequences];

    uint64_t *hashList = d_hashList;

    // cudaError_t err;

    // printf("Total Hashes: %d", numSequences*1000);

    for (int i = 0; i < totalNumSequences * 1000; ++i)
        h_hashList[i] = hashList[i];

    // err = cudaMemcpy(h_hashList, hashList, totalNumSequences*1000*sizeof(uint64_t), cudaMemcpyDeviceToHost);
    // if (err != cudaSuccess) {
    //     fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
    //     exit(1);
    // }

    // printf("i\thashList[i] (%zu)\n");
    for (int j = 0; j < totalNumSequences; j++)
    {
        fprintf(stderr, "Sequence (%d)\n", j);
        for (int i = 0; i < numValues; i++)
        {
            fprintf(stderr, "%i\t%lu\n", i, h_hashList[i * totalNumSequences + j]);
        }
    }
}