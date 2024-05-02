#ifndef MASH_CUH
#define MASH_CUH

#include <stdint.h>
#include <iostream>
#include <vector>

// typedef uint64_t hash_t;


namespace GpuSketch
{
    struct Param
    {
        uint64_t kmerSize;
        uint64_t sketchSize;
        uint64_t threshold;
        int numBlocks;
        int blockSize;

        Param(uint64_t t_kmerSize, uint64_t t_sketchSize, uint64_t t_threshold, int t_numBlocks, int t_blockSize)
        {
            kmerSize = t_kmerSize; sketchSize = t_sketchSize; threshold = t_threshold; numBlocks = t_numBlocks; blockSize = t_blockSize;
        };
    };

    struct DeviceArrays
    {
        char * d_seqs;
        uint64_t * d_aggrSeqsLen;
        uint64_t * d_seqsLen;
        size_t d_numSequences;
        uint64_t * d_hashList;
        uint64_t * d_hashListPruned;
        float * d_mashDist;

        void allocateDeviceArrays (char ** h_seqs, uint64_t ** h_seqsLen, size_t numSequences, Param& params);
        void printSketchValues(uint64_t numValues);
        void printMashDist(uint64_t h_numSequences, std::vector<std::string>seqs);
        void deallocateDeviceArrays ();
    };

    static DeviceArrays deviceArrays;

    

    void sketchConstructionOnGpu
    (
        char * d_seqs,
        uint64_t * d_aggrSeqsLen,
        uint64_t * d_seqsLen,
        size_t d_numSequences,
        uint64_t * d_hashList,
        uint64_t * d_hashListPruned,
        Param& params,
        uint64_t seed
    );

};

#endif