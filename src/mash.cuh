#ifndef MASH_CUH
#define MASH_CUH

#include <stdint.h>
#include <iostream>
#include <vector>

// typedef uint32_t hash_t;


namespace GpuSketch
{
    struct Param
    {
        uint32_t kmerSize;
        uint32_t sketchSize;
        uint32_t threshold;
        int numBlocks;
        int blockSize;

        Param(uint32_t t_kmerSize, uint32_t t_sketchSize, uint32_t t_threshold, int t_numBlocks, int t_blockSize)
        {
            kmerSize = t_kmerSize; sketchSize = t_sketchSize; threshold = t_threshold; numBlocks = t_numBlocks; blockSize = t_blockSize;
        };
    };

    struct DeviceArrays
    {
        uint32_t * d_compressedSeqs;
        uint32_t * d_aggseqLengths;
        uint32_t * d_seqLengths;
        size_t d_numSequences;
        uint32_t * d_hashList;
        uint32_t * d_hashListPruned;
        float * d_mashDist;

        void allocateDeviceArrays (uint32_t ** h_compressedSeqs, uint32_t * h_seqLengths, size_t numSequences, Param& params);
        void printSketchValues(int numValues, uint32_t * h_seqLengths);
        void printMashDist(uint32_t h_numSequences);
        void deallocateDeviceArrays ();
    };

    static DeviceArrays deviceArrays;

    

    void sketchConstructionOnGpu
    (
        uint32_t * d_compressedSeqs,
        uint32_t * d_aggseqLengths,
        uint32_t * d_seqLengths,
        size_t d_numSequences,
        uint32_t * d_hashList,
        uint32_t * d_hashListPruned,
        uint32_t * h_seqLengths,
        Param& params
    );

    void mashDistConstructionOnGpu
    (
        uint32_t * d_hashList,
        uint32_t * d_seqLengths,
        size_t d_numSequences,
        float * d_mashDist,
        uint32_t * h_seqLengths,
        Param& params
    );

};

#endif