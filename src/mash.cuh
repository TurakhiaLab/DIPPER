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
        uint64_t * d_compressedSeqs;
        uint64_t * d_aggseqLengths;
        uint64_t * d_seqLengths;
        size_t d_numSequences;
        uint64_t * d_hashList;
        float * d_mashDist;

        void allocateDeviceArrays (uint64_t ** h_compressedSeqs, uint64_t * h_seqLengths, size_t numSequences, Param& params);
        void printSketchValues(int numValues);
        void printMashDist(uint64_t h_numSequences);
        void deallocateDeviceArrays ();
    };

    static DeviceArrays deviceArrays;

    

    void sketchConstructionOnGpu
    (
        uint64_t * d_compressedSeqs,
        uint64_t * d_aggseqLengths,
        uint64_t * d_seqLengths,
        size_t d_numSequences,
        uint64_t * d_hashList,
        uint64_t * h_seqLengths,
        Param& params
    );

    void mashDistConstructionOnGpu
    (
        uint64_t * d_hashList,
        uint64_t * d_seqLengths,
        size_t d_numSequences,
        float * d_mashDist,
        uint64_t * h_seqLengths,
        Param& params
    );

};

#endif