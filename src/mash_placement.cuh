#ifndef MASHPL_CUH
#define MASHPL_CUH

#include <stdint.h>
#include <iostream>
#include <vector>

// typedef uint64_t hash_t;


namespace MashPlacement
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
        uint64_t * d_prefixCompressed;
        uint64_t * d_aggseqLengths;
        uint64_t * d_seqLengths;
        size_t d_numSequences;
        uint64_t * d_hashList;
        double * d_mashDist;
        int numSequences;
        double * d_dist;
        int * d_head;
        int * d_e;
        int * d_nxt;
        int idx, bd;
        double * d_len;
        int * d_belong;
        double * d_closest_dis;
        int * d_closest_id;

        void allocateDeviceArrays (uint64_t ** h_compressedSeqs, uint64_t * h_seqLengths, size_t num, Param& params);
        void printSketchValues(int numValues);
        void printMashDist(uint64_t h_numSequences, std::vector<std::string> seqs);
        void deallocateDeviceArrays ();
    };

    static DeviceArrays deviceArrays;

    

    void sketchConstructionOnGpu
    (
        uint64_t * d_compressedSeqs,
        uint64_t * d_prefixCompressed,
        uint64_t * d_aggseqLengths,
        uint64_t * d_seqLengths,
        size_t d_numSequences,
        uint64_t * d_hashList,
        uint64_t * h_seqLengths,
        Param& params
    );

    void findPlacementTree(
        int numSequences,
        int bd, // id to place
        int idx, // id of linked-list position
        double * d_dist,
        int * d_head,
        int * d_e,
        double * d_len,
        int * d_nxt,
        int * d_belong,
        double * d_closest_dis,
        int * d_closest_id,
        std::vector <std::string> name,
        uint64_t * d_hashList,
        uint64_t * d_seqLengths,
        Param& params
    );


};

#endif