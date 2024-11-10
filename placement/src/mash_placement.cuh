#ifndef MASHPL_CUH
#define MASHPL_CUH

#include <stdint.h>
#include <iostream>
#include <vector>
#include <cstdio>
#include <string>

// typedef uint64_t hash_t;

void clusterKernelWrapper(int *clusterMap, int numSequences, uint64_t ** twoBitCompressedSeqs, int MAX_LEVELS, int  *d_stopFlag, int *d_sharedCount,
    int *d_clusterMap,int *d_cInstr,int nodesInThisLevel,int *stopFlag);

namespace MashPlacement
{
    struct Param
    {
        uint64_t kmerSize;
        uint64_t sketchSize;
        uint64_t threshold;
        uint64_t distanceType;
        std::string in;
        std::string out;

        Param(uint64_t t_kmerSize, uint64_t t_sketchSize, uint64_t t_threshold, uint64_t t_distanceType, std::string t_in, std::string t_out)
        {
            kmerSize = t_kmerSize; sketchSize = t_sketchSize; threshold = t_threshold, distanceType=t_distanceType;
            in = t_in, out = t_out;
        };
    };
    
    struct MashDeviceArrays{
        uint64_t * d_compressedSeqs;
        uint64_t * d_prefixCompressed;
        uint64_t * d_aggseqLengths;
        uint64_t * d_seqLengths;
        size_t numSequences;
        uint64_t * d_hashList;

        uint64_t * h_hashList;

        void allocateDeviceArrays (uint64_t ** h_compressedSeqs, uint64_t * h_seqLengths, size_t num, Param& params);
        void printSketchValues(int numValues);
        void sketchConstructionOnGpu (Param& params);
        void deallocateDeviceArrays ();
        void distConstructionOnGpu(Param& params, int rowId, double* d_mashDist) const;

    };
    static MashDeviceArrays mashDeviceArrays;


    struct MSADeviceArrays{
        uint64_t * d_compressedSeqs;
        uint64_t * d_seqLengths;
        size_t numSequences;
        int seqLen;

        void allocateDeviceArrays (uint64_t ** h_compressedSeqs, uint64_t * h_seqLengths, size_t num, Param& params);
        void deallocateDeviceArrays ();
        void distConstructionOnGpu(Param& params, int rowId, double* d_mashDist) const;

    };
    static MSADeviceArrays msaDeviceArrays;


    struct MatrixReader{
        int numSequences;
        double * h_dist;
        std::vector <std::string> name;
        FILE* filePtr;
        char * buffer;
        void allocateDeviceArrays (int num, FILE* fPtr);
        void distConstructionOnGpu(Param& params, int rowId, double* d_dist);
    };
    static MatrixReader matrixReader;


    struct PlacementDeviceArrays{
        int idx, bd;
        int numSequences;
        int * d_head;
        int * d_e;
        int * d_nxt;
        int * d_belong;
        int * d_bfsorder;
        int * d_dfsorder;
        int * d_dep;
        int * d_dfsrk;
        int * d_levelst;
        int * d_leveled;
        int * d_rev;
        double * d_dist;
        double * d_len;
        double * d_lim;

        void allocateDeviceArrays (size_t num);
        void deallocateDeviceArrays ();
        void findPlacementTree(
            Param& params,
            const MashDeviceArrays& mashDeviceArrays,
            MatrixReader& matrixReader,
            const MSADeviceArrays& msaDeviceArrays
        );
        void printTree(std::vector <std::string> name);
    };
    static PlacementDeviceArrays placementDeviceArrays;

};

#endif