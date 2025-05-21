#ifndef MASHPL_CUH
#define MASHPL_CUH

#include <stdint.h>
#include <iostream>
#include <vector>
#include <cstdio>
#include <string>

// typedef uint64_t hash_t;


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
        uint64_t batchSize;
        uint64_t backboneSize;

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
        uint64_t * d_hashList;
        uint64_t * d_hashListConst;
        uint64_t * d_hashListBackbone;
        
        uint64_t * h_hashList;

        size_t      totalNumSequences;
        size_t      backboneSize;

        int * d_leafID_map;
        
        void allocateDeviceArrays (uint64_t ** h_compressedSeqs, uint64_t * h_seqLengths, size_t num, Param& params);
        void printSketchValues(int numValues);
        void sketchConstructionOnGpu (Param& params, uint64_t** twoBitCompressedSeqs, uint64_t * h_seqLengths, uint64_t numSequences);
        void deallocateDeviceArrays ();
        void distConstructionOnGpu(Param& params, int rowId, double* d_mashDist) const;
        void distConstructionOnCpu(Param& params, int rowId, double* d_mashDist) const;
        void distConstructionOnGpuForBackbone(Param& params, int rowId, double* d_mashDist) const;
        void distRangeConstructionOnGpu(Param& params, int rowId, double* d_mashDist, int l, int r, bool clustering = false) const;
        void distRangeConstructionOnCpu(Param& params, int rowId, double* d_mashDist, int l, int r) const;
        void distSpecialIDConstructionOnGpu(Param& params, 
                                            int rowId, 
                                            double* d_mashDist, 
                                            int numToConstruct, 
                                            int * d_id,
                                            int * d_leafMap) const;
        void distSpecialIDConstructionOnCpu(Param& params, int rowId, std::vector<double> & h_mashDist, int numToConstruct, std::vector<int> & h_id) const;

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

    struct KPlacementDeviceArraysHost{
        int idx, bd;
        int numSequences;
        int totalNumSequences;
        int * h_head;
        int * h_e;
        int * h_nxt;
        int * h_belong;
        int * h_closest_id;
        int * h_closest_id_cluster;
        double * h_dist;
        double * h_len;
        double * h_closest_dis;
        double * h_closest_dis_cluster;
        int * clusterID;

        void allocateHostArrays (size_t num, size_t totalNum);
        void deallocateHostArrays ();
        void findClusterTree(
            Param& params,
            const MashDeviceArrays& mashDeviceArrays,
            MatrixReader& matrixReader,
            const MSADeviceArrays& msaDeviceArrays
        );
        void printTreeCpu(std::vector <std::string> name);
    };
    static KPlacementDeviceArraysHost kplacementDeviceArraysHost;

    struct NJDeviceArrays
    {
        int d_numSequences;
        double * d_mashDist;
        double * d_U;
        void deallocateDeviceArrays ();
        void getDismatrix(
            int numSequences,
            Param& params,
            const MashDeviceArrays& mashDeviceArrays,
            MatrixReader& matrixReader,
            const MSADeviceArrays& msaDeviceArrays
        );
        void findNeighbourJoiningTree(std::vector <std::string> &name);
    };
    static NJDeviceArrays njDeviceArrays;

    struct KPlacementDeviceArrays{
        int idx, bd;
        int numSequences;
        int totalNumSequences;
        int * d_head;
        int * d_e;
        int * d_nxt;
        int * d_belong;
        int * d_closest_id;
        int * d_closest_id_cluster;
        double * d_dist;
        double * d_len;
        double * d_closest_dis;
        double * d_closest_dis_cluster;

        void allocateDeviceArrays (size_t num, size_t totalNum);
        void deallocateDeviceArrays ();
        
        void findBackboneTree(
            Param& params,
            const MashDeviceArrays& mashDeviceArrays,
            MatrixReader& matrixReader,
            const MSADeviceArrays& msaDeviceArrays,
            const KPlacementDeviceArraysHost& kplacementDeviceArraysHost
        );

        void findClusters(
            Param& params,
            const MashDeviceArrays& mashDeviceArrays,
            MatrixReader& matrixReader,
            const MSADeviceArrays& msaDeviceArrays,
            KPlacementDeviceArraysHost& kplacementDeviceArraysHost
        );

        void findClusterTree(
            Param& params,
            MashDeviceArrays& mashDeviceArrays,
            MatrixReader& matrixReader,
            const MSADeviceArrays& msaDeviceArrays,
            KPlacementDeviceArraysHost& kplacementDeviceArraysHost
        );
        void printTree(std::vector <std::string> name);
    };
    static KPlacementDeviceArrays kplacementDeviceArrays;

};

#endif