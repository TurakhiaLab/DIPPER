/**
 * mash_placement.cuh
 *
 * Declarations for DIPPER MASH/placement pipeline: Param, MashDeviceArrays,
 * MSADeviceArrays, MatrixReader, PlacementDeviceArrays (exact), KPlacementDeviceArrays
 * (k-closest), NJDeviceArrays, and divide-and-conquer variants (DC structs).
 */

#ifndef MASHPL_CUH
#define MASHPL_CUH

#include <stdint.h>
#include <cmath>
#include <functional>
#include <iostream>
#include <vector>
#include <cstdio>
#include <string>
#include "tree.hpp"

namespace MashPlacement
{
    /** Default false: print multifurcating Newick by collapsing ~zero-length internal branches when present.
     *  Set true (e.g. --print-binary-newick) to print the fully resolved binary tree with explicit :0 edges. */
    inline bool g_printBinaryNewick = false;

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
        uint64_t totalNumSeqs;
        std::pair<int, int> range;
        bool isProtein=false;

        Param(uint64_t t_kmerSize, uint64_t t_sketchSize, uint64_t t_threshold, uint64_t t_distanceType, std::string t_in, std::string t_out)
        {
            kmerSize = t_kmerSize; sketchSize = t_sketchSize; threshold = t_threshold, distanceType=t_distanceType;
            in = t_in, out = t_out; 
        };
    };

    /** Device arrays for MASH: compressed seqs, prefix/agg lengths, hash list, sketch construction. */
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

    /** MASH device arrays for divide-and-conquer (backbone vs full, leaf map, etc.). */
    struct MashDeviceArraysDC{
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
        
        void allocateDeviceArraysDC (uint64_t ** h_compressedSeqs, uint64_t * h_seqLengths, size_t num, Param& params);
        void printSketchValuesDC(int numValues);
        void sketchConstructionOnGpuDC (Param& params, uint64_t** twoBitCompressedSeqs, uint64_t * h_seqLengths, uint64_t numSequences);
        void deallocateDeviceArraysDC ();
        void distConstructionOnGpuDC(Param& params, int rowId, double* d_mashDist) const;
        void distConstructionOnCpu(Param& params, int rowId, double* d_mashDist) const;
        void distConstructionOnGpuForBackboneDC(Param& params, int rowId, double* d_mashDist) const;
        void distRangeConstructionOnGpuDC(Param& params, int rowId, double* d_mashDist, int l, int r, bool clustering = false) const;
        void distRangeConstructionOnCpu(Param& params, int rowId, double* d_mashDist, int l, int r) const;
        void distSpecialIDConstructionOnGpuDC(Param& params, 
                                            int rowId, 
                                            double* d_mashDist, 
                                            int numToConstruct, 
                                            int * d_id,
                                            int * d_leafMap) const;
        void distSpecialIDConstructionOnCpu(Param& params, int rowId, std::vector<double> & h_mashDist, int numToConstruct, std::vector<int> & h_id) const;

    };
    static MashDeviceArraysDC mashDeviceArraysDC;

    /** Device arrays for aligned sequences: 4-bit compressed MSA, seq lengths, distance construction. */
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

    /** MSA device arrays for DC (backbone vs const, seq len). */
    struct MSADeviceArraysDC{
        uint64_t * d_compressedSeqsBackBone;
        uint64_t * d_compressedSeqsConst;
        // uint64_t * d_seqLengths;
        // size_t numSequences;
        int d_seqLen;

        uint64_t * h_compressedSeqs=nullptr;

        size_t      totalNumSequences;
        size_t      backboneSize;

        void allocateDeviceArraysDC (uint64_t ** h_compressedSeqs, uint64_t * h_seqLengths, size_t num, Param& params, int gpuNum=0);
        void transferToDeviceArraysDC (uint64_t ** h_compressedSeqs, uint64_t * h_seqLengths, size_t num, int gpuCluster, Param& params);
        void distConstructionOnGpuDC(Param& params, int rowId, double* d_mashDist) const;
        void distConstructionOnGpuForBackboneDC(Param& params, int rowId, double* d_mashDist) const;
        void distRangeConstructionOnGpuDC(Param& params, int rowId, double* d_mashDist, int l, int r, bool clustering = false) const;
        void distSpecialIDConstructionOnGpuDC(Param& params, int rowId, double* d_mashDist, int numToConstruct, int* d_id, int * d_leafMap) const;
        void deallocateDeviceArraysDC ();

    };
    static MSADeviceArraysDC msaDeviceArraysDC;

    /** PHYLIP distance matrix reader and device copy. */
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

    /** Exact placement: adjacency, BFS/DFS/depth/levelst/leveled, lim; no k-closest. */
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
        void printTree(std::vector <std::string> name, std::ofstream& output_);
    };
    static PlacementDeviceArrays placementDeviceArrays;

    /** K-closest placement: closest_id/closest_dis per edge; build tree or add query to backbone. */
    struct KPlacementDeviceArrays{
        int idx, bd;
        int numSequences;
        int backboneSize;
        int * d_head;
        int * d_e;
        int * d_nxt;
        int * d_belong;
        int * d_closest_id;
        double * d_dist;
        double * d_len;
        double * d_closest_dis;

        void allocateDeviceArrays (size_t num, int backboneSize=-1);
        void deallocateDeviceArrays ();
        void initializeDeviceArrays(Tree* t);
        void findPlacementTree(
            Param& params,
            const MashDeviceArrays& mashDeviceArrays,
            MatrixReader& matrixReader,
            const MSADeviceArrays& msaDeviceArrays
        );

        void addQuery(
            Param& params,
            const MashDeviceArrays& mashDeviceArrays,
            MatrixReader& matrixReader,
            const MSADeviceArrays& msaDeviceArrays
        );
        void updateBranchLengthsFromTopology(
            Tree* t,
            double* h_mashDist,
            size_t numSeqs,
            uint64_t seed,
            Param& params,
            const MashDeviceArrays& mashDeviceArrays,
            MatrixReader& matrixReader,
            const MSADeviceArrays& msaDeviceArrays
        );
        void printTree(std::vector <std::string> name, std::ofstream& output_);
        void printTree(std::vector <std::string> name, std::ofstream& output_, UnrootedTree* t, std::vector<int>& edgeIdsMappingToNewTreeEdges);

        void estimateBranchLengthsFromTopology(
            Param& params,
            const MashDeviceArrays& mashDeviceArrays,
            MatrixReader& matrixReader,
            const MSADeviceArrays& msaDeviceArrays,
            std::vector<int>& edgeIdsMappingToNewTreeEdges,
            UnrootedTree* t,
            std::vector<std::string>& names
        );
    };
    static KPlacementDeviceArrays kplacementDeviceArrays;

    /** Conventional neighbor-joining: distance matrix and U matrix. */
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
        void findNeighbourJoiningTree(std::vector <std::string> &name, std::ofstream& output_);
    };
    static NJDeviceArrays njDeviceArrays;

    /** Host-side DC arrays for cluster trees and k-closest. */
    struct KPlacementDeviceArraysHostDC{
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

        void allocateHostArraysDC (size_t num, size_t totalNum);
        void deallocateHostArraysDC ();
        void findClusterTreeDC(
            Param& params,
            const MashDeviceArraysDC& mashDeviceArraysDC,
            MatrixReader& matrixReader,
            const MSADeviceArraysDC& msaDeviceArraysDC
        );
        void printTreeCpuDC(std::vector <std::string> name);
    };
    static KPlacementDeviceArraysHostDC kplacementDeviceArraysHostDC;

    /** Device-side DC k-closest placement: backbone, clusters, per-cluster trees. */
    struct KPlacementDeviceArraysDC{
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

        void allocateDeviceArraysDC (size_t num, size_t totalNum, int gpuNum=0);
        void deallocateDeviceArraysDC ();
        
        void findBackboneTreeDC(
            Param& params,
            const MashDeviceArraysDC& mashDeviceArrays,
            MatrixReader& matrixReader,
            const MSADeviceArraysDC& msaDeviceArrays,
            const KPlacementDeviceArraysHostDC& kplacementDeviceArraysHostDC,
            int gpuNum=0
        );

        void findClustersDC(
            Param& params,
            const MashDeviceArraysDC& mashDeviceArrays,
            MatrixReader& matrixReader,
            const MSADeviceArraysDC& msaDeviceArrays,
            KPlacementDeviceArraysHostDC& kplacementDeviceArraysHostDC
        );
        
        void findClustersDC_batch(
            Param& params,
            const MashDeviceArraysDC& mashDeviceArrays,
            MatrixReader& matrixReader,
            const MSADeviceArraysDC& msaDeviceArrays,
            KPlacementDeviceArraysHostDC& kplacementDeviceArraysHostDC,
            const int clusteringBatchIdx
        );

        void findClusterTreeDC(
            Param& params,
            MashDeviceArraysDC& mashDeviceArrays,
            MatrixReader& matrixReader,
            MSADeviceArraysDC& msaDeviceArrays,
            KPlacementDeviceArraysHostDC& kplacementDeviceArraysHostDC,
            std::vector<int>& largeClustersIdx
        );
        void findClusterTreeDC_batch(
            Param& params,
            MashDeviceArraysDC& mashDeviceArrays,
            MatrixReader& matrixReader,
            MSADeviceArraysDC& msaDeviceArrays,
            KPlacementDeviceArraysHostDC& kplacementDeviceArraysHostDC,
            const std::string dir,
            std::vector<bool>& isCluster
        );

        void findBackboneTreeDCRecursive(
            Param& params,
            MashDeviceArraysDC& mashDeviceArrays,
            MatrixReader& matrixReader,
            MSADeviceArraysDC& msaDeviceArrays,
            const KPlacementDeviceArraysHostDC& kplacementDeviceArraysHostDC,
            int clusterIdx
        );
        

        void printTreeDC(std::vector <std::string> name, std::ofstream& output_);
    };
    static KPlacementDeviceArraysDC kplacementDeviceArraysDC;

/** Write rooted Newick from host-side adjacency (h_head/h_e/h_nxt/h_len) used by GPU placement trees. */
inline void printNewickFromHostAdjacency(std::ostream& out, const std::vector<std::string>& name, int* h_head,
                                         int* h_e, int* h_nxt, double* h_len, int rootNode, bool binaryNewick) {
    const double kLenEps = 1e-12;
    auto hasKids = [h_head, h_e, h_nxt](int node, int from) -> bool {
        for (int i = h_head[node]; i != -1; i = h_nxt[i]) {
            if (h_e[i] != from) {
                return true;
            }
        }
        return false;
    };
    std::function<void(int, int)> dfs;
    dfs = [&](int node, int from) {
        if (!hasKids(node, from)) {
            out << name[node];
            return;
        }
        if (binaryNewick) {
            out << "(";
            std::vector<int> pos;
            for (int i = h_head[node]; i != -1; i = h_nxt[i]) {
                if (h_e[i] != from) {
                    pos.push_back(i);
                }
            }
            for (size_t i = 0; i < pos.size(); ++i) {
                dfs(h_e[pos[i]], node);
                out << ":" << h_len[pos[i]] << (i + 1 == pos.size() ? ")" : ",");
            }
            return;
        }
        std::vector<std::pair<int, double>> kids;
        std::function<void(int, int)> addCollapsed;
        addCollapsed = [&](int n, int fr) {
            for (int i = h_head[n]; i != -1; i = h_nxt[i]) {
                int ch = h_e[i];
                if (ch == fr) {
                    continue;
                }
                double L = h_len[i];
                if (std::fabs(L) <= kLenEps && hasKids(ch, n)) {
                    addCollapsed(ch, n);
                } else {
                    kids.push_back({ch, L});
                }
            }
        };
        addCollapsed(node, from);
        out << "(";
        for (size_t i = 0; i < kids.size(); ++i) {
            dfs(kids[i].first, node);
            out << ":" << kids[i].second << (i + 1 == kids.size() ? ")" : ",");
        }
    };
    dfs(rootNode, -1);
}

};



#endif