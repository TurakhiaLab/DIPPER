/**
 * cpu/divide_and_conquer/mash.cpp
 *
 * CPU MASH for DC: distRangeConstructionOnCpu (MASH distances for rows [l,r]),
 * distSpecialIDConstructionOnCpu (distances for a set of IDs). Used by
 * placement_close_k DC flow.
 */

#include "../mash_placement.hpp"

#include <stdio.h>
#include <queue>
#include <chrono>
#include <iostream>
#include <cmath>
#include <tbb/parallel_for.h>

/** MASH distances from row rowId to rows [st,ed] (TBB over range). */
void mashDistConstructionRangeCpu(
    int rowId,
    uint64_t *h_hashList,
    double *h_mashDist,
    uint64_t kmerSize,
    uint64_t sketchSize,
    int numSequences,
    int st,
    int ed)
{
    tbb::parallel_for(tbb::blocked_range<int>(st, ed + 1), [&](tbb::blocked_range<int> range)
                      {
    for (int idx = range.begin(); idx < range.end(); ++idx) {
        int uni = 0, bPos = rowId*sketchSize, inter = 0;
        uint64_t aval, bval;
        for(int i=idx*sketchSize; uni < sketchSize; i++, uni++){
            aval = h_hashList[i];
            while(uni < sketchSize && bPos < rowId*sketchSize + sketchSize){
                bval = h_hashList[bPos];
                // printf("%ull %ull\n",aval,bval);
                if(bval > aval) break;
                if(bval < aval) uni++;
                else inter++;
                bPos ++;
            }
            if(uni >= sketchSize) break;
        }
        double jaccardEstimate = std::max(double(inter),1.0)/uni;
        double temp = fabs(std::log(2.0*jaccardEstimate/(1.0+jaccardEstimate))/kmerSize);
        h_mashDist[idx] = std::min(1.0, temp);
        // printf("idx: %d, mashDist: %lf\n", idx, h_mashDist[idx]);
    } });
}

/** MASH distances from row rowId to rows in h_id; output into h_mashDist. */
void mashDistConstructionSpecialIDCpu(
    int rowId,
    uint64_t *h_hashList,
    // double * h_mashDist,
    std::vector<double> &h_mashDist,
    uint64_t kmerSize,
    uint64_t sketchSize,
    int totalNumSequences,
    int numToConstruct,
    std::vector<int> &h_id
    // int * h_id
)
{
    for (int idx_ = 0; idx_ < numToConstruct; idx_++)
    {
        // tbb::parallel_for(tbb::blocked_range<int>(0, numToConstruct), [&](tbb::blocked_range<int> range){
        // for (int idx_ = range.begin(); idx_ < range.end(); ++idx_) {
        int idx = h_id[idx_];
        if (idx == -1)
            continue;
        int uni = 0, bPos = rowId * sketchSize, inter = 0;
        uint64_t aval, bval;

        for (int i = idx * sketchSize; uni < sketchSize; i++, uni++)
        {
            aval = h_hashList[i];
            while (uni < sketchSize && bPos < rowId * sketchSize + sketchSize)
            {
                bval = h_hashList[bPos];
                if (bval > aval)
                    break;
                if (bval < aval)
                    uni++;
                else
                    inter++;
                bPos++;
            }
            if (uni >= sketchSize)
                break;
        }

        double jaccardEstimate = std::max(double(inter), 1.0) / uni;
        double temp = fabs(std::log(2.0 * jaccardEstimate / (1.0 + jaccardEstimate)) / kmerSize);
        h_mashDist[idx] = std::min(1.0, temp);
        // }});
    }
}

/** Dispatch mashDistConstructionRangeCpu for [l,r]. */
void MashPlacement::MashDeviceArraysDC::distRangeConstructionOnCpu(Param &params, int rowId, double *h_mashDist, int l, int r) const
{
    mashDistConstructionRangeCpu(
        rowId,
        h_hashList,
        h_mashDist,
        params.kmerSize,
        params.sketchSize,
        this->totalNumSequences,
        l,
        r);
}

/** Dispatch mashDistConstructionSpecialIDCpu for h_id. */
void MashPlacement::MashDeviceArraysDC::distSpecialIDConstructionOnCpu(Param &params, int rowId, std::vector<double> &h_mashDist, int numToConstruct, std::vector<int> &h_id) const
{

    mashDistConstructionSpecialIDCpu(
        rowId,
        h_hashList,
        h_mashDist,
        params.kmerSize,
        params.sketchSize,
        this->totalNumSequences,
        numToConstruct,
        h_id);
}