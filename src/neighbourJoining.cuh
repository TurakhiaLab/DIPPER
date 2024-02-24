#ifndef NJ_CUH
#define NJ_CUH

#include <stdint.h>
#include <iostream>
#include <vector>


#ifndef MASH_CUH
#include "../src/mash.cuh"
#endif

// typedef uint32_t hash_t;

namespace neighbourJoining
{
    
    struct DeviceArrays
    {
        //char **seqName;
        uint32_t d_numSequences;
        double * d_mashDist;
        double * d_U;
        uint32_t * d_flag;
        double * d_oriMashDist;
        uint32_t * d_rowID;
        void allocateDeviceArrays (uint32_t numSequences, double *mashDist);
        void deallocateDeviceArrays ();
        void inputDismatrix(uint32_t numSequences, double *mashDist);
    };

    static DeviceArrays deviceArrays;
    void findNeighbourJoiningTree(uint32_t d_numSequences, double *d_mashDist, double *d_U, uint32_t *d_flag, double *d_oriMashDist, uint32_t *d_rowID);
};

#endif
