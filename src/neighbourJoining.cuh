#ifndef NJ_CUH
#define NJ_CUH

#include <stdint.h>
#include <iostream>
#include <vector>
#include <string>

// typedef uint32_t hash_t;

namespace neighbourJoining
{
    
    struct DeviceArrays
    {
        //char **seqName;
        uint32_t d_numSequences;
        double * d_mashDist;
        double * d_U;
        double * d_oriMashDist;
        void deallocateDeviceArrays ();
        void inputDismatrix(uint32_t numSequences, double *mashDist);// mashDist on CPU
        void allocateDeviceArrays(uint32_t numSequences, double *mashDist);// mashDist on GPU
    };

    static DeviceArrays deviceArrays;
    void findNeighbourJoiningTree(uint32_t d_numSequences, double *d_mashDist, double *d_U, double *d_oriMashDist, std::vector <std::string> &name);
};

#endif