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
        char **seqName;
        size_t d_numSequences;
        float * d_mashDist;

        void allocateDeviceArrays (size_t , GpuSketch::DeviceArrays);
        void deallocateDeviceArrays ();
    };

    static DeviceArrays deviceArrays;


};

#endif