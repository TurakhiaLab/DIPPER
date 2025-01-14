#ifndef RNJ_CUH
#define RNJ_CUH

#include <stdint.h>
#include <iostream>
#include <vector>
#include <string>

// typedef uint32_t hash_t;

namespace rapidNJ
{
    
    struct DeviceArrays
    {
        //char **seqName;
        uint32_t d_numSequences;
        double * d_mashDist;
        uint32_t * d_id;
        double * d_U;
        double * d_oriMashDist;
        double * d_globalMin;
        uint32_t * d_timeStamp;
        uint32_t * global_id;
        uint32_t * d_rowId;
        double * d_preMaxU;
        void deallocateDeviceArrays ();
        void inputDismatrix(uint32_t numSequences, double *mashDist);// mashDist on CPU
        void allocateDeviceArrays(uint32_t numSequences, double *mashDist);// mashDist on GPU
    };
    
    static DeviceArrays deviceArrays;
    void findNeighbourJoiningTree(uint32_t d_numSequences, double *d_mashDist, double *d_U, double *d_oriMashDist, uint32_t *d_id, std::vector <std::string> &name, double * globalMin, uint32_t *d_timeStamp, uint32_t * global_id, uint32_t * d_rowId, double * d_preMaxU);
};


#endif