#ifndef PL_CUH
#define PL_CUH

#include <stdint.h>
#include <iostream>
#include <vector>
#include <string>

namespace placement
{
    struct DeviceArrays
    {
        int numSequences;
        int bd;
        int tk;
        double * d_dist;
        int * d_head;
        int * d_e;
        int * d_nxt;
        int idx;
        double * d_len;
        int * d_belong;
        double * d_closest_dis;
        int * d_closest_id;
        void deallocateDeviceArrays ();
        void allocateDeviceArrays(
            int k,
            int num, 
            double *dist, 
            int *head, 
            int *e,
            int *nxt,
            int cnt,
            int cur,
            double * len,
            int * belong
        );
    };
    
    static DeviceArrays deviceArrays;
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
        int tk
    );
};


#endif