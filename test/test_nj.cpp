#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <boost/program_options.hpp> 
#include "../src/kseq.h"
#include "zlib.h"


#ifndef MURMURHASH3_HPP
#include "../src/murmurHash3.hpp"
#endif

#ifndef TWOBITCOMPRESSOR_HPP
#include "../src/twoBitCompressor.hpp"
#endif

#ifndef MASH_CUH
#include "../src/mash.cuh"
#endif

#ifndef MASHDISTANCE_HPP
#include "../src/mashDistance.hpp"
#endif

#ifndef DISTANCEMATRIX_HPP
#include "../src/distanceMatrix.hpp"
#endif

#ifndef NJ_CUH
#include "../src/neighbourJoining.cuh"
#endif


namespace po = boost::program_options;

KSEQ_INIT2(, gzFile, gzread)

int main() {
    uint32_t numSequences;
    std::cin>>numSequences;
    double *dismatrix = new double[numSequences*(numSequences-1)/2];
    std::string seqname;
    for(int i=0,id;i<numSequences;i++){
        std::cin>>seqname;
        for(int j=0;j<numSequences;j++){
            double temp;
            scanf("%lf",&temp);
            if(j>i) dismatrix[id++]=temp;
        }
    }
    //for(int i=0;i<numSequences*(numSequences-1)/2;i++) std::cout<<dismatrix[i]<<"\n";
    auto njStart = std::chrono::high_resolution_clock::now();
    neighbourJoining::deviceArrays.inputDismatrix(numSequences, dismatrix);
    auto njEnd = std:: chrono::high_resolution_clock::now();
    std::chrono::nanoseconds njTime = njEnd - njStart;
    std::cout << "Allocate memory in: "<< njTime.count() <<" ns\n";
    neighbourJoining::findNeighbourJoiningTree(
        neighbourJoining::deviceArrays.d_numSequences,
        neighbourJoining::deviceArrays.d_mashDist,
        neighbourJoining::deviceArrays.d_U,
        neighbourJoining::deviceArrays.d_flag,
        neighbourJoining::deviceArrays.d_oriMashDist,
        neighbourJoining::deviceArrays.d_rowID
    );
    njEnd = std:: chrono::high_resolution_clock::now();
    njTime = njEnd - njStart;
    std::cout << "Neighbor Joining in: " <<  njTime.count() << " ns\n";
    //GpuSketch::deviceArrays.deallocateDeviceArrays();
    return 0;
}

