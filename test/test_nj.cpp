#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <vector>
#include <boost/program_options.hpp> 

#ifndef NJ_CUH
#include "../src/neighbourJoining.cuh"
#endif


namespace po = boost::program_options;
int main(int argc, char** argv) {
    std::string matrixFileName;
    int verbosity = 0;
    int upperTriangle = 0;
    po::options_description desc{"Options"};
    desc.add_options()
    ("matrixFileName,i", po::value<std::string>(&matrixFileName)->required(), "Input PHYLIP format distance matrix file name [REQUIRED]")
    ("verbose,v", po::value<int>(&verbosity)->implicit_value(1)->default_value(0), "Enable verbose mode.")
    ("help,h", "Print help messages");
    //Takes a PHYLIP format full distance matrix as input
    po::options_description allOptions;
    allOptions.add(desc);
    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).options(allOptions).run(), vm);
        po::notify(vm);
    } catch(std::exception &e) {
        std::cerr << desc << std::endl;
        exit(1);
    }
    freopen(matrixFileName.c_str(),"r",stdin);
    uint32_t numSequences;
    std::cin>>numSequences;
    auto ImportStart = std::chrono::high_resolution_clock::now();
    std::vector <std::string> name(numSequences);
    double *dismatrix = new double[numSequences*(numSequences-1)/2];
    for(int i=0,id;i<numSequences;i++){
        std::cin>>name[i];
        for(int j=0;j<numSequences;j++){
            double temp;
            scanf("%lf",&temp);
            if(j>i) dismatrix[id++]=temp;
            // In order to match with other applications,
            // only upper-triangle matrix is passed to the neighbourJoining namespace
        }
    }
    auto ImportEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds ImportTime = ImportEnd - ImportStart;
    if(verbosity) std::cerr << "Import distance matrix in: "<< ImportTime.count() <<" ns\n";
    // Output the time on loading distance matrix from file
    auto njStart = std::chrono::high_resolution_clock::now();
    neighbourJoining::deviceArrays.inputDismatrix(numSequences, dismatrix);
    //If having dismatrix available on GPU, use allocateDeviceArrays instead of inputDismatrix
    //Note that dismatrix has to be upper-triangle matrix
    auto njEnd = std:: chrono::high_resolution_clock::now();
    std::chrono::nanoseconds njTime = njEnd - njStart;
    if(verbosity) std::cerr << "Allocate memory in: "<< njTime.count() <<" ns\n";
    // Output the time on allocating memory
    neighbourJoining::findNeighbourJoiningTree(
        neighbourJoining::deviceArrays.d_numSequences,
        neighbourJoining::deviceArrays.d_mashDist,
        neighbourJoining::deviceArrays.d_U,
        neighbourJoining::deviceArrays.d_oriMashDist,
        name
    );
    njEnd = std:: chrono::high_resolution_clock::now();
    njTime = njEnd - njStart;
    if(verbosity) std::cerr << "Neighbor Joining in: " <<  njTime.count() << " ns\n";
    // Output the time on calculating Neighbour Joining tree
    neighbourJoining::deviceArrays.deallocateDeviceArrays();
    return 0;
}
