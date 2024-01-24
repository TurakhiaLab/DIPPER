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

#ifndef MSA_HPP
#include "../src/msa.hpp"
#endif

#ifndef CREATESKETCH_HPP
#include "../src/createSketch.hpp"
#endif

#ifndef MASHDISTANCE_HPP
#include "../src/mashDistance.hpp"
#endif

#ifndef DISTANCEMATRIX_HPP
#include "../src/distanceMatrix.hpp"
#endif

#ifndef TWOBITCOMPRESSOR_HPP
#include "../src/twoBitCompressor.hpp"
#endif

#include "../src/neighbourJoining.cpp"

namespace po = boost::program_options;

KSEQ_INIT2(, gzFile, gzread)

po::options_description mainDesc("MSA Command Line Arguments");


void parseArguments(int argc, char** argv)
{
    // Setup boost::program_options
    mainDesc.add_options()
        ("tree,t", po::value<std::string>()->required(), "Initial Tree - Newick format (required)")
        ("sequences,f", po::value<std::string>()->required(), "Tip sequences - Fasta format (required)")
        ("kmer-size,k", po::value<std::string>(), "Kmer-size (Default = 16)")
        ("sketch-size,s", po::value<std::string>(), "Sketch-size (Default = 10000)")
        ("threshold,r", po::value<std::string>(), "Erroneous k-mer threshold (Default = 1)")
        ("help,h", "Print help messages");

}


void readSequences(po::variables_map& vm, msa::utility* util)
{
    auto seqReadStart = std::chrono::high_resolution_clock::now();
    std::string seqFileName = vm["sequences"].as<std::string>();

    gzFile f_rd = gzopen(seqFileName.c_str(), "r");
    if (!f_rd) {
        fprintf(stderr, "ERROR: cant open file: %s\n", seqFileName.c_str());
        exit(1);
    }

    kseq_t* kseq_rd = kseq_init(f_rd);

    while (kseq_read(kseq_rd) >= 0) {
        size_t seqLen = kseq_rd->seq.l;
        util->seqs[kseq_rd->name.s] = std::string(kseq_rd->seq.s, seqLen);
    }
    
    auto seqReadEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds seqReadTime = seqReadEnd - seqReadStart;
    std::cout << "Sequences read in: " <<  seqReadTime.count() << " ns\n";
}

void checkAlignment(std::vector<std::string>& ref)
{
    size_t len = 0;
    bool set = false;

    for (auto &r: ref)
    {
        if (!set) len = r.size();
        if (r.size() != len)
        {
            fprintf(stderr, "Error: Alignment Size do not match\n");
        }
    }
}



int main(int argc, char** argv) {

    std::string seqFileName;
    

    parseArguments(argc, argv);

    po::variables_map vm;
    try{
        po::store(po::command_line_parser(argc, argv).options(mainDesc).run(), vm);
        po::notify(vm);
    }
    catch(std::exception &e){
        std::cerr << mainDesc << std::endl;
        // Return with error code 1 unless the user specifies help
        if(vm.count("help"))
            return 0;
        else
            return 1;
    }

    // Kmer Size
    int k = 16;
    try {k= std::stoi(vm["kmer-size"].as<std::string>());}
    catch(std::exception &e){}

    // Sketch Size
    size_t sketchSize = 1000;
    try {sketchSize= (size_t)std::stoi(vm["sketch-size"].as<std::string>());}
    catch(std::exception &e){}

    // Erroneous k-mer thresold
    int threshold = 1;
    try {threshold= std::stoi(vm["threshold"].as<std::string>());}
    catch(std::exception &e){}

    std::cout << "kmer-size: " << k << "\nSketch-size: " << sketchSize << "\nErroneous k-mer threshold: " << threshold << "\n";

    // Define MSA utility
    msa::utility* util = new msa::utility;

    // Read Input Sequences (Fasta format)
    readSequences(vm, util);

    sketch* sketchLocal  = new sketch(util->seqs, k, sketchSize, threshold);

    std::vector<float> distMatrix;
    std::vector<std::string> distMatrixSeqOrder;

    distanceMatrix(sketchLocal->sketchMap, k, sketchSize, distMatrix, distMatrixSeqOrder);


    // int count = 0;
    // for (size_t i=0; i<sketchLocal->sketchMap.size(); i++)
    // {
    //     for (size_t j=0; j<sketchLocal->sketchMap.size(); j++)
    //     {
    //         float d = 0;
    //         if (j>i) {d = distMatrix[count];count++;}
    //         std::cout << d << "\t";
    //     }
    //     std::cout << "\n";
    // }

    findNeighbourJoiningTree(distMatrix, sketchLocal->sketchMap.size(), distMatrixSeqOrder);

    // for (auto& s: sketchLocal->sketchMap)
    // {
    //     std::cout<<s.first<<"\t";
    //     while(!s.second.empty())
    //     {
    //         std::cout <<s.second.top()<<"\t";
    //         s.second.pop();
    //     }
    //     std::cout<<"\n";
    // }
    
    // const char* in = "0";
    // std::vector<std::string> in_string = {"A", "AA", "AAA", "AAAAAAAAAAAAAAAA"};
    // int len = std::strlen(in);
    // char out[4];
    // uint32_t out_help;
    // uint32_t seed = 54;
    // uint32_t * in_help;

    // MurmurHash3_x86_32  ( in, len, seed, out);
    // std::cout << in << ": " << *((uint32_t *)out) << "\n";

    // for (auto &s: in_string)
    // {
    //     uint32_t compLen = (s.size()+15)/16;
    //     uint32_t * twoBitCompressed = new uint32_t[compLen];
    //     twoBitCompressor(s, len, twoBitCompressed);
    //     out_help = MurmurHash3_x86_32_help  ( twoBitCompressed, s.size(), seed);
    //     std::cout << s << ": " << (out_help) << "\n";
    // }

    return 0;
}

