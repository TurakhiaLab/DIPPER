#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <bits/stdc++.h>
#include <boost/program_options.hpp> 
#include "../src/kseq.h"
#include "zlib.h"


#ifndef MURMURHASH3_HPP
#include "../src/murmurHash3.hpp"
#endif

#ifndef TWOBITCOMPRESSOR_HPP
#include "../src/twoBitCompressor.hpp"
#endif

#ifndef MASHPL_CUH
#include "../src/mash_placement.cuh"
#endif

#ifndef MASHDISTANCE_HPP
#include "../src/mashDistance.hpp"
#endif

#ifndef DISTANCEMATRIX_HPP
#include "../src/distanceMatrix.hpp"
#endif

namespace po = boost::program_options;

KSEQ_INIT2(, gzFile, gzread)

po::options_description mainDesc("MSA Command Line Arguments");


void parseArguments(int argc, char** argv)
{
    // Setup boost::program_options
    mainDesc.add_options()
        // ("tree,t", po::value<std::string>()->required(), "Initial Tree - Newick format (required)")
        ("sequences,f", po::value<std::string>()->required(), "Tip sequences - Fasta format (required)")
        ("kmer-size,k", po::value<std::string>(), "Kmer-size (Valid: 2-15, Default = 15)")
        ("sketch-size,s", po::value<std::string>(), "Sketch-size (Default = 10000)")
        ("threshold,r", po::value<std::string>(), "Erroneous k-mer threshold (Default = 1)")
        ("numBlocks,b", po::value<std::string>(), "Number of cuda blocks (Default = 1)")
        ("blockSize,b", po::value<std::string>(), "Number of cuda threads per block (Default = 1)")
        ("help,h", "Print help messages");

}


void readSequences(po::variables_map& vm, std::vector<std::string>& seqs, std::vector<std::string>& names)
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
        seqs.push_back(std::string(kseq_rd->seq.s, seqLen));
        names.push_back(std::string(kseq_rd->name.s, kseq_rd->name.l));
    }
    
    auto seqReadEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds seqReadTime = seqReadEnd - seqReadStart;
    // std::cout << "Sequences read in: " <<  seqReadTime.count() << " ns\n";
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
    auto inputStart = std::chrono::high_resolution_clock::now();

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
    uint64_t k = 15;
    try {k= (uint64_t)std::stoi(vm["kmer-size"].as<std::string>());}
    catch(std::exception &e){}

    // Sketch Size
    uint64_t sketchSize = 1000;
    try {sketchSize= (uint64_t)std::stoi(vm["sketch-size"].as<std::string>());}
    catch(std::exception &e){}

    // Erroneous k-mer thresold
    uint64_t threshold = 1;
    try {threshold= (uint64_t)std::stoi(vm["threshold"].as<std::string>());}
    catch(std::exception &e){}

    // Number of cuda Blocks
    int numBlocks = 1;
    try {numBlocks= std::stoi(vm["numBlocks"].as<std::string>());}
    catch(std::exception &e){}

    // Number of threads per cuda block
    int blockSize = 1;
    try {blockSize= std::stoi(vm["blockSize"].as<std::string>());}
    catch(std::exception &e){}

    // std::cout << "kmer-size: " << k << 
    // "\nSketch-size: " << sketchSize << 
    // "\nErroneous k-mer threshold: " << threshold << 
    // "\nNo. of cuda blocks: " << numBlocks <<
    // "\nNo. of thread per cuda block: " << blockSize << 
    // "\n";

    MashPlacement::Param params(k, sketchSize, threshold, numBlocks, blockSize);

    std::vector<std::string> seqs,names;

    // Read Input Sequences (Fasta format)
    readSequences(vm, seqs, names);
    size_t numSequences = seqs.size();
    std::vector<int> ids(numSequences);
    std::vector<std::string> temp1(numSequences),temp2(numSequences);
    for(int i=0;i<numSequences;i++) ids[i]=i;
    std::mt19937 rnd(time(NULL));
    std::shuffle(ids.begin(),ids.end(),rnd);
    for(int i=0;i<numSequences;i++){
        temp1[i]=seqs[ids[i]];
        temp2[i]=names[ids[i]];
    }
    seqs=temp1,names=temp2;
    // Compress Sequences (2-bit compressor)
    auto compressStart = std::chrono::high_resolution_clock::now();
    // fprintf(stdout, "Compressing input sequence using two-bit encoding.\n");
    uint64_t ** twoBitCompressedSeqs = new uint64_t*[numSequences];
    uint64_t * seqLengths = new uint64_t[numSequences];
    for (size_t i=0; i<numSequences; i++)
    {   
        uint64_t twoBitCompressedSize = (seqs[i].size()+15)/16;
        uint64_t * twoBitCompressed = new uint64_t[twoBitCompressedSize];
        twoBitCompressor(seqs[i], seqs[i].size(), twoBitCompressed);

        seqLengths[i] = seqs[i].size();
        twoBitCompressedSeqs[i] = twoBitCompressed;

    }
    auto compressEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds compressTime = compressEnd - compressStart;
    // std::cout << "Compressed in: " <<  compressTime.count() << " ns\n";
    auto inputEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds inputTime = inputEnd - inputStart; 
    std::cerr << "Input in: " <<  inputTime.count()/1000000 << " ms\n";

    // Create arrays
    auto createArrayStart = std::chrono::high_resolution_clock::now();
    // fprintf(stdout, "\nAllocating Gpu device arrays.\n");
    MashPlacement::deviceArrays.allocateDeviceArrays(twoBitCompressedSeqs, seqLengths, numSequences, params);
    auto createArrayEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds createArrayTime = createArrayEnd - createArrayStart; 
    std::cerr << "Allocated in: " <<  createArrayTime.count()/1000000 << " ms\n";

    // Build sketch on Gpu
    auto createSketchStart = std::chrono::high_resolution_clock::now();
    MashPlacement::sketchConstructionOnGpu
    (
        MashPlacement::deviceArrays.d_compressedSeqs,
        MashPlacement::deviceArrays.d_prefixCompressed,
        MashPlacement::deviceArrays.d_aggseqLengths,
        MashPlacement::deviceArrays.d_seqLengths,
        MashPlacement::deviceArrays.d_numSequences,
        MashPlacement::deviceArrays.d_hashList,
        seqLengths,
        params
    );
    auto createSketchEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds createSketchTime = createSketchEnd - createSketchStart; 
    std::cerr << "Sketch Created in: " <<  createSketchTime.count()/1000000 << " ms\n";

    // Print first 10 hash values corresponding to each sequence
    // MashPlacement::deviceArrays.printSketchValues(10);

    // Build distance matrix
    auto createTreeStart = std::chrono::high_resolution_clock::now();
    MashPlacement::findPlacementTree(
        MashPlacement::deviceArrays.numSequences,
        MashPlacement::deviceArrays.bd, // id to place
        MashPlacement::deviceArrays.idx, // id of linked-list position
        MashPlacement::deviceArrays.d_dist,
        MashPlacement::deviceArrays.d_head,
        MashPlacement::deviceArrays.d_e,
        MashPlacement::deviceArrays.d_len,
        MashPlacement::deviceArrays.d_nxt,
        MashPlacement::deviceArrays.d_belong,
        MashPlacement::deviceArrays.d_closest_dis,
        MashPlacement::deviceArrays.d_closest_id,
        names,
        MashPlacement::deviceArrays.d_hashList,
        MashPlacement::deviceArrays.d_seqLengths,
        params
    );
    auto createTreeEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds createTreeTime = createTreeEnd - createTreeStart; 
    std::cerr << "Tree Created in: " <<  createTreeTime.count()/1000000 << " ms\n";

    // MashPlacement::deviceArrays.printMashDist(numSequences, names);

    /*Allocate NJ device arrays before deallocating MashPlacement device arrays*/
    
    // neighbourJoining::deviceArrays.allocateDeviceArrays(numSequences, MashPlacement::deviceArrays.d_mashDist);
    
    MashPlacement::deviceArrays.deallocateDeviceArrays();





    return 0;
}
