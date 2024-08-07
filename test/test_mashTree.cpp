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

po::options_description mainDesc("MSA Command Line Arguments");


void parseArguments(int argc, char** argv)
{
    // Setup boost::program_options
    mainDesc.add_options()
            // ("tree,t", po::value<std::string>()->required(), "Initial Tree - Newick format (required)")
            ("sequences,f", po::value<std::string>()->required(), "Tip sequences - Fasta format (required)")
            ("seed,e", po::value<std::string>(), "Seed for the hash function (Default = 50)")
            ("kmer-size,k", po::value<std::string>(), "Kmer-size (Valid: 2-15, Default = 15)")
            ("sketch-size,s", po::value<std::string>(), "Sketch-size (Default = 1000)")
            ("threshold,r", po::value<std::string>(), "Erroneous k-mer threshold (Default = 1)")
            ("numBlocks,b", po::value<std::string>(), "Number of cuda blocks (Default = 128)")
            ("blockSize,bs", po::value<std::string>(), "Number of cuda threads per block (Default = 160)")
            ("verbose,v", po::value<int>(&verbosity)->implicit_value(0)->default_value(1), "Enable verbose mode.")
            ("help,h", "Print help messages");

}


void readSequences(po::variables_map& vm, std::vector<std::string>& seqs, std::vector<std::string>& seqNames)
{
    auto seqReadStart = std::chrono::high_resolution_clock::now();
    std::string seqFileName = vm["sequences"].as<std::string>();
    // std::cout << "Sequence file name " <<  seqFileName << "\n";

    gzFile f_rd = gzopen(seqFileName.c_str(), "r");
    if (!f_rd) {
        fprintf(stderr, "ERROR: cant open file: %s\n", seqFileName.c_str());
        exit(1);
    }

    kseq_t* kseq_rd = kseq_init(f_rd);

    while (kseq_read(kseq_rd) >= 0) {
        size_t seqLen = kseq_rd->seq.l;
        seqs.push_back(std::string(kseq_rd->seq.s, seqLen));
        // std::cout << kseq_rd->name.s << std::endl;
        seqNames.push_back(kseq_rd->name.s);
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

    // Seed for hash function
    uint32_t seed = 50;
    try {seed = (uint32_t)std::stoi(vm["seed"].as<std::string>());}
    catch(std::exception &e){}

    // Kmer Size
    uint32_t k = 15;
    try {k= (uint32_t)std::stoi(vm["kmer-size"].as<std::string>());}
    catch(std::exception &e){}

    // Sketch Size
    uint32_t sketchSize = 1000;
    try {sketchSize= (uint32_t)std::stoi(vm["sketch-size"].as<std::string>());}
    catch(std::exception &e){}

    // Erroneous k-mer thresold
    uint32_t threshold = 1;
    try {threshold= (uint32_t)std::stoi(vm["threshold"].as<std::string>());}
    catch(std::exception &e){}

    // Number of cuda Blocks
    int numBlocks = 128;
    try {numBlocks= std::stoi(vm["numBlocks"].as<std::string>());}
    catch(std::exception &e){}

    // Number of threads per cuda block
    int blockSize = 160;
    try {blockSize= std::stoi(vm["blockSize"].as<std::string>());}
    catch(std::exception &e){}

    if (verbosity)
    {

        std::cout << "kmer-size: " << k <<
                  "\nSketch-size: " << sketchSize <<
                  "\nErroneous k-mer threshold: " << threshold <<
                  "\nNo. of cuda blocks: " << numBlocks <<
                  "\nNo. of thread per cuda block: " << blockSize <<
                  "\n";

    }

    GpuSketch::Param params(k, sketchSize, threshold, numBlocks, blockSize);

    std::vector<std::string> seqs;
    std::vector<std::string> seqNames;

    // Read Input Sequences (Fasta format)
    readSequences(vm, seqs, seqNames);

    size_t numSequences = seqs.size();

    // Compress Sequences (2-bit compressor)
    auto compressStart = std::chrono::high_resolution_clock::now();
    if (verbose) fprintf(stdout, "Compressing input sequence using two-bit encoding.\n");
    uint32_t ** twoBitCompressedSeqs = new uint32_t*[numSequences];
    uint32_t * seqLengths = new uint32_t[numSequences];
    for (size_t i=0; i<numSequences; i++)
    {
        uint32_t twoBitCompressedSize = (seqs[i].size()+15)/16;
        uint32_t * twoBitCompressed = new uint32_t[twoBitCompressedSize];
        twoBitCompressor(seqs[i], seqs[i].size(), twoBitCompressed);

        seqLengths[i] = seqs[i].size();
        twoBitCompressedSeqs[i] = twoBitCompressed;

    }
    auto compressEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds compressTime = compressEnd - compressStart;
    if (verbose) std::cout << "Compressed in: " <<  compressTime.count() << " ns\n";


    // Create arrays
    auto createArrayStart = std::chrono::high_resolution_clock::now();
    if (verbose)  fprintf(stdout, "\nAllocating Gpu device arrays.\n");
    GpuSketch::deviceArrays.allocateDeviceArrays(twoBitCompressedSeqs, seqLengths, numSequences, params);
    auto createArrayEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds createArrayTime = createArrayEnd - createArrayStart;
    if (verbose)  std::cout << "Allocated in: " <<  createArrayTime.count() << " ns\n";

    // Build sketch on Gpu
    auto createSketchStart = std::chrono::high_resolution_clock::now();
    GpuSketch::sketchConstructionOnGpu
            (
                    GpuSketch::deviceArrays.d_compressedSeqs,
                    GpuSketch::deviceArrays.d_aggseqLengths,
                    GpuSketch::deviceArrays.d_seqLengths,
                    GpuSketch::deviceArrays.d_numSequences,
                    GpuSketch::deviceArrays.d_hashList,
                    GpuSketch::deviceArrays.d_hashListPruned,
                    seqLengths,
                    params,
                    seed
            );
    auto createSketchEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds createSketchTime = createSketchEnd - createSketchStart;
    //std::cout << "Sketch Created in: " <<  createSketchTime.count() << " ns\n";

    // Print first 10 hash values corresponding to each sequence
    //GpuSketch::deviceArrays.printSketchValues(10, seqLengths, params);

    // Build distance matrix
    auto createDistMatStart = std::chrono::high_resolution_clock::now();
    GpuSketch::mashDistConstructionOnGpu
            (
                    GpuSketch::deviceArrays.d_hashListPruned,
                    GpuSketch::deviceArrays.d_seqLengths,
                    GpuSketch::deviceArrays.d_numSequences,
                    GpuSketch::deviceArrays.d_mashDist,
                    seqLengths,
                    params
            );
    auto createDistMatEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds createDistMatTime = createDistMatEnd - createDistMatStart;
    if (verbose)  std::cout << "Distance Matrix Created in: " <<  createDistMatTime.count() << " ns\n";
    GpuSketch::deviceArrays.printMashDist(numSequences, seqNames);




    GpuSketch::deviceArrays.deallocateDeviceArrays();





    return 0;
}
