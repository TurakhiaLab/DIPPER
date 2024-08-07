#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <boost/program_options.hpp> 
#include "../src/kseq.h"
#include "zlib.h"





#ifndef MASH_CUH
#include "../src/mash_sumit.cuh"
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
        ("help,h", "Print help messages");

}


void readSequences(po::variables_map& vm, char ** seqsIn, uint64_t ** seqsLenIn, std::vector<std::string>& seqName)
{
    std::vector<std::string> seqVector;

    auto seqReadStart = std::chrono::high_resolution_clock::now();
    std::string seqFileName = vm["sequences"].as<std::string>();
    // std::cout << "Sequence file name " <<  seqFileName << "\n";

    gzFile f_rd = gzopen(seqFileName.c_str(), "r");
    if (!f_rd) {
        fprintf(stderr, "ERROR: cant open file: %s\n", seqFileName.c_str());
        exit(1);
    }

    kseq_t* kseq_rd = kseq_init(f_rd);

    uint64_t totalLen = 0;
    while (kseq_read(kseq_rd) >= 0) {
        size_t seqLen = kseq_rd->seq.l;
        seqVector.push_back(std::string(kseq_rd->seq.s, seqLen));
        seqName.push_back(kseq_rd->name.s);
        totalLen += seqLen;
    }
    
    char * seqs = new char [totalLen];
    uint64_t * seqsLen = new uint64_t [seqVector.size()];


    uint64_t ptr = 0;
    for (auto&s: seqVector)
    {
        std::strcpy(seqs, s.c_str());
        seqsLen[ptr] = s.size();
        seqs += s.size();
        ptr++;
    }

    seqs -= totalLen;

    auto seqReadEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds seqReadTime = seqReadEnd - seqReadStart;

    *seqsIn = seqs;
    *seqsLenIn = seqsLen;
}

void printSeqs(char ** seqsIn, uint64_t ** seqsLenIn, std::vector<std::string>& seqName)
{
    char * seqs = *seqsIn;
    uint64_t * seqsLen = *seqsLenIn;

    for (int i = 0; i<seqName.size(); i++)
    {
        std::cout << seqName[i]  << ":" <<  seqsLen[i] << std::endl;
        std::string s (seqs, seqsLen[i]);
        std::cout << s << std::endl;
        seqs += seqsLen[i];
    }
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
    uint64_t k = 15;
    try {k= (uint64_t)std::stoi(vm["kmer-size"].as<std::string>());}
    catch(std::exception &e){}

    // Sketch Size
    uint64_t sketchSize = 1280;
    try {sketchSize= (uint64_t)std::stoi(vm["sketch-size"].as<std::string>());}
    catch(std::exception &e){}

    // Seed Size
    uint64_t seed = 50;
    try {seed= (uint64_t)std::stoi(vm["seed"].as<std::string>());}
    catch(std::exception &e){}

    // Erroneous k-mer thresold
    uint64_t threshold = 1;
    try {threshold= (uint64_t)std::stoi(vm["threshold"].as<std::string>());}
    catch(std::exception &e){}

    // Number of cuda Blocks
    int numBlocks = 64;
    try {numBlocks= std::stoi(vm["numBlocks"].as<std::string>());}
    catch(std::exception &e){}

    // Number of threads per cuda block
    int blockSize = 128;
    try {blockSize= std::stoi(vm["blockSize"].as<std::string>());}
    catch(std::exception &e){}

    std::cout << "kmer-size: " << k << 
    "\nSketch-size: " << sketchSize << 
    "\nErroneous k-mer threshold: " << threshold << 
    "\nNo. of cuda blocks: " << numBlocks <<
    "\nNo. of thread per cuda block: " << blockSize << 
    "\n";

    assert((sketchSize%blockSize==0));

    GpuSketch::Param params(k, sketchSize, threshold, numBlocks, blockSize);

    char * seqs;
    uint64_t * seqsLen;
    std::vector<std::string> seqName;

    readSequences(vm, &seqs, &seqsLen, seqName);     // Read Input Sequences (Fasta format)
    uint64_t numSequences = seqName.size();

    // printSeqs(&seqs, &seqsLen, seqName);

    // Create arrays
    auto createArrayStart = std::chrono::high_resolution_clock::now();
    fprintf(stdout, "\nAllocating Gpu device arrays.\n");
    GpuSketch::deviceArrays.allocateDeviceArrays(&seqs, &seqsLen, numSequences, params);
    auto createArrayEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds createArrayTime = createArrayEnd - createArrayStart; 
    std::cout << "Allocated in: " <<  createArrayTime.count() << " ns\n";

    // Build sketch on Gpu
    auto createSketchStart = std::chrono::high_resolution_clock::now();
    GpuSketch::sketchConstructionOnGpu
    (
        GpuSketch::deviceArrays.d_seqs,
        GpuSketch::deviceArrays.d_aggrSeqsLen,
        GpuSketch::deviceArrays.d_seqsLen,
        GpuSketch::deviceArrays.d_numSequences,
        GpuSketch::deviceArrays.d_hashList,
        GpuSketch::deviceArrays.d_hashListPruned,
        params,
        seed
    );
    auto createSketchEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds createSketchTime = createSketchEnd - createSketchStart; 
    std::cout << "Sketch Created in: " <<  createSketchTime.count() << " ns\n";

    return 0;
}

