#include <cstdint>
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <bits/stdc++.h>
#include <boost/program_options.hpp> 
#include "../src/kseq.h"
#include "zlib.h"


#ifndef TWOBITCOMPRESSOR_HPP
#include "../src/twoBitCompressor.hpp"
#endif

#ifndef FOURBITCOMPRESSOR_HPP
#include "../src/fourBitCompressor.hpp"
#endif

#ifndef MASHPL_CUH
#include "../src/mash_placement.cuh"
#endif

namespace po = boost::program_options;

KSEQ_INIT2(, gzFile, gzread)

po::options_description mainDesc("Placement Command Line Arguments");

#define CHECK_CUDA_ERROR(error) checkCudaError(error, __FILE__, __LINE__)

#define THREADS_PER_BLOCK 256

void parseArguments(int argc, char** argv)
{
    // Setup boost::program_options
    mainDesc.add_options()
        // ("tree,t", po::value<std::string>()->required(), "Initial Tree - Newick format (required)")
        ("input-file,f", po::value<std::string>()->required(), "Input format, phylip format for distance matrix, fasta format for MSA or raw sequences, required")
        ("kmer-size,k", po::value<std::string>(), "Kmer-size (Valid: 2-15, Default = 15)")
        ("sketch-size,s", po::value<std::string>(), "Sketch-size (Default = 10000)")
        ("threshold,r", po::value<std::string>(), "Erroneous k-mer threshold (Default = 1)")
        ("distance-type,t", po::value<std::string>(), "Distance type to be calculated, range from 1 to 5, please try 1 and 2 for now")
        ("input-format,i", po::value<std::string>()->required(), "Input format (d - distance matrix, r - raw sequences, a - MSA), required")
        ("output-format,o", po::value<std::string>()->required(), "Output format (d - distance matrix, t - phylogenetic tree), required")
        ("help,h", "Print help messages");

}


void readSequences(po::variables_map& vm, std::vector<std::string>& seqs, std::vector<std::string>& names)
{
    auto seqReadStart = std::chrono::high_resolution_clock::now();
    std::string seqFileName = vm["input-file"].as<std::string>();

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

// struct treeNode {
//     int nodeNum;
//     int nodechild1;
//     int nodechild2;
// };



void getTwoRandomIndices(int *clusterMap, int numSequences, int searchIndex, MashPlacement::treeNode *node)
{
    std::unordered_set<int> uniqueIndices;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, numSequences - 1); // Range from 0 to clusterSize - 1

    int randomNumber = dis(gen);
    for (int i = 0; i <numSequences; i++)
    {
        int indx = ( randomNumber+ i) % numSequences;
        if (clusterMap[indx] == searchIndex && uniqueIndices.find(indx) == uniqueIndices.end())
        {
            uniqueIndices.insert(indx);
            if (uniqueIndices.size() == 1)
                randomNumber = dis(gen);
            else if (uniqueIndices.size() == 2)
                break;
        }
    }

    if (uniqueIndices.size() >= 2)
    {

        auto it = uniqueIndices.begin();
        node->nodechild1 = *it++;
        node->nodechild2 = *it;
        node->nodeNum = searchIndex;
    }
}


 std::unordered_map<int, std::vector<int>> findUniqueNumbersAndIndices(const std::vector<int>& arr) {
            std::unordered_map<int, std::vector<int>> numberIndices;

            for (int i = 0; i < arr.size(); ++i) {
                numberIndices[arr[i]].push_back(i);
            }

            return numberIndices;
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

    uint64_t distanceType = 1;
    try {threshold= (uint64_t)std::stoi(vm["distance-type"].as<std::string>());}
    catch(std::exception &e){}

    std::string in = "r";
    try {in = vm["input-format"].as<std::string>();}
    catch(std::exception &e){}

    std::string out = "t";
    try {out = vm["output-format"].as<std::string>();}
    catch(std::exception &e){}


    MashPlacement::Param params(k, sketchSize, threshold, distanceType, in, out);
    if (in == "m" && out == "t"){
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
        uint64_t ** fourBitCompressedSeqs = new uint64_t*[numSequences];
        uint64_t * seqLengths = new uint64_t[numSequences];
        for (size_t i=0; i<numSequences; i++)
        {   
            uint64_t fourBitCompressedSize = (seqs[i].size()+15)/16;
            uint64_t * fourBitCompressed = new uint64_t[fourBitCompressedSize];
            fourBitCompressor(seqs[i], seqs[i].size(), fourBitCompressed);

            seqLengths[i] = seqs[i].size();
            fourBitCompressedSeqs[i] = fourBitCompressed;

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
        // std::cerr<<"########\n";
        MashPlacement::msaDeviceArrays.allocateDeviceArrays(fourBitCompressedSeqs, seqLengths, numSequences, params);
        // std::cerr<<"########\n";
        MashPlacement::placementDeviceArrays.allocateDeviceArrays(numSequences);
        auto createArrayEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds createArrayTime = createArrayEnd - createArrayStart; 
        std::cerr << "Allocated in: " <<  createArrayTime.count()/1000000 << " ms\n";


        //Build Tree on Gpu
        auto createTreeStart = std::chrono::high_resolution_clock::now();
        MashPlacement::placementDeviceArrays.findPlacementTree(params, MashPlacement::mashDeviceArrays, MashPlacement::matrixReader, MashPlacement::msaDeviceArrays);
        auto createTreeEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds createTreeTime = createTreeEnd - createTreeStart; 
        MashPlacement::placementDeviceArrays.printTree(names);
        std::cerr << "Tree Created in: " <<  createTreeTime.count()/1000000 << " ms\n";

        // Print first 10 hash values corresponding to each sequence
        // MashPlacement::mashDeviceArrays.printSketchValues(10);


        MashPlacement::msaDeviceArrays.deallocateDeviceArrays();
        MashPlacement::placementDeviceArrays.deallocateDeviceArrays();
    }
    else if (in == "r" && out == "t"){
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
            uint64_t twoBitCompressedSize = (seqs[i].size()+31)/32;
            uint64_t * twoBitCompressed = new uint64_t[twoBitCompressedSize];
            twoBitCompressor(seqs[i], seqs[i].size(), twoBitCompressed);

            seqLengths[i] = seqs[i].size();
            twoBitCompressedSeqs[i] = twoBitCompressed;
        }
 
        std::cout<<std::endl;
        auto compressEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds compressTime = compressEnd - compressStart;
        // std::cout << "Compressed in: " <<  compressTime.count() << " ns\n";
        auto inputEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds inputTime = inputEnd - inputStart; 
        std::cerr << "Input in: " <<  inputTime.count()/1000000 << " ms\n";

        // Create arrays
        auto createArrayStart = std::chrono::high_resolution_clock::now();
        // fprintf(stdout, "\nAllocating Gpu device arrays.\n");
        MashPlacement::mashDeviceArrays.allocateDeviceArrays(twoBitCompressedSeqs, seqLengths, numSequences, params);
       
        auto createArrayEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds createArrayTime = createArrayEnd - createArrayStart; 
        std::cerr << "Allocated in: " <<  createArrayTime.count()/1000000 << " ms\n";

        // Build sketch on Gpu
        auto createSketchStart = std::chrono::high_resolution_clock::now();
        MashPlacement::mashDeviceArrays.sketchConstructionOnGpu(params);
        auto createSketchEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds createSketchTime = createSketchEnd - createSketchStart; 
        std::cerr << "Sketch Created in: " <<  createSketchTime.count()/1000000 << " ms\n";

        // Clustering operation
        int *clusterMap = new int[numSequences](); // Allocate memory for cluster mapping
        int MAX_LEVELS = 0;
        int tempNumSequences = numSequences;
        // Calculate the maximum number of levels in the tree
        while (tempNumSequences >>= 1) ++MAX_LEVELS; // Determine tree depth based on sequence count
        // Create an array of pointers to treeNode objects
        MashPlacement::treeNode **nodes = new MashPlacement::treeNode*[1 << MAX_LEVELS]; // Allocate memory for tree nodes

        // Initialize each node in the array
        for (int i = 0; i < 1 << MAX_LEVELS; i++)  nodes[i] = new MashPlacement::treeNode(); // Create treeNode objects

        // Process cluster levels using the MashPlacement::mashDeviceArrays object
        MashPlacement::mashDeviceArrays.processClusterLevels(clusterMap, numSequences, nodes, MAX_LEVELS, params); // Perform clustering

        // Find unique numbers and their indices in the clusterMap
        auto uniqueNumberIndices = findUniqueNumbersAndIndices(std::vector<int>(clusterMap, clusterMap + numSequences)); // Identify unique clusters

        // Clean up dynamically allocated memory
        delete[] nodes; // Free memory for tree nodes
        delete[] clusterMap; // Free memory for cluster mapping

        // Print Mash distances and deallocate device arrays
        // MashPlacement::mashDeviceArrays.printMashDist(params, numSequences, names); // Output Mash distances
        MashPlacement::mashDeviceArrays.deallocateDeviceArrays(); // Free GPU memory

        // Print unique numbers and their indices
        std::cout << "Unique numbers and their indices:" << std::endl;
        for (const auto& pair : uniqueNumberIndices) {
            std::cout << "Number " << pair.first << " occurs at indices: ";
            for (int index : pair.second) {
                std::cout << index << " ";
            }
            std::cout << std::endl;
        }

        // Process each unique cluster
        for (const auto& pair : uniqueNumberIndices) {
            std::cout << std::endl << std::endl << std::endl << std::endl;
            std::vector<std::string> namesCluster;

            // Prepare data for the current cluster
            uint64_t ** twoBitCompressedSeqsCluster;
            int numSequencesCluster = pair.second.size();
            uint64_t * seqLengthsCluster;

            twoBitCompressedSeqsCluster = new uint64_t*[numSequencesCluster]; // Allocate memory for compressed sequences
            seqLengthsCluster = new uint64_t[numSequencesCluster]; // Allocate memory for sequence lengths

            // Populate cluster data
            for (int k = 0; k < numSequencesCluster; k++) {
                int index = pair.second[k]; 
                std::cout << names[index] << " ";

                seqLengthsCluster[k] = seqLengths[index]; // Store sequence length
                twoBitCompressedSeqsCluster[k] = twoBitCompressedSeqs[index]; // Store compressed sequence
                namesCluster.push_back(names[index]); // Store sequence name
            }

            // Allocate device arrays for placement
            MashPlacement::placementDeviceArrays.allocateDeviceArrays(numSequencesCluster); // Prepare GPU memory

            // Build Tree on GPU for the current cluster
            auto createTreeStart = std::chrono::high_resolution_clock::now(); // Start timing
            MashPlacement::placementDeviceArrays.findPlacementTree(params, MashPlacement::mashDeviceArrays, MashPlacement::matrixReader, MashPlacement::msaDeviceArrays); // Construct tree
            auto createTreeEnd = std::chrono::high_resolution_clock::now(); // End timing
            std::chrono::nanoseconds createTreeTime = createTreeEnd - createTreeStart; 

            // Print the tree and creation time
            MashPlacement::placementDeviceArrays.printTree(namesCluster); // Output tree structure
            std::cerr << "Tree Created in: " <<  createTreeTime.count()/1000000 << " ms\n"; // Report tree creation time

            // Deallocate device arrays
            MashPlacement::placementDeviceArrays.deallocateDeviceArrays(); // Free GPU memory
        }       
    }
    else if(in == "d" && out == "t") {
        std::string fileName = vm["input-file"].as<std::string>();
        FILE* filePtr = fopen(fileName.c_str(), "r");
        if (filePtr == nullptr){
            std::cerr << "Cannot open file: " << fileName << std::endl;
            return 1;
        }
        const size_t bufferSize = 64 * 1024 * 1024; 
        char* buffer = new char[bufferSize];
        if (setvbuf(filePtr, buffer, _IOFBF, bufferSize) != 0) {
            std::cerr << "Failed in setting buffer" << std::endl;
            delete[] buffer;
            fclose(filePtr);
            return 1;
        }
        char *temp = new char[20];
        int numSequences;
        fscanf(filePtr, "%d", &numSequences);
        fgets(temp, 20, filePtr);
        MashPlacement::matrixReader.allocateDeviceArrays(numSequences, filePtr);
        MashPlacement::placementDeviceArrays.allocateDeviceArrays(numSequences);
        MashPlacement::placementDeviceArrays.findPlacementTree(params, MashPlacement::mashDeviceArrays, MashPlacement::matrixReader, MashPlacement::msaDeviceArrays);
        MashPlacement::placementDeviceArrays.printTree(MashPlacement::matrixReader.name);
        fclose(filePtr);
    }
    else{
        printf("Invalid input-output combinations!!!!!\n");
        exit(1);
    }
    return 0;
}


