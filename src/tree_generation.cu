#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <chrono>
#include <cassert>
#include <bits/stdc++.h>
#include <boost/program_options.hpp> 
#include "../src/kseq.h"
#include "zlib.h"
#include <cuda_runtime.h>
#include <tbb/tbb.h>
#include <tbb/parallel_for.h>
#include "version.hpp"

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

po::options_description mainDesc("DIPPER Command Line Arguments");


void parseArguments(int argc, char** argv)
{
    // Setup boost::program_options
    po::options_description requiredDesc("Required Options");
    requiredDesc.add_options()
        ("input-format,i",     po::value<std::string>()->required(),
        "Input format:\n"
        "  d - distance matrix in PHYLIP format\n"
        "  r - unaligned sequences in FASTA format\n"
        "  m - aligned sequences in FASTA format")

        ("input-file,I",       po::value<std::string>()->required(),
        "Input file path:\n"
        "  PHYLIP format for distance matrix\n"
        "  FASTA format for aligned or unaligned sequences")

        ("output-file,O",      po::value<std::string>()->required(),
        "Output file path");


    po::options_description optionalDesc("Optional Options");
    optionalDesc.add_options()
        ("output-format,o",    po::value<std::string>(),
        "Output format:\n"
        "  t - phylogenetic tree in Newick format (default)\n"
        "  d - distance matrix in PHYLIP format (coming soon)")

        ("algorithm,m",        po::value<std::string>(),
        "Algorithm selection:\n"
        "  0 - default mode\n"
        "  1 - force placement\n"
        "  2 - force conventional NJ\n"
        "  3 - force divide-and-conquer")

        ("K-closest,K",   po::value<std::string>(),
        "Placement mode:\n"
        "  -1 - exact mode\n"
        "  10 - default")

        ("kmer-size,k",        po::value<std::string>(),
        "K-mer size:\n"
        "  Valid range: 2-15 (default: 15)")

        ("sketch-size,s",      po::value<std::string>(),
        "Sketch size (default: 1000)")

        ("distance-type,d",    po::value<std::string>(),
        "Distance type to calculate:\n"
        "  1 - uncorrected\n"
        "  2 - JC (default)\n"
        "  3 - Tajima-Nei\n"
        "  4 - K2P\n"
        "  5 - Tamura\n"
        "  6 - Jinnei")

        ("add,a",
        "Add query to backbone using k-closest placement")

        ("threads,T", po::value<int>(),
        "Number of CPU threads (default: all available threads)")

        ("gap-col,g", po::value<double>()->default_value(1.0),
        "Threshold for gappy columns (default: 1)\n Columns with a gap ratio exceeding this value will be removed to speed up computation.")

        ("nongap-row,n", po::value<double>()->default_value(0),
        "Threshold for nongap characters (default: 0)\nSequences with fewer nongap characters than the median count multiplied by this value will be removed to reduce noise.")

        ("gpu-index", po::value<int>()->default_value(0),
        "GPU device index to use (default: 0)")

        ("input-tree,t",       po::value<std::string>(),
        "Input backbone tree (Newick format), required with --add option")

        ("help,h",
        "Print this help message")
        
        ("version,v", "Print DIPPER version");

    mainDesc.add(requiredDesc).add(optionalDesc);

}

std::vector<std::string> filter_msa(std::vector<std::string>& sequences, std::vector<std::string> names, po::variables_map vm)
{
    double colTH = vm["gap-col"].as<double>(), rowTH = vm["nongap-row"].as<double>();

    if (sequences.empty()) {
        return sequences;
    }

    const size_t num_sequences = sequences.size();
    const size_t num_columns = sequences[0].length();
    std::vector<std::string> filtered_sequences (num_sequences);
    
    if (colTH < 1) {
        auto rmColStart = std::chrono::high_resolution_clock::now();
        // Calculate the maximum number of allowed gaps for a column to be kept.
        const size_t max_allowed_gaps = static_cast<size_t>(
            std::floor(colTH * num_sequences)
        );
        std::cerr << "Maximum allowed gaps: " << max_allowed_gaps << '\n';
        // This boolean vector will store the result of our parallel analysis:
        // true if the column should be KEPT, false if it should be removed.
        std::vector<bool> keep_column(num_columns);

        // --- Step 1: Parallel Gap Counting using TBB ---
        // The parallel_for iterates over all columns (0 to num_columns-1).
        tbb::parallel_for( tbb::blocked_range<size_t>(0, num_columns), [&](const tbb::blocked_range<size_t>& r) {
        for (size_t col_idx = r.begin(); col_idx != r.end(); ++col_idx) {
                size_t gap_count = 0;
                // Iterate over all sequences for the current column
                for (size_t seq_idx = 0; seq_idx < num_sequences; ++seq_idx) {
                    if (sequences[seq_idx][col_idx] == '-') {
                        gap_count++;
                    }
                }
                // Store the decision
                keep_column[col_idx] = (gap_count <= max_allowed_gaps);
            }
        }
        );

        // --- Step 2: Sequential Filtering ---
        // Calculate the length of the new sequences and pre-size the result strings
        size_t new_length = 0;
        for (bool keep : keep_column) {
            if (keep) {
                new_length++;
            }
        }
        
        for (size_t i = 0; i < num_sequences; ++i) {
            filtered_sequences.reserve(new_length); // Initialize string with new length
        }
        // Iterate through the original sequences, copying only the KEPT columns
        tbb::parallel_for( tbb::blocked_range<size_t>(0, num_sequences), [&](const tbb::blocked_range<size_t>& r) {
        for (size_t seq_idx = r.begin(); seq_idx != r.end(); ++seq_idx) {
            size_t new_col_idx = 0;
            for (size_t old_col_idx = 0; old_col_idx < num_columns; ++old_col_idx) {
                if (keep_column[old_col_idx]) {
                    filtered_sequences[seq_idx] += sequences[seq_idx][old_col_idx];
                    new_col_idx++;
                }
            }
        }
        }
        );
        std::cerr << "Original length: " << num_columns 
                  << ", length after removing gappy columns: " << new_length << '\n';
        auto rmColEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds rmColTime = rmColEnd - rmColStart;
        std::cerr << "Remove gappy columns in " << rmColTime.count() / 1000000 << " ms\n";
    }
    
    if (colTH == 1) filtered_sequences = std::move(sequences);
    if (rowTH > 0) {
        auto rmRowStart = std::chrono::high_resolution_clock::now();
        int numCols = filtered_sequences[0].size();
        std::vector<int> nongapCount (num_sequences, 0);
        // Count nongap characters of all sequences
        tbb::parallel_for( tbb::blocked_range<size_t>(0, num_sequences), [&](const tbb::blocked_range<size_t>& r) {
        for (size_t row_idx = r.begin(); row_idx != r.end(); ++row_idx) {
                int nongap_count = 0;
                for (size_t col_idx = 0; col_idx < numCols; ++col_idx) {
                    if (filtered_sequences[row_idx][col_idx] != '-') nongap_count++;
                }
                nongapCount[row_idx] = nongap_count;
            }
        }
        );
        std::vector<int> nongapCount_copy = nongapCount;
        std::sort(nongapCount_copy.begin(), nongapCount_copy.end());
        int min_nongap = static_cast<int>( std::floor(nongapCount_copy[num_sequences/2] * rowTH));
        nongapCount_copy.clear();

        std::cerr << "Minimum required non-gaps: " << min_nongap << '\n';

        std::vector<std::string> preserve_seqs;
        std::vector<std::string> preserve_names;
        for (int sIdx = 0; sIdx < num_sequences; ++sIdx) {
            if (nongapCount[sIdx] >= min_nongap) {
                preserve_seqs.push_back(filtered_sequences[sIdx]);
                preserve_names.push_back(names[sIdx]);
            }
        }

        std::cerr << "Original: " << filtered_sequences.size()
                  << ", sequence number after removing gappy sequences: " << preserve_seqs.size() << '\n';
        auto rmRowEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds rmRowTime = rmRowEnd - rmRowStart;
        std::cerr << "Remove gappy sequences in " << rmRowTime.count() / 1000000 << " ms\n";
        filtered_sequences = std::move(preserve_seqs);
        names = std::move(preserve_names);
        assert(filtered_sequences.size() == names.size());
    }
    return filtered_sequences;
}

void readAllSequences(po::variables_map& vm, std::vector<std::string>& seqs, std::vector<std::string>& names, std::unordered_map<std::string, int>& nameToIdx)
{
    auto seqReadStart = std::chrono::high_resolution_clock::now();
    std::string seqFileName = vm["input-file"].as<std::string>();

    gzFile f_rd = gzopen(seqFileName.c_str(), "r");
    if (!f_rd) {
        fprintf(stderr, "ERROR: cant open file: %s\n", seqFileName.c_str());
        exit(1);
    }

    kseq_t* kseq_rd = kseq_init(f_rd);

    seqs.resize(names.size());

    while (kseq_read(kseq_rd) >= 0) {
        size_t seqLen = kseq_rd->seq.l;
        if (nameToIdx.find(std::string(kseq_rd->name.s, kseq_rd->name.l)) == nameToIdx.end()) {
            seqs.push_back(std::string(kseq_rd->seq.s, seqLen));
            names.push_back(std::string(kseq_rd->name.s, kseq_rd->name.l));
        } else {
            int id = nameToIdx[std::string(kseq_rd->name.s, kseq_rd->name.l)];
            seqs[id] = std::string(kseq_rd->seq.s, seqLen);
        }
    }

    auto seqReadEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds seqReadTime = seqReadEnd - seqReadStart;
    // std::cout << "Sequences read in: " <<  seqReadTime.count() << " ns\n";
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

void writeAlignment(std::string fileName, std::vector<std::string> seqs, std::vector<std::string> names) {
    std::ofstream outFile;
    std::vector<std::pair<std::string,std::string>> seqPairs;
    for (int i = 0; i < seqs.size(); ++i) {
        seqPairs.push_back({names[i], seqs[i]});
    } 
    
    std::sort(seqPairs.begin(), seqPairs.end(),[](const auto& a, const auto& b) { return a.second < b.second; });
    
    outFile.open(fileName);
    if (!outFile) {
        fprintf(stderr, "ERROR: cant open file: %s\n", fileName.c_str());
        exit(1);
    }
    for (int i = 0; i < seqPairs.size(); ++i) {
        outFile << ('>' + seqPairs[i].first + '\n');
        outFile << (seqPairs[i].second + '\n');
    }
    outFile.close();

    auto current_second = seqPairs[0].second;
    bool first_in_line = true;

    std::vector<std::string> ids;
    for (const auto& p : seqPairs) {
        if (p.second != current_second) {
            // 換新的一行（新的 second）
            if (ids.size() > 30) {
                std::cout << ids.size() << ": \n";
                for (auto id: ids) std::cout << id << ',';
                std::cout << "\n";
            }
            current_second = p.second;
            ids.clear();
            ids.push_back(p.first);
        }
        else {
            ids.push_back(p.first);
        }
    }
    return;
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
        if(vm.count("help")) {
            std::cerr << mainDesc << std::endl;
            return 0;
        }
        else if (vm.count("version")) {
            std::cout << "DIPPER Version " << PROJECT_VERSION << std::endl;
            return;
        } 
        std::cerr << "\033[31m" << e.what() << "\033[0m"  << std::endl;
        std::cerr << mainDesc << std::endl;
        return 1;
    }

    

    
    if (vm.count("add")) {
        if (!vm.count("input-tree")) {
            std::cerr << "\033[31m" << "Backbone tree (--input-tree/-t) is required with --add option" << "\033[0m" << std::endl;
            std::cerr << mainDesc << std::endl;
            return 1;
        }
    }
    

    // Kmer Size
    uint64_t k = 15;
    try {k= (uint64_t)std::stoi(vm["kmer-size"].as<std::string>());}
    catch(std::exception &e){}

    // Sketch Size
    uint64_t sketchSize = 1000;
    // try {sketchSize= (uint64_t)std::stoi(vm["sketch-size"].as<std::string>());}
    // catch(std::exception &e){}

    // Erroneous k-mer thresold
    uint64_t threshold = 1;
    // try {threshold= (uint64_t)std::stoi(vm["threshold"].as<std::string>());}
    // catch(std::exception &e){}

    uint64_t distanceType = 1;
    try {distanceType= (uint64_t)std::stoi(vm["distance-type"].as<std::string>());}
    catch(std::exception &e){}

    std::string in = "r";
    try {in = vm["input-format"].as<std::string>();}
    catch(std::exception &e){}

    std::string out = "t";
    try {out = vm["output-format"].as<std::string>();}
    catch(std::exception &e){}

    std::string algo = "0";
    try {algo = vm["algorithm"].as<std::string>();}
    catch(std::exception &e){}

    std::string placemode = "10";
    try {placemode = vm["K-closest"].as<std::string>();}
    catch(std::exception &e){}

    bool add = false;
    if (vm.count("add")) add = true;

    std::string treeFile = "";
    try {treeFile = vm["input-tree"].as<std::string>();}
    catch(std::exception &e){}
    if (add && treeFile == "") {
        std::cerr << "ERROR: Input tree file is required for adding query to a backbone tree.\n";
        return 1;
    }

    std::string outputFile = vm["output-file"].as<std::string>();
    std::ofstream output_(outputFile.c_str());

    // int device_id = 1;  // Use GPU 1 (second GPU)
    // cudaError_t err = cudaSetDevice(device_id);
    // if (err != cudaSuccess) {
    //     std::cerr << "Failed to set CUDA device: " << cudaGetErrorString(err) << std::endl;
    //     return -1;
    // }

    int placement_thr = 30000; 
    int dc_thr = 1000000; 

    cudaSetDevice(vm["gpu-index"].as<int>());
    int cpuThreads = (vm.count("threads")) ? vm["threads"].as<int>() : tbb::this_task_arena::max_concurrency();
    tbb::global_control init(tbb::global_control::max_allowed_parallelism, cpuThreads);
    fprintf(stderr, "Using GPU #%d.\n", vm["gpu-index"].as<int>());
    fprintf(stderr, "Maximum available CPU cores: %d. Using %d CPU cores.\n", tbb::this_task_arena::max_concurrency(), cpuThreads);


    MashPlacement::Param params(k, sketchSize, threshold, distanceType, in, out);

    if (add) {
        // Load the tree from the file
        std::ifstream treeFileStream(treeFile);
        if (!treeFileStream) {
            std::cerr << "ERROR: Unable to open input tree file: " << treeFile << "\n";
            return 1;
        }
        std::vector<std::string> seqs, names, namesDump;
        readSequences(vm, seqs, namesDump);
        std::cerr << "Read " << seqs.size() << " sequences from input file.\n";
        assert(seqs.size() > 0 && "No sequences found in the input file.");

        std::string newickTree;
        std::getline(treeFileStream, newickTree);
        Tree *t = new Tree(newickTree, namesDump.size());
        std::cerr << "Tree loaded successfully with "<< t->allNodes.size()<<" nodes and root " << t->root->name << ".\n";
        size_t backboneSize = t->m_numLeaves;
        size_t numSequences = seqs.size();

        std::unordered_map<int, int> idMap;

        names.resize(backboneSize);
        for (int i=0; i<numSequences;i++){
            if (t->allNodes.find(namesDump[i]) == t->allNodes.end()) {
                names.push_back(namesDump[i]);
                idMap[i] = names.size()-1;
            } else {
                names[t->allNodes[namesDump[i]]->idx] = namesDump[i];
                idMap[i]=t->allNodes[namesDump[i]]->idx;
            }
        }

        if (in == "r" && out == "t") {
            uint64_t ** twoBitCompressedSeqs = new uint64_t*[numSequences];
            uint64_t * seqLengths = new uint64_t[numSequences];
            tbb::parallel_for(tbb::blocked_range<int>(0, numSequences), [&](tbb::blocked_range<int> range){
            for (int idx_= range.begin(); idx_ < range.end(); ++idx_) {
                uint64_t i = static_cast<uint64_t>(idx_);
                uint64_t twoBitCompressedSize = (seqs[i].size()+31)/32;
                uint64_t * twoBitCompressed = new uint64_t[twoBitCompressedSize];
                twoBitCompressor(seqs[i], seqs[i].size(), twoBitCompressed);

                int newId = idMap[i];
                seqLengths[newId] = seqs[i].size();
                twoBitCompressedSeqs[newId] = twoBitCompressed;
            }});
            std::cerr << "Allocating Mash Device Arrays" << std::endl;
            MashPlacement::mashDeviceArrays.allocateDeviceArrays(twoBitCompressedSeqs, seqLengths, numSequences, params);
            
            std::cerr << "Sketch Construction in Progress" << std::endl;
            MashPlacement::mashDeviceArrays.sketchConstructionOnGpu(params);

            MashPlacement::kplacementDeviceArrays.allocateDeviceArrays(numSequences, backboneSize);
            MashPlacement::kplacementDeviceArrays.initializeDeviceArrays(t);
            MashPlacement::kplacementDeviceArrays.addQuery(params, MashPlacement::mashDeviceArrays, MashPlacement::matrixReader, MashPlacement::msaDeviceArrays);
            MashPlacement::kplacementDeviceArrays.printTree(names, output_);
        } else if (in == "m" && out == "t") {
            uint64_t ** fourBitCompressedSeqs = new uint64_t*[numSequences];
            uint64_t * seqLengths = new uint64_t[numSequences];
            tbb::parallel_for(tbb::blocked_range<int>(0, numSequences), [&](tbb::blocked_range<int> range){
            for (int idx_= range.begin(); idx_ < range.end(); ++idx_) {
                uint64_t i = static_cast<uint64_t>(idx_);
                uint64_t fourBitCompressedSize = (seqs[i].size()+15)/16;
                uint64_t * fourBitCompressed = new uint64_t[fourBitCompressedSize];
                fourBitCompressor(seqs[i], seqs[i].size(), fourBitCompressed);

                int newId = idMap[i];
                seqLengths[newId] = seqs[i].size();
                fourBitCompressedSeqs[newId] = fourBitCompressed;
            }});
            MashPlacement::msaDeviceArrays.allocateDeviceArrays(fourBitCompressedSeqs, seqLengths, numSequences, params);
            MashPlacement::kplacementDeviceArrays.allocateDeviceArrays(numSequences, backboneSize);
            MashPlacement::kplacementDeviceArrays.initializeDeviceArrays(t);
            MashPlacement::kplacementDeviceArrays.addQuery(params, MashPlacement::mashDeviceArrays, MashPlacement::matrixReader, MashPlacement::msaDeviceArrays);
            MashPlacement::kplacementDeviceArrays.printTree(names, output_);
        } else {
            std::cerr << "Adding new sequnces only supported with input aligned and unaligned sequences\n";
            exit(1);
        }
        return;
    }

    if (in == "m" && out == "t"){
        std::vector<std::string> seqs,names_, names;

        // Read Input Sequences (Fasta format)
        readSequences(vm, seqs, names_);
        seqs = filter_msa(seqs, names_, vm);
        std::string fileName = vm["input-file"].as<std::string>()+".masked.aln";
        // writeAlignment(fileName, seqs, names_);
        size_t numSequences = seqs.size();
        names.resize(numSequences);
        std::vector<int> ids(numSequences);
        for(int i=0;i<numSequences;i++) ids[i]=i;
        std::mt19937 rnd(time(NULL));
        std::shuffle(ids.begin(),ids.end(),rnd);

    
        // Compress Sequences (2-bit compressor)
        auto compressStart = std::chrono::high_resolution_clock::now();
        // fprintf(stdout, "Compressing input sequence using two-bit encoding.\n");
        uint64_t ** fourBitCompressedSeqs = new uint64_t*[numSequences];
        uint64_t * seqLengths = new uint64_t[numSequences];
        tbb::parallel_for(tbb::blocked_range<int>(0, numSequences), [&](tbb::blocked_range<int> range){
        for (int idx_= range.begin(); idx_ < range.end(); ++idx_) {
            uint64_t i = static_cast<uint64_t>(idx_);
            uint64_t fourBitCompressedSize = (seqs[i].size()+15)/16;
            uint64_t * fourBitCompressed = new uint64_t[fourBitCompressedSize];
            fourBitCompressor(seqs[i], seqs[i].size(), fourBitCompressed);

            seqLengths[ids[i]] = seqs[i].size();
            fourBitCompressedSeqs[ids[i]] = fourBitCompressed;
            names[ids[i]] = names_[i];
        }});
        
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
        // std::cerr<<"########\n";
        MashPlacement::msaDeviceArrays.allocateDeviceArrays(fourBitCompressedSeqs, seqLengths, numSequences, params);
        if(algo=="1"||algo=="0"&&numSequences>=placement_thr&&numSequences<dc_thr){
            std::cerr<<"Using ";
            if(placemode=="-1"){
                std::cerr<<" exact placement mode\n";
                MashPlacement::placementDeviceArrays.allocateDeviceArrays(numSequences);
                auto createArrayEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds createArrayTime = createArrayEnd - createArrayStart; 
                std::cerr << "Allocated in: " <<  createArrayTime.count()/1000000 << " ms\n";


                //Build Tree on Gpu
                auto createTreeStart = std::chrono::high_resolution_clock::now();
                MashPlacement::placementDeviceArrays.findPlacementTree(params, MashPlacement::mashDeviceArrays, MashPlacement::matrixReader, MashPlacement::msaDeviceArrays);
                auto createTreeEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds createTreeTime = createTreeEnd - createTreeStart; 
                MashPlacement::placementDeviceArrays.printTree(names, output_);
                std::cerr << "Tree Created in: " <<  createTreeTime.count()/1000000 << " ms\n";

                // Print first 10 hash values corresponding to each sequence
                // MashPlacement::mashDeviceArrays.printSketchValues(10);
                MashPlacement::msaDeviceArrays.deallocateDeviceArrays();
                MashPlacement::placementDeviceArrays.deallocateDeviceArrays();
            }
            else{
                std::cerr<<"k-closest placement mode\n";
                MashPlacement::kplacementDeviceArrays.allocateDeviceArrays(numSequences);
                auto createArrayEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds createArrayTime = createArrayEnd - createArrayStart; 
                std::cerr << "Allocated in: " <<  createArrayTime.count()/1000000 << " ms\n";


                //Build Tree on Gpu
                auto createTreeStart = std::chrono::high_resolution_clock::now();
                MashPlacement::kplacementDeviceArrays.findPlacementTree(params, MashPlacement::mashDeviceArrays, MashPlacement::matrixReader, MashPlacement::msaDeviceArrays);
                auto createTreeEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds createTreeTime = createTreeEnd - createTreeStart; 
                MashPlacement::kplacementDeviceArrays.printTree(names, output_);
                std::cerr << "Tree Created in: " <<  createTreeTime.count()/1000000 << " ms\n";

                // Print first 10 hash values corresponding to each sequence
                // MashPlacement::mashDeviceArrays.printSketchValues(10);
                MashPlacement::msaDeviceArrays.deallocateDeviceArrays();
                MashPlacement::kplacementDeviceArrays.deallocateDeviceArrays();
            }
        }
        else if (algo=="3"|| algo=="0"&&numSequences>=dc_thr){
            std::cerr<<"Using divide-and-conquer mode\n";
            int totalNumSequences = numSequences;
            int backboneSize = numSequences/20;
            params.batchSize = backboneSize;
            params.backboneSize = backboneSize;
            MashPlacement::msaDeviceArraysDC.allocateDeviceArraysDC(fourBitCompressedSeqs, seqLengths, numSequences, params);
            MashPlacement::kplacementDeviceArraysDC.allocateDeviceArraysDC(backboneSize, totalNumSequences);
            auto createArrayEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds createArrayTime = createArrayEnd - createArrayStart; 
            std::cerr << "Allocated in: " <<  createArrayTime.count()/1000000 << " ms\n";

            //Build Tree on Gpu
            auto createTreeStart = std::chrono::high_resolution_clock::now();
            MashPlacement::kplacementDeviceArraysDC.findBackboneTreeDC(params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, MashPlacement::msaDeviceArraysDC, MashPlacement::kplacementDeviceArraysHostDC);
            MashPlacement::kplacementDeviceArraysDC.findClustersDC(params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, MashPlacement::msaDeviceArraysDC, MashPlacement::kplacementDeviceArraysHostDC);
            MashPlacement::kplacementDeviceArraysDC.findClusterTreeDC(params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, MashPlacement::msaDeviceArraysDC, MashPlacement::kplacementDeviceArraysHostDC);

            auto createTreeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds createTreeTime = createTreeEnd - createTreeStart; 
            MashPlacement::kplacementDeviceArraysDC.printTreeDC(names, output_);
            std::cerr << "Tree Created in: " <<  createTreeTime.count()/1000000 << " ms\n";

            // Print first 10 hash values corresponding to each sequence
            // MashPlacement::mashDeviceArrays.printSketchValues(10);
            MashPlacement::msaDeviceArraysDC.deallocateDeviceArraysDC();
            MashPlacement::kplacementDeviceArraysDC.deallocateDeviceArraysDC();
        }
        else{
            std::cerr<<"Using conventional NJ\n";
            if(numSequences>=40000){
                std::cerr<<"Warning: forcing conventional NJ on large datasets might result in unexpected behavior\n";
            }
            MashPlacement::njDeviceArrays.getDismatrix(
                numSequences,params, MashPlacement::mashDeviceArrays, MashPlacement::matrixReader, MashPlacement::msaDeviceArrays
            );
            MashPlacement::njDeviceArrays.findNeighbourJoiningTree(names, output_);
            MashPlacement::msaDeviceArrays.deallocateDeviceArrays();
            MashPlacement::njDeviceArrays.deallocateDeviceArrays();
        }
    }
    else if (in == "r" && out == "t"){
        std::vector<std::string> seqs,names_, names;

        // Read Input Sequences (Fasta format)
        readSequences(vm, seqs, names_);
        size_t numSequences = seqs.size();
        names.resize(numSequences);
        std::vector<int> ids(numSequences);
        for(int i=0;i<numSequences;i++) ids[i]=i;
        std::mt19937 rnd(time(NULL));
        std::shuffle(ids.begin(),ids.end(),rnd);
        
        // Compress Sequences (2-bit compressor)
        auto compressStart = std::chrono::high_resolution_clock::now();
        // fprintf(stdout, "Compressing input sequence using two-bit encoding.\n");
        uint64_t ** twoBitCompressedSeqs = new uint64_t*[numSequences];
        uint64_t * seqLengths = new uint64_t[numSequences];
        tbb::parallel_for(tbb::blocked_range<int>(0, numSequences), [&](tbb::blocked_range<int> range){
        for (int idx_= range.begin(); idx_ < range.end(); ++idx_) {
            uint64_t i = static_cast<uint64_t>(idx_);
            uint64_t twoBitCompressedSize = (seqs[i].size()+31)/32;
            uint64_t * twoBitCompressed = new uint64_t[twoBitCompressedSize];
            twoBitCompressor(seqs[i], seqs[i].size(), twoBitCompressed);

            seqLengths[ids[i]] = seqs[i].size();
            twoBitCompressedSeqs[ids[i]] = twoBitCompressed;
            names[ids[i]] = names_[i];
        }});
        
        auto compressEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds compressTime = compressEnd - compressStart;
        // std::cout << "Compressed in: " <<  compressTime.count() << " ns\n";
        auto inputEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds inputTime = inputEnd - inputStart; 
        std::cerr << "Input in: " <<  inputTime.count()/1000000 << " ms\n";


        //Build Tree on Gpu
        if(algo=="1"||algo=="0"&&numSequences>=placement_thr&&numSequences<dc_thr){
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
            if(placemode=="-1"){
                std::cerr<<"Using exact placement mode\n";
                MashPlacement::placementDeviceArrays.allocateDeviceArrays(numSequences);
                auto createTreeStart = std::chrono::high_resolution_clock::now();
                MashPlacement::placementDeviceArrays.findPlacementTree(params, MashPlacement::mashDeviceArrays, MashPlacement::matrixReader, MashPlacement::msaDeviceArrays);
                auto createTreeEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds createTreeTime = createTreeEnd - createTreeStart; 
                MashPlacement::placementDeviceArrays.printTree(names, output_);
                std::cerr << "Tree Created in: " <<  createTreeTime.count()/1000000 << " ms\n";
                MashPlacement::mashDeviceArrays.deallocateDeviceArrays();
                MashPlacement::placementDeviceArrays.deallocateDeviceArrays();
            }
            else{
                std::cerr<<"Using k-closest placement mode\n";
                MashPlacement::kplacementDeviceArrays.allocateDeviceArrays(numSequences);
                auto createTreeStart = std::chrono::high_resolution_clock::now();
                MashPlacement::kplacementDeviceArrays.findPlacementTree(params, MashPlacement::mashDeviceArrays, MashPlacement::matrixReader, MashPlacement::msaDeviceArrays);
                auto createTreeEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds createTreeTime = createTreeEnd - createTreeStart; 
                MashPlacement::kplacementDeviceArrays.printTree(names, output_);
                std::cerr << "Tree Created in: " <<  createTreeTime.count()/1000000 << " ms\n";
                MashPlacement::mashDeviceArrays.deallocateDeviceArrays();
                MashPlacement::kplacementDeviceArrays.deallocateDeviceArrays();
            }
        }
        else if (algo=="3"||algo=="0"&&numSequences>=dc_thr){
            std::cerr<<"Using divide-and-conquer mode\n";
            
            int totalNumSequences = numSequences;
            int backboneSize = numSequences/100;
            params.batchSize = backboneSize;
            params.backboneSize = backboneSize;

            auto createArrayStart = std::chrono::high_resolution_clock::now();
            MashPlacement::mashDeviceArraysDC.allocateDeviceArraysDC(twoBitCompressedSeqs, seqLengths, numSequences, params);
            auto createArrayEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds createArrayTime = createArrayEnd - createArrayStart; 
            std::cerr << "Allocated in: " <<  createArrayTime.count()/1000000 << " ms\n";

            auto createSketchStart = std::chrono::high_resolution_clock::now();
            MashPlacement::mashDeviceArraysDC.sketchConstructionOnGpuDC(params, twoBitCompressedSeqs, seqLengths, numSequences);
            auto createSketchEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds createSketchTime = createSketchEnd - createSketchStart; 
            std::cerr << "Sketch Created in: " <<  createSketchTime.count()/1000000 << " ms\n";
            
            MashPlacement::kplacementDeviceArraysDC.allocateDeviceArraysDC(backboneSize, totalNumSequences);
            MashPlacement::kplacementDeviceArraysHostDC.allocateHostArraysDC(backboneSize, totalNumSequences);
            auto createTreeStart = std::chrono::high_resolution_clock::now();
            
            MashPlacement::kplacementDeviceArraysDC.findBackboneTreeDC(params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, MashPlacement::msaDeviceArraysDC, MashPlacement::kplacementDeviceArraysHostDC);
            MashPlacement::kplacementDeviceArraysDC.findClustersDC(params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, MashPlacement::msaDeviceArraysDC, MashPlacement::kplacementDeviceArraysHostDC);
            auto createTreeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds createTreeTime = createTreeEnd - createTreeStart;

            MashPlacement::kplacementDeviceArraysDC.findClusterTreeDC(params, MashPlacement::mashDeviceArraysDC, MashPlacement::matrixReader, MashPlacement::msaDeviceArraysDC, MashPlacement::kplacementDeviceArraysHostDC);
            MashPlacement::kplacementDeviceArraysDC.printTreeDC(names, output_);
            MashPlacement::kplacementDeviceArraysDC.deallocateDeviceArraysDC();
            MashPlacement::mashDeviceArraysDC.deallocateDeviceArraysDC();
            std::cerr << "Tree Created in: " <<  createTreeTime.count()/1000000 << " ms\n";
        }
        else{
            std::cerr<<"Using conventional NJ\n";
            if(numSequences>=40000){
                std::cerr<<"Warning: forcing conventional NJ on large datasets might result in unexpected behavior\n";
            }
            // Create arrays
            auto createArrayStart = std::chrono::high_resolution_clock::now();
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

            MashPlacement::njDeviceArrays.getDismatrix(
                numSequences,params, MashPlacement::mashDeviceArrays, MashPlacement::matrixReader, MashPlacement::msaDeviceArrays
            );
            MashPlacement::njDeviceArrays.findNeighbourJoiningTree(names, output_);
            MashPlacement::mashDeviceArrays.deallocateDeviceArrays();
            MashPlacement::njDeviceArrays.deallocateDeviceArrays();
        }

        // Print first 10 hash values corresponding to each sequence
        // MashPlacement::mashDeviceArrays.printSketchValues(10);

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
        if(algo=="1"||algo=="0"&&numSequences>=placement_thr&&numSequences<dc_thr){
            if(placemode=="-1"){
                std::cerr<<"Using exact placement mode\n";
                MashPlacement::placementDeviceArrays.allocateDeviceArrays(numSequences);
                MashPlacement::placementDeviceArrays.findPlacementTree(params, MashPlacement::mashDeviceArrays, MashPlacement::matrixReader, MashPlacement::msaDeviceArrays);
                MashPlacement::placementDeviceArrays.printTree(MashPlacement::matrixReader.name, output_);
            }
            else{
                std::cerr<<"Using k-closest placement mode\n";
                MashPlacement::kplacementDeviceArrays.allocateDeviceArrays(numSequences);
                MashPlacement::kplacementDeviceArrays.findPlacementTree(params, MashPlacement::mashDeviceArrays, MashPlacement::matrixReader, MashPlacement::msaDeviceArrays);
                MashPlacement::kplacementDeviceArrays.printTree(MashPlacement::matrixReader.name, output_);
            }
        } else if (algo=="3"|| algo=="0"&&numSequences>=dc_thr){
            std::cerr<<"Divide-and-conquer mode not supported with input matrix\n";
            exit(1);
        }
        else{
            std::cerr<<"Using conventional NJ\n";
            if(numSequences>=40000){
                std::cerr<<"Warning: forcing conventional NJ on large datasets might result in unexpected behavior\n";
            }
            MashPlacement::njDeviceArrays.getDismatrix(
                numSequences,params, MashPlacement::mashDeviceArrays, MashPlacement::matrixReader, MashPlacement::msaDeviceArrays
            );
            MashPlacement::njDeviceArrays.findNeighbourJoiningTree(MashPlacement::matrixReader.name, output_);
            MashPlacement::njDeviceArrays.deallocateDeviceArrays();
        }
        fclose(filePtr);
    }
    else{
        printf("Invalid input-output combinations!!!!!\n");
        exit(1);
    }
    return 0;
}