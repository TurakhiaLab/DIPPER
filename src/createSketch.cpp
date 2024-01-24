#ifndef CREATESKETCH_HPP
#include "createSketch.hpp"
#endif

#ifndef MURMURHASH3_HPP
#include "../src/murmurHash3.hpp"
#endif

#include <limits>

void sketchHelper(const char* seq, int& k, size_t& sketchSize, int& threshold, std::priority_queue<hash_t>& sketch)
{
    std::unordered_map<hash_t, int> hashCount;
    hash_t maxHashValue=std::numeric_limits<uint32_t>::max();

    size_t strLen=std::strlen(seq);

    for (size_t i=0; i<strLen; i+=k)
    {
        int len = i+k<strLen? k: strLen-i;

        // for (int z=0; z<len;z++) std::cout << (seq+i+z)[0]; std::cout << "\n";
        
        char out[4]; // Only valid for uint32_t
        hash_t seed = 54;

        MurmurHash3_x86_32  (seq+i, len, seed, out);

        hash_t hash = *((uint32_t *)out);

        if (hash>=maxHashValue) continue;
        else if (hashCount.find(hash) != hashCount.end()) hashCount[hash] += 1;
        else hashCount[hash] = 1;

        if (hashCount[hash] >= threshold)
        {
            if (sketch.size() == sketchSize) sketch.pop();

            sketch.push(hash);
            hashCount.erase(hash);
            if (sketch.size() == sketchSize)
                maxHashValue=sketch.top();
        }

    }

    // If no of hashes in sketch is less than sketchSize, insert errorneous k-mers?

}

sketch::sketch(std::unordered_map<std::string, std::string>& seqs, int& k, size_t& sketchSize, int& threshold)
{

    for (auto& seqData: seqs)
    {
        if (PRINT_HASH_STATUS) std::cout<<"Building sketch for " << seqData.first;
        std::priority_queue<hash_t> sketch;
        const char* seq = (seqData.second).c_str();
        sketchHelper(seq, k, sketchSize, threshold, sketch);
        if (PRINT_HASH_STATUS) std::cout << "-- Sketch size: " << sketch.size() << "\n";
        sketchMap[seqData.first] = sketch;
    }

}


