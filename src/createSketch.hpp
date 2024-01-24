#ifndef CREATESKETCH_HPP
#define CREATESKETCH_HPP

#include <queue>
#include <string>
#include <cstring>
#include <vector>
#include <iostream>
#include <fstream>
#include <unordered_map>

typedef uint32_t hash_t; 
#define INF 10000000
#define PRINT_HASH_STATUS true


class sketch
{
public:
    std::unordered_map<std::string, std::priority_queue<hash_t>> sketchMap;
    sketch(std::unordered_map<std::string, std::string>& seqs, int& k, size_t& sketchSize, int& threshold);
};


#endif