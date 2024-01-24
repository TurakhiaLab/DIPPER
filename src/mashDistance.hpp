#ifndef MASHDISTANCE_HPP
#define MASHDISTANCE_HPP

#ifndef CREATESKETCH_HPP
#include "createSketch.hpp"
#endif

#include <bits/stdc++.h>

float jaccardEstimate(std::priority_queue<hash_t>A, std::priority_queue<hash_t>B);
float mashDistance(std::priority_queue<hash_t>&A, std::priority_queue<hash_t>&B, int k, size_t sketchSize);

#endif