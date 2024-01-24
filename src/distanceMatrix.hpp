#ifndef DISTANCEMATRIX_HPP
#define DISTANCEMATRIX_HPP

#ifndef MASHDISTANCE_HPP
#include "mashDistance.hpp"
#endif


void distanceMatrix(std::unordered_map<std::string, std::priority_queue<hash_t>>& sketch, int k, size_t sketchSize, std::vector<float>& distMatrix, std::vector<std::string>& distMatrixSeqOrder);


#endif