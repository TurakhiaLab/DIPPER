#ifndef DISTANCEMATRIX_HPP
#include "distanceMatrix.hpp"
#endif

void distanceMatrix(std::unordered_map<std::string, std::priority_queue<hash_t>>& sketch, int k, size_t sketchSize, std::vector<float>& distMatrix, std::vector<std::string>& distMatrixSeqOrder)
{
    size_t sequenceCount = sketch.size();
    size_t matrixSize = sequenceCount*(sequenceCount - 1)/2;
    distMatrix.resize(matrixSize);
    distMatrixSeqOrder.resize(sketch.size());

    std::vector<std::priority_queue<hash_t>> sketchLocal;
    int count = 0;
    for (auto& s: sketch) 
    {
        // std::cout << count++ << "\t" << s.first << "\n";
        distMatrixSeqOrder[count++] = s.first;
        sketchLocal.push_back(s.second);
    }
    int matCount = 0;
    for (size_t i=0; i<sequenceCount; i++)
    {
        for (size_t j=i+1; j<sequenceCount; j++)
        {
            // std::cout << i << "\t" << j << "\t";
            distMatrix[matCount] = mashDistance(sketchLocal[i], sketchLocal[j], k, sketchSize);
            matCount++;
        }
    }

}
