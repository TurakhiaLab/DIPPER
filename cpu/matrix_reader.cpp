/**
 * cpu/matrix_reader.cpp
 *
 * PHYLIP distance matrix reader. Allocates host buffer, reads names and
 * distances per row; distConstructionOnGpu fills d_mashDist for a given row.
 */

#include "mash_placement.hpp"

#include <stdio.h>
#include <cstdio>
#include <queue>
#include <chrono>
#include <iostream>

/** Allocate host arrays and attach FILE for PHYLIP matrix. */
void MashPlacement::MatrixReader::allocateDeviceArrays(int num, FILE *fPtr)
{
    numSequences = num;
    h_dist = new double[num];
    name.resize(num);
    buffer = new char[num * 20];
    filePtr = fPtr;
}

/** Read row rowId (name + distances), copy into d_mashDist. */
void MashPlacement::MatrixReader::distConstructionOnGpu(Param &params, int rowId, double *d_mashDist)
{
    fgets(buffer, numSequences * 20, filePtr);
    char *p = buffer;
    name[rowId] = "";
    while (1)
    {
        char c = *p;
        p++;
        if (c == '\t' || c == '\n' || c == ' ')
            break;
        name[rowId] += c;
    }
    if (rowId == 0)
        return;
    for (int i = 0; i < rowId; i++)
    {
        std::string num = "";
        while (1)
        {
            char c = *p;
            p++;
            if (c == '\t' || c == '\n' || c == ' ')
                break;
            num += c;
        }
        h_dist[i] = stof(num);
    }
    for (size_t i = 0; i < rowId; ++i)
        d_mashDist[i] = h_dist[i];
    // cudaMemcpy(d_mashDist, h_dist, rowId*sizeof(double), cudaMemcpyHostToDevice);
}