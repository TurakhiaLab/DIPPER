/**
 * cpu/MSA.cpp
 *
 * CPU aligned-sequence arrays and distance construction. 4-bit compressed MSA,
 * seq lengths; pairwise distances (uncorrected, JC, Tajima-Nei, K2P, Tamura,
 * Jin-Nei) via TBB over target rows.
 */

#include "mash_placement.hpp"

#include <stdio.h>
#include <queue>
#include <chrono>
#include <cmath>
#include <iostream>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

/** Allocate host MSA arrays, flatten 4-bit compressed data. */
void MashPlacement::MSADeviceArrays::allocateDeviceArrays(uint64_t ** h_compressedSeqs, uint64_t * h_seqLengths, size_t num, Param& params)
{
    // cudaError_t err;

    numSequences = int(num);
    seqLen = h_seqLengths[0];
    // Allocate memory
    d_seqLengths = new uint64_t [1ll*numSequences];
    /*
    err = cudaMalloc(&d_seqLengths, 1ll*numSequences*sizeof(uint64_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    */
    // std::cerr<<"????????\n";
    /* Flatten data */
    uint64_t flatStringLength=0;
    for (size_t i =0; i<numSequences; i++) flatStringLength+= params.isProtein ? (h_seqLengths[i]+7)/8 : (h_seqLengths[i]+15)/16;
    uint64_t * h_flattenCompressSeqs = new uint64_t[flatStringLength];
    flatStringLength=0;
    for (size_t i =0; i<numSequences; i++) 
    {
        uint64_t flatStringLengthLocal = params.isProtein ? (h_seqLengths[i]+7)/8 : (h_seqLengths[i]+15)/16;
        flatStringLength+=flatStringLengthLocal;
        for (size_t j=0; j<flatStringLengthLocal;j++)  
        {
            h_flattenCompressSeqs[j] = h_compressedSeqs[i][j];
            // if (i==9) printf("%u\n",h_flattenCompressSeqs[j]); 
        }
        h_flattenCompressSeqs += flatStringLengthLocal;
    }
    h_flattenCompressSeqs -= flatStringLength;
    params.flatStringLength = flatStringLength;
    // std::cerr<<"?????????\n";

    // err = cudaMalloc(&d_compressedSeqs, 1ll*flatStringLength*sizeof(uint64_t));
    // if (err != cudaSuccess)
    // {
    //     fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
    //     exit(1);
    // }

    d_compressedSeqs = new uint64_t[1ll*flatStringLength];

    for (size_t i = 0; i < numSequences; i++) {
        d_seqLengths[i] = h_seqLengths[i];
    }
    for (size_t i = 0; i < flatStringLength; i++) {
        d_compressedSeqs[i] = h_flattenCompressSeqs[i];
    }


    // Transfer data
    /*
    err = cudaMemcpy(d_seqLengths, h_seqLengths, 1ll*numSequences*sizeof(uint64_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) 
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
        exit(1);
    }

    err = cudaMemcpy(d_compressedSeqs, h_flattenCompressSeqs, 1ll*flatStringLength*sizeof(uint64_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) 
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
        exit(1);
    }
    cudaDeviceSynchronize();
    */
}

void MashPlacement::MSADeviceArrays::deallocateDeviceArrays(){
    delete [] d_compressedSeqs;
    delete [] d_seqLengths;
    // cudaFree(d_compressedSeqs);
    // cudaFree(d_seqLengths);
}

#define DIST_UNCORRECTED 1
#define DIST_JUKESCANTOR 2
#define DIST_TAJIMANEI 3
#define DIST_KIMURA2P 4
#define DIST_TAMURA 5
#define DIST_JINNEI 6

/** Count useful (non-gap) sites and matches between two compressed rows. */
void calculateParams(int tarRowId, int curRowId, int seqLen, uint64_t * compressedSeqs, int & useful, int & match){
    int compLen=(seqLen+15)/16;
    long long px=1ll*curRowId*compLen, py=1ll*tarRowId*compLen;
    for(int i=0;i<compLen;i++){
        long long vt=compressedSeqs[px+i], vc=compressedSeqs[py+i];
        for(int j=0;j<16&&i*16+j<seqLen;j++){
            int et=(vt>>(j*4))&15, ec=(vc>>(j*4))&15;
            if(et<4||ec<4) useful++;
            if(et<4&&et==ec) match++;
        }
    }
}


void calculateParamsParallel(int tarRowId, int curRowId, int seqLen, uint64_t * compressedSeqs, int & useful, int & match, bool isProtein=false){
    int compLen = isProtein ? (seqLen+7)/8 : (seqLen+15)/16;
    long long px=1ll*curRowId*compLen, py=1ll*tarRowId*compLen;
    for (int i = 0; i < compLen; ++i) {
        long long vt=compressedSeqs[px+i], vc=compressedSeqs[py+i];
        if (isProtein) {
            for(int j=0;j<8&&i*8+j<seqLen;j++){
                int et=(vt>>(j*8))&255, ec=(vc>>(j*8))&255;
                if(et<20||ec<20) useful++;
                if(et<20&&et==ec) match++;
            }
        } else {
            for(int j=0;j<16&&i*16+j<seqLen;j++){
                int et=(vt>>(j*4))&15, ec=(vc>>(j*4))&15;
                if(et<4||ec<4) useful++;
                if(et<4&&et==ec) match++;
            }
        }
    }
}


void calculateParams_TJ(int tarRowId, int curRowId, int seqLen, uint64_t * compressedSeqs, int * frac, int &tot, int &match, int * pr){
    int compLen=(seqLen+15)/16;
    long long px=1ll*curRowId*compLen, py=1ll*tarRowId*compLen;
    for(int i=0;i<compLen;i++){
        long long vt=compressedSeqs[px+i], vc=compressedSeqs[py+i];
        for(int j=0;j<16&&i*16+j<seqLen;j++){
            int et=(vt>>(j*4))&15, ec=(vc>>(j*4))&15;
            if(et>=4||ec>=4) continue;
            frac[ec]++, frac[et]++, tot++;
            if(ec>et){
                int temp=ec;
                ec=et,et=temp;
            }
            if(ec==et) match++;
            if(ec==0&&et==2) pr[0]++;
            else if(ec==0&&et==3) pr[1]++;
            else if(ec==1&&et==2) pr[2]++;
            else if(ec==1&&et==3) pr[3]++;
        }
    }
}

void calculateParamsParallel_TJ(int tarRowId, int curRowId, int seqLen, uint64_t * compressedSeqs, int * frac, int &tot, int &match, int * pr){
    int compLen=(seqLen+15)/16;
    long long px=1ll*curRowId*compLen, py=1ll*tarRowId*compLen;

    for (int i=0; i<compLen; i++) {
        long long vt=compressedSeqs[px+i], vc=compressedSeqs[py+i];
        for(int j=0;j<16&&i*16+j<seqLen;j++){
            int et=(vt>>(j*4))&15, ec=(vc>>(j*4))&15;
            if(et>=4||ec>=4) continue;
            if (ec == 0) frac[0]++;
            else if (ec == 1) frac[1]++;
            else if (ec == 2) frac[2]++;
            else if (ec == 3) frac[3]++;
            if (et == 0) frac[0]++;
            else if (et == 1) frac[1]++;
            else if (et == 2) frac[2]++;
            else if (et == 3) frac[3]++;
            tot++;
            if(ec>et){
                int temp=ec;
                ec=et,et=temp;
            }
            if(ec==et) match++;
            if(ec==0&&et==2) pr[0]++;
            else if(ec==0&&et==3) pr[1]++;
            else if(ec==1&&et==2) pr[2]++;
            else if(ec==1&&et==3) pr[3]++;
        }
    }
}

void calculateParams_K2P(int tarRowId, int curRowId, int seqLen, uint64_t * compressedSeqs, int &p, int &q, int &tot){
    int compLen=(seqLen+15)/16;
    long long px=1ll*curRowId*compLen, py=1ll*tarRowId*compLen;
    for(int i=0;i<compLen;i++){
        long long vt=compressedSeqs[px+i], vc=compressedSeqs[py+i];
        for(int j=0;j<16&&i*16+j<seqLen;j++){
            int et=(vt>>(j*4))&15, ec=(vc>>(j*4))&15;
            if(et>=4||ec>=4) continue;
            tot++;
            if(et==ec) continue;
            if(et%2==ec%2) p++;
            else q++;
        }
    }
}

void calculateParamsParallel_K2P(int tarRowId, int curRowId, int seqLen, uint64_t * compressedSeqs, int &p, int &q, int &tot){
    int compLen=(seqLen+15)/16;
    long long px=1ll*curRowId*compLen, py=1ll*tarRowId*compLen;
    for (int i=0; i<compLen; i++) {
        long long vt=compressedSeqs[px+i], vc=compressedSeqs[py+i];
        for(int j=0;j<16&&i*16+j<seqLen;j++){
            int et=(vt>>(j*4))&15, ec=(vc>>(j*4))&15;
            if(et>=4||ec>=4) continue;
            tot++;
            if(et==ec) continue;
            if(et%2==ec%2) p++;
            else q++;
        }
    }
}

void calculateParams_TAMURA(int tarRowId, int curRowId, int seqLen, uint64_t * compressedSeqs, int &p, int &q, int &tot, int &gc1, int &gc2){
    int compLen=(seqLen+15)/16;
    long long px=1ll*curRowId*compLen, py=1ll*tarRowId*compLen;
    for(int i=0;i<compLen;i++){
        long long vt=compressedSeqs[px+i], vc=compressedSeqs[py+i];
        for(int j=0;j<16&&i*16+j<seqLen;j++){
            int et=(vt>>(j*4))&15, ec=(vc>>(j*4))&15;
            if(et>=4||ec>=4) continue;
            tot++;
            if(et==ec) continue;
            if(et%2==ec%2) p++;
            else q++;
            if(ec==1||ec==2) gc1++;
            if(et==1||et==2) gc2++;
        }
    }
}

void calculateParamsParallel_TAMURA(int tarRowId, int curRowId, int seqLen, uint64_t * compressedSeqs, int &p, int &q, int &tot, int &gc1, int &gc2){
    int compLen=(seqLen+15)/16;
    long long px=1ll*curRowId*compLen, py=1ll*tarRowId*compLen;
    for (int i=0; i<compLen; i++) {
        long long vt=compressedSeqs[px+i], vc=compressedSeqs[py+i];
        for(int j=0;j<16&&i*16+j<seqLen;j++){
            int et=(vt>>(j*4))&15, ec=(vc>>(j*4))&15;
            if(et>=4||ec>=4) continue;
            tot++;
            if(et==ec) continue;
            if(et%2==ec%2) p++;
            else q++;
            if(ec==1||ec==2) gc1++;
            if(et==1||et==2) gc2++;
        }
    }
}

/** Build distance row rowId; dispatch by distanceType (TBB over target rows). */
void MSADistConstruction(
    int rowId,
    uint64_t * compressedSeqs,
    double * dist,
    int seqLen,
    int numSequences,
    int distanceType,
    bool isProtein = false
){
    int threadNum = 1024;
    if(distanceType==DIST_UNCORRECTED||distanceType==DIST_JUKESCANTOR) {
        tbb::parallel_for(tbb::blocked_range<int>(0, rowId), [&](const tbb::blocked_range<int>& range) {
        for (int blockID = range.begin(); blockID < range.end(); ++blockID) {
        // for (int blockID = 0; blockID < rowId; ++blockID) {
            int useful=0, match=0;
            calculateParamsParallel(rowId, blockID, seqLen, compressedSeqs, useful, match, isProtein);
            double uncor = (useful > 0) ? (1.0 - double(match) / useful) : 1.0;
            if (distanceType == DIST_UNCORRECTED) {
                dist[blockID] = uncor;
            } else {
                if (isProtein) {
                    double denom = 1.0 - uncor;
                    if (denom <= 1e-6) denom = 1e-6;
                    dist[blockID] = 2.0 * (1.0 / sqrt(denom)) - 1.0;
                } else {
                    double arg = 1.0 - uncor / 0.75;
                    if (arg <= 1e-6) arg = 1e-6;
                    dist[blockID] = -0.75 * log(arg);
                }
            }
        }
        });
    }
    else if(distanceType==DIST_TAJIMANEI){
        tbb::parallel_for(tbb::blocked_range<int>(0, rowId), [&](const tbb::blocked_range<int>& range) {
        for (int blockID = range.begin(); blockID < range.end(); ++blockID) {
            int frac[4]={},pr[4]={},tot=0,match=0;
            double fr[4]={};
            calculateParamsParallel_TJ(rowId, blockID, seqLen, compressedSeqs, frac, tot, match, pr);
            if (tot == 0) {
                dist[blockID] = 1.0;
                continue;
            }
            for(int i=0;i<4;i++) fr[i]=double(frac[i])/tot/2.0;
            double h=0;
            h+=0.5*pr[0]*fr[0]*fr[2];
            h+=0.5*pr[1]*fr[0]*fr[3];
            h+=0.5*pr[2]*fr[1]*fr[2];
            h+=0.5*pr[3]*fr[1]*fr[3];
            double D=double(tot-match)/tot;
            double b=0.5*(1.0-fr[0]*fr[0]-fr[2]*fr[2]+D*D/h);
            double arg = (b > 1e-9) ? (1.0 - D / b) : 1e-6;
            if (arg <= 1e-6) arg = 1e-6;
            dist[blockID] = -b * log(arg);
        }
        });
    }
    else if(distanceType==DIST_KIMURA2P||distanceType==DIST_JINNEI){
        tbb::parallel_for(tbb::blocked_range<int>(0, rowId), [&](const tbb::blocked_range<int>& range) {
        for (int blockID = range.begin(); blockID < range.end(); ++blockID) {
            int p=0,q=0,tot=0;
            calculateParamsParallel_K2P(rowId, blockID, seqLen, compressedSeqs, p, q, tot);
            if (tot == 0) {
                dist[blockID] = 1.0;
                continue;
            }
            double pp=double(p)/tot,qq=double(q)/tot;
            if(distanceType==DIST_KIMURA2P) {
                double a1 = 1.0 - 2.0*pp - qq;
                double a2 = 1.0 - 2.0*qq;
                if (a1 <= 1e-6) a1 = 1e-6;
                if (a2 <= 1e-6) a2 = 1e-6;
                dist[blockID] = -0.5 * log(a1 * sqrt(a2));
            } else {
                double d1 = 1.0 - 2.0*pp - qq;
                double d2 = 1.0 - 2.0*qq;
                if (d1 <= 1e-6) d1 = 1e-6;
                if (d2 <= 1e-6) d2 = 1e-6;
                dist[blockID] = 0.5 * (1.0 / d1 + 0.5 / d2 - 1.5);
            }
        }
        });
    }
    else if(distanceType==DIST_TAMURA){
        tbb::parallel_for(tbb::blocked_range<int>(0, rowId), [&](const tbb::blocked_range<int>& range) {
        for (int blockID = range.begin(); blockID < range.end(); ++blockID) {
            int p=0,q=0,tot=0,gc1=0,gc2=0;
            calculateParamsParallel_TAMURA(rowId, blockID, seqLen, compressedSeqs, p, q, tot, gc1, gc2);
            if (tot == 0) {
                dist[blockID] = 1.0;
                continue;
            }
            double pp=double(p)/tot,qq=double(q)/tot, c=double(gc1)/tot+double(gc2)/tot-2*double(gc1)*double(gc2)/tot/tot;
            if (c <= 1e-9) c = 1e-9;
            double a1 = 1.0 - pp/c - qq;
            double a2 = 1.0 - 2.0*qq;
            if (a1 <= 1e-6) a1 = 1e-6;
            if (a2 <= 1e-6) a2 = 1e-6;
            dist[blockID] = -c * log(a1) - 0.5 * (1.0 - c) * log(a2);
        }
        });
    }
    else {
        tbb::parallel_for(tbb::blocked_range<int>(0, rowId), [&](const tbb::blocked_range<int>& range) {
        for (int blockID = range.begin(); blockID < range.end(); ++blockID) {
            dist[blockID] = 0.0;
        }
        });
    }
}

/** Launch MSADistConstruction for row rowId. */
void MashPlacement::MSADeviceArrays::distConstructionOnGpu(Param& params, int rowId, double* d_mashDist) const{
    MSADistConstruction (
        rowId, 
        d_compressedSeqs, 
        d_mashDist, 
        seqLen,
        numSequences,
        params.distanceType,
        params.isProtein
    );
}


