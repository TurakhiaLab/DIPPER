#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <vector>
#include <bits/stdc++.h>
#include <boost/program_options.hpp> 

#ifndef PL_CUH
#include "../src/placement.cuh"
#endif


namespace po = boost::program_options;
int main(int argc, char** argv) {
    std::string matrixFileName;
    char *buffer = new char[20000000];
    int numSequences, bd=2, tk=5;
    std::cin>>numSequences;
    fgets(buffer,20000000,stdin);
    std::vector <std::string> name(numSequences);
    double *dismatrix = new double[1LL*numSequences*numSequences];
    double *len = new double[numSequences*8]();
    double *u = new double[numSequences*2]();
    int *head = new int[numSequences*2]();
    int *e = new int[numSequences*8]();
    int *nxt = new int[numSequences*8]();
    int *id = new int[numSequences];
    int *f = new int[numSequences]();
    int *belong = new int[numSequences*8];
    memset(head, -1, numSequences*2*sizeof(int));
    memset(e,-1,numSequences*8*sizeof(int));
    memset(nxt,-1,numSequences*8*sizeof(int));
    memset(belong,-1,numSequences*8*sizeof(int));
    int idx = 0;
    std::mt19937 rnd(time(NULL));
    auto inputstart = std:: chrono::high_resolution_clock::now();
    std::string num;
    long long cc=0;
    for(int i=0;i<numSequences;i++){
        fgets(buffer,20000000,stdin);
        char *p=buffer;
        // while(*p!='\n') std::cout<<"*p";std::cout<<'\n';
        while(1){
            char c=*p;
            p++;
            if(c=='\t'||c=='\n'||c==' ') break;
            name[i]+=c;
        }
        // name[i]="Name_i";
        id[i]=i;
        for(int j=0;j<numSequences;j++){
            num="";
            while(1){
                char c=*p;
                p++;
                if(c=='\t'||c=='\n'||c==' ') break;
                num+=c;
            }
            // temp = double(rnd())/UINT_MAX;
            // if(i==j) temp=0;
            dismatrix[cc++]=stof(num);
        }
    }
    auto inputend = std:: chrono::high_resolution_clock::now();
    std::chrono::nanoseconds inputTime = inputend - inputstart;
    std::cerr << "Input in: "<< inputTime.count()/1000000 <<" ms\n";
    std::function<void(int,int)>  print=[&](int node,int from){
        // printf("####%d####",node);
        if(head[node]!=-1&&nxt[head[node]]!=-1){
            std::vector <int> v;
            for(int i=head[node];i!=-1;i=nxt[i]) if(e[i]!=from) v.push_back(i);
            printf("(");
            for(int i=0;i<v.size();i++){
                // printf("????%d????",v[i]);
                print(e[v[i]],node);
                printf(":");
                printf("%.5g%c",len[v[i]],i+1==v.size()?')':',');
            }
        }
        else std::cout<<name[node];
    };
    int cnt = numSequences;
    for (int i = 0; i < bd - 1; i++) {
        if (i == bd - 2) {
            std::vector<int> vv;
            for (int j = 0; j < bd; j++)
                if (!f[j]) vv.push_back(j);
            double dd = dismatrix[vv[0]*numSequences+vv[1]];
            e[idx] = id[vv[0]], len[idx] = dd/2, nxt[idx] = head[cnt], head[cnt] = idx, belong[idx] = cnt, idx++;
            e[idx] = id[vv[1]], len[idx] = dd/2, nxt[idx] = head[cnt], head[cnt] = idx, belong[idx] = cnt, idx++;
            e[idx] = cnt, len[idx] = dd/2, nxt[idx] = head[id[vv[0]]], head[id[vv[0]]] = idx, belong[idx] = id[vv[0]], idx++;            
            e[idx] = cnt, len[idx] = dd/2, nxt[idx] = head[id[vv[1]]], head[id[vv[1]]] = idx, belong[idx] = id[vv[1]], idx++; 
            cnt++;
            continue;
        }
        for (int j = 0; j < bd; j++) {
            u[j] = 0;
            for (int k = 0; k < bd; k++)
                if (f[j] == 0 && f[k] == 0 && j != k)
                    u[j] += dismatrix[j*numSequences+k];
            u[j] /= (bd - i - 2);
        }
        double mn = 10000;
        int x = -1, y = -1;
        for (int j = 0; j < bd; j++) {
            if (f[j]) continue;
            for (int k = 0; k < bd; k++) {
                if (j == k || f[k]) continue;
                double val = dismatrix[j*numSequences+k] - u[j] - u[k];
                if (val < mn) mn = val, x = j, y = k;
            }
        }
        double blx = (dismatrix[x*numSequences+y] + u[x] - u[y]) * 0.5;
        double bly = dismatrix[x*numSequences+y] - blx;
        if (blx < 0) bly += blx, blx = 0;
        if (bly < 0) blx += bly, bly = 0;
        e[idx] = id[x], len[idx] = blx, nxt[idx] = head[cnt], head[cnt] = idx, belong[idx] = cnt, idx++;
        e[idx] = id[y], len[idx] = bly, nxt[idx] = head[cnt], head[cnt] = idx, belong[idx] = cnt, idx++;
        e[idx] = cnt, len[idx] = blx, nxt[idx] = head[id[x]], head[id[x]] = idx, belong[idx] = id[x], idx++;            
        e[idx] = cnt, len[idx] = bly, nxt[idx] = head[id[y]], head[id[y]] = idx, belong[idx] = id[y], idx++; 
        // std::cerr << cnt<<" "<<id[x] << " " << id[y] << " " <<mn<<" "<<idx<<"\n";
        //printf("%.10lf %.10lf\n",mn, dis[x][y]);
        f[y] = 1, id[x] = cnt++;
        for (int j = 0; j < bd; j++)
            if (j != x && j != y && !f[j]) {
                double val = (dismatrix[x*numSequences+j] + dismatrix[y*numSequences+j] - dismatrix[x*numSequences+y]) / 2.0;
                //if(val<0) val=0;
                dismatrix[x*numSequences+j] = dismatrix[j*numSequences+x] = val;
            }
        // print(cnt-1, -1),printf(";\n");
    }
    // Root of the tree has an index of d_numSequneces*2-2
    // print(cnt-1, -1),printf(";\n");
    auto start = std:: chrono::high_resolution_clock::now();
    placement::deviceArrays.allocateDeviceArrays(
        tk,
        numSequences, 
        dismatrix, 
        head, 
        e,
        nxt,
        idx,
        bd,
        len,
        belong
    );
    placement::findPlacementTree(
        placement::deviceArrays.numSequences,
        placement::deviceArrays.bd,
        placement::deviceArrays.idx,
        placement::deviceArrays.d_dist,
        placement::deviceArrays.d_head,
        placement::deviceArrays.d_e,
        placement::deviceArrays.d_len,
        placement::deviceArrays.d_nxt,
        placement::deviceArrays.d_belong,
        placement::deviceArrays.d_closest_dis,
        placement::deviceArrays.d_closest_id,
        name,
        tk
    );
    auto end = std:: chrono::high_resolution_clock::now();
    placement::deviceArrays.deallocateDeviceArrays();
    std::chrono::nanoseconds placementTime = end - start;
    std::cerr << "Calculation in: "<< placementTime.count()/1000000 <<" ms\n";
    return 0;
}