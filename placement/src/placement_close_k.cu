#include "mash_placement.cuh"

#include <stdio.h>
#include <queue>
#include <thrust/sort.h>
#include <thrust/scan.h>
#include <thrust/binary_search.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <chrono>
#include <iostream>
#include <cub/cub.cuh>

void MashPlacement::KPlacementDeviceArrays::allocateDeviceArrays(size_t num){
    cudaError_t err;
    numSequences = int(num);
    bd = 2, idx = 0;
    err = cudaMalloc(&d_dist, numSequences*sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_head, numSequences*2*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_e, numSequences*8*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_len, numSequences*8*sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_nxt, numSequences*8*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_belong, numSequences*8*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_closest_dis, numSequences*20*sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_closest_id, numSequences*20*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
}

__global__ void initialize(
    int lim,
    int nodes,
    double * d_closest_dis,
    int * d_closest_id,
    int * head,
    int * nxt,
    int * belong,
    int * e
){
    int tx=threadIdx.x,bs=blockDim.x;
    int bx=blockIdx.x,gs=gridDim.x;
    int idx=tx+bs*bx;
    if(idx<lim){
        for(int i=0;i<5;i++){
            d_closest_dis[idx*5+i]=2;
            d_closest_id[idx*5+i]=-1;
        }
        nxt[idx] = -1;
        e[idx] = -1;
        belong[idx] = -1;
    }
    if(idx<nodes) head[idx] = -1;
}

struct compare_tuple {
  __host__ __device__
  bool operator()(thrust::tuple<int,double,double> lhs, thrust::tuple<int,double,double> rhs)
  {
    return thrust::get<2>(lhs) < thrust::get<2>(rhs);
    //Always find the tuple whose third value (the criteria we want to minimize) is minimized
  }
};
/*
Three variables in tuple:
ID of branch in linked list,
distance to new node inserted on branch from starting vertex (belong[id]),
distance from new node inserted on branch to new node inserted outside branch
*/

__global__ void calculateBranchLength(
    int num, // should be bd, not numSequences 
    int * head,
    int * nxt,
    double * dis, 
    int * e, 
    double * len, 
    int * belong,
    thrust::tuple<int,double,double> * minPos,
    int lim,
    double * closest_dis,
    int * closest_id,
    int totSeqNum
){
    int tx=threadIdx.x,bs=blockDim.x,bx=blockIdx.x,gs=gridDim.x;
    int idx=tx+bs*bx;
    if(idx>=lim) return;
    if(idx>=num*4-4||belong[idx]<e[idx]){
        thrust::tuple <int,double,double> minTuple(0,0,2);
        minPos[bx*bs+tx]=minTuple;
        return;
    }
    int x=belong[idx],oth=e[idx];
    int eid=idx,otheid;
    double dis1=0, dis2=0, val;
    for(int i=0;i<5;i++)
        if(closest_id[eid*5+i]!=-1){
            val = dis[closest_id[eid*5+i]]-closest_dis[eid*5+i];
            if(val>dis1) dis1=val;
        }
    otheid=head[oth];
    while(e[otheid]!=x) assert(otheid!=-1),otheid=nxt[otheid];
    for(int i=0;i<5;i++)
        if(closest_id[otheid*5+i]!=-1){
            val = dis[closest_id[otheid*5+i]]-closest_dis[otheid*5+i];
            if(val>dis2) dis2=val;
        }
    double additional_dis=(dis1+dis2-len[eid])/2;
    if(additional_dis<0) additional_dis=0;
    dis1-=additional_dis,dis2-=additional_dis;
    if(dis1<0) dis1=0;
    if(dis2<0) dis2=0;
    if(dis1>len[eid]) additional_dis+=dis1-len[eid],dis1=len[eid];
    if(dis2>len[eid]) additional_dis+=dis2-len[eid],dis2=len[eid];
    // assert(dis1+dis2-1e-6<=len[eid]);
    double rest=len[eid]-dis1-dis2;
    dis1+=rest/2,dis2+=rest/2;
    thrust::tuple <int,double,double> minTuple(eid,dis1,additional_dis);
    minPos[bx*bs+tx]=minTuple;
}

__global__ void updateClosestNodes(
    int * head,
    int * nxt,
    int * e,
    double * len,
    double * closest_dis,
    int * closest_id,
    int x,
    int * id,
    int * from,
    double * dis
){
    int l=0,r=-1;
    id[++r]=x,dis[x]=0,from[x]=-1;
    while(l<=r){
        int node=id[l],fb=from[l];
        double d=dis[l];
        l++;
        for(int i=head[node];i!=-1;i=nxt[i]){
            if(e[i]==fb) continue;
            for(int j=0;j<5;j++){
                double nowd=closest_dis[i*5+j];
                if(nowd>d){
                    for(int k=4;k>j;k--){
                        closest_dis[i*5+k]=closest_dis[i*5+k-1];
                        closest_id[i*5+k]=closest_id[i*5+k-1];
                    }
                    closest_dis[i*5+j]=d;
                    closest_id[i*5+j]=x;
                    id[++r]=e[i],dis[r]=d+len[i],from[r]=node;
                    break;
                }
            }
        }
    }
}

__global__ void updateTreeStructure(
    int * head,
    int * nxt,
    int * e,
    double * len,
    double * closest_dis,
    int * closest_id,
    int * belong,
    int eid,
    double fracLen,
    double addLen,
    int placeId, // Id of the newly placed node
    int edgeCount, // Position to insert a new edge in linked list
    int numSequences
){
    int middle=placeId+numSequences-1, outside=placeId;
    int x=belong[eid],y=e[eid];
    double originalDis=len[eid];
    int xe,ye;
    for(int i=head[x];i!=-1;i=nxt[i])
        if(e[i]==y){
            e[i]=middle,len[i]=fracLen,xe=i;
            break;
        }
    for(int i=head[y];i!=-1;i=nxt[i])
        if(e[i]==x){
            e[i]=middle,len[i]-=fracLen,ye=i;
            break;
        }
    /*
    Need to update:
    e, len, nxt, head, belong, closest_dis, closest_id
    */
    //middle -> x
    e[edgeCount]=x,len[edgeCount]=fracLen,nxt[edgeCount]=head[middle],head[middle]=edgeCount,belong[edgeCount]=middle;
    for(int i=0;i<5;i++)
        if(closest_id[ye*5+i]!=-1){
            closest_id[edgeCount*5+i]=closest_id[ye*5+i];
            closest_dis[edgeCount*5+i]=closest_dis[ye*5+i]+originalDis-fracLen;
        }
    edgeCount++;
    //middle -> y
    e[edgeCount]=y,len[edgeCount]=originalDis-fracLen,nxt[edgeCount]=head[middle],head[middle]=edgeCount,belong[edgeCount]=middle;
    for(int i=0;i<5;i++)
        if(closest_id[xe*5+i]!=-1){
            closest_id[edgeCount*5+i]=closest_id[xe*5+i];
            closest_dis[edgeCount*5+i]=closest_dis[xe*5+i]+fracLen;
        }
    edgeCount++;
    //outside -> middle
    e[edgeCount]=middle,len[edgeCount]=addLen,nxt[edgeCount]=head[outside],head[outside]=edgeCount,belong[edgeCount]=outside;
    edgeCount++;
    //middle -> outside
    e[edgeCount]=outside,len[edgeCount]=addLen,nxt[edgeCount]=head[middle],head[middle]=edgeCount,belong[edgeCount]=middle;
    int e1=edgeCount-2, e2=edgeCount-3;
    for(int i=0;i<5;i++){
        if(closest_id[e1*5+i]==-1) break;
        for(int j=0;j<5;j++)
            if(closest_dis[edgeCount*5+j]>closest_dis[e1*5+i]){
                for(int k=4;k>j;k--){
                    closest_dis[edgeCount*5+k]=closest_dis[edgeCount*5+k-1];
                    closest_id[edgeCount*5+k]=closest_id[edgeCount*5+k-1];
                }
                closest_dis[edgeCount*5+j]=closest_dis[e1*5+i];
                closest_id[edgeCount*5+j]=closest_id[e1*5+i];
                break;
            }
    }
    for(int i=0;i<5;i++){
        if(closest_id[e2*5+i]==-1) break;
        for(int j=0;j<5;j++)
            if(closest_dis[edgeCount*5+j]>closest_dis[e2*5+i]){
                for(int k=4;k>j;k--){
                    closest_dis[edgeCount*5+k]=closest_dis[edgeCount*5+k-1];
                    closest_id[edgeCount*5+k]=closest_id[edgeCount*5+k-1];
                }
                closest_dis[edgeCount*5+j]=closest_dis[e2*5+i];
                closest_id[edgeCount*5+j]=closest_id[e2*5+i];
                break;
            }
    }
    edgeCount++;
}

__global__ void buildInitialTree(
    int numSequences,
    int * head,
    int * e,
    double * len,
    int * nxt,
    int * belong,
    double * dis,
    int edgeCount
){
    int nv = numSequences;
    double d = dis[0];
    // 0 -> nv
    e[edgeCount]=nv,len[edgeCount]=d/2,nxt[edgeCount]=head[0],head[0]=edgeCount,belong[edgeCount]=0;
    edgeCount++;
    // 1 -> nv
    e[edgeCount]=nv,len[edgeCount]=d/2,nxt[edgeCount]=head[1],head[1]=edgeCount,belong[edgeCount]=1;
    edgeCount++;
    // nv -> 0
    e[edgeCount]=0,len[edgeCount]=d/2,nxt[edgeCount]=head[nv],head[nv]=edgeCount,belong[edgeCount]=nv;
    edgeCount++;
    // nv -> 0
    e[edgeCount]=1,len[edgeCount]=d/2,nxt[edgeCount]=head[nv],head[nv]=edgeCount,belong[edgeCount]=nv;
    edgeCount++;
}

void MashPlacement::KPlacementDeviceArrays::deallocateDeviceArrays(){
    cudaFree(d_head);
    cudaFree(d_e);
    cudaFree(d_nxt);
    cudaFree(d_belong);
    cudaFree(d_closest_id);
    cudaFree(d_dist);
    cudaFree(d_len);
    cudaFree(d_closest_dis);
}


void MashPlacement::KPlacementDeviceArrays::printTree(std::vector <std::string> name){
    int * h_head = new int[numSequences*2];
    int * h_e = new int[numSequences*8];
    int * h_nxt = new int[numSequences*8];
    double * h_len = new double[numSequences*8];
    double * h_closest_dis = new double[numSequences*20];
    int * h_closest_id = new int[numSequences*20];
    std::function<void(int,int)>  print=[&](int node, int from){
        if(h_nxt[h_head[node]]!=-1){
            printf("(");
            std::vector <int> pos;
            for(int i=h_head[node];i!=-1;i=h_nxt[i])
                if(h_e[i]!=from)
                    pos.push_back(i);
            for(size_t i=0;i<pos.size();i++){
                print(h_e[pos[i]],node);
                printf(":");
                printf("%.5g%c",h_len[pos[i]],i+1==pos.size()?')':',');
            }
        }
        else std::cout<<name[node];
    };
    auto err = cudaMemcpy(h_head, d_head, numSequences*2*sizeof(int),cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        exit(1);
    }
    err = cudaMemcpy(h_e, d_e, numSequences*8*sizeof(int),cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        exit(1);
    }
    err = cudaMemcpy(h_nxt, d_nxt, numSequences*8*sizeof(int),cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        exit(1);
    }
    err = cudaMemcpy(h_len, d_len, numSequences*8*sizeof(double),cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        exit(1);
    }
    print(numSequences+bd-2,-1);
    std::cout<<";\n";
}


void MashPlacement::KPlacementDeviceArrays::findPlacementTree(
    Param& params,
    const MashDeviceArrays& mashDeviceArrays,
    MatrixReader& matrixReader,
    const MSADeviceArrays& msaDeviceArrays
){ 
    if(params.in == "d"){
        matrixReader.distConstructionOnGpu(params, 0, d_dist);
    }
    int * d_id;
    auto err = cudaMalloc(&d_id, numSequences*2*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    int * d_from;
    err = cudaMalloc(&d_from, numSequences*2*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    double * d_dis;
    err = cudaMalloc(&d_dis, numSequences*2*sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    
    /*
    Initialize closest nodes by inifinite
    */
    int threadNum = 256, blockNum = (numSequences*4-4+threadNum-1)/threadNum;
    initialize <<<blockNum, threadNum>>> (
        numSequences*4-4,
        numSequences*2,
        d_closest_dis,
        d_closest_id,
        d_head,
        d_nxt,
        d_belong,
        d_e
    );
    /*
    Build Initial Tree
    */
    if(params.in == "r"){
        mashDeviceArrays.distConstructionOnGpu(
            params,
            1,
            d_dist
        );
    }
    else if(params.in == "d"){
        matrixReader.distConstructionOnGpu(
            params,
            1,
            d_dist
        );
    }
    else if(params.in == "m"){
        msaDeviceArrays.distConstructionOnGpu(
            params,
            1,
            d_dist
        );
    }
    // double * h_dis = new double[numSequences];
    // cudaMemcpy(h_dis,d_dist,numSequences*sizeof(double),cudaMemcpyDeviceToHost);
    // fprintf(stderr, "%d\n",1);
    // for(int j=0;j<1;j++) fprintf(stderr,"%.8lf ",h_dis[j]);std::cerr<<'\n';
    buildInitialTree <<<1,1>>> (
        numSequences,
        d_head,
        d_e,
        d_len,
        d_nxt,
        d_belong,
        d_dist,
        idx
    );
    idx += 4;
    /*
    Initialize closest nodes by inital tree
    */
    for(int i=0;i<bd;i++){
        updateClosestNodes <<<1,1>>> (
            d_head,
            d_nxt,
            d_e,
            d_len,
            d_closest_dis,
            d_closest_id,
            i,
            d_id,
            d_from,
            d_dis
        );
    }
 
    thrust::device_vector <thrust::tuple<int,double,double>> minPos(numSequences*4-4);
    // std::cout<<"FFF\n";
    std::chrono::nanoseconds disTime(0), treeTime(0);
    for(int i=bd;i<numSequences;i++){
        auto disStart = std::chrono::high_resolution_clock::now();
        blockNum = (i + 255) / 256;
        if(params.in == "r"){
            mashDeviceArrays.distConstructionOnGpu(
                params,
                i,
                d_dist
            );
        }
        else if(params.in == "d"){
            matrixReader.distConstructionOnGpu(
                params,
                i,
                d_dist
            );
        }
        else if(params.in == "m"){
            msaDeviceArrays.distConstructionOnGpu(
                params,
                i,
                d_dist
            );
        }
        cudaDeviceSynchronize();

        // double * h_dis = new double[numSequences];
        // cudaMemcpy(h_dis,d_dist,numSequences*sizeof(double),cudaMemcpyDeviceToHost);
        // fprintf(stderr, "%d\n",i);
        // for(int j=0;j<i;j++) std::cerr<<h_dis[j]<<" ";std::cerr<<'\n';


        auto disEnd = std::chrono::high_resolution_clock::now();
        auto treeStart = std::chrono::high_resolution_clock::now();
        blockNum = (numSequences*4-4 + 255) / 256;
        calculateBranchLength <<<blockNum,threadNum>>> (
            i,
            d_head,
            d_nxt,
            d_dist,
            d_e,
            d_len,
            d_belong,
            thrust::raw_pointer_cast(minPos.data()),
            numSequences*4-4,
            d_closest_dis,
            d_closest_id,
            numSequences
        );
        auto iter=thrust::min_element(minPos.begin(),minPos.end(),compare_tuple());
        thrust::tuple<int,double,double> smallest=*iter;
        /*
        Update Tree (and assign closest nodes to newly added nodes)
        */
        int eid=thrust::get<0>(smallest);
        double fracLen=thrust::get<1>(smallest),addLen=thrust::get<2>(smallest);
        updateTreeStructure <<<1,1>>>(
            d_head,
            d_nxt,
            d_e,
            d_len,
            d_closest_dis,
            d_closest_id,
            d_belong,
            eid,
            fracLen,
            addLen,
            i,
            idx,
            numSequences
        );
        idx+=4;
        /*
        Update closest nodes
        */
        updateClosestNodes <<<1,1>>> (
            d_head,
            d_nxt,
            d_e,
            d_len,
            d_closest_dis,
            d_closest_id,
            i,
            d_id,
            d_from,
            d_dis
        );
        cudaDeviceSynchronize();
        auto treeEnd = std::chrono::high_resolution_clock::now();
        disTime += disEnd - disStart;
        treeTime += treeEnd - treeStart;
    }
    std::cerr << "Distance Operation Time " <<  disTime.count()/1000000 << " ms\n";
    std::cerr << "Tree Operation Time " <<  treeTime.count()/1000000 << " ms\n";
}
