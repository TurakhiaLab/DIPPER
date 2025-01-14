#ifndef PL_CUH
#include "placement.cuh"
#endif

#include <stdio.h>
#include <queue>
#include <vector>
#include <thrust/sort.h>
#include <thrust/scan.h>
#include <thrust/binary_search.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/tuple.h>
#include <cuda_runtime.h>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>
#include <cublas_v2.h>


void placement::DeviceArrays::allocateDeviceArrays(
    int k,
    int num, 
    double *dist, 
    int *head, 
    int *e,
    int *nxt,
    int cnt,
    int cur,
    double * len,
    int * belong
){
    bd = cur, idx = cnt, numSequences = num, tk=k;
    cudaError_t err;
    err = cudaMalloc(&d_dist, 1LL*numSequences*numSequences*sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMemcpy(d_dist, dist, 1LL*numSequences*numSequences*sizeof(double),cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_head, numSequences*2*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMemcpy(d_head, head, numSequences*2*sizeof(int),cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_e, numSequences*8*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMemcpy(d_e, e, numSequences*8*sizeof(int),cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_len, numSequences*8*sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMemcpy(d_len, len, numSequences*8*sizeof(double),cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_nxt, numSequences*8*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMemcpy(d_nxt, nxt, numSequences*8*sizeof(int),cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_belong, numSequences*8*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMemcpy(d_belong, belong, numSequences*8*sizeof(int),cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_closest_dis, numSequences*4*tk*sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_closest_id, numSequences*4*tk*sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
}

void placement::DeviceArrays::deallocateDeviceArrays(){

}


__global__ void initialize(
    int lim,
    double * d_closest_dis,
    int * d_closest_id,
    int tk
){
    int tx=threadIdx.x,bs=blockDim.x;
    int bx=blockIdx.x,gs=gridDim.x;
    int idx=tx+bs*bx;
    if(idx<lim){
        for(int i=0;i<tk;i++){
            d_closest_dis[idx*tk+i]=2;
            d_closest_id[idx*tk+i]=-1;
        }
    }
}

struct compare_tuple
{
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
    int totSeqNum,
    int tk
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
    for(int i=0;i<tk;i++)
        if(closest_id[eid*tk+i]!=-1){
            val = dis[1ll*num*totSeqNum+closest_id[eid*tk+i]]-closest_dis[eid*tk+i];
            if(val>dis1) dis1=val;
        }
    otheid=head[oth];
    while(e[otheid]!=x) assert(otheid!=-1),otheid=nxt[otheid];
    for(int i=0;i<tk;i++)
        if(closest_id[otheid*tk+i]!=-1){
            val = dis[1ll*num*totSeqNum+closest_id[otheid*tk+i]]-closest_dis[otheid*tk+i];
            if(val>dis2) dis2=val;
        }
    double additional_dis=(dis1+dis2-len[eid])/2;
    if(additional_dis<0) additional_dis=0;
    dis1-=additional_dis,dis2-=additional_dis;
    if(dis1<0) dis1=0;
    if(dis2<0) dis2=0;
    if(dis1>len[eid]) additional_dis+=dis1-len[eid],dis1=len[eid];
    if(dis2>len[eid]) additional_dis+=dis2-len[eid],dis2=len[eid];
    assert(dis1+dis2-1e-6<=len[eid]);
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
    double * dis,
    int tk
){
    int l=0,r=-1;
    id[++r]=x,dis[x]=0,from[x]=-1;
    while(l<=r){
        int node=id[l],fb=from[l];
        double d=dis[l];
        l++;
        for(int i=head[node];i!=-1;i=nxt[i]){
            if(e[i]==fb) continue;
            for(int j=0;j<tk;j++){
                double nowd=closest_dis[i*tk+j];
                if(nowd>d){
                    for(int k=4;k>j;k--){
                        closest_dis[i*tk+k]=closest_dis[i*tk+k-1];
                        closest_id[i*tk+k]=closest_id[i*tk+k-1];
                    }
                    closest_dis[i*tk+j]=d;
                    closest_id[i*tk+j]=x;
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
    int numSequences,
    int tk
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
    for(int i=0;i<tk;i++)
        if(closest_id[ye*tk+i]!=-1){
            closest_id[edgeCount*tk+i]=closest_id[ye*tk+i];
            closest_dis[edgeCount*tk+i]=closest_dis[ye*tk+i]+originalDis-fracLen;
        }
    edgeCount++;
    //middle -> y
    e[edgeCount]=y,len[edgeCount]=originalDis-fracLen,nxt[edgeCount]=head[middle],head[middle]=edgeCount,belong[edgeCount]=middle;
    for(int i=0;i<tk;i++)
        if(closest_id[xe*tk+i]!=-1){
            closest_id[edgeCount*tk+i]=closest_id[xe*tk+i];
            closest_dis[edgeCount*tk+i]=closest_dis[xe*tk+i]+fracLen;
        }
    edgeCount++;
    //outside -> middle
    e[edgeCount]=middle,len[edgeCount]=addLen,nxt[edgeCount]=head[outside],head[outside]=edgeCount,belong[edgeCount]=outside;
    edgeCount++;
    //middle -> outside
    e[edgeCount]=outside,len[edgeCount]=addLen,nxt[edgeCount]=head[middle],head[middle]=edgeCount,belong[edgeCount]=middle;
    int e1=edgeCount-2, e2=edgeCount-3;
    for(int i=0;i<tk;i++){
        if(closest_id[e1*tk+i]==-1) break;
        for(int j=0;j<tk;j++)
            if(closest_dis[edgeCount*tk+j]>closest_dis[e1*tk+i]){
                for(int k=4;k>j;k--){
                    closest_dis[edgeCount*tk+k]=closest_dis[edgeCount*tk+k-1];
                    closest_id[edgeCount*tk+k]=closest_id[edgeCount*tk+k-1];
                }
                closest_dis[edgeCount*tk+j]=closest_dis[e1*tk+i];
                closest_id[edgeCount*tk+j]=closest_id[e1*tk+i];
                break;
            }
    }
    for(int i=0;i<tk;i++){
        if(closest_id[e2*tk+i]==-1) break;
        for(int j=0;j<tk;j++)
            if(closest_dis[edgeCount*tk+j]>closest_dis[e2*tk+i]){
                for(int k=4;k>j;k--){
                    closest_dis[edgeCount*tk+k]=closest_dis[edgeCount*tk+k-1];
                    closest_id[edgeCount*tk+k]=closest_id[edgeCount*tk+k-1];
                }
                closest_dis[edgeCount*tk+j]=closest_dis[e2*tk+i];
                closest_id[edgeCount*tk+j]=closest_id[e2*tk+i];
                break;
            }
    }
    edgeCount++;
}


void placement::findPlacementTree(
    int numSequences,
    int bd, // id to place
    int idx, // id of linked-list position
    double * d_dist,
    int * d_head,
    int * d_e,
    double * d_len,
    int * d_nxt,
    int * d_belong,
    double * d_closest_dis,
    int * d_closest_id,
    std::vector<std::string> name,
    int tk
){
    int * h_head = new int[numSequences*2];
    int * h_e = new int[numSequences*8];
    int * h_nxt = new int[numSequences*8];
    double * h_len = new double[numSequences*8];
    double * h_closest_dis = new double[numSequences*20];
    int * h_closest_id = new int[numSequences*20];
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
    int threadNum = 256, blockNum = (numSequences*4-4+threadNum-1)/threadNum;
    initialize <<<blockNum, threadNum>>> (
        numSequences*4-4,
        d_closest_dis,
        d_closest_id,
        tk
    );
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
            d_dis,
            tk
        );
    }
    thrust::device_vector <thrust::tuple<int,double,double>> minPos(numSequences*4-4);
    for(int i=bd;i<numSequences;i++){
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
            numSequences,
            tk
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
            numSequences,
            tk
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
            d_dis,
            tk
        );

        /*
        ----------------------Testing-------------------------
        */
        // std::cout<<i<<" "<<eid<<" "<<fracLen<<" "<<addLen<<'\n';
        // cudaError_t err;
        // err = cudaMemcpy(h_head, d_head, numSequences*2*sizeof(int),cudaMemcpyDeviceToHost);
        // if (err != cudaSuccess)
        // {
        //     fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        //     exit(1);
        // }
        // err = cudaMemcpy(h_e, d_e, numSequences*8*sizeof(int),cudaMemcpyDeviceToHost);
        // if (err != cudaSuccess)
        // {
        //     fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        //     exit(1);
        // }
        // err = cudaMemcpy(h_nxt, d_nxt, numSequences*8*sizeof(int),cudaMemcpyDeviceToHost);
        // if (err != cudaSuccess)
        // {
        //     fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        //     exit(1);
        // }
        // err = cudaMemcpy(h_len, d_len, numSequences*8*sizeof(double),cudaMemcpyDeviceToHost);
        // if (err != cudaSuccess)
        // {
        //     fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        //     exit(1);
        // }
        // err = cudaMemcpy(h_closest_dis, d_closest_dis, numSequences*20*sizeof(double),cudaMemcpyDeviceToHost);
        // if (err != cudaSuccess)
        // {
        //     fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        //     exit(1);
        // }
        // err = cudaMemcpy(h_closest_id, d_closest_id, numSequences*20*sizeof(int),cudaMemcpyDeviceToHost);
        // if (err != cudaSuccess)
        // {
        //     fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        //     exit(1);
        // }
        // print(numSequences+bd-2,-1);
        // std::cout<<";\n";
        // for(int j=0;j<i*4;j++){
        //     std::cout<<"Edge "<<j<<" to "<<h_e[j]<<":\n";
        //     for(int k=0;k<5;k++) printf("%.5g %d\n",h_closest_dis[j*5+k],h_closest_id[j*5+k]);
        // }
        /*
        ----------------Testing Ends-----------------------------
        */
    }
    /*
    Copy the tree back to CPU and output
    */
    // cudaError_t err;
    err = cudaMemcpy(h_head, d_head, numSequences*2*sizeof(int),cudaMemcpyDeviceToHost);
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