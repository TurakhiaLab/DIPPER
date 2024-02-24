#ifndef NJ_CUH
#include "neighbourJoining.cuh"
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


void neighbourJoining::DeviceArrays::allocateDeviceArrays(uint32_t numSequences, double *mashDist){
    d_oriMashDist = mashDist;
    d_numSequences = numSequences;
    cudaError_t err;
    err = cudaMalloc(&d_flag, numSequences*sizeof(uint32_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!111\n");
        exit(1);
    }
    err = cudaMalloc(&d_U, numSequences*sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!333\n");
        exit(1);
    }
    err = cudaMalloc(&d_mashDist, numSequences*numSequences*sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!333\n");
        exit(1);
    }
    err = cudaMalloc(&d_rowID, numSequences*numSequences*sizeof(uint32_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!333\n");
        exit(1);
    }
    //d_tree.resize(d_numSequences*2);
    cudaDeviceSynchronize();
    err = cudaMemset(d_flag, 0, numSequences*sizeof(uint32_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!222\n");
        exit(1);
    }
    cudaDeviceSynchronize();
}

void neighbourJoining::DeviceArrays::inputDismatrix(uint32_t numSequences, double *mashDist){
    auto err = cudaMalloc(&d_oriMashDist, numSequences*(numSequences-1)/2*sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!111\n");
        exit(1);
    }
    err = cudaMemcpy(d_oriMashDist, mashDist, numSequences*(numSequences-1)/2*sizeof(double),cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        exit(1);
    }
    d_numSequences = numSequences;
    err = cudaMalloc(&d_flag, numSequences*sizeof(uint32_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!111\n");
        exit(1);
    }
    err = cudaMalloc(&d_U, numSequences*sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!333\n");
        exit(1);
    }
    err = cudaMalloc(&d_mashDist, numSequences*numSequences*sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!333\n");
        exit(1);
    }
    err = cudaMalloc(&d_rowID, numSequences*numSequences*sizeof(uint32_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!333\n");
        exit(1);
    }
    //d_tree.resize(d_numSequences*2);
    err = cudaMemset(d_flag, 0, numSequences*sizeof(uint32_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!222\n");
        exit(1);
    }
    cudaDeviceSynchronize();
}

void neighbourJoining::DeviceArrays::deallocateDeviceArrays(){
    cudaFree(d_mashDist);
    cudaFree(d_flag);
    cudaFree(d_U);
    cudaDeviceSynchronize();
}

__global__ void updateUMatrix(uint32_t d_numSequences, double *d_mashDist, double *d_U, uint32_t *d_flag){
    uint32_t numSequences=d_numSequences;
    __shared__ double temp_U[128];
    for(auto i=0;i<128;i++)
        if(threadIdx.x==0)
            temp_U[i]=0;
    for(auto i=blockIdx.x;i<numSequences;i+=gridDim.x){
        if(d_flag[i]) continue;
        double temp=0;
        for(auto j=threadIdx.x;j<numSequences;j+=blockDim.x)
            if(d_flag[j]==0&&i!=j) temp += d_mashDist[i*numSequences+j];
        atomicAdd(&temp_U[i/gridDim.x],temp);
    }
    __syncthreads();
    if(threadIdx.x==0){
        for(auto i=blockIdx.x;i<numSequences;i+=gridDim.x) d_U[i]=temp_U[i/gridDim.x];
    }
}

__global__ void buildDist(uint32_t d_numSequences, double *d_mashDist, double *d_oriMashDist, uint32_t * d_rowID){
    int tx=threadIdx.x, bx=blockIdx.x;
    int bs=blockDim.x, gs=gridDim.x;
    int st=bx, step=gs;
    uint32_t numSequences=d_numSequences,pos=(numSequences*2-1-st)*st/2;
    for(auto i=st;i<numSequences;pos+=(numSequences*2-i*2-1-step)*step/2,i+=step){
        for(auto j=tx;j<numSequences;j+=bs){
            d_rowID[i*numSequences+j]=i;
            if(i==j) d_mashDist[i*numSequences+j]=0;
            if(i>=j) continue;
            double val=d_oriMashDist[pos+j-i-1];
            d_mashDist[i*numSequences+j]=val;
            d_mashDist[j*numSequences+i]=val;
        }
    }
}

__global__ void findMinDist(uint32_t d_numSequences, double *d_mashDist, double *d_U, thrust::tuple<uint32_t,uint32_t,double> *minPos, uint32_t len){
    __shared__ double U[128];//length=d_numSequences/gridDim.x+1;
    int tx=threadIdx.x, bx=blockIdx.x;
    int bs=blockDim.x, gs=gridDim.x;
    uint32_t numSequences=d_numSequences;
    uint32_t sz=numSequences/gs,st=sz*bx;
    if(numSequences%gs>bx) sz++;
    st+=min(bx,numSequences%gs);
    uint32_t ed=st+sz;
    for(auto i=st+tx;i<ed;i+=bs)
        U[i-st]=d_U[i]/(numSequences-2);
    __syncthreads();
    double colU;
    double minD=10000;
    uint32_t x=0,y=0;
    for(auto j=tx;j<numSequences;j+=bs){
        colU=d_U[j]/(numSequences-2);
        for(auto i=st;i<ed;i++){
            double temp=d_mashDist[i*len+j]-U[i-st]-colU;
            if(i!=j&&temp<minD)
                minD=temp,x=i,y=j;
        }
    }

    thrust::tuple <uint32_t,uint32_t,double> minTuple(x,y,minD);
    minPos[bx*bs+tx]=minTuple;
}

struct compare_tuple
{
  __host__ __device__
  bool operator()(thrust::tuple<uint32_t,uint32_t,double> lhs, thrust::tuple<uint32_t,uint32_t,double> rhs)
  {
    return thrust::get<2>(lhs) < thrust::get<2>(rhs);
  }
};


__global__ void updateDisMatrix(uint32_t d_numSequences, double *d_mashDist, double disxy, uint32_t x,uint32_t y, double* d_U, uint32_t total){
    int tx=threadIdx.x, bx=blockIdx.x;
    int bs=blockDim.x, gs=gridDim.x;
    int start=bx*bs+tx,step=bs*gs,numSequences=total;
    if(tx==0&&bx==0) d_U[x]=0;
    for(int i=start;i<d_numSequences-1;i+=step)
        if(i!=x&&i!=y){
            double val=(d_mashDist[x*numSequences+i]+d_mashDist[y*numSequences+i]-disxy)*0.5;
            double fardis=d_mashDist[(d_numSequences-1)*numSequences+i];
            d_U[i]+=-d_mashDist[x*numSequences+i]-d_mashDist[y*numSequences+i]+val;
            atomicAdd(&d_U[x],val);
            d_mashDist[x*numSequences+i]=val;
            d_mashDist[i*numSequences+x]=val;
            d_mashDist[y*numSequences+i]=fardis;
            d_mashDist[i*numSequences+y]=fardis;
        }
    if(tx==0&&bx==0){
        d_U[y]=d_U[d_numSequences-1];
        uint32_t i=d_numSequences-1;
        double val=(d_mashDist[x*numSequences+i]+d_mashDist[y*numSequences+i]-disxy)*0.5;
        d_U[y]+=-d_mashDist[x*numSequences+i]-d_mashDist[y*numSequences+i]+val;
        atomicAdd(&d_U[x],val);
        d_mashDist[x*numSequences+y]=val;
        d_mashDist[y*numSequences+x]=val;
    }
}


void neighbourJoining::findNeighbourJoiningTree(uint32_t d_numSequences, double *d_mashDist, double *d_U, uint32_t *d_flag, double *d_oriMashDist, uint32_t * d_rowID){
    std::vector <std::vector<std::pair<uint32_t,double>>> tree(d_numSequences*2);

   
    int cntThread=256,cntBlock=256;
    int realID[d_numSequences];
    for(int i=0;i<d_numSequences;i++) realID[i]=i;
    buildDist <<<cntBlock,cntThread>>> (d_numSequences, d_mashDist, d_oriMashDist, d_rowID);
    // cudaDeviceSynchronize();
    // double * h_mashDist = new double[d_numSequences*d_numSequences];
    // auto err = cudaMemcpy(h_mashDist, d_mashDist, (d_numSequences*d_numSequences)*sizeof(double), cudaMemcpyDeviceToHost);
    // if (err != cudaSuccess) {
    //     fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
    //     exit(1);
    // }
    // for(int i=0;i<d_numSequences*d_numSequences;i++) printf("%f ",h_mashDist[i]);
    // printf("\n");
    thrust::device_vector <thrust::tuple<uint32_t,uint32_t,double>> minPos(cntThread*cntBlock);
    updateUMatrix <<<cntBlock,cntThread>>> (d_numSequences, d_mashDist, d_U, d_flag);
    //thrust::tuple<uint32_t,uint32_t,double> minPos[cntThread*cntBlock];
    int ID=d_numSequences;
    for(int i=0;i<d_numSequences-2;i++){
        //thrust::reduce_by_key(thrust::device, d_rowID, d_rowID+d_numSequences*d_numSequences,d_mashDist, tempMemory, d_U);
        // cudaDeviceSynchronize();
        // double * h_U = new double[d_numSequences];
        // auto err = cudaMemcpy(h_U, d_U, (d_numSequences)*sizeof(double), cudaMemcpyDeviceToHost);
        // if (err != cudaSuccess) {
        //     fprintf(stderr, "Gpu_ERROR: cudaMemCpy failed!\n");
        //     exit(1);
        // }
        // for(int j=0;j<d_numSequences;j++) printf("%f ",h_U[j]);
        // printf("\n");
        
        findMinDist <<<cntBlock,cntThread>>> (d_numSequences-i,d_mashDist,d_U,thrust::raw_pointer_cast(minPos.data()),d_numSequences);
        // cudaDeviceSynchronize();
        // printf("%d----------\n",i);
        // auto iter=minPos.begin();
        // while(iter!=minPos.end()){
        //     thrust::tuple<uint32_t,uint32_t,double> smallest=*iter;
        //     uint32_t x=thrust::get<0>(smallest),y=thrust::get<1>(smallest);
        //     double d=thrust::get<2>(smallest);
        //     printf("%u %u %f\n", x,y,d);
        //     iter++;
        // }

        auto iter=thrust::min_element(minPos.begin(),minPos.end(),compare_tuple());
        thrust::tuple<uint32_t,uint32_t,double> smallest=*iter;
        uint32_t x=thrust::get<0>(smallest),y=thrust::get<1>(smallest);
        double modified_dis=thrust::get<2>(smallest);
        
        double* dis=new double;
        double *ux=new double;
        double *uy=new double;
        cudaMemcpy(dis,d_mashDist+x*d_numSequences+y,sizeof(double),cudaMemcpyDeviceToHost);
        cudaMemcpy(ux,d_U+x,sizeof(double),cudaMemcpyDeviceToHost);
        cudaMemcpy(uy,d_U+y,sizeof(double),cudaMemcpyDeviceToHost);
        double blX = (*dis+(*ux)/(d_numSequences-i-2)-(*uy)/(d_numSequences-i-2))*0.5;
        double blY = *dis - blX;
        if(blX<0) blY+=blX, blX=0;
        if(blY<0) blX+=blY, blY=0;
        //std::cout << blX <<" "<<blY<<'\n';
        tree[ID].push_back({realID[x],blX});
        tree[ID].push_back({realID[y],blY});
        //std::cout<<ID<<" "<<realID[x]<<" "<<realID[y]<<"\n";
        realID[x]=ID++, realID[y]=realID[d_numSequences-i-1];
        
        updateDisMatrix <<<cntBlock,cntThread>>> (d_numSequences-i,d_mashDist,*dis,x,y,d_U, d_numSequences);
        cudaDeviceSynchronize();
        
        //printf("%d %u %u %f\n", i, x,y,modified_dis);

    }
    uint32_t x=0,y=1;
    double dis=0;
    cudaMemcpy(&dis,d_mashDist+x*d_numSequences+y,sizeof(double),cudaMemcpyDeviceToHost);
    tree[d_numSequences*2-2].push_back({realID[x],dis*0.5});
    tree[d_numSequences*2-2].push_back({realID[y],dis*0.5});
    cudaDeviceSynchronize();
    std::function<void(int)>  print=[&](int node){
        if(tree[node].size()){
            printf("(");
            for(size_t i=0;i<tree[node].size();i++){
                print(tree[node][i].first);
                printf(":");
                printf("%.5f%c",tree[node][i].second,i+1==tree[node].size()?')':',');
            }
        }
        else printf("Node_%d",node);
    };
    print(d_numSequences*2-2);
    std::cout<<"\n\n";
}


