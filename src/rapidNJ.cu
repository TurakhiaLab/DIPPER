#ifndef RNJ_CUH
#include "rapidNJ.cuh"
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
#include <thrust/gather.h>



void rapidNJ::DeviceArrays::allocateDeviceArrays(uint32_t numSequences, double *mashDist){
    //Allocate memories on GPU 
    //mashDist should be memory on GPU, should contain upper triangle distance matrix, row by row
    //distances on diagonal (i=j) should not be contained in mashDist
    d_oriMashDist = mashDist;
    d_numSequences = numSequences;
    cudaError_t err;
    err = cudaMalloc(&d_U, numSequences*sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_timeStamp, numSequences*sizeof(uint32_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_mashDist, numSequences*numSequences*sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_id, numSequences*numSequences*sizeof(uint32_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_globalMin, sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    cudaDeviceSynchronize();
}

void rapidNJ::DeviceArrays::inputDismatrix(uint32_t numSequences, double *mashDist){
    //Allocate memories on GPU and copy data to GPU
    //mashDist should be memory on CPU, should contain upper triangle distance matrix, row by row
    //distances on diagonal (i=j) should not be contained in mashDist
    d_numSequences = numSequences;
    auto err = cudaMalloc(&d_oriMashDist, numSequences*(numSequences-1)/2*sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    cudaDeviceSynchronize();
    err = cudaMemcpy(d_oriMashDist, mashDist, numSequences*(numSequences-1)/2*sizeof(double),cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_U, numSequences*sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_timeStamp, numSequences*sizeof(uint32_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_mashDist, numSequences*numSequences*sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_id, numSequences*numSequences*sizeof(uint32_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_globalMin, sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&global_id, sizeof(uint32_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_rowId, numSequences*sizeof(uint32_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_preMaxU, numSequences*sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    cudaDeviceSynchronize();
}

void rapidNJ::DeviceArrays::deallocateDeviceArrays(){
    //Free space on GPU at the end
    cudaFree(d_mashDist);
    cudaFree(d_oriMashDist);
    cudaFree(d_U);
    cudaDeviceSynchronize();
}

__global__ void calculateU(uint32_t d_numSequences, double *d_mashDist, double *d_U){
    //Only executes once, calculate the sum of distances in each row
    uint32_t numSequences=d_numSequences;
    __shared__ double temp_U[128]; // Make sure to set this value larger than numSequences/gridSize
    for(auto i=0;i<128;i++)
        if(threadIdx.x==0)
            temp_U[i]=0;
    __syncthreads();
    for(auto i=blockIdx.x;i<numSequences;i+=gridDim.x){
        double temp=0;
        for(auto j=threadIdx.x;j<numSequences;j+=blockDim.x)
            if(i!=j) temp += d_mashDist[i*numSequences+j];
        atomicAdd(&temp_U[i/gridDim.x],temp);
        // For each block, calculate the sum of distances in some rows
        // For each thread, use atomicAdd to update the sum to shared memory
    }
    __syncthreads();
    if(threadIdx.x==0){
        for(auto i=blockIdx.x;i<numSequences;i+=gridDim.x) d_U[i]=temp_U[i/gridDim.x];
        // For each block, update the sum to global memory
    }
}

__global__ void buildDist(uint32_t d_numSequences, double *d_mashDist, double *d_oriMashDist){
    // Transform the upper triangle matrix to full matrix
    int tx=threadIdx.x, bx=blockIdx.x;
    int bs=blockDim.x, gs=gridDim.x;
    int st=bx, step=gs;
    uint32_t numSequences=d_numSequences,pos=(numSequences*2-1-st)*st/2;
    for(uint32_t i=st;i<numSequences;pos+=(numSequences*2-i*2-1-step)*step/2,i+=step){
        for(uint32_t j=tx;j<numSequences;j+=bs){
            if(i==j) d_mashDist[i*numSequences+j]=100; // Distances on diagnal is set to 0
            if(i>=j) continue; // Only do the following process when we are in upper triangle
            double val=d_oriMashDist[pos+j-i-1];
            d_mashDist[i*numSequences+j]=val;
            d_mashDist[j*numSequences+i]=val;
        }
    }
}

struct compare_tuple
{
  __host__ __device__
  bool operator()(thrust::tuple<uint32_t,uint32_t,double, double> lhs, thrust::tuple<uint32_t,uint32_t,double,double> rhs)
  {
    return thrust::get<2>(lhs) < thrust::get<2>(rhs);
    //Always find the tuple whose third value (the criteria we want to minimize) is minimized
  }
};



__global__ void buildID(uint32_t d_numSequences, uint32_t *d_id, uint32_t *d_timeStamp, uint32_t * d_rowId){
    int tx=threadIdx.x, bx=blockIdx.x;
    int bs=blockDim.x, gs=gridDim.x;
    for(int i=bx;i<d_numSequences;i+=gs){
        for(int j=tx;j<d_numSequences;j+=bs)
            d_id[i*d_numSequences+j]=j;
        if(tx==0) d_timeStamp[i]=i;
        if(tx==0) d_rowId[i]=i;
    }
}
__device__ double atomicMin(double* address, double val)
{
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        if(__longlong_as_double(assumed)<val) break;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(fmin(val, __longlong_as_double(assumed))));
    } while (assumed != old);

    return __longlong_as_double(old);
}
__global__ void findMinDist(uint32_t curNum, double *globalMin, double *d_mashDist, double *d_U, uint32_t *d_id, thrust::tuple<uint32_t,uint32_t,double,double> *minPos, uint32_t d_numSequences, uint32_t *d_timeStamp,uint32_t* global_id, uint32_t * d_rowId, double * d_preMaxU, uint32_t lim){
    int tx=threadIdx.x, bx=blockIdx.x;
    int bs=blockDim.x, gs=gridDim.x;
    double minD=100,oriD;
    uint32_t x=0,y=0;
    uint32_t timeStamp, id,i, p;
    double curU, dis, mU;
    while(1){
        i=atomicAdd(global_id,1);
        if(i>=lim) break;
        mU=d_preMaxU[i]/(curNum-2);
        // __syncthreads();
        p=i;
        i=d_rowId[i];
        d_rowId[p]=0;
        timeStamp=d_timeStamp[i];
        if(timeStamp>d_numSequences*2) continue;
        curU=d_U[i]/(curNum-2);
        for(int j=0;j<d_numSequences;j++){
            id=d_id[i*d_numSequences+j];
            dis=d_mashDist[i*d_numSequences+id];
            //atomicAdd(&d_rowId[p],1);
            if(dis-mU-curU>*globalMin) break;
            if(id==i) continue;
            // assert(d_U[i]<=mU);
            double val=dis-curU-d_U[id]/(curNum-2);
            if(val<minD&&d_timeStamp[id]<=timeStamp){
                minD=val,x=i,y=id,oriD=dis;
                if(val<*globalMin) atomicMin(globalMin,val);
            }
        }
    }
    thrust::tuple <uint32_t,uint32_t,double,double> minTuple(x,y,minD,oriD);
    minPos[bx*bs+tx]=minTuple;
}


__global__ void updateDisMatrix(uint32_t d_numSequences, double * d_mashDist, uint32_t *d_id, double * d_U, uint32_t * d_timeStamp,uint32_t x,uint32_t y,double disxy, uint32_t * d_rowId){
    int tx=threadIdx.x, bx=blockIdx.x;
    int bs=blockDim.x, gs=gridDim.x;
    __shared__ double indexSum[256];
    indexSum[tx]=0;
    for(int i=tx+bx*bs;i<d_numSequences;i+=bs*gs){
        d_rowId[i]=i;
        d_id[x*d_numSequences+i]=i;
        if(i==x||i==y||d_timeStamp[i]>d_numSequences*2){
            d_mashDist[x*d_numSequences+i]=100;
            continue;
        }
        double dx=d_timeStamp[x]>=d_timeStamp[i]?d_mashDist[x*d_numSequences+i]:d_mashDist[i*d_numSequences+x];
        double dy=d_timeStamp[y]>=d_timeStamp[i]?d_mashDist[y*d_numSequences+i]:d_mashDist[i*d_numSequences+y];
        double newdis=(dx+dy-disxy)*0.5;
        indexSum[tx]+=newdis;
        d_U[i]-=dx+dy-newdis;
        d_mashDist[x*d_numSequences+i]=newdis;
    }
    __syncthreads();
    for(int len=1;len<256;len<<=1){
        if(tx%(len*2)==0&&tx+len<256) indexSum[tx]+=indexSum[tx+len];
        __syncthreads();
    }
    if(tx==0) atomicAdd(&d_U[x],indexSum[0]);
}

struct my_policy : thrust::device_execution_policy<my_policy> {
        my_policy(int gs, int bs) : grid_size(gs), block_size(bs) {}
        int grid_size;
        int block_size;
};

void rapidNJ::findNeighbourJoiningTree(uint32_t d_numSequences, double *d_mashDist, double *d_U, double *d_oriMashDist, uint32_t *d_id, std::vector <std::string> &name, double * d_globalMin, uint32_t *d_timeStamp,uint32_t *global_id, uint32_t * d_rowId, double * d_preMaxU){
    //"name" should be a vector contain the name of each sample/sequence
    //see test_nj.cpp for more information
    auto superStart = std::chrono::high_resolution_clock::now();
    std::vector <std::vector<std::pair<uint32_t,double>>> tree(d_numSequences*2);
    // Store the tree in a vector of vectors, while each small vector consists of several pairs
    // First element in the pair is the index of its child
    // Second element in the pair is the distance to its child
    int cntThread=256,cntBlock=256; // Could be adjusted based on experiment
    int realID[d_numSequences]; // Actual sequence id corresponding to i-th row/column in the matrix
    for(int i=0;i<d_numSequences;i++) realID[i]=i;
    buildDist <<<1024,64>>> (d_numSequences, d_mashDist, d_oriMashDist);
    cudaDeviceSynchronize();
    // Build the full distance matrix
    thrust::device_vector <thrust::tuple<uint32_t,uint32_t,double,double>> minPos(cntThread*cntBlock);
    // Device_vector that stores minimum value found by each thread
    calculateU <<<1024,64>>> (d_numSequences, d_mashDist, d_U);
    cudaDeviceSynchronize();
    buildID <<<1024,64>>> (d_numSequences, d_id, d_timeStamp, d_rowId);
    cudaDeviceSynchronize();
    my_policy policy_reduce((d_numSequences+31)/32, 32);
    my_policy policy_sort((d_numSequences+31)/32,32);
    thrust::device_vector<double> tempDis(d_numSequences);
    thrust::device_vector<uint32_t> tempId(d_numSequences);
    for(int i=0;i<d_numSequences;i++){
        thrust::sort_by_key(thrust::device, d_mashDist+i*d_numSequences,d_mashDist+(i+1)*d_numSequences,d_id+i*d_numSequences);
        thrust::scatter(thrust::device, d_mashDist+i*d_numSequences,d_mashDist+(i+1)*d_numSequences,d_id+i*d_numSequences,tempDis.begin());
        thrust::copy(thrust::device, tempDis.begin(), tempDis.begin()+d_numSequences, thrust::device_pointer_cast(d_mashDist+i*d_numSequences));
    }
    auto preProcess = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds preProcessTime = preProcess - superStart;
    std::cout<<"Preprocess in "<<preProcessTime.count()/1000000<<" ms.\n";
    // double *h_dis=new double[d_numSequences*d_numSequences];
    // cudaMemcpy(h_dis, d_mashDist, d_numSequences*d_numSequences*sizeof(double),cudaMemcpyDeviceToHost);
    // for(int i=0;i<d_numSequences;i++)
    //     for(int j=0;j<d_numSequences;j++)
    //         printf("%.5g%c",h_dis[i*d_numSequences+j],j+1==d_numSequences?'\n':' ');
    // uint32_t *h_id=new uint32_t[d_numSequences*d_numSequences];
    // cudaMemcpy(h_id, d_id, d_numSequences*d_numSequences*sizeof(uint32_t),cudaMemcpyDeviceToHost);
    // for(int i=0;i<d_numSequences;i++)
    //     for(int j=0;j<d_numSequences;j++)
    //         printf("%u%c",h_id[i*d_numSequences+j],j+1==d_numSequences?'\n':' ');
    long long t1=0,t2=0,t3=0;
    int ID=d_numSequences;
    for(int i=0;i<d_numSequences-2;i++){
        double *tp=new double;
        uint32_t *temp=new uint32_t;
        *tp=100.0;
        *temp=0;
        cudaMemcpy(global_id,temp,sizeof(uint32_t),cudaMemcpyHostToDevice);
        cudaMemcpy(d_globalMin,tp,sizeof(double),cudaMemcpyHostToDevice);
        thrust::sort_by_key(thrust::device, d_timeStamp,d_timeStamp+d_numSequences,d_rowId, thrust::less<uint32_t>());
        thrust::scatter(thrust::device,d_timeStamp,d_timeStamp+d_numSequences,d_rowId,tempId.begin());
        thrust::copy(thrust::device, tempId.begin(), tempId.begin()+d_numSequences, thrust::device_pointer_cast(d_timeStamp));
        thrust::gather(thrust::device, d_rowId,d_rowId+d_numSequences,d_U,d_preMaxU);
        //thrust::reverse(thrust::device, d_preMaxU,d_preMaxU+d_numSequences);
        thrust::exclusive_scan(thrust::device, d_preMaxU, d_preMaxU+d_numSequences, d_preMaxU, 0.0, thrust::maximum<double>());
        //thrust::reverse(thrust::device, d_preMaxU,d_preMaxU+d_numSequences);
        //double mU=thrust::reduce(thrust::device,d_U,d_U+d_numSequences,0.0,thrust::maximum<double>());
        findMinDist <<<cntBlock,cntThread>>> (d_numSequences-i,d_globalMin,d_mashDist,d_U,d_id,thrust::raw_pointer_cast(minPos.data()),d_numSequences,d_timeStamp,global_id, d_rowId, d_preMaxU, d_numSequences-i);
        
        // if(i==0){
        //     uint32_t *h_dis=new uint32_t[d_numSequences];
        //     cudaMemcpy(h_dis, d_rowId, d_numSequences*sizeof(uint32_t),cudaMemcpyDeviceToHost);
        //     for(int j=0;j<d_numSequences;j++)
        //         printf("%d: %d times\n",j,h_dis[j]);
        // }

        //cudaDeviceSynchronize();
        auto njEnd = std::chrono::high_resolution_clock::now();
        auto iter=thrust::min_element(thrust::device, minPos.begin(),minPos.end(),compare_tuple());
        // Find the global minimum
        thrust::tuple<uint32_t,uint32_t,double,double> smallest=*iter;
        uint32_t x=thrust::get<0>(smallest),y=thrust::get<1>(smallest);
        //if(x<y) std::swap(x,y);
        double modified_dis=thrust::get<2>(smallest),dis=thrust::get<3>(smallest);
        double *ux=new double;
        double *uy=new double;
        //auto beforeCopy = std::chrono::high_resolution_clock::now();
        cudaMemcpy(ux,d_U+x,sizeof(double),cudaMemcpyDeviceToHost);
        cudaMemcpy(uy,d_U+y,sizeof(double),cudaMemcpyDeviceToHost);
        //auto beforeUpdate = std::chrono::high_resolution_clock::now();
        //auto mid=beforeUpdate-beforeCopy;
        //std::cout<<"Three cudaMemcpy takes "<<mid.count()<<" ns.\n";
        double blX = (dis+(*ux)/(d_numSequences-i-2)-(*uy)/(d_numSequences-i-2))*0.5;
        double blY = dis - blX;
        if(blX<0) blY+=blX, blX=0;
        if(blY<0) blX+=blY, blY=0;
        // Avoid negative branches
        tree[ID].push_back({realID[x],blX});
        tree[ID].push_back({realID[y],blY});
        // Insert the two edges into the tree
        realID[x]=ID++;
        *tp=0;
        cudaMemcpy(d_U+x,tp,sizeof(double),cudaMemcpyHostToDevice);
        *tp=0;
        cudaMemcpy(d_U+y,tp,sizeof(double),cudaMemcpyHostToDevice);
        updateDisMatrix <<<32,256>>>(d_numSequences,d_mashDist,d_id, d_U,d_timeStamp,x,y,dis,d_rowId);
        thrust::sort_by_key(thrust::device, d_mashDist+x*d_numSequences,d_mashDist+(x+1)*d_numSequences,d_id+x*d_numSequences);
        thrust::scatter(thrust::device, d_mashDist+x*d_numSequences,d_mashDist+(x+1)*d_numSequences,d_id+x*d_numSequences,tempDis.begin());
        //assert(tempDis.size()==d_numSequences);
        thrust::copy(thrust::device, tempDis.begin(), tempDis.begin()+d_numSequences, thrust::device_pointer_cast(d_mashDist+x*d_numSequences));
        *temp=1<<25;
        cudaMemcpy(d_timeStamp+y,temp,sizeof(uint32_t),cudaMemcpyHostToDevice);
        *temp=d_numSequences+i;
        cudaMemcpy(d_timeStamp+x,temp,sizeof(uint32_t),cudaMemcpyHostToDevice);
        printf("%d %d %d %.8lf %.8lf\n",i, realID[x],realID[y],modified_dis,dis);
        // double h_min;
        // cudaMemcpy(&h_min,d_globalMin,sizeof(double),cudaMemcpyDeviceToHost);
        // std::chrono::nanoseconds njTime = njEnd - njStart;
        // cudaDeviceSynchronize();
        // auto iterEnd = std::chrono::high_resolution_clock::now();
        // std::chrono::nanoseconds otherIterTime = iterEnd - beforeUpdate;
        // std::chrono::nanoseconds middleIterTime = beforeUpdate - njEnd;
        // //printf("%.8lf\n",h_min);
        // t1+=njTime.count(),t2+=middleIterTime.count(),t3+=otherIterTime.count();
    }
    cudaDeviceSynchronize();
    // auto allProcess = std::chrono::high_resolution_clock::now();
    // std::chrono::nanoseconds forLoopTime = allProcess - preProcess;
    // std::cout<<"For loop in "<<forLoopTime.count()/1000000<<" ms.\n";
    // std::cout<<"Find Minimum in "<<t1/1000000<<" ms.\n";
    // std::cout<<"Other parts before update in "<<t2/1000000<<" ms.\n";
    // std::cout<<"Other parts of iteration in "<<t3/1000000<<" ms.\n";

    // // If only two nodes are left, their distance is fixed
    // uint32_t x=0,y=1;
    // double dis=0;
    // cudaMemcpy(&dis,d_mashDist+x*d_numSequences+y,sizeof(double),cudaMemcpyDeviceToHost);
    // tree[d_numSequences*2-2].push_back({realID[x],dis*0.5});
    // tree[d_numSequences*2-2].push_back({realID[y],dis*0.5});
    // cudaDeviceSynchronize();
    // // Output the tree recursively
    // std::function<void(int)>  print=[&](int node){
    //     if(tree[node].size()){
    //         printf("(");
    //         for(size_t i=0;i<tree[node].size();i++){
    //             print(tree[node][i].first);
    //             printf(":");
    //             printf("%.5g%c",tree[node][i].second,i+1==tree[node].size()?')':',');
    //         }
    //     }
    //     else std::cout<<"\'"<<name[node]<<"\':";
    // };
    // // Root of the tree has an index of d_numSequneces*2-2
    // print(d_numSequences*2-2);
    // std::cout<<";\n\n";
}