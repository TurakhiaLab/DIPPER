#include "mash_placement.cuh"

#include <stdio.h>
#include <queue>
#include <unordered_map>
#include <thrust/sort.h>
#include <thrust/scan.h>
#include <thrust/binary_search.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <chrono>
#include <iostream>
#include <fstream>
#include <cub/cub.cuh>

void MashPlacement::KPlacementDeviceArrays::allocateDeviceArrays(size_t num, int backboneSize)
{
    cudaError_t err;
    numSequences = int(num);
    bd = 2, idx = 0;
    this->backboneSize = backboneSize;
    err = cudaMalloc(&d_dist, numSequences * sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_head, numSequences * 2 * sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_e, numSequences * 8 * sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_len, numSequences * 8 * sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_nxt, numSequences * 8 * sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_belong, numSequences * 8 * sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_closest_dis, numSequences * 20 * sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    err = cudaMalloc(&d_closest_id, numSequences * 20 * sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
}

__global__ void initializeID(
    int lim,
    double *d_closest_dis,
    int *d_closest_id)
{
    int tx = threadIdx.x, bs = blockDim.x;
    int bx = blockIdx.x, gs = gridDim.x;
    int idx = tx + bs * bx;
    if (idx < lim)
    {
        for (int i = 0; i < 5; i++)
        {
            d_closest_dis[idx * 5 + i] = 2;
            d_closest_id[idx * 5 + i] = -1;
        }
    }
}

__global__ void updateClosestNodes(
    int *head,
    int *nxt,
    int *e,
    double *len,
    double *closest_dis,
    int *closest_id,
    int x,
    int *id,
    int *from,
    double *dis)
{
    int l = 0, r = -1;
    id[++r] = x, dis[x] = 0, from[x] = -1;
    while (l <= r)
    {
        int node = id[l], fb = from[l];
        double d = dis[l];
        l++;
        for (int i = head[node]; i != -1; i = nxt[i])
        {
            // printf("%d %d: \n", node, head[node]);
            if (e[i] == fb)
                continue;
            for (int j = 0; j < 5; j++)
            {
                double nowd = closest_dis[i * 5 + j];
                if (nowd > d)
                {
                    for (int k = 4; k > j; k--)
                    {
                        closest_dis[i * 5 + k] = closest_dis[i * 5 + k - 1];
                        closest_id[i * 5 + k] = closest_id[i * 5 + k - 1];
                    }
                    // printf("%d: (%d %lf)\t", i*5+j, x, d);
                    closest_dis[i * 5 + j] = d;
                    closest_id[i * 5 + j] = x;
                    id[++r] = e[i], dis[r] = d + len[i], from[r] = node;
                    break;
                }
            }
            // printf("\n");
        }
    }
}

void MashPlacement::KPlacementDeviceArrays::initializeDeviceArrays(Tree *t)
{
    size_t numSequences = this->numSequences;
    size_t totalNodes = t->allNodes.size();

    Node *root = t->root;
    if (root == nullptr)
    {
        fprintf(stderr, "Error: Tree root is null.\n");
        return;
    }
    if (root->children.empty())
    {
        fprintf(stderr, "Error: Tree has no children.\n");
        return;
    }
    if (totalNodes < 2)
    {
        fprintf(stderr, "Error: Tree has less than two nodes.\n");
        return;
    }

    int *h_head = new int[numSequences * 2];
    int *h_e = new int[numSequences * 8];
    int *h_nxt = new int[numSequences * 8];
    int *h_belong = new int[numSequences * 8];
    double *h_len = new double[numSequences * 8];

    // initialize arrays
    for (int i = 0; i < numSequences * 2; i++)
        h_head[i] = -1;
    for (int i = 0; i < numSequences * 8; i++)
    {
        h_e[i] = -1;
        h_nxt[i] = -1;
        h_len[i] = 2;
        h_belong[i] = -1;
    }

    // dfs postorder traversal to initialize the tree structure
    std::function<void(Node *, size_t &)> dfs = [&](Node *node, size_t &edgeCount)
    {
        if (node == nullptr)
            return;
        for (Node *child : node->children)
        {
            dfs(child, edgeCount);
        }
        if (node->parent == nullptr)
            return;
        // if (node->parent->parent == nullptr) return; // skip root's children
        int x = node->idx;
        int y = node->parent->idx;
        // child->parent edge (x->y)
        h_e[edgeCount] = y;
        h_len[edgeCount] = node->bl;
        h_belong[edgeCount] = x;
        h_nxt[edgeCount] = h_head[x];
        h_head[x] = edgeCount;
        edgeCount++;

        // parent->child edge (y->x)
        h_e[edgeCount] = x;
        h_len[edgeCount] = node->bl;
        h_belong[edgeCount] = y;
        h_nxt[edgeCount] = h_head[y];
        h_head[y] = edgeCount;
        edgeCount++;
        // }
    };
    size_t edgeCount = 0;
    dfs(root, edgeCount);

    // transfer data to device
    cudaError_t err;
    err = cudaMemcpy(d_head, h_head, numSequences * 2 * sizeof(int), cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed for d_head!\n");
        exit(1);
    }
    err = cudaMemcpy(d_e, h_e, numSequences * 8 * sizeof(int), cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed for d_e!\n");
        exit(1);
    }
    err = cudaMemcpy(d_len, h_len, numSequences * 8 * sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed for d_len!\n");
        exit(1);
    }
    err = cudaMemcpy(d_nxt, h_nxt, numSequences * 8 * sizeof(int), cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed for d_nxt!\n");
        exit(1);
    }
    err = cudaMemcpy(d_belong, h_belong, numSequences * 8 * sizeof(int), cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed for d_belong!\n");
        exit(1);
    }

    // initialize closest_dis and closest_id
    int *d_id;
    err = cudaMalloc(&d_id, numSequences * 2 * sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    int *d_from;
    err = cudaMalloc(&d_from, numSequences * 2 * sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    double *d_dis;
    err = cudaMalloc(&d_dis, numSequences * 2 * sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    initializeID<<<1024, 1024>>>(
        numSequences * 4 - 4,
        d_closest_dis,
        d_closest_id);

    for (int i = 0; i < this->backboneSize; i++)
    {
        updateClosestNodes<<<1, 1>>>(
            d_head,
            d_nxt,
            d_e,
            d_len,
            d_closest_dis,
            d_closest_id,
            i,
            d_id,
            d_from,
            d_dis);
    }

    return;
}

__global__ void initialize(
    int lim,
    int nodes,
    double *d_closest_dis,
    int *d_closest_id,
    int *head,
    int *nxt,
    int *belong,
    int *e)
{
    int tx = threadIdx.x, bs = blockDim.x;
    int bx = blockIdx.x, gs = gridDim.x;
    int idx = tx + bs * bx;
    for (int t = idx; t < lim; t += bs * gs)
    {
        if (t < lim)
        {
            for (int i = 0; i < 5; i++)
            {
                d_closest_dis[t * 5 + i] = 2;
                d_closest_id[t * 5 + i] = -1;
            }
            nxt[t] = -1;
            e[t] = -1;
            belong[t] = -1;
        }
        if (t < nodes)
            head[t] = -1;
    }
}




struct compare_tuple {
  __host__ __device__
  bool operator()(thrust::tuple<int,double,double> lhs, thrust::tuple<int,double,double> rhs)
  {
    return thrust::get<2>(lhs) < thrust::get<2>(rhs);
    //Always find the tuple whose third value (the criteria we want to minimize) is minimized
  }
};


struct find_if_tuple_bidir {
    int bidirNewedgeIdToAttach;
    __host__ __device__
    find_if_tuple_bidir(int id) : bidirNewedgeIdToAttach(id) {}
    __host__ __device__
    bool operator()(thrust::tuple<int,double,double> tup) const
    {
      return thrust::get<0>(tup) == bidirNewedgeIdToAttach;
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
    int *head,
    int *nxt,
    double *dis,
    int *e,
    double *len,
    int *belong,
    thrust::tuple<int, double, double> *minPos,
    int lim,
    double * closest_dis,
    int * closest_id,
    int totSeqNum
){
    int tx=threadIdx.x,bs=blockDim.x,bx=blockIdx.x,gs=gridDim.x;
    int idx_=tx+bs*bx;
    for (int idx=idx_; idx<lim; idx+=bs*gs){
        if(idx>=lim) return;
        if(idx>=num*4-4||belong[idx]<e[idx]){
            thrust::tuple <int,double,double> minTuple(0,0,2);
            minPos[idx]=minTuple;
            continue;
        }
        int x = belong[idx], oth = e[idx];
        int eid = idx, otheid;
        double dis1 = 0, dis2 = 0, val;
        for (int i = 0; i < 5; i++)
            if (closest_id[eid * 5 + i] != -1)
            {
                val = dis[closest_id[eid * 5 + i]] - closest_dis[eid * 5 + i];
                if (val > dis1)
                    dis1 = val;
            }
        otheid=head[oth];
        /* modified for fixed topology */
        while(otheid!=-1 && e[otheid]!=x) otheid=nxt[otheid];
        if(otheid==-1){
            thrust::tuple <int,double,double> minTuple(0,0,2);
            minPos[idx]=minTuple;
            continue;
        }
        for(int i=0;i<5;i++)
            if(closest_id[otheid*5+i]!=-1){
                val = dis[closest_id[otheid*5+i]]-closest_dis[otheid*5+i];
                if(val>dis2) dis2=val;
            }
        double additional_dis = (dis1 + dis2 - len[eid]) / 2;
        if (additional_dis < 0)
            additional_dis = 0;
        dis1 -= additional_dis, dis2 -= additional_dis;
        if (dis1 < 0)
            dis1 = 0;
        if (dis2 < 0)
            dis2 = 0;
        if (dis1 > len[eid])
            additional_dis += dis1 - len[eid], dis1 = len[eid];
        if (dis2 > len[eid])
            additional_dis += dis2 - len[eid], dis2 = len[eid];
        // assert(dis1+dis2-1e-6<=len[eid]);
        double rest=len[eid]-dis1-dis2;
        dis1+=rest/2,dis2+=rest/2;
        // print eif, dis1, and additional_dis
        // printf("eid %d (nodes %d-%d): dis1=%.8lf, dis2=%.8lf, len=%.8lf, add=%.8lf\n", eid, x, oth, dis1, dis2, len[eid], additional_dis);
        thrust::tuple <int,double,double> minTuple(eid,dis1,additional_dis);
        minPos[idx]=minTuple;
    }
}

__global__ void updateTreeStructuretoAddQuery(
    int *head,
    int *nxt,
    int *e,
    double *len,
    double *closest_dis,
    int *closest_id,
    int *belong,
    int eid,
    double fracLen,
    double addLen,
    int placeId,   // Id of the newly placed node
    int edgeCount, // Position to insert a new edge in linked list
    int numSequences)
{
    int middle = placeId + numSequences - 1, outside = placeId;
    int x = belong[eid], y = e[eid];
    double originalDis = len[eid];
    int xe, ye;
    for (int i = head[x]; i != -1; i = nxt[i])
        if (e[i] == y)
        {
            e[i] = middle, len[i] = fracLen, xe = i;
            break;
        }
    for (int i = head[y]; i != -1; i = nxt[i])
        if (e[i] == x)
        {
            e[i] = middle, len[i] -= fracLen, ye = i;
            break;
        }
    /*
    Need to update:
    e, len, nxt, head, belong, closest_dis, closest_id
    */
    // middle -> x
    e[edgeCount] = x, len[edgeCount] = fracLen, nxt[edgeCount] = head[middle], head[middle] = edgeCount, belong[edgeCount] = middle;
    for (int i = 0; i < 5; i++)
        if (closest_id[ye * 5 + i] != -1)
        {
            closest_id[edgeCount * 5 + i] = closest_id[ye * 5 + i];
            closest_dis[edgeCount * 5 + i] = closest_dis[ye * 5 + i] + originalDis - fracLen;
        }
    edgeCount++;
    // middle -> y
    e[edgeCount] = y, len[edgeCount] = originalDis - fracLen, nxt[edgeCount] = head[middle], head[middle] = edgeCount, belong[edgeCount] = middle;
    for (int i = 0; i < 5; i++)
        if (closest_id[xe * 5 + i] != -1)
        {
            closest_id[edgeCount * 5 + i] = closest_id[xe * 5 + i];
            closest_dis[edgeCount * 5 + i] = closest_dis[xe * 5 + i] + fracLen;
        }
    edgeCount++;
    // outside -> middle
    e[edgeCount] = middle, len[edgeCount] = addLen, nxt[edgeCount] = head[outside], head[outside] = edgeCount, belong[edgeCount] = outside;
    edgeCount++;
    // middle -> outside
    e[edgeCount] = outside, len[edgeCount] = addLen, nxt[edgeCount] = head[middle], head[middle] = edgeCount, belong[edgeCount] = middle;
    int e1 = edgeCount - 2, e2 = edgeCount - 3;
    for (int i = 0; i < 5; i++)
    {
        if (closest_id[e1 * 5 + i] == -1)
            break;
        for (int j = 0; j < 5; j++)
            if (closest_dis[edgeCount * 5 + j] > closest_dis[e1 * 5 + i])
            {
                for (int k = 4; k > j; k--)
                {
                    closest_dis[edgeCount * 5 + k] = closest_dis[edgeCount * 5 + k - 1];
                    closest_id[edgeCount * 5 + k] = closest_id[edgeCount * 5 + k - 1];
                }
                closest_dis[edgeCount * 5 + j] = closest_dis[e1 * 5 + i];
                closest_id[edgeCount * 5 + j] = closest_id[e1 * 5 + i];
                break;
            }
    }
    for (int i = 0; i < 5; i++)
    {
        if (closest_id[e2 * 5 + i] == -1)
            break;
        for (int j = 0; j < 5; j++)
            if (closest_dis[edgeCount * 5 + j] > closest_dis[e2 * 5 + i])
            {
                for (int k = 4; k > j; k--)
                {
                    closest_dis[edgeCount * 5 + k] = closest_dis[edgeCount * 5 + k - 1];
                    closest_id[edgeCount * 5 + k] = closest_id[edgeCount * 5 + k - 1];
                }
                closest_dis[edgeCount * 5 + j] = closest_dis[e2 * 5 + i];
                closest_id[edgeCount * 5 + j] = closest_id[e2 * 5 + i];
                break;
            }
    }
    edgeCount++;
}

__global__ void updateTreeStructure(
    int *head,
    int *nxt,
    int *e,
    double *len,
    double *closest_dis,
    int *closest_id,
    int *belong,
    int eid,
    double fracLen,
    double addLen,
    int placeId,   // Id of the newly placed node
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
            // printf("x -> middle: eid: %d\n", i);
            break;
        }
    for(int i=head[y];i!=-1;i=nxt[i])
        if(e[i]==x){
            e[i]=middle,len[i]-=fracLen,ye=i;
            // printf("y -> middle: eid: %d\n", i);
            break;
        }
    /*
    Need to update:
    e, len, nxt, head, belong, closest_dis, closest_id
    */
    //middle -> x
    e[edgeCount]=x,len[edgeCount]=fracLen,nxt[edgeCount]=head[middle],head[middle]=edgeCount,belong[edgeCount]=middle;
    // printf("middle -> x: eid: %d\n", edgeCount);
    for(int i=0;i<5;i++)
        if(closest_id[ye*5+i]!=-1){
            closest_id[edgeCount*5+i]=closest_id[ye*5+i];
            closest_dis[edgeCount*5+i]=closest_dis[ye*5+i]+originalDis-fracLen;
        }
    edgeCount++;
    //middle -> y
    e[edgeCount]=y,len[edgeCount]=originalDis-fracLen,nxt[edgeCount]=head[middle],head[middle]=edgeCount,belong[edgeCount]=middle;
    // printf("middle -> y: eid: %d\n", edgeCount);
    for(int i=0;i<5;i++)
        if(closest_id[xe*5+i]!=-1){
            closest_id[edgeCount*5+i]=closest_id[xe*5+i];
            closest_dis[edgeCount*5+i]=closest_dis[xe*5+i]+fracLen;
        }
    edgeCount++;
    //outside -> middle
    e[edgeCount]=middle,len[edgeCount]=addLen,nxt[edgeCount]=head[outside],head[outside]=edgeCount,belong[edgeCount]=outside;
    // printf("outside -> middle: eid: %d\n", edgeCount);
    edgeCount++;
    //middle -> outside
    e[edgeCount]=outside,len[edgeCount]=addLen,nxt[edgeCount]=head[middle],head[middle]=edgeCount,belong[edgeCount]=middle;
    // printf("middle -> outside: eid: %d\n", edgeCount);
    int e1=edgeCount-2, e2=edgeCount-3;
    for(int i=0;i<5;i++){
        if(closest_id[e1*5+i]==-1) break;
        for(int j=0;j<5;j++)
            if(closest_dis[edgeCount*5+j]>closest_dis[e1*5+i]){
                for(int k=4;k>j;k--){
                    closest_dis[edgeCount*5+k]=closest_dis[edgeCount*5+k-1];
                    closest_id[edgeCount*5+k]=closest_id[edgeCount*5+k-1];
                }
                closest_dis[edgeCount * 5 + j] = closest_dis[e1 * 5 + i];
                closest_id[edgeCount * 5 + j] = closest_id[e1 * 5 + i];
                break;
            }
    }
    for (int i = 0; i < 5; i++)
    {
        if (closest_id[e2 * 5 + i] == -1)
            break;
        for (int j = 0; j < 5; j++)
            if (closest_dis[edgeCount * 5 + j] > closest_dis[e2 * 5 + i])
            {
                for (int k = 4; k > j; k--)
                {
                    closest_dis[edgeCount * 5 + k] = closest_dis[edgeCount * 5 + k - 1];
                    closest_id[edgeCount * 5 + k] = closest_id[edgeCount * 5 + k - 1];
                }
                closest_dis[edgeCount * 5 + j] = closest_dis[e2 * 5 + i];
                closest_id[edgeCount * 5 + j] = closest_id[e2 * 5 + i];
                break;
            }
    }
    edgeCount++;
}

__global__ void buildInitialTree(
    int numSequences,
    int *head,
    int *e,
    double *len,
    int *nxt,
    int *belong,
    double *dis,
    int edgeCount)
{
    int nv = numSequences;
    double d = dis[0];
    // 0 -> nv
    e[edgeCount]=nv,len[edgeCount]=d/2,nxt[edgeCount]=head[0],head[0]=edgeCount,belong[edgeCount]=0;
    // printf("0 -> nv: eid: %d, from: %d, to: %d\n", edgeCount);
    edgeCount++;
    // 1 -> nv
    e[edgeCount]=nv,len[edgeCount]=d/2,nxt[edgeCount]=head[1],head[1]=edgeCount,belong[edgeCount]=1;
    // printf("1 -> nv: eid: %d\n", edgeCount);
    edgeCount++;
    // nv -> 0
    e[edgeCount] = 0, len[edgeCount] = d / 2, nxt[edgeCount] = head[nv], head[nv] = edgeCount, belong[edgeCount] = nv;
    edgeCount++;
    // nv -> 0
    e[edgeCount] = 1, len[edgeCount] = d / 2, nxt[edgeCount] = head[nv], head[nv] = edgeCount, belong[edgeCount] = nv;
    edgeCount++;
}

void MashPlacement::KPlacementDeviceArrays::deallocateDeviceArrays()
{
    cudaFree(d_head);
    cudaFree(d_e);
    cudaFree(d_nxt);
    cudaFree(d_belong);
    cudaFree(d_closest_id);
    cudaFree(d_dist);
    cudaFree(d_len);
    cudaFree(d_closest_dis);
}


void MashPlacement::KPlacementDeviceArrays::printTree(std::vector <std::string> name, std::ofstream& output_){
    int * h_head = new int[numSequences*2];
    int * h_e = new int[numSequences*8];
    int * h_nxt = new int[numSequences*8];
    double * h_len = new double[numSequences*8];
    
    double * h_closest_dis = new double[numSequences*20];
    int * h_closest_id = new int[numSequences*20];
    auto err = cudaMemcpy(h_head, d_head, numSequences*2*sizeof(int),cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        exit(1);
    }
    err = cudaMemcpy(h_e, d_e, numSequences * 8 * sizeof(int), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        exit(1);
    }
    err = cudaMemcpy(h_nxt, d_nxt, numSequences * 8 * sizeof(int), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        exit(1);
    }
    err = cudaMemcpy(h_len, d_len, numSequences * 8 * sizeof(double), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        exit(1);
    }


    // print h_head, h_e, h_nxt, h_len
    // std::cerr<<"Tree: ";
    // for(int i=0;i<numSequences*2;i++){
    //     std::cerr << h_head[i] << "\t";
    // }
    // std::cerr<<"\n";
    // for(int i=0;i<numSequences*8;i++){
    //     std::cerr << h_e[i] << "\t";
    // }
    // std::cerr<<"\n";
    // for(int i=0;i<numSequences*8;i++){
    //     std::cerr << h_nxt[i] << "\t";
    // }
    // std::cerr<<"\n";
    // for(int i=0;i<numSequences*8;i++){
    //     std::cerr << h_len[i] << "\t";
    // }
    // std::cerr<<"\n";


    printNewickFromHostAdjacency(output_, name, h_head, h_e, h_nxt, h_len, numSequences + bd - 2,
                                 g_printBinaryNewick);
    output_ << ";\n";
}

void MashPlacement::KPlacementDeviceArrays::printTree(std::vector <std::string> name, std::ofstream& output_, UnrootedTree* t, std::vector<int>& edgeIdsMappingToNewTreeEdges){
    double * h_len = new double[numSequences*8];
    auto err = cudaMemcpy(h_len, d_len, numSequences*8*sizeof(double),cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Gpu_ERROR: cudaMemcpy failed!\n");
        exit(1);
    }

    // Iterate over all edges in the tree and update the edge length based on edgeIdsMappingToNewTreeEdges
    for (auto &nodePair : t->getNodes()) {
        UnrootedNode* node = nodePair.second;
        for (auto &edge : node->neighbors) {
            int edge_id = edge.edge_id;
            assert(edge_id >= 0 && edge_id < numSequences * 8 && "Edge ID out of bounds");
            edge.length = h_len[edge_id];
        }
    }

    if (!g_printBinaryNewick) {
        t->collapseZeroLengthInternalEdges();
    }
    output_ << t->toNewick() << std::endl;
    return;

}

void MashPlacement::KPlacementDeviceArrays::findPlacementTree(
    Param &params,
    const MashDeviceArrays &mashDeviceArrays,
    MatrixReader &matrixReader,
    const MSADeviceArrays &msaDeviceArrays)
{
    if (params.in == "d")
    {
        matrixReader.distConstructionOnGpu(params, 0, d_dist);
    }
    int *d_id;
    auto err = cudaMalloc(&d_id, numSequences * 2 * sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    int *d_from;
    err = cudaMalloc(&d_from, numSequences * 2 * sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    double *d_dis;
    err = cudaMalloc(&d_dis, numSequences * 2 * sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }

    /*
    Initialize closest nodes by inifinite
    */
    // int threadNum = 256, blockNum = (numSequences*4-4+threadNum-1)/threadNum;
    int threadNum = 1024, blockNum = 1024;
    initialize<<<blockNum, threadNum>>>(
        numSequences * 4 - 4,
        numSequences * 2,
        d_closest_dis,
        d_closest_id,
        d_head,
        d_nxt,
        d_belong,
        d_e);
    /*
    Build Initial Tree
    */
    if (params.in == "r")
    {
        mashDeviceArrays.distConstructionOnGpu(
            params,
            1,
            d_dist);
    }
    else if (params.in == "d")
    {
        matrixReader.distConstructionOnGpu(
            params,
            1,
            d_dist);
    }
    else if (params.in == "m")
    {
        msaDeviceArrays.distConstructionOnGpu(
            params,
            1,
            d_dist);
    }
    // cudaDeviceSynchronize();

    // return;
    // double * h_dis = new double[numSequences];
    // cudaMemcpy(h_dis,d_dist,numSequences*sizeof(double),cudaMemcpyDeviceToHost);
    // for(int j=0;j<1;j++) fprintf(stderr,"%.8lf ",h_dis[j]);std::cerr<<'\n';

    buildInitialTree <<<1,1>>> (
        numSequences,
        d_head,
        d_e,
        d_len,
        d_nxt,
        d_belong,
        d_dist,
        idx);
    idx += 4;
    /*
    Initialize closest nodes by inital tree
    */
    for (int i = 0; i < bd; i++)
    {
        updateClosestNodes<<<1, 1>>>(
            d_head,
            d_nxt,
            d_e,
            d_len,
            d_closest_dis,
            d_closest_id,
            i,
            d_id,
            d_from,
            d_dis);
    }

    thrust::device_vector<thrust::tuple<int, double, double>> minPos(numSequences * 4 - 4);
    // std::cout<<"FFF\n";
    std::chrono::nanoseconds disTime(0), treeTime(0);
    for (int i = bd; i < numSequences; i++)
    {
        auto disStart = std::chrono::high_resolution_clock::now();
        // blockNum = (i + 255) / 256;
        // blockNum = 1024;
        if (params.in == "r")
        {
            mashDeviceArrays.distConstructionOnGpu(
                params,
                i,
                d_dist);
        }
        else if (params.in == "d")
        {
            matrixReader.distConstructionOnGpu(
                params,
                i,
                d_dist);
        }
        else if (params.in == "m")
        {
            msaDeviceArrays.distConstructionOnGpu(
                params,
                i,
                d_dist);
        }
        cudaDeviceSynchronize();

        // double * h_dis = new double[numSequences];
        // cudaMemcpy(h_dis,d_dist,numSequences*sizeof(double),cudaMemcpyDeviceToHost);
        // fprintf(stderr, "%d\n",i);
        // for(int j=0;j<i;j++) std::cerr<<h_dis[j]<<" ";std::cerr<<'\n';

        auto disEnd = std::chrono::high_resolution_clock::now();
        auto treeStart = std::chrono::high_resolution_clock::now();
        
        calculateBranchLength <<<blockNum,threadNum>>> (
            i,
            d_head,
            d_nxt,
            d_dist,
            d_e,
            d_len,
            d_belong,
            thrust::raw_pointer_cast(minPos.data()),
            numSequences * 4 - 4,
            d_closest_dis,
            d_closest_id,
            numSequences
        );
        
        auto iter=thrust::min_element(minPos.begin(),minPos.end(),compare_tuple());
        thrust::tuple<int,double,double> smallest=*iter;
        // /* print top 5 sorted elements */
        // for(int j=0;j<5;j++){
        //     thrust::tuple<int,double,double> s=*iter;
        //     std::cerr<<thrust::get<0>(s)<<" "<<thrust::get<1>(s)<<" "<<thrust::get<2>(s)<<'\n';
        //     iter++;
        // } 
        
        /*
        Update Tree (and assign closest nodes to newly added nodes)
        */
        int eid = thrust::get<0>(smallest);
        double fracLen = thrust::get<1>(smallest), addLen = thrust::get<2>(smallest);
        updateTreeStructure<<<1, 1>>>(
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
            numSequences);
        idx += 4;
        /*
        Update closest nodes
        */
        updateClosestNodes<<<1, 1>>>(
            d_head,
            d_nxt,
            d_e,
            d_len,
            d_closest_dis,
            d_closest_id,
            i,
            d_id,
            d_from,
            d_dis);
        // cudaDeviceSynchronize();
        auto treeEnd = std::chrono::high_resolution_clock::now();
        disTime += disEnd - disStart;
        treeTime += treeEnd - treeStart;
        // std::cerr << "Seq " << i << " at " << eid << " with fracLen " << fracLen 
        //           << " and addLen " << addLen << std::endl;
    }
    std::cerr << "Distance Operation Time " << disTime.count() / 1000000 << " ms\n";
    std::cerr << "Tree Operation Time " << treeTime.count() / 1000000 << " ms\n";
}





void MashPlacement::KPlacementDeviceArrays::addQuery(
    Param &params,
    const MashDeviceArrays &mashDeviceArrays,
    MatrixReader &matrixReader,
    const MSADeviceArrays &msaDeviceArrays)
{

    int *d_id;
    auto err = cudaMalloc(&d_id, numSequences * 2 * sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    int *d_from;
    err = cudaMalloc(&d_from, numSequences * 2 * sizeof(int));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    double *d_dis;
    err = cudaMalloc(&d_dis, numSequences * 2 * sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Gpu_ERROR: cudaMalloc failed!\n");
        exit(1);
    }
    int threadNum = 256, blockNum = (numSequences * 4 - 4 + threadNum - 1) / threadNum;
    idx += 4 * this->backboneSize - 4; // Adjust idx to account for the backbone tree size

    thrust::device_vector<thrust::tuple<int, double, double>> minPos(numSequences * 4 - 4);
    std::chrono::nanoseconds disTime(0), treeTime(0);
    for (int i = this->backboneSize; i < numSequences; i++)
    {
        auto disStart = std::chrono::high_resolution_clock::now();
        blockNum = (i + 255) / 256;
        // blockNum = 1024;
        if (params.in == "r")
        {
            mashDeviceArrays.distConstructionOnGpu(
                params,
                i,
                d_dist);
        }
        else if (params.in == "d")
        {
            matrixReader.distConstructionOnGpu(
                params,
                i,
                d_dist);
        }
        else if (params.in == "m")
        {
            msaDeviceArrays.distConstructionOnGpu(
                params,
                i,
                d_dist);
        }
        cudaDeviceSynchronize();

        // double * h_dis = new double[numSequences];
        // cudaMemcpy(h_dis,d_dist,numSequences*sizeof(double),cudaMemcpyDeviceToHost);
        // fprintf(stderr, "%d\n",i);
        // for(int j=0;j<i;j++) std::cerr<<h_dis[j]<<" ";std::cerr<<'\n';

        auto disEnd = std::chrono::high_resolution_clock::now();
        auto treeStart = std::chrono::high_resolution_clock::now();
        blockNum = (numSequences * 4 - 4 + 255) / 256;
        // blockNum = 1024;
        calculateBranchLength<<<blockNum, threadNum>>>(
            i,
            d_head,
            d_nxt,
            d_dist,
            d_e,
            d_len,
            d_belong,
            thrust::raw_pointer_cast(minPos.data()),
            numSequences * 4 - 4,
            d_closest_dis,
            d_closest_id,
            numSequences);
        auto iter = thrust::min_element(minPos.begin(), minPos.end(), compare_tuple());
        thrust::tuple<int, double, double> smallest = *iter;
        /*
        Update Tree (and assign closest nodes to newly added nodes)
        */
        int eid = thrust::get<0>(smallest);
        double fracLen = thrust::get<1>(smallest), addLen = thrust::get<2>(smallest);
        updateTreeStructure<<<1, 1>>>(
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
            numSequences);
        idx += 4;
        /*
        Update closest nodes
        */
        updateClosestNodes<<<1, 1>>>(
            d_head,
            d_nxt,
            d_e,
            d_len,
            d_closest_dis,
            d_closest_id,
            i,
            d_id,
            d_from,
            d_dis);
        // cudaDeviceSynchronize();
        auto treeEnd = std::chrono::high_resolution_clock::now();
        disTime += disEnd - disStart;
        treeTime += treeEnd - treeStart;
    }
    std::cerr << "Distance Operation Time " <<  disTime.count()/1000000 << " ms\n";
    std::cerr << "Tree Operation Time " <<  treeTime.count()/1000000 << " ms\n";

}

void MashPlacement::KPlacementDeviceArrays::estimateBranchLengthsFromTopology(
    Param& params,
    const MashDeviceArrays& mashDeviceArrays,
    MatrixReader& matrixReader,
    const MSADeviceArrays& msaDeviceArrays,
    std::vector<int>& edgeMapOldToNew,
    UnrootedTree* t,
    std::vector<std::string>& names
){

    std::unordered_map<int, int> bidirEdgeMap;
    std::unordered_map<int, std::vector<int>> edgeMapNewToOld;
    std::unordered_map<int, bool> visistedEdges;
    std::unordered_map<int, bool> visistedNodes;

    
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

    
    int threadNum = 1024, blockNum = 1024;
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
    // cudaDeviceSynchronize();

    // return;
    // double * h_dis = new double[numSequences];
    // cudaMemcpy(h_dis,d_dist,numSequences*sizeof(double),cudaMemcpyDeviceToHost);
    // for(int j=0;j<1;j++) fprintf(stderr,"%.8lf ",h_dis[j]);std::cerr<<'\n';

    // std::cout << "==================== Building initial tree ====================" << std::endl;
    // std::cout << names[0] << " " << names[1] << std::endl;

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
    cudaDeviceSynchronize();
    idx += 4;

    /* edges 0 and 2 are same edge in different directions */
    std::vector<UnrootedEdge> edges = t->edgesBetween(names[0], t->getRoot());
    edgeMapNewToOld[0]=std::vector<int>();
    for(auto &e: edges){
        visistedEdges[e.edge_id]=true;
        edgeMapOldToNew[e.edge_id]=0;
        edgeMapNewToOld[0].push_back(e.edge_id);
    }
    edges = t->edgesBetween(names[1], t->getRoot());
    edgeMapNewToOld[1]=std::vector<int>();
    for(auto &e: edges){
        visistedEdges[e.edge_id]=true;
        edgeMapOldToNew[e.edge_id]=1;
        edgeMapNewToOld[1].push_back(e.edge_id);
    }
    
    /* edges 0 and 2 are same edge in different directions (based on buildInitialTree) */
    bidirEdgeMap[0]=2;
    bidirEdgeMap[2]=0;
    /* edges 1 and 3 are same edge in different directions (based on buildInitialTree) */
    bidirEdgeMap[1]=3;
    bidirEdgeMap[3]=1;

    // std::cout << "bidirEdgeMap: " << std::endl;
    // for(auto &e: bidirEdgeMap){
    //     std::cout << "edge " << e.first << " -> " << e.second << std::endl;
    // }

    // std::cout << "edgeMapNewToOld: " << std::endl;
    // for(auto &e: edgeMapNewToOld){
    //     std::cout << "edge " << e.first << " -> ";
    //     for(auto &f: e.second){
    //         std::cout << f << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // std::cout << "edgeMapOldToNew: " << std::endl;
    // for (int i=0;i<edgeMapOldToNew.size();i++){
    //     std::cout << "edge " << i << " -> " << edgeMapOldToNew[i] << std::endl;
    // }
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

    for(int i=bd;i<numSequences;i++){
        std::string queryName = names[i];
        std::vector<UnrootedEdge> edges = t->edgesBetween(names[i], names[i-1]);

        int oldedgeIdToAttach = -1;
        int newedgeIdToAttach = -1;
        for (auto &e: edges){
            if (visistedEdges.find(e.edge_id) != visistedEdges.end() && visistedEdges[e.edge_id]) {
                oldedgeIdToAttach = e.edge_id;
                break;
            }
        }
        if(oldedgeIdToAttach == -1){
            std::cerr << "Error: Edge not found" << std::endl;
            exit(1);
        }
        // std::cout << "oldedgeIdToAttach: " << oldedgeIdToAttach << std::endl;

        /* update edgeMapNewToOld */
        int bidirNewedgeIdToAttach=-1;
        assert(edgeMapOldToNew[oldedgeIdToAttach] != -1 && "Edge not found in edgeMapOldToNew");
        newedgeIdToAttach=edgeMapOldToNew[oldedgeIdToAttach];
        if (edgeMapNewToOld.find(newedgeIdToAttach) == edgeMapNewToOld.end()){
            newedgeIdToAttach = bidirEdgeMap[newedgeIdToAttach];
        }
        bidirNewedgeIdToAttach = bidirEdgeMap[newedgeIdToAttach];
        
        edgeMapNewToOld[bidirNewedgeIdToAttach]=std::vector<int>();
        // std::cout << "moving edges: ";
        for (auto &e: edges){
            if (visistedEdges[e.edge_id] && edgeMapOldToNew[e.edge_id] == newedgeIdToAttach) {
                // std::cout << e.edge_id << " ";
                edgeMapNewToOld[bidirNewedgeIdToAttach].push_back(e.edge_id);
                // Remove e.edge_id from edgeMapNewToOld[newedgeIdToAttach] vector
                auto& vec = edgeMapNewToOld[newedgeIdToAttach];
                vec.erase(std::remove(vec.begin(), vec.end(), e.edge_id), vec.end());
                edgeMapOldToNew[e.edge_id]=bidirNewedgeIdToAttach;
            }
        }
        // std::cout << std::endl;

        /* update bidirEdgeMap */
        bidirEdgeMap[newedgeIdToAttach]=idx;
        bidirEdgeMap[idx]=newedgeIdToAttach;
        bidirEdgeMap[bidirNewedgeIdToAttach]=idx+1;
        bidirEdgeMap[idx+1]=bidirNewedgeIdToAttach;

        /* add new tip information */
        edgeMapNewToOld[idx+2]=std::vector<int>();
        for (auto &e: edges){
            if (!visistedEdges[e.edge_id]) {
                edgeMapNewToOld[idx+2].push_back(e.edge_id);
                edgeMapOldToNew[e.edge_id]=idx+2;
                visistedEdges[e.edge_id]=true;
            }
        }
        bidirEdgeMap[idx+2]=idx+3;
        bidirEdgeMap[idx+3]=idx+2;

        // std::cout << "bidirEdgeMap: " << std::endl;
        // for(auto &e: bidirEdgeMap){
        //     std::cout << "edge " << e.first << " -> " << e.second << std::endl;
        // }

        // std::cout << "edgeMapNewToOld: " << std::endl;
        // for(auto &e: edgeMapNewToOld){
        //     std::cout << "edge " << e.first << " -> ";
        //     for(auto &f: e.second){
        //         std::cout << f << " ";
        //     }
        //     std::cout << std::endl;
        // }

        // std::cout << "edgeMapOldToNew: " << std::endl;
        // for (int i=0;i<edgeMapOldToNew.size();i++){
        //     std::cout << "edge " << i << " -> " << edgeMapOldToNew[i] << std::endl;
        // }
        
        if(params.in == "r"){
            mashDeviceArrays.distConstructionOnGpu(
                params,
                i,
                d_dist
            );
        } else if(params.in == "m"){
            msaDeviceArrays.distConstructionOnGpu(
                params,
                i,
                d_dist
            );
        }
        cudaDeviceSynchronize();

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
        
        auto iter = thrust::find_if(minPos.begin(), minPos.end(), find_if_tuple_bidir(newedgeIdToAttach));
        if (iter == minPos.end()) {
            iter = thrust::find_if(minPos.begin(), minPos.end(), find_if_tuple_bidir(bidirNewedgeIdToAttach));
        }
        if (iter == minPos.end()) {
            std::cerr << "Error: Edge not found in minPos" << std::endl;
            exit(1);
        }
        
        thrust::tuple<int,double,double> smallest=*iter;
        
        /*
        Update Tree (and assign closest nodes to newly added nodes)
        */
        int eid=newedgeIdToAttach;
        double fracLen=thrust::get<1>(smallest),addLen=thrust::get<2>(smallest);
        // std::cout << "eid: " << eid << " fracLen: " << fracLen << " addLen: " << addLen << std::endl;
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
        cudaDeviceSynchronize();
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
    }
}