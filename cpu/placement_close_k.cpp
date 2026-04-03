/* cpu/placement_close_k.cpp: CPU k-closest placement. Host arrays, backbone init, findPlacementTree/addQuery. */

#include "mash_placement.hpp"

#include <stdio.h>
#include <queue>
#include <tbb/parallel_for.h>
#include <chrono>
#include <iostream>
#include <fstream>
#include <cassert>
#include <functional>

/** Allocate host arrays for k-closest placement; backboneSize for add-query mode. */
void MashPlacement::KPlacementDeviceArrays::allocateDeviceArrays(size_t num, int backboneSize)
{

    numSequences = int(num);
    bd = 2, idx = 0;
    this->backboneSize = backboneSize;
    d_dist = new double[numSequences];
    d_head = new int[numSequences * 2];
    d_e = new int[numSequences * 8];
    d_len = new double[numSequences * 8];
    d_nxt = new int[numSequences * 8];
    d_belong = new int[numSequences * 8];
    d_closest_id = new int[numSequences * 20];
    d_closest_dis = new double[numSequences * 20];
    for (int i = 0; i < numSequences * 20; ++i)
    {
        d_closest_dis[i] = 0;
        d_closest_id[i] = 0;
        if (i < numSequences * 8)
        {
            d_e[i] = 0;
            d_len[i] = 0;
            d_nxt[i] = 0;
            d_belong[i] = 0;
        }
        if (i < numSequences * 2)
        {
            d_head[i] = 0;
        }
        if (i < numSequences)
        {
            d_dist[i] = 0;
        }
    }
}

void initializeID(
    int lim,
    double *d_closest_dis,
    int *d_closest_id)
{
    tbb::parallel_for(tbb::blocked_range<int>(0, lim), [&](tbb::blocked_range<int> range)
                      { 
    for (int idx = range.begin(); idx < range.end(); ++idx) {
        for (int i = 0; i < 5; i++) {
            d_closest_dis[idx * 5 + i] = 2;
            d_closest_id[idx * 5 + i] = -1;
        }
    } });
}

void updateClosestNodes(
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

/** Build adjacency from Tree t via DFS; init closest from backbone. */
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
    std::memcpy(d_head, h_head, numSequences * 2 * sizeof(int));
    std::memcpy(d_e, h_e, numSequences * 8 * sizeof(int));
    std::memcpy(d_len, h_len, numSequences * 8 * sizeof(double));
    std::memcpy(d_nxt, h_nxt, numSequences * 8 * sizeof(int));
    std::memcpy(d_belong, h_belong, numSequences * 8 * sizeof(int));

    // initialize closest_dis and closest_id
    int *d_id = new int[numSequences * 2];
    int *d_from = new int[numSequences * 2];
    double *d_dis = new double[numSequences * 2];

    initializeID(
        numSequences * 4 - 4,
        d_closest_dis,
        d_closest_id);

    for (int i = 0; i < this->backboneSize; i++)
    {
        updateClosestNodes(
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

/** Reset adjacency and closest_id/closest_dis to sentinels. */
void initialize(
    int lim,
    int nodes,
    double *d_closest_dis,
    int *d_closest_id,
    int *head,
    int *nxt,
    int *belong,
    int *e)
{
    tbb::parallel_for(tbb::blocked_range<int>(0, lim), [&](tbb::blocked_range<int> range)
                      { 
    for (int idx = range.begin(); idx < range.end(); ++idx) {
    // for (int idx = 0; idx < lim; ++idx) {
        for(int i=0;i<5;i++){
            d_closest_dis[idx*5+i]=2;
            d_closest_id[idx*5+i]=-1;
        }
        nxt[idx] = -1;
        e[idx] = -1;
        belong[idx] = -1;
        if(idx<nodes) head[idx] = -1;
    } });
}

/** Compare placement tuples by addLen (third element). */
bool compare_tuple(std::tuple<int, double, double> lhs, std::tuple<int, double, double> rhs)
{
    return std::get<2>(lhs) < std::get<2>(rhs);
}

/** For each candidate edge, compute (eid, fracLen, addLen) from closest_id/closest_dis. */
void calculateBranchLength(
    int num, // should be bd, not numSequences
    int *head,
    int *nxt,
    double *dis,
    int *e,
    double *len,
    int *belong,
    std::tuple<int, double, double> *minPos,
    int lim,
    double *closest_dis,
    int *closest_id,
    int totSeqNum)
{

    tbb::parallel_for(tbb::blocked_range<int>(0, lim), [&](tbb::blocked_range<int> range)
                      { 
    for (int idx = range.begin(); idx < range.end(); ++idx) {
    // for (int idx = 0; idx < lim; ++idx) {
        if(idx>=num*4-4||belong[idx]<e[idx]){
            std::tuple <int,double,double> minTuple(0,0,2);
            minPos[idx]=minTuple;
            continue;
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
        std::tuple <int,double,double> minTuple(eid,dis1,additional_dis);
        // printf("%d: <%d,%f,%f>\n", idx, std::get<0>(minTuple),std::get<1>(minTuple),std::get<2>(minTuple));
        minPos[idx]=minTuple;
    } });
}

/*
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
{s
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
    // !!!
    // Need to update:
    // e, len, nxt, head, belong, closest_dis, closest_id
    // !!!

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
*/

/** Insert new leaf on edge eid; update closest_id/closest_dis for new edges. */
void updateTreeStructure(
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

/** Build 3-node tree (0, 1, nv) with two edges of length dis[0]/2. */
void buildInitialTree(
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
    e[edgeCount] = nv, len[edgeCount] = d / 2, nxt[edgeCount] = head[0], head[0] = edgeCount, belong[edgeCount] = 0;
    edgeCount++;
    // 1 -> nv
    e[edgeCount] = nv, len[edgeCount] = d / 2, nxt[edgeCount] = head[1], head[1] = edgeCount, belong[edgeCount] = 1;
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
    delete[] d_head;
    delete[] d_e;
    delete[] d_nxt;
    delete[] d_belong;
    delete[] d_closest_id;
    delete[] d_dist;
    delete[] d_len;
    delete[] d_closest_dis;
}

/** Copy tree to local buffers and emit Newick via recursive DFS. */
void MashPlacement::KPlacementDeviceArrays::printTree(std::vector<std::string> name, std::ofstream &output_)
{
    int *h_head = new int[numSequences * 2];
    int *h_e = new int[numSequences * 8];
    int *h_nxt = new int[numSequences * 8];
    double *h_len = new double[numSequences * 8];
    double *h_closest_dis = new double[numSequences * 20];
    int *h_closest_id = new int[numSequences * 20];
    std::function<void(int,int)> print = [&](int node, int from) {
        if (h_nxt[h_head[node]] != -1) {
            output_ << "(";
            std::vector<std::pair<int,int>> pos; // {edge_index, parent_node}
            for (int i = h_head[node]; i != -1; i = h_nxt[i]) {
                if (h_e[i] != from) {
                    if (h_len[i] == 0) {
                        int collapsed = h_e[i];
                        if (h_head[collapsed] != -1) {
                            for (int j = h_head[collapsed]; j != -1; j = h_nxt[j]) {
                                if (h_e[j] != node) {
                                    pos.push_back({j, collapsed});
                                }
                            }
                        }
                    } else {
                        pos.push_back({i, node});
                    }
                }
            }
            for (size_t i = 0; i < pos.size(); i++) {
                auto [edgeIdx, parent] = pos[i];
                print(h_e[edgeIdx], parent);
                output_ << ":";
                output_ << h_len[edgeIdx] << (i+1 == pos.size() ? ')' : ',');
            }
        } else {
            output_ << name[node];
        }
    };
    
    for (int i = 0; i < numSequences * 2; ++i)
    {
        h_head[i] = d_head[i];
    }
    for (int i = 0; i < numSequences * 8; ++i)
    {
        h_e[i] = d_e[i];
        h_nxt[i] = d_nxt[i];
        h_len[i] = d_len[i];
    }

    // auto err = cudaMemcpy(h_head, d_head, numSequences*2*sizeof(int),cudaMemcpyDeviceToHost);
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

    print(numSequences + bd - 2, -1);
    // std::cout<<";\n";
    output_ << ";\n";
}

/** K-closest placement: init, build initial tree, updateClosestNodes for initial nodes;
 * then for each new taxon: dist -> calculateBranchLength -> updateTree -> updateClosestNodes. */
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
    int *d_from;
    double *d_dis;
    d_id = new int[numSequences * 2];
    d_from = new int[numSequences * 2];
    d_dis = new double[numSequences * 2];

    /*
    Initialize closest nodes by inifinite
    */

    initialize(
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

    // return;
    // double * h_dis = new double[numSequences];
    // cudaMemcpy(h_dis,d_dist,numSequences*sizeof(double),cudaMemcpyDeviceToHost);
    // fprintf(stderr, "%d\n",1);
    // for(int j=0;j<1;j++) fprintf(stderr,"%.8lf ",d_dist[j]);std::cerr<<'\n';
    buildInitialTree(
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
        updateClosestNodes(
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

    std::vector<std::tuple<int, double, double>> minPos(numSequences * 4 - 4);
    // std::cout<<"FFF\n";
    std::chrono::nanoseconds disTime(0), treeTime(0);
    for (int i = bd; i < numSequences; i++)
    {
        auto disStart = std::chrono::high_resolution_clock::now();
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
        // cudaDeviceSynchronize();

        // double * h_dis = new double[numSequences];
        // cudaMemcpy(h_dis,d_dist,numSequences*sizeof(double),cudaMemcpyDeviceToHost);
        // fprintf(stderr, "%d\n",i);
        // for(int j=0;j<i;j++) std::cerr<<h_dis[j]<<" ";std::cerr<<'\n';

        auto disEnd = std::chrono::high_resolution_clock::now();
        auto treeStart = std::chrono::high_resolution_clock::now();
        // blockNum = (numSequences*4-4 + 255) / 256;
        // blockNum = 1024;
        calculateBranchLength(
            i,
            d_head,
            d_nxt,
            d_dist,
            d_e,
            d_len,
            d_belong,
            minPos.data(),
            numSequences * 4 - 4,
            d_closest_dis,
            d_closest_id,
            numSequences);
        auto iter = std::min_element(minPos.begin(), minPos.end(), compare_tuple);
        std::tuple<int, double, double> smallest = *iter;
        /*
        Update Tree (and assign closest nodes to newly added nodes)
        */
        int eid = std::get<0>(smallest);
        double fracLen = std::get<1>(smallest), addLen = std::get<2>(smallest);
        updateTreeStructure(
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
        updateClosestNodes(
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
        // std::cerr << "Added seq " << i << "at " << eid << " with fracLen " << fracLen
        //           << " and addLen " << addLen << "\n";
    }
    std::cerr << "Distance Operation Time " << disTime.count() / 1000000 << " ms\n";
    std::cerr << "Tree Operation Time " << treeTime.count() / 1000000 << " ms\n";
}

/** Add queries to backbone: place each new sequence via k-closest, updateTree, updateClosestNodes. */
void MashPlacement::KPlacementDeviceArrays::addQuery(
    Param &params,
    const MashDeviceArrays &mashDeviceArrays,
    MatrixReader &matrixReader,
    const MSADeviceArrays &msaDeviceArrays)
{

    int *d_id;
    int *d_from;
    double *d_dis;
    d_id = new int[numSequences * 2];
    d_from = new int[numSequences * 2];
    d_dis = new double[numSequences * 2];

    int threadNum = 256, blockNum = (numSequences * 4 - 4 + threadNum - 1) / threadNum;
    idx += 4 * this->backboneSize - 4; // Adjust idx to account for the backbone tree size

    std::vector<std::tuple<int, double, double>> minPos(numSequences * 4 - 4);
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
        // cudaDeviceSynchronize();

        // double * h_dis = new double[numSequences];
        // cudaMemcpy(h_dis,d_dist,numSequences*sizeof(double),cudaMemcpyDeviceToHost);
        // fprintf(stderr, "%d\n",i);
        // for(int j=0;j<i;j++) std::cerr<<h_dis[j]<<" ";std::cerr<<'\n';

        auto disEnd = std::chrono::high_resolution_clock::now();
        auto treeStart = std::chrono::high_resolution_clock::now();
        blockNum = (numSequences * 4 - 4 + 255) / 256;
        // blockNum = 1024;
        calculateBranchLength(
            i,
            d_head,
            d_nxt,
            d_dist,
            d_e,
            d_len,
            d_belong,
            minPos.data(),
            numSequences * 4 - 4,
            d_closest_dis,
            d_closest_id,
            numSequences);
        auto iter = std::min_element(minPos.begin(), minPos.end(), compare_tuple);
        std::tuple<int, double, double> smallest = *iter;
        /*
        Update Tree (and assign closest nodes to newly added nodes)
        */
        int eid = std::get<0>(smallest);
        double fracLen = std::get<1>(smallest), addLen = std::get<2>(smallest);
        updateTreeStructure(
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
        updateClosestNodes(
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
    std::cerr << "Distance Operation Time " << disTime.count() / 1000000 << " ms\n";
    std::cerr << "Tree Operation Time " << treeTime.count() / 1000000 << " ms\n";
}