/**
 * cpu/divide_and_conquer/placement_close_k_cu.cpp
 *
 * CPU DC k-closest placement implementation. Allocates host DC arrays,
 * findBackboneTreeDC / findClustersDC / findClusterTreeDC (and _batch
 * variants), printTreeDC. Uses TBB over edges/levels; mirrors GPU
 * placement_close_k.cu logic on host.
 */

#include "../mash_placement.hpp"

#include <stdio.h>
#include <zlib.h>
#include <queue>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <chrono>
#include <iostream>
#include <tbb/parallel_for.h>
#include <cassert>
#include <functional>

/** Allocate host DC placement arrays (adjacency, closest_id/dis, cluster copies). */
void MashPlacement::KPlacementDeviceArraysDC::allocateDeviceArraysDC(size_t num, size_t totalNum, int gpuNum){
    // cudaError_t err;

    numSequences = int(num);
    totalNumSequences = int(totalNum);
    bd = 2, idx = 0;
    d_dist = new double [totalNumSequences];
    d_head = new int [totalNumSequences*2];
    d_e = new int [totalNumSequences*8];
    d_len = new double [totalNumSequences*8];
    d_nxt = new int [totalNumSequences*8];
    d_belong = new int [totalNumSequences*8];
    d_closest_dis = new double [totalNumSequences*20];
    d_closest_id = new int [totalNumSequences*20];
    d_closest_dis_cluster = new double [totalNumSequences*20];
    d_closest_id_cluster = new int [totalNumSequences*20];
}

/** Reset adjacency and closest_id/closest_dis to sentinels (TBB over edges). */
void initializeDC(
    int lim,
    int nodes,
    double * d_closest_dis,
    int * d_closest_id,
    int * head,
    int * nxt,
    int * belong,
    int * e
){
    // int tx=threadIdx.x,bs=blockDim.x;
    // int bx=blockIdx.x,gs=gridDim.x;
    // int idx=tx+bs*bx;
    tbb::parallel_for(tbb::blocked_range<int>(0, lim), [&](tbb::blocked_range<int> range){
    for (int idx_ = range.begin(); idx_ < range.end(); ++idx_) {
    // for (int idx_=idx; idx_<lim; idx_+=gs*bs) {
        // if(idx_>=lim) return;
        for(int i=0;i<5;i++){
            d_closest_dis[idx_*5+i]=2;
            d_closest_id[idx_*5+i]=-1;
        }
        nxt[idx_] = -1;
        e[idx_] = -1;
        belong[idx_] = -1;
        if(idx_<nodes) head[idx_] = -1;
    }
    });
}

/** Compare placement tuples by addLen (third element). */
struct compare_tupleDC {
  bool operator()(std::tuple<int,double,double> lhs, std::tuple<int,double,double> rhs)
  {
    return std::get<2>(lhs) < std::get<2>(rhs);
  }
};

/** For each candidate edge, compute (eid, fracLen, addLen) from closest_id/dis (TBB). */
void calculateBranchLengthDC(
    int num, // should be bd, not numSequences 
    int * head,
    int * nxt,
    double * dis, 
    int * e, 
    double * len, 
    int * belong,
    std::tuple<int,double,double> * minPos,
    int lim,
    double * closest_dis,
    int * closest_id
){
    // int tx=threadIdx.x,bs=blockDim.x,bx=blockIdx.x,gs=gridDim.x;
    // int idx=tx+bs*bx;
    // if(idx>=lim) return;
    tbb::parallel_for(tbb::blocked_range<int>(0, lim), [&](tbb::blocked_range<int> range){
    for (int idx_ = range.begin(); idx_ < range.end(); ++idx_) {
    // for (int idx_=idx; idx_<lim; idx_+=gs*bs) {
        if(idx_>=num*4-4||belong[idx_]<e[idx_]){
            std::tuple <int,double,double> minTuple(0,0,2);
            minPos[idx_]=minTuple;
            continue;
        }
        int x=belong[idx_],oth=e[idx_];
        int eid=idx_,otheid;
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
        minPos[idx_]=minTuple;
    }
    });
}


void calculateBranchLengthSpecialIDDC(
    int num, // useless here
    int * head,
    int * nxt,
    double * dis, 
    int * e, 
    double * len, 
    int * belong,
    std::tuple<int,double,double> * minPos,
    int lim,
    double * closest_dis,
    int * closest_id,
    int numToCalculate,
    int * d_edgeMask 
){
    // int tx=threadIdx.x,bs=blockDim.x,bx=blockIdx.x,gs=gridDim.x;
    // int idx=tx+bs*bx;
    // if(idx>=numToCalculate) return;
    tbb::parallel_for(tbb::blocked_range<int>(0, lim), [&](tbb::blocked_range<int> range){
    for (int idx_ = range.begin(); idx_ < range.end(); ++idx_) {
    // for (int idx_=idx; idx_<numToCalculate; idx_+=gs*bs) {
        if(idx_>=numToCalculate) return;
        idx_ = d_edgeMask[idx_];
        if(belong[idx_]<e[idx_]){
            std::tuple <int,double,double> minTuple(0,0,2);
            minPos[idx_]=minTuple;
            continue;
        }
        int x=belong[idx_],oth=e[idx_];
        int eid=idx_,otheid;
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
        minPos[idx_]=minTuple;
    }
    });
}

void updateClosestNodesDC(
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

void updateClosestNodesClusterDC(
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
    int * d_edgeMaskIndex
){
    int l=0,r=-1;
    id[++r]=x,dis[x]=0,from[x]=-1;
    while(l<=r){
        int node=id[l],fb=from[l];
        double d=dis[l];
        l++;
        for(int i=head[node];i!=-1&&d_edgeMaskIndex[head[node]]==head[node];i=nxt[i]){
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

void updateClosestNodesInClusterDC(
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
    int cluster_eid,
    int * belong,
    int * d_edgeMaskIndex
){
    int l=0,r=-1;
    id[++r]=x,dis[x]=0,from[x]=-1;
    int ed1=e[cluster_eid], ed2=belong[cluster_eid];
    while(l<=r){
        int node=id[l],fb=from[l];
        double d=dis[l];
        l++;
        if(node==ed1||node==ed2) continue;
        
        for(int i=head[node];i!=-1;i=nxt[i]){
            // printf("node: %d, head[node]: %d, belong[node]: %d i : %d \n", node, head[node], belong[node], i);
            if (d_edgeMaskIndex[i]!=i) continue;
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

void updateTreeStructureDC(
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
    int totalNumSequences
){
    int middle=placeId+totalNumSequences-1, outside=placeId;
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

void updateTreeStructureInClusterDC(
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
    int totalNumSequences,
    int placeCount // this is the placeCount-th leave
){
    int middle=placeCount+totalNumSequences-1, outside=placeId;
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

void buildInitialTreeDC(
    int totalNumSequences,
    int * head,
    int * e,
    double * len,
    int * nxt,
    int * belong,
    double * dis,
    int edgeCount
){
    int nv = totalNumSequences;
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
    // nv -> 1
    e[edgeCount]=1,len[edgeCount]=d/2,nxt[edgeCount]=head[nv],head[nv]=edgeCount,belong[edgeCount]=nv;
    edgeCount++;
}

void updateClusterInfoDC (
    int leafID,
    int edgeidx,
    int * d_leafMask,
    int * d_edgeMask,
    int * d_edgeMaskIndex,
    int edgeCount,
    int leafCount,
    int * d_leafMap,
    int leaf_idx_in_cluster
){

    d_leafMap[leafCount]=leaf_idx_in_cluster;
    d_leafMask[leafCount++]=leafID;
    for(int i=1;i<=4;i++) {
        d_edgeMask[edgeCount++]=edgeidx-i;
        d_edgeMaskIndex[edgeidx-i]=edgeidx-i;
    }

}

void copyClosestIdsDC(
    int * d_closest_id,
    int * d_closest_id_cluster,
    double * d_closest_dis,
    double * d_closest_dis_cluster,
    int edgeCount
){
    for(int i=0;i<edgeCount;i++){
        for(int j=0;j<5;j++){
            d_closest_id_cluster[i*5+j]=d_closest_id[i*5+j];
            d_closest_dis_cluster[i*5+j]=d_closest_dis[i*5+j];
        }
    }
}

void copyBackClosestIdsDC(
    int * d_closest_id,
    int * d_closest_id_cluster,
    double * d_closest_dis,
    double * d_closest_dis_cluster,
    int edgeCount
){
    for(int i=0;i<edgeCount;i++){
        for(int j=0;j<5;j++){
            d_closest_id[i*5+j]=d_closest_id_cluster[i*5+j];
            d_closest_dis[i*5+j]=d_closest_dis_cluster[i*5+j];
        }
    }
}

void initializeClusterDC (
    int eid,
    int * e,
    int * belong,
    int * head,
    int * nxt,
    int * closest_id,
    int * edgeMask,
    int * leafMask,
    int * edgeMaskIndex,
    int * d_leafMap
){
    int x=belong[eid],y=e[eid];
    int otheid=head[y];
    while(e[otheid]!=x) assert(e[otheid]!=-1),otheid=nxt[otheid];
    
    int leafCount=0;
    for(int i=0;i<5;i++) leafMask[leafCount++]=closest_id[eid*5+i];
    for(int i=0;i<5;i++) d_leafMap[i]=leafMask[i];
    for(int i=0;i<5;i++) leafMask[leafCount++]=closest_id[otheid*5+i];
    for(int i=0;i<5;i++) d_leafMap[i+5]=leafMask[i+5];

    int edgeCount=0;
    edgeMask[edgeCount++]=eid, edgeMask[edgeCount++]=otheid;
    edgeMaskIndex[eid]=eid, edgeMaskIndex[otheid]=otheid;

    // printf("closest ids:\n");
    // for(int i=0;i<5;i++) printf("%d\t", closest_id[624*5+i]);
    // for(int i=0;i<5;i++) printf("%d\t", closest_id[541*5+i]);
    // printf("\n");
    
}



void MashPlacement::KPlacementDeviceArraysDC::deallocateDeviceArraysDC(){
    delete [] d_head;
    delete [] d_e;
    delete [] d_nxt;
    delete [] d_belong;
    delete [] d_closest_id;
    delete [] d_dist;
    delete [] d_len;
    delete [] d_closest_dis;

    // cudaFree(d_head);
    // cudaFree(d_e);
    // cudaFree(d_nxt);
    // cudaFree(d_belong);
    // cudaFree(d_closest_id);
    // cudaFree(d_dist);
    // cudaFree(d_len);
    // cudaFree(d_closest_dis);
}
/** Copy tree to local buffers and emit Newick via recursive DFS. */
void MashPlacement::KPlacementDeviceArraysDC::printTreeDC(std::vector <std::string> name, std::ofstream& output_){
    int * h_head = new int[totalNumSequences*2];
    int * h_e = new int[totalNumSequences*8];
    int * h_nxt = new int[totalNumSequences*8];
    double * h_len = new double[totalNumSequences*8];
    double * h_closest_dis = new double[totalNumSequences*20];
    int * h_closest_id = new int[totalNumSequences*20];
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
    std::memcpy(h_head, d_head, totalNumSequences*2 * sizeof(int));
    std::memcpy(h_e, d_e, totalNumSequences*8*sizeof(int));
    std::memcpy(h_nxt, d_nxt, totalNumSequences*8*sizeof(int));
    std::memcpy(h_len, d_len, totalNumSequences*8*sizeof(double));
    print(totalNumSequences+bd-2,-1);
    output_<<";\n";
}

/*
d_head: index of each node in the arrays ()

*/

void printSeqsDC(
    uint64_t * d_compressedSeqs,
    int num,
    int size
){
    for (int i=490;i<490+num;i++){
        for (int j=0;j<size;j++){
            printf("%ld\n",(d_compressedSeqs[i*size+j]));
        }
    }
}

/** Build backbone tree on first numSequences taxa: init, distances, initial tree,
 * updateClosestNodes, then sequential placement (dist, calculateBranchLength, updateTree, updateClosest). */
void MashPlacement::KPlacementDeviceArraysDC::findBackboneTreeDC(
    Param& params,
    const MashDeviceArraysDC& mashDeviceArrays,
    MatrixReader& matrixReader,
    const MSADeviceArraysDC& msaDeviceArrays,
    const KPlacementDeviceArraysHostDC& kplacementDeviceArraysHost,
    int gpuNum
){ 
    // cudaError_t err;
    if(params.in == "d"){
        matrixReader.distConstructionOnGpu(params, 0, d_dist);
    }
    int * d_id = new int [totalNumSequences*2];
    int * d_from = new int [totalNumSequences*2];
    double * d_dis = new double [totalNumSequences*2];

    /*
    Initialize closest nodes by inifinite
    */
    // int threadNum = 1024, blockNum = (totalNumSequences*4-4+threadNum-1)/threadNum;
    int threadNum = 1024, blockNum = 1024;
    initializeDC (
        totalNumSequences*4-4,
        totalNumSequences*2,
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
        mashDeviceArrays.distConstructionOnGpuForBackboneDC(
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
        msaDeviceArrays.distConstructionOnGpuForBackboneDC(
            params,
            1,
            d_dist
        );
    }
    
    buildInitialTreeDC (
        totalNumSequences,
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
        updateClosestNodesDC (
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
 
    std::cerr << "Finished initial tree construction" << std::endl;
    auto backboneStart = std::chrono::high_resolution_clock::now();
    std::vector <std::tuple<int,double,double>> minPos(numSequences*4-4);


    std::chrono::nanoseconds disTime(0), treeTime(0);
    for(int i=bd;i<numSequences;i++){
        auto disStart = std::chrono::high_resolution_clock::now();
        if(params.in == "r"){
            mashDeviceArrays.distRangeConstructionOnGpuDC(
                params,
                i,
                d_dist,
                0,
                i-1
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
            msaDeviceArrays.distRangeConstructionOnGpuDC(
                params,
                i,
                d_dist,
                0,
                i-1
            );
        }

        auto disEnd = std::chrono::high_resolution_clock::now();
        auto treeStart = std::chrono::high_resolution_clock::now();
        calculateBranchLengthDC (
            i,
            d_head,
            d_nxt,
            d_dist,
            d_e,
            d_len,
            d_belong,
            minPos.data(),
            numSequences*4-4,
            d_closest_dis,
            d_closest_id
        );

        auto iter=std::min_element(minPos.begin(),minPos.begin()+numSequences*4-4,compare_tupleDC());
        std::tuple<int,double,double> smallest=*iter;
        /*
        Update Tree (and assign closest nodes to newly added nodes)
        */
        int eid=std::get<0>(smallest);
        double fracLen=std::get<1>(smallest),addLen=std::get<2>(smallest);
        updateTreeStructureDC (
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
            totalNumSequences
        );
        idx+=4;

        /*
        Update closest nodes
        */
        updateClosestNodesDC (
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
        auto treeEnd = std::chrono::high_resolution_clock::now();
        disTime += disEnd - disStart;
        treeTime += treeEnd - treeStart;
    }
    // cudaDeviceSynchronize();
    auto backboneEnd = std::chrono::high_resolution_clock::now();
    auto backboneTime = backboneEnd - backboneStart;
    std::cerr << "Finished backbone construction in: "<< backboneTime.count()/1000000 << " ms\n";
    return;
}

/** Assign each non-backbone sequence to k-closest edge (clusterID); batch dist + calculateBranchLength. */
void MashPlacement::KPlacementDeviceArraysDC::findClustersDC(
    Param& params,
    const MashDeviceArraysDC& mashDeviceArrays,
    MatrixReader& matrixReader,
    const MSADeviceArraysDC& msaDeviceArrays,
    KPlacementDeviceArraysHostDC& kplacementDeviceArraysHost
){ 
    // cudaError err;
    int idx=params.backboneSize*4-4;
    
    kplacementDeviceArraysHost.clusterID = new int[totalNumSequences];
    std::vector <std::tuple<int,double,double>> minPos(totalNumSequences*4-4);
    int threadNum = 1024, blockNum = 1024;
    auto clusterStart = std::chrono::high_resolution_clock::now();
    uint64_t localBatchSize = params.batchSize;
    for(int i=numSequences;i<totalNumSequences;i+=localBatchSize){
        if (i+localBatchSize>totalNumSequences) localBatchSize=totalNumSequences-i;
        // std::cerr << "Processing batch: "<< i << " to " << i+localBatchSize << "\n";

        /* copy data to d_hashListConst or d_compressedSeqsConst */
        if (params.in == "r") 
            std::memcpy(mashDeviceArrays.d_hashListConst, mashDeviceArrays.h_hashList+i*params.sketchSize, localBatchSize*params.sketchSize*sizeof(uint64_t));
            // err = cudaMemcpy(mashDeviceArrays.d_hashListConst, mashDeviceArrays.h_hashList+i*params.sketchSize, localBatchSize*params.sketchSize*sizeof(uint64_t),cudaMemcpyHostToDevice);
        else if (params.in == "m"){
            size_t maxLengthCompressed = (msaDeviceArrays.d_seqLen + 15) / 16;
            std::memcpy(msaDeviceArrays.d_compressedSeqsConst, msaDeviceArrays.h_compressedSeqs+i*maxLengthCompressed, 1ll*localBatchSize*maxLengthCompressed*sizeof(uint64_t));
            // err = cudaMemcpy(msaDeviceArrays.d_compressedSeqsConst, msaDeviceArrays.h_compressedSeqs+i*maxLengthCompressed, 1ll*localBatchSize*maxLengthCompressed*sizeof(uint64_t),cudaMemcpyHostToDevice);
        }
        else {
            std::cerr << "Error: Input type must be unaligned or aligned for clustering based approach\n";
            exit(1);
        }
        
        for (int j=i;j<i+localBatchSize;j++){
            auto disStart = std::chrono::high_resolution_clock::now();
            if(params.in == "r"){
                mashDeviceArrays.distRangeConstructionOnGpuDC(
                    params,
                    j-i,
                    d_dist,
                    0,
                    numSequences-1,
                    true
                );
            }
            else if(params.in == "d"){
                matrixReader.distConstructionOnGpu(
                    params,
                    j,
                    d_dist
                );
            }
            else if(params.in == "m"){
                msaDeviceArrays.distRangeConstructionOnGpuDC(
                    params,
                    j-i,
                    d_dist,
                    0,
                    numSequences-1,
                    true
                );
            }

            // double * h_dis = new double[numSequences];
            // cudaMemcpy(h_dis,d_dist,numSequences*sizeof(double),cudaMemcpyDeviceToHost);
            // fprintf(stderr, "%d\n",i);
            // for(int j=0;j<i;j++) std::cerr<<h_dis[j]<<" ";std::cerr<<'\n';


            auto disEnd = std::chrono::high_resolution_clock::now();
            auto treeStart = std::chrono::high_resolution_clock::now();
            
            calculateBranchLengthDC (
                j,
                d_head,
                d_nxt,
                d_dist,
                d_e,
                d_len,
                d_belong,
                minPos.data(),
                numSequences*4-4,
                d_closest_dis,
                d_closest_id
            );
            auto iter=std::min_element(minPos.begin(),minPos.begin()+numSequences*4-4,compare_tupleDC());
            std::tuple<int,double,double> smallest=*iter;
            kplacementDeviceArraysHost.clusterID[j] = std::get<0>(smallest);
            // std::cerr << "Sequence " << j << " assigned to cluster " << kplacementDeviceArraysHost.clusterID[j] << "\n";
        }
    }
    auto clusterEnd = std::chrono::high_resolution_clock::now();
    auto clusterTime = clusterEnd - clusterStart;

    std::cerr << "Finished clustering in: "<< clusterTime.count()/1000000 << " ms\n";

    /* Copy data from device to host */
    // err = cudaMemcpy(kplacementDeviceArraysHost.h_dist, d_dist, totalNumSequences*sizeof(double),cudaMemcpyDeviceToHost);
    // if (err != cudaSuccess)
    // {
    //     fprintf(stderr, "Gpu_ERROR: d_dist cudaMemcpy failed!\n");
    //     exit(1);
    // }

    // err = cudaMemcpy(kplacementDeviceArraysHost.h_head, d_head, totalNumSequences*2*sizeof(int),cudaMemcpyDeviceToHost);
    // if (err != cudaSuccess)
    // {
    //     fprintf(stderr, "Gpu_ERROR: d_head cudaMemcpy failed!\n");
    //     exit(1);
    // }

    // err = cudaMemcpy(kplacementDeviceArraysHost.h_e, d_e, totalNumSequences*8*sizeof(int),cudaMemcpyDeviceToHost);
    // if (err != cudaSuccess)
    // {
    //     fprintf(stderr, "Gpu_ERROR: d_e cudaMemcpy failed!\n");
    //     exit(1);
    // }

    // err = cudaMemcpy(kplacementDeviceArraysHost.h_len, d_len, totalNumSequences*8*sizeof(double),cudaMemcpyDeviceToHost);
    // if (err != cudaSuccess)
    // {
    //     fprintf(stderr, "Gpu_ERROR: d_len cudaMemcpy failed!\n");
    //     exit(1);
    // }
    // err = cudaMemcpy(kplacementDeviceArraysHost.h_nxt, d_nxt, totalNumSequences*8*sizeof(int),cudaMemcpyDeviceToHost);
    // if (err != cudaSuccess)
    // {
    //     fprintf(stderr, "Gpu_ERROR: d_nxt cudaMemcpy failed!\n");
    //     exit(1);
    // }

    // err = cudaMemcpy(kplacementDeviceArraysHost.h_belong, d_belong, totalNumSequences*8*sizeof(int),cudaMemcpyDeviceToHost);
    // if (err != cudaSuccess)
    // {
    //     fprintf(stderr, "Gpu_ERROR: d_belong cudaMemcpy failed!\n");
    //     exit(1);
    // }

    // err = cudaMemcpy(kplacementDeviceArraysHost.h_closest_dis, d_closest_dis, totalNumSequences*20*sizeof(double),cudaMemcpyDeviceToHost);
    // if (err != cudaSuccess)
    // {
    //     fprintf(stderr, "Gpu_ERROR: d_closest_dis cudaMemcpy failed!\n");
    //     exit(1);
    // }

    // err = cudaMemcpy(kplacementDeviceArraysHost.h_closest_id, d_closest_id, totalNumSequences*20*sizeof(int),cudaMemcpyDeviceToHost);
    // if (err != cudaSuccess)
    // {
    //     fprintf(stderr, "Gpu_ERROR: d_closest_id cudaMemcpy failed!\n");
    //     exit(1);
    // }

    // err = cudaMemcpy(kplacementDeviceArraysHost.h_closest_dis_cluster, d_closest_dis_cluster, totalNumSequences*20*sizeof(double),cudaMemcpyDeviceToHost);
    // if (err != cudaSuccess)
    // {
    //     fprintf(stderr, "Gpu_ERROR: d_closest_dis_cluster cudaMemcpy failed!\n");
    //     exit(1);
    // }

    // err = cudaMemcpy(kplacementDeviceArraysHost.h_closest_id_cluster, d_closest_id_cluster, totalNumSequences*20*sizeof(int),cudaMemcpyDeviceToHost);
    // if (err != cudaSuccess)
    // {
    //     fprintf(stderr, "Gpu_ERROR: d_closest_id_cluster cudaMemcpy failed!\n");
    //     exit(1);
    // }

    // cudaDeviceSynchronize();

    std::cerr << "Finished data transfer\n";
    return;
}

/** Batch variant of findClustersDC for a single clustering batch index. */
void MashPlacement::KPlacementDeviceArraysDC::findClustersDC_batch(
    Param& params,
    const MashDeviceArraysDC& mashDeviceArrays,
    MatrixReader& matrixReader,
    const MSADeviceArraysDC& msaDeviceArrays,
    KPlacementDeviceArraysHostDC& kplacementDeviceArraysHost,
    const int clusteringBatchIdx
){ 
    // cudaError err;
    std::vector <std::tuple<int,double,double>> minPos(totalNumSequences*4-4);
    int threadNum = 1024, blockNum = 1024;

    auto clusterStart = std::chrono::high_resolution_clock::now();

    for(int i=0;i<params.batchSize && i<params.totalNumSeqs-clusteringBatchIdx*params.batchSize;i++){
        
        if(params.in == "r"){
            mashDeviceArrays.distRangeConstructionOnGpuDC(
                params,
                i,
                d_dist,
                0,
                numSequences-1,
                true
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
            msaDeviceArrays.distRangeConstructionOnGpuDC(
                params,
                i,
                d_dist,
                0,
                numSequences-1,
                true
            );
        }


        auto disEnd = std::chrono::high_resolution_clock::now();
        auto treeStart = std::chrono::high_resolution_clock::now();
        
        calculateBranchLengthDC (
            clusteringBatchIdx*params.batchSize + i,
            d_head,
            d_nxt,
            d_dist,
            d_e,
            d_len,
            d_belong,
            minPos.data(),
            numSequences*4-4,
            d_closest_dis,
            d_closest_id
        );
        auto iter=std::min_element(minPos.begin(),minPos.begin()+numSequences*4-4,compare_tupleDC());
        std::tuple<int,double,double> smallest=*iter;
        kplacementDeviceArraysHost.clusterID[clusteringBatchIdx*params.batchSize + i] = std::get<0>(smallest);

    }
    // cudaDeviceSynchronize();

    return;
}

void rearrangeHashListInClusterDC(
    int numSequences,
    int sketchSize,
    uint64_t * original,
    uint64_t * target
){
    // int tx = threadIdx.x, bx = blockIdx.x;
    // int bs = blockDim.x;
    // int idx = tx+bs*bx;
    tbb::parallel_for(tbb::blocked_range<int>(0, numSequences), [&](tbb::blocked_range<int> range){ 
    for (int idx_ = range.begin(); idx_ < range.end(); ++idx_) {
    // if(idx>=numSequences) return;
    // for (int idx_=idx; idx_<numSequences; idx_+=bs*gridDim.x){
        // if (idx_ >= numSequences) return;
        for(int i=0;i<sketchSize;i++){
            target[i*numSequences+idx_] = original[idx_*sketchSize+i];
        }
    }
    });
}

void transferMashClusterInfoDC(
    MashPlacement::MashDeviceArraysDC& mashDeviceArrays,
    std::vector<int> leafList,
    MashPlacement::Param& params
){
    uint64_t * hashListLocal = new uint64_t[params.backboneSize*params.sketchSize];
    uint64_t * hashList = mashDeviceArrays.h_hashList;
    int l=0;
    for (auto &leaf: leafList){
        std::memcpy(hashListLocal+l*params.sketchSize, hashList+leaf*params.sketchSize, params.sketchSize*sizeof(uint64_t));
        l++;
    }
    std::memcpy(mashDeviceArrays.d_hashListConst, hashListLocal, params.backboneSize*params.sketchSize*sizeof(uint64_t));

    /* Rearrange only for backbone tree */
    uint64_t * temp_hashList = new uint64_t [params.sketchSize*params.backboneSize];
    
    int threadsPerBlock = 1024;
    int blocksPerGrid = 1024;
    rearrangeHashListInClusterDC (
        params.backboneSize,
        int(params.sketchSize),
        mashDeviceArrays.d_hashListConst,
        temp_hashList
    );
    std::swap(mashDeviceArrays.d_hashListConst, temp_hashList);
    delete [] temp_hashList;
    delete [] hashListLocal;
    return;
}


void transferMsaClusterInfoDC(
    MashPlacement::MSADeviceArraysDC& mashDeviceArrays,
    std::vector<int> leafList,
    MashPlacement::Param& params
){
    size_t maxLengthCompressed = (mashDeviceArrays.d_seqLen + 15) / 16;
    uint64_t * compressedSeqs_local = new uint64_t[params.backboneSize*maxLengthCompressed];
    uint64_t * compressedSeqs = mashDeviceArrays.h_compressedSeqs;

    int l=0;
    for (auto &leaf: leafList){
        memcpy(compressedSeqs_local+1ll*l*maxLengthCompressed, compressedSeqs+1ll*leaf*maxLengthCompressed, 1ll*maxLengthCompressed*sizeof(uint64_t));
        l++;
    }
    std::memcpy(mashDeviceArrays.d_compressedSeqsConst, compressedSeqs_local, 1ll*params.backboneSize*maxLengthCompressed*sizeof(uint64_t));
    delete[] compressedSeqs_local;
    return;
}




void print_d_hashListConstDC(
    uint64_t * d_hashListConst,
    int sketchSize,
    int batchSize
){
    for (int i=0;i<batchSize;i++){
        for (int j=0;j<sketchSize;j++){
            printf("%llu ", d_hashListConst[i*sketchSize+j]);
        }
        printf("\n");
    }
}

void resetEdgeMaskIndexDC(int * d_edgeMaskIndex, int size){
    // int tx = threadIdx.x, bx = blockIdx.x;
    // int bs = blockDim.x, gs = gridDim.x;
    // int idx = tx+bs*bx;
    // if(idx>=size) return;
    tbb::parallel_for(tbb::blocked_range<int>(0, size), [&](tbb::blocked_range<int> range){ 
    for (int idx_ = range.begin(); idx_ < range.end(); ++idx_) {
    // for (int idx_=idx; idx_<size; idx_+=bs*gs){
        // if (idx_ >= size) return;
        d_edgeMaskIndex[idx_] = -1;
    }
    });
    // d_edgeMaskIndex[idx] = -1;
    
}

void resetIdFromDisDC(int * id, int * from, double * dis, int size){
    // int tx = threadIdx.x, bx = blockIdx.x;
    // int bs = blockDim.x;
    // int idx = tx+bs*bx;
    // if(idx>=size) return;
    tbb::parallel_for(tbb::blocked_range<int>(0, size), [&](tbb::blocked_range<int> range){ 
    for (int idx_ = range.begin(); idx_ < range.end(); ++idx_) {
        id[idx_] = -1;
        from[idx_] = -1;
        dis[idx_] = 2.0;   
    }
    });
    
}

/** For each cluster: transfer MASH/MSA, init, place each leaf; emit Newick. */
void MashPlacement::KPlacementDeviceArraysDC::findClusterTreeDC(
    Param& params,
    MashDeviceArraysDC& mashDeviceArrays,
    MatrixReader& matrixReader,
    MSADeviceArraysDC& msaDeviceArrays,
    KPlacementDeviceArraysHostDC& kplacementDeviceArraysHost
){ 
    int idx=params.backboneSize*4-4;
    int threadNum = 1024, blockNum = 1024;
    int * d_id;
    int * d_from;
    double * d_dis;
    d_id = new int [totalNumSequences*2];
    d_from = new int [totalNumSequences*2];
    d_dis = new double [totalNumSequences*2];
    
    auto * cluster_id = kplacementDeviceArraysHost.clusterID;
    std::vector<std::vector <int>> contains(numSequences*4-4);

    for(int i=numSequences;i<totalNumSequences;i++) contains[cluster_id[i]].push_back(i);
    
    int * d_edgeMask;
    d_edgeMask = new int [totalNumSequences*4];

    int * d_edgeMaskIndex;
    d_edgeMaskIndex = new int [totalNumSequences*4];
    
    int * d_leafMask;
    d_leafMask = new int [totalNumSequences*2];

    int * d_leafMap;
    d_leafMap = new int [params.backboneSize];

    int insertLeafCount=numSequences;
    std::vector <std::tuple<int,double,double>> minPos(totalNumSequences*4-4);
    
    std::vector<int> leafList (params.backboneSize);
    std::unordered_map <int,int> h_leafMap;
    
    // for(int i=0;i<numSequences*4-4;i++){ 
    int actual_placement=0;
    int i=0;
    while (i<numSequences*4-4) {
        std::cerr << "Processing batch "<< i << " out of " << numSequences*4-4 <<  std::endl;
        h_leafMap.clear();
        int startClusterID = i;
        int batchSize=0;
        if (params.backboneSize<contains[i].size()) {
            std::cerr << "Cluster " << i << " size (" << contains[i].size() <<") is larger than backbone size\n";
            exit(1);
        }
        while (batchSize<params.backboneSize && i<numSequences*4-4) {
            if (contains[i].size() + batchSize >= params.backboneSize) break;
            for (auto &leaf: contains[i]) { 
                h_leafMap[leaf] = batchSize;
                leafList[batchSize++] = leaf; 
            }
            i++;
        }

        if(params.in == "r")
            transferMashClusterInfoDC(mashDeviceArrays, leafList, params);
        else if (params.in == "m")
            transferMsaClusterInfoDC(msaDeviceArrays, leafList, params);
        else {
            std::cerr << "Error: Input type must be unaligned or aligned for clustering based approach\n";
            exit(1);
        }

        int localCount_=0;
        for (int j=startClusterID;j<i;j++) {

            if (contains[j].size() == 0) continue;
            resetEdgeMaskIndexDC (d_edgeMaskIndex, totalNumSequences*4);
            initializeClusterDC (
                j,
                d_e,
                d_belong,
                d_head,
                d_nxt,
                d_closest_id,
                d_edgeMask,
                d_leafMask,
                d_edgeMaskIndex,
                d_leafMap
            );
            int edgeCount=2, leafCount=10;
            for(auto &leaf:contains[j]) {
                actual_placement++;
                // std::cerr << "Processing leaf: "<< leaf << std::endl;
                // copy d_leafMask to host and print
                // if (j == 544) {

                //     int * h_leafMask = new int[totalNumSequences*2];
                //     err = cudaMemcpy(h_leafMask, d_leafMask, totalNumSequences*2*sizeof(int),cudaMemcpyDeviceToHost);
                //     if (err != cudaSuccess)
                //     {
                //         fprintf(stderr, "Gpu_ERROR: d_leafMask cudaMemcpy failed!\n");
                //         exit(1);
                //     }
                //     for (int z=0;z<leafCount;z++){
                //         std::cerr << h_leafMask[z] << "\t";
                //     }
                //     std::cerr << "\n";

                //     int * h_edgeMask = new int[totalNumSequences*4];
                //     err = cudaMemcpy(h_edgeMask, d_edgeMask, totalNumSequences*4*sizeof(int),cudaMemcpyDeviceToHost);
                //     if (err != cudaSuccess)
                //     {
                //         fprintf(stderr, "Gpu_ERROR: d_edgeMask cudaMemcpy failed!\n");
                //         exit(1);
                //     }
                //     for (int z=0;z<edgeCount;z++){
                //         std::cerr << h_edgeMask[z] << "\t";
                //     }
                //     std::cerr << "\n";

                //     int * h_edgeMaskIndex = new int[totalNumSequences*4];
                //     err = cudaMemcpy(h_edgeMaskIndex, d_edgeMaskIndex, totalNumSequences*4*sizeof(int),cudaMemcpyDeviceToHost);
                //     if (err != cudaSuccess)
                //     {
                //         fprintf(stderr, "Gpu_ERROR: d_edgeMaskIndex cudaMemcpy failed!\n");
                //         exit(1);
                //     }
                //     for (int z=0;z<edgeCount;z++){
                //         std::cerr << h_edgeMaskIndex[h_edgeMask[z]] << "\t";
                //     }
                //     std::cerr << "\n";
                // }



                int leaf_idx_in_cluster = h_leafMap[leaf];
                if(params.in == "r"){
                    mashDeviceArrays.distSpecialIDConstructionOnGpuDC(
                        params,
                        localCount_,
                        d_dist,
                        leafCount,
                        d_leafMask,
                        d_leafMap
                    );

                } else if(params.in == "m"){
                    msaDeviceArrays.distSpecialIDConstructionOnGpuDC(
                        params,
                        localCount_,
                        d_dist,
                        leafCount,
                        d_leafMask,
                        d_leafMap
                    );
                } 
                localCount_++;    

                calculateBranchLengthSpecialIDDC (
                    j,
                    d_head,
                    d_nxt,
                    d_dist,
                    d_e,
                    d_len,
                    d_belong,
                    minPos.data(),
                    numSequences*4-4,
                    d_closest_dis,
                    d_closest_id,
                    edgeCount,
                    d_edgeMask
                );
                auto iter=std::min_element(minPos.begin(),minPos.begin()+edgeCount,compare_tupleDC());
                std::tuple<int,double,double> smallest=*iter;

                /*
                Update Tree Structure
                */

                int eid=std::get<0>(smallest);
                double fracLen=std::get<1>(smallest),addLen=std::get<2>(smallest);
                // std::cerr << "eid: " << eid << " Cluster ID: " << j << " dist " << addLen << " dist2 " << fracLen << "\n";
                updateTreeStructureInClusterDC (
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
                    leaf,
                    idx,
                    totalNumSequences,
                    insertLeafCount
                );
                idx+=4, insertLeafCount++;

                /*
                Update edgeMask and leafMask
                */

                updateClusterInfoDC (
                    leaf,
                    idx,
                    d_leafMask,
                    d_edgeMask,
                    d_edgeMaskIndex,
                    edgeCount,
                    leafCount,
                    d_leafMap,
                    leaf_idx_in_cluster
                );

                /*
                Update closest nodes
                */

                updateClosestNodesInClusterDC (
                    d_head,
                    d_nxt,
                    d_e,
                    d_len,
                    d_closest_dis,
                    d_closest_id,
                    leaf,
                    d_id,
                    d_from,
                    d_dis,
                    j,
                    d_belong,
                    d_edgeMaskIndex
                );

                edgeCount+=4, leafCount++;
                // std::cerr << "leaf: " << leaf << " Cluster ID: " << j << " eid " << eid << " dist " << addLen << " dist2 " << fracLen << "\n";
                // if (leaf == 879) exit(0);
                // exit(0);
            }
        }
        // exit(0);
        assert(localCount_==batchSize);
        
    }
    // cudaDeviceSynchronize();

}

bool read_line(gzFile file, std::string& line, int maxLength) {
    char buffer[maxLength];
    if (gzgets(file, buffer, maxLength) == Z_NULL) {
        return false;
    }
    line = buffer;
    if (!line.empty() && line.back() == '\n') line.pop_back();
    return true;
}

bool read_binary_blob(gzFile file, size_t n, uint64_t* out_data) {
    size_t total_bytes = n * sizeof(uint64_t);
    char* raw_ptr = reinterpret_cast<char*>(out_data);

    size_t bytes_read = 0;
    while (bytes_read < total_bytes) {
        int chunk = gzread(file, raw_ptr + bytes_read, total_bytes - bytes_read);
        if (chunk <= 0) {
            std::cerr << "Error: Unexpected end of file while reading binary data." << std::endl;
            delete[] out_data;
            out_data = nullptr;
            return false;
        }
        bytes_read += chunk;
    }
    return true;
}

/** Batch variant of findClusterTreeDC (cluster data from dir, isCluster filter). */
void MashPlacement::KPlacementDeviceArraysDC::findClusterTreeDC_batch(
    Param& params,
    MashDeviceArraysDC& mashDeviceArrays,
    MatrixReader& matrixReader,
    MSADeviceArraysDC& msaDeviceArrays,
    KPlacementDeviceArraysHostDC& kplacementDeviceArraysHost,
    const std::string dir,
    std::vector<bool>& isCluster
){ 
    int idx=params.backboneSize*4-4;
    int threadNum = 1024, blockNum = 1024;
    int * d_id = new int [totalNumSequences*2];
    int * d_from = new int [totalNumSequences*2];
    double * d_dis = new double [totalNumSequences*2];
    auto * cluster_id = kplacementDeviceArraysHost.clusterID;
    std::vector<std::vector <int>> contains(numSequences*4-4);

    for(int i=numSequences;i<totalNumSequences;i++) contains[cluster_id[i]].push_back(i);
    
   
    int * d_edgeMask = new int [totalNumSequences*4];
    int * d_edgeMaskIndex = new int [totalNumSequences*4];
    int * d_leafMask = new int [totalNumSequences*2];
    int * d_leafMap = new int [params.backboneSize];
    
    int insertLeafCount=numSequences;
    std::vector <std::tuple<int,double,double>> minPos(totalNumSequences*4-4);
    std::vector<int> leafList (params.backboneSize);
    std::unordered_map <int,int> h_leafMap;
    
    int test_var = 0;
    int actual_placement = 0;
    int clustersPerBatchFile = 1000;
    int clusterFiles = (numSequences*4-4)/clustersPerBatchFile + ((numSequences*4-4)%clustersPerBatchFile != 0);
    
    for (int bc=0; bc<clusterFiles;bc++){
        /* Find number of sequences in the file*/
        int seqInFile=0;
        std::unordered_map<int, std::vector<int>> clusterToCompressSeqIdxMap;
        std::unordered_map<int, int> localIdxToOriginalIdxMap;
        for (int i=bc*clustersPerBatchFile;i<(bc+1)*clustersPerBatchFile && i<numSequences*4-4;i++){
            if (!isCluster[i]) continue;
            clusterToCompressSeqIdxMap[i] = std::vector<int>();
            seqInFile+=contains[i].size();
        }

        if (seqInFile == 0) continue;


        /* Read all sequences in the file*/
        std::string path = dir + "/" + std::to_string(bc) + ".gz";
        std::cerr << "Reading file: " << path << " with " << seqInFile << " sequences\n";
        size_t maxLengthCompressed = (msaDeviceArrays.d_seqLen + 15) / 16;
        uint64_t * compressedSeqs_local = new uint64_t[seqInFile*maxLengthCompressed];

        gzFile file = gzopen(path.c_str(), "rb");
        if (!file) {
            std::cerr << "gzopen failed: " << path << " " << strerror(errno) << "\n";
            exit(0);
        }
        std::string line;
        uint64_t counter=0;
        while (read_line(file, line, maxLengthCompressed)) {
            if (line.empty() || line[0] != '>') continue;

            std::string name, id;
            std::istringstream ss(line.substr(1));
            std::getline(ss, name, '\t');
            std::getline(ss, id, '\t');
            localIdxToOriginalIdxMap[counter] = std::stoi(id);

            if (!read_binary_blob(file, maxLengthCompressed, compressedSeqs_local+1ll*counter*maxLengthCompressed)) {
                std::cerr << "Failed to read binary data for sequence: " << name << std::endl;
                break;
            }

            if (clusterToCompressSeqIdxMap.find(cluster_id[std::stoi(id)]) != clusterToCompressSeqIdxMap.end()) {
                clusterToCompressSeqIdxMap[cluster_id[std::stoi(id)]].push_back(counter);
            }
            counter++;
        }
        gzclose(file);


        for(int i=bc*clustersPerBatchFile;i<(bc+1)*clustersPerBatchFile && i<numSequences*4-4;i++){
            if (!isCluster[i]) continue;
            std::cerr << "Handling cluster " << i << " with size " << contains[i].size() << "\n";
            
            // open file
            uint64_t * compressedSeqs_local_per_cluster = new uint64_t[params.backboneSize*maxLengthCompressed];

            int counter=0;
            for (auto &compressIdx: clusterToCompressSeqIdxMap[i]){
                if (counter + 1> numSequences) {
                    std::cerr << "Cluster " << i << " size is larger than backbone size\n";
                    exit(0);
                }
                
                memcpy(compressedSeqs_local_per_cluster+1ll*counter*maxLengthCompressed, compressedSeqs_local+1ll*compressIdx*maxLengthCompressed, 1ll*maxLengthCompressed*sizeof(uint64_t));
                
                h_leafMap[localIdxToOriginalIdxMap[compressIdx]] = counter;
                leafList[counter] = localIdxToOriginalIdxMap[compressIdx];
                counter++;
                
            }

            if(params.in == "m") {
                memcpy(msaDeviceArrays.d_compressedSeqsConst, compressedSeqs_local_per_cluster, 1ll*params.backboneSize*maxLengthCompressed*sizeof(uint64_t));
                // auto err = cudaMemcpy(msaDeviceArrays.d_compressedSeqsConst, compressedSeqs_local_per_cluster, 1ll*params.backboneSize*maxLengthCompressed*sizeof(uint64_t),cudaMemcpyHostToDevice);
            }
            

            if (contains[i].size() == 0) continue;
            int localCount_=0;

            resetEdgeMaskIndexDC (d_edgeMaskIndex, totalNumSequences*4);
            initializeClusterDC (
                i,
                d_e,
                d_belong,
                d_head,
                d_nxt,
                d_closest_id,
                d_edgeMask,
                d_leafMask,
                d_edgeMaskIndex,
                d_leafMap
            );
            int edgeCount=2, leafCount=10;
            for(auto &leaf:contains[i]) {
                actual_placement++;
                int leaf_idx_in_cluster = h_leafMap[leaf];
                if(params.in == "m"){
                    msaDeviceArrays.distSpecialIDConstructionOnGpuDC(
                        params,
                        localCount_,
                        d_dist,
                        leafCount,
                        d_leafMask,
                        d_leafMap
                    );
                } 
                localCount_++;

                calculateBranchLengthSpecialIDDC (
                    i,
                    d_head,
                    d_nxt,
                    d_dist,
                    d_e,
                    d_len,
                    d_belong,
                    minPos.data(),
                    numSequences*4-4,
                    d_closest_dis,
                    d_closest_id,
                    edgeCount,
                    d_edgeMask
                );
                auto iter=std::min_element(minPos.begin(),minPos.begin()+edgeCount,compare_tupleDC());
                std::tuple<int,double,double> smallest=*iter;

                /*
                Update Tree Structure
                */

                int eid=std::get<0>(smallest);
                double fracLen=std::get<1>(smallest),addLen=std::get<2>(smallest);
                // std::cerr << "eid: " << eid << " Cluster ID: " << j << " dist " << addLen << " dist2 " << fracLen << "\n";
                updateTreeStructureInClusterDC (
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
                    leaf,
                    idx,
                    totalNumSequences,
                    insertLeafCount
                );
                idx+=4, insertLeafCount++;

                /*
                Update edgeMask and leafMask
                */

                updateClusterInfoDC (
                    leaf,
                    idx,
                    d_leafMask,
                    d_edgeMask,
                    d_edgeMaskIndex,
                    edgeCount,
                    leafCount,
                    d_leafMap,
                    leaf_idx_in_cluster
                );

                /*
                Update closest nodes
                */

                updateClosestNodesInClusterDC (
                    d_head,
                    d_nxt,
                    d_e,
                    d_len,
                    d_closest_dis,
                    d_closest_id,
                    leaf,
                    d_id,
                    d_from,
                    d_dis,
                    i,
                    d_belong,
                    d_edgeMaskIndex
                );

                edgeCount+=4, leafCount++;
                // std::cerr << "leaf: " << leaf << " Cluster ID: " << j << " eid " << eid << " dist " << addLen << " dist2 " << fracLen << "\n";
                // if (leaf == 879) exit(0);
                // exit(0);
            }
            delete[] compressedSeqs_local_per_cluster;

        }

        delete[] compressedSeqs_local;
    }
    // cudaDeviceSynchronize();


    return;
}

