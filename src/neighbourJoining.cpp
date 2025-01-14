#ifndef TREE_HPP
#include "tree.hpp"
#endif

#include <utility>
#include <limits>
#include <bits/stdc++.h>

class neighbourJoining {
    public:
        size_t numNodes;
        //size_t totalNumNodes;
        std::vector<float> U;
        std::vector<std::vector<float>> distMatrix;
        std::vector<std::string> nodesIdentifier;
        std::unordered_map<std::string, Node*> nodesInfo;

        void updateUMat();
        std::pair<size_t, size_t> findMinDist();
        void updateNodes(size_t, size_t);
        void updateDistMatrix(size_t, size_t);
};

void neighbourJoining::updateUMat(){
    //std::cout << "Updating UMAT" << std::endl;
    /* Initialize U vector*/
    for (size_t i=0; i<numNodes; i++){
        float Ui = 0.0;
        for (size_t j=0; j<numNodes; j++){
            if (j>i) Ui += distMatrix[i][j];
            else Ui += distMatrix[j][i];
        }
        U[i] = Ui/(numNodes-2);
    }
    //std::cout << "Done" << std::endl;
}

std::pair<size_t, size_t> neighbourJoining::findMinDist(){
    //std::cout << "Finding Minimum Distance" << std::endl;
    size_t mergeIndexI=0, mergeIndexJ=0;
    float minDist = std::numeric_limits<float>::max();

    for (size_t i=0; i<numNodes; i++) {
        for(size_t j=i+1; j<numNodes; j++) {
            float modifiedDist = distMatrix[i][j] - U[i] - U[j];
            if (modifiedDist<minDist) {
                minDist = modifiedDist;
                mergeIndexI = i;
                mergeIndexJ = j;
            }
        }
    }
    //std::cout << "Done" << std::endl;
    return std::make_pair(mergeIndexI, mergeIndexJ);
}

void neighbourJoining::updateNodes(size_t minI, size_t minJ) {
    //std::cout << "Updating Nodes" << std::endl;

    float blI = (((minI<minJ)?distMatrix[minI][minJ]:distMatrix[minJ][minI]) + U[minI] - U[minJ])*0.5;
    float blJ = ((minI<minJ)?distMatrix[minI][minJ]:distMatrix[minJ][minI]) - blI;
    //std::cout << blI <<" "<<blJ<<'\n';
    /* Create New Node*/
    std::string newSeqName = "node_" + std::to_string(numNodes++);
    Node * newNode = new Node(newSeqName, 0.0);
    nodesInfo[newSeqName] = newNode;
    nodesIdentifier.push_back(newSeqName);

    /* Update both nodes*/
    std::string seqName;
    Node* node;
    
    seqName = nodesIdentifier[minI];
    node = nodesInfo[seqName];
    node->branchLength = blI;
    //std::cout<<"######"<<blI<<'\n';
    node->parent = newNode;
    newNode->children.push_back(node);

    seqName = nodesIdentifier[minJ];
    node = nodesInfo[seqName];
    node->branchLength = blJ;
    node->parent = newNode;
    newNode->children.push_back(node);
    // std::cout << "Done" << std::endl;

}


void neighbourJoining::updateDistMatrix(size_t minI, size_t minJ) {
    // std::cout << "Updating Distance Matrix" << std::endl;

    float distIJ = (minI<minJ) ? distMatrix[minI][minJ]:distMatrix[minJ][minI];
    for (size_t i=0; i<numNodes - 1; i++) { /*Additional node is already added but no need to compute distance from it*/
        if (i!=minI && i!=minJ){
            float newdist =( ((i<minI) ? distMatrix[i][minI]:distMatrix[minI][i]) + 
                             ((i<minJ) ? distMatrix[i][minJ]:distMatrix[minJ][i]) -
                            distIJ ) * (0.5);
            distMatrix[i].push_back(newdist);
        }
    }

    distMatrix.push_back(std::vector<float>(numNodes, 0.0));

    for (size_t i=0; i<numNodes; i++) { /*Additional node is already added but that node is not removed*/
        std::vector<float>::iterator iter;
        
        iter = distMatrix[i].begin() + minI;
        distMatrix[i].erase(iter);
        
        /* Check if minJ is smaller or greater than minJ*/
        iter = distMatrix[i].begin() + minJ;
        if (minJ>minI) iter--;
        distMatrix[i].erase(iter);
    }

    /* Remove minI and minJ complete vectors */
    std::vector<std::vector<float>>::iterator iter;
    
    iter = distMatrix.begin() + minI;
    distMatrix.erase(iter); 
    
        /* Check if minJ is smaller or greater than minJ*/
    iter = distMatrix.begin() + minJ;
    if (minJ>minI) iter--;
    distMatrix.erase(iter);

    /* update nodesIdentifier */
    std::vector<std::string>::iterator sIter;
    sIter = nodesIdentifier.begin() + minI;
    nodesIdentifier.erase(sIter);

        /* Check if minJ is smaller or greater than minJ*/
    sIter = nodesIdentifier.begin() + minJ;
    if (minJ>minI) sIter--;
    nodesIdentifier.erase(sIter);

    /*Reduce number of nodes*/
    numNodes--;
    numNodes--;
    // std::cout << "Done" << std::endl;

}
void print(Node* node){
    if(node->children.size()){
        printf("(");
        for(size_t i=0;i<node->children.size();i++){
            print(node->children[i]);
            printf(":");
            printf("%.12f%c",node->children[i]->branchLength,i+1==node->children.size()?')':',');
        }
    }
    else std::cout<<node->identifier;
}
void findNeighbourJoiningTree(std::vector<float>& distMatrix, size_t numTips, std::vector<std::string> distMatrixSeqOrder)
{
    neighbourJoining NJ;
    NJ.numNodes = numTips;
    NJ.distMatrix.resize(numTips);
    for (auto &d: NJ.distMatrix) d.resize(numTips);
    NJ.U.resize(numTips);

    /*Create Nodes for every Tip*/
    for (auto seqName: distMatrixSeqOrder)
    {
        Node* node = new Node(seqName, 0.0);
        //std::cout<<seqName<<" "<<'\n';
        NJ.nodesIdentifier.push_back(seqName);
        NJ.nodesInfo[seqName] = node;
    }

    /* Initialize distance Matrix*/
    size_t index = 0;
    for (size_t i=0; i<numTips; i++){
        for (size_t j=i; j<numTips; j++){
            if (i==j) NJ.distMatrix[i][j] = 0;
            else{
                NJ.distMatrix[i][j] = distMatrix[index];
                index++;
            }
        }
    }

    //std::vector <std::vector<std::pair<string,float>>> treeStructure;
    //treeStructure.resize(NJ.numNodes);
    size_t numNodesToMerge = numTips;
    while ( NJ.numNodes > 1)
    {
        /* Checking sizes after each iteration */
        //std::cout << "Nodes count: " << NJ.numNodes << "\nDistance Matrix size: " << NJ.distMatrix.size() << std::endl;
        //for (size_t i=0; i<NJ.distMatrix.size();i++)  std::cout << "distMatrix["<<i<<"]" << "size: " << NJ.distMatrix[i].size() << std::endl;
    

        NJ.updateUMat();
        std::pair<size_t, size_t> position = NJ.findMinDist();
        size_t minI = position.first, minJ = position.second;
        NJ.updateNodes(minI, minJ);
        NJ.updateDistMatrix(minI, minJ);
        //std::cout << minI << "," << minJ << "\t" <<  NJ.nodesIdentifier[NJ.nodesIdentifier.size() - 1] << std::endl;
    }

    /* if numNodes is not equal to 1, something is wrong */
    assert(NJ.numNodes == 1);

    Node* root = NJ.nodesInfo[NJ.nodesIdentifier[0]];

    //std::cout << root->identifier << std::endl;
    print(root);
    //for(auto &t:NJ.nodesIdentifier) std::cout<<t<<'\n';
}

