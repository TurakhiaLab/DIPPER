#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_map>
#include <algorithm>

class Node 
{
public:
    Node(const std::string& id, double len);
    Node(const std::string& id, Node* par, double len);
    size_t getNumLeaves();
    size_t getNumNodes();
    bool is_leaf() {return !(name.substr(0,4) == "node");}
    void setNumleaves() {numLeaves = getNumLeaves();};

    std::string name;
    int idx;
    double bl;
    size_t level;

    Node* parent;
    std::vector< Node* > children;
    size_t numLeaves = {0};
    float weight = {0};
};

class Tree
{
public:

    size_t m_currInternalNode{ 0 };
    size_t m_maxDepth{ 0 };
    size_t m_numLeaves{ 0 };
    size_t m_numLeafID{ 0 };
    float m_meanDepth{ 0 };
    std::string newInternalNodeId() { return "node_" + std::to_string(++m_currInternalNode);}

    Node* root;
    std::unordered_map< std::string, Node* > allNodes;
    Tree(std::string newick, size_t totalLeaves);
    Tree(Node* node);
    Tree() {root = nullptr;}
    ~Tree();
    std::string getNewickString(Node* node = nullptr);
    void dfsExpansion(Node* node, std::vector< Node* >& vec);
};