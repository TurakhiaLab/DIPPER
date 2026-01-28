/**
 * tree.hpp
 *
 * In-memory phylogenetic tree: Node (name, branch length, parent, children, level)
 * and Tree (root, allNodes, Newick parse/serialization, DFS, leaf/internal counts).
 */

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
    Node(const std::string& id, double len);       /** Root or standalone node. */
    Node(const std::string& id, Node* par, double len);  /** Child of par. */
    size_t getNumLeaves();
    size_t getNumNodes();
    bool is_leaf() {return !(name.substr(0,4) == "node");}
    bool is_leaf_new() {return children.size() == 0;}
    void setNumleaves() {numLeaves = getNumLeaves();};

    std::string name;
    int idx;
    double bl;
    size_t level;
    int lca=0;

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
    Tree(std::string newick, size_t totalLeaves=0);  /** Parse Newick and build tree. */
    Tree(Node* node);
    Tree() {root = nullptr;}
    ~Tree();
    std::string getNewickString(Node* node = nullptr);
    void updateNewick();
    void dfsExpansion(Node* node, std::vector< Node* >& vec);
};