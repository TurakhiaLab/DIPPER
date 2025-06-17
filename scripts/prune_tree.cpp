#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <random>
#include <sstream>
#include <algorithm>

using namespace std;

struct TreeNode {
    string name;
    vector<TreeNode*> children;
    TreeNode* parent = nullptr;
};

// Simple Newick parser (assumes valid Newick without branch lengths)
TreeNode* parseNewick(const string& input) {
    vector<TreeNode*> stack;
    TreeNode* current = new TreeNode();
    TreeNode* root = current;

    string token;
    for (char c : input) {
        if (c == '(') {
            TreeNode* newNode = new TreeNode();
            newNode->parent = current;
            current->children.push_back(newNode);
            stack.push_back(current);
            current = newNode;
        } else if (c == ',' || c == ')') {
            if (!token.empty()) {
                TreeNode* leaf = new TreeNode();
                leaf->name = token;
                leaf->parent = current;
                current->children.push_back(leaf);
                token.clear();
            }
            if (c == ')') {
                current = stack.back();
                stack.pop_back();
            }
        } else if (c == ';') {
            if (!token.empty()) {
                current->name = token;
                token.clear();
            }
        } else {
            token += c;
        }
    }
    return root;
}

// Collect all leaf nodes
void collectLeaves(TreeNode* node, vector<TreeNode*>& leaves) {
    if (node->children.empty()) {
        leaves.push_back(node);
    } else {
        for (auto child : node->children) {
            collectLeaves(child, leaves);
        }
    }
}

// Prune leaves not in the set
bool pruneTree(TreeNode* node, const unordered_set<string>& keepTips) {
    if (node->children.empty()) {
        return keepTips.count(node->name) > 0;
    }

    vector<TreeNode*> newChildren;
    for (auto child : node->children) {
        if (pruneTree(child, keepTips)) {
            newChildren.push_back(child);
        } else {
            delete child;
        }
    }

    node->children = newChildren;
    return !node->children.empty();
}

// Convert back to Newick format
string toNewick(TreeNode* node) {
    stringstream ss;
    if (!node->children.empty()) {
        ss << "(";
        for (size_t i = 0; i < node->children.size(); ++i) {
            ss << toNewick(node->children[i]);
            if (i != node->children.size() - 1)
                ss << ",";
        }
        ss << ")";
    }
    if (!node->name.empty())
        ss << node->name;
    return ss.str();
}

void deleteTree(TreeNode* node) {
    for (auto child : node->children)
        deleteTree(child);
    delete node;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " tree.nwk N" << endl;
        return 1;
    }

    string newickFile = argv[1];
    int keepN = stoi(argv[2]);

    ifstream in(newickFile);
    if (!in) {
        cerr << "Error: cannot open Newick file " << newickFile << endl;
        return 1;
    }

    stringstream buffer;
    buffer << in.rdbuf();
    string newick = buffer.str();

    TreeNode* tree = parseNewick(newick);

    vector<TreeNode*> allLeaves;
    collectLeaves(tree, allLeaves);

    if (keepN > allLeaves.size()) {
        cerr << "Error: Requested more tips than available (" << allLeaves.size() << ")." << endl;
        deleteTree(tree);
        return 1;
    }
    std::cerr << "Read input tree with " << allLeaves.size() << " tips." << endl;

    // Randomly select N tips
    shuffle(allLeaves.begin(), allLeaves.end(), mt19937(random_device{}()));
    unordered_set<string> keepTips;
    for (int i = 0; i < keepN; ++i) {
        keepTips.insert(allLeaves[i]->name);
    }
    std::cerr << "Keeping " << keepN << " tips." << endl;

    pruneTree(tree, keepTips);

    cout << toNewick(tree) << ";" << endl;

    deleteTree(tree);
    return 0;
}
