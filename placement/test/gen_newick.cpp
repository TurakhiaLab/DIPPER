#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <fstream>

class Node {
    public:
        Node(double edge_length, int id) : id(id), e(edge_length), child(nullptr), next(nullptr), parent(nullptr) {};
	    int id;
	    double e;
	    Node* child;
	    Node* next;
	    Node* parent;
};


void printTree(Node* root, int* idx, std::string& outputString) {
	Node* n;
	if(root->child!=nullptr) {
		outputString += "(";
		for(n=root->child;n!=nullptr;n=n->next) {
			if(n!=root->child) outputString += ",";
			printTree(n, idx, outputString);
            outputString += ":1";
		}
		outputString += ")";
	}
	else {
		outputString += "id_"; // leaves name
		outputString += std::to_string(root->id);
    }
}

int main(int argc,char** argv) {
	int i,n=atoi(argv[1]);
    srand(1);
    std::unordered_map<int, Node*> allNodes;
	std::vector<std::pair<int, Node*>> noParent;
    // Get all leaves
	for(i=0;i<n;++i) {
		// double edge_len = (double)rand()/(double)RAND_MAX + 0.1; // random branch length
		double edge_len = 0.0; // 0 branch length
        Node* node = new Node(edge_len, i);
        allNodes[i] = node;
		noParent.push_back(std::make_pair(i, node));
	}

    int totalNodes = n;
	std::vector<int> selected;
	std::vector<int> selectedKeys;
    while (!noParent.empty()) {
		int childrenNum = 2; // Binary tree
        // int childrenNum = 2 + rand() % 4; // Non-binary tree
		if (childrenNum > noParent.size()) childrenNum = noParent.size();
        while (selected.size() < childrenNum) {
            int select = rand() % noParent.size();
            if (std::find(selected.begin(), selected.end(), select) == selected.end()) {
                selected.push_back(select);
            }
        }
        // double edge_len = (double)rand()/(double)RAND_MAX + 0.1;
		double edge_len = 0.0;

        Node* node = new Node(edge_len, totalNodes);
        
        auto iter = noParent.begin();
		for (auto id: selected) selectedKeys.push_back(noParent[id].first);
		
        node->child = allNodes[selectedKeys[0]];
        for (int i = 1; i < selectedKeys.size(); ++i) {
            allNodes[selectedKeys[i-1]]->next = allNodes[selectedKeys[i]];
        }
		for (int i = 1; i < selectedKeys.size(); ++i) {
            allNodes[selectedKeys[i]]->parent = node;
        }
        allNodes[totalNodes] = node;
		std::sort(selected.rbegin(), selected.rend()); 
		for (auto id: selected) noParent.erase(noParent.begin()+id);
		noParent.push_back(std::make_pair(totalNodes, node));
        totalNodes += 1;

		selected.clear();
		selectedKeys.clear();
        if (noParent.size() == 1) break;
    }

    int start = 0;
	std::string out_str = "";
	printTree(allNodes[totalNodes-1], &start, out_str);
	out_str += ";\n";
	std::cout << out_str;
	return 0;
}