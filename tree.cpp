#include "tree.hpp"

Node::Node(const std::string& id, double len){
    this->name = id;
    level = 1;
    bl = len;
    parent = nullptr;
}


Node::Node(const std::string& id, Node* par, double len){
    this->name = id;
    bl = len;
    parent = par;
    level = par->level + 1;
    par->children.push_back(this);
}

size_t Node::getNumLeaves(){
    size_t num_leaves = 0;
    if (children.size() == 0) return num_leaves;
    for (auto ch: children){
        if (ch->is_leaf()) num_leaves += 1;
        else num_leaves += ch->getNumLeaves();
    }
    return num_leaves;
}

size_t Node::getNumNodes(){
    size_t num_nodes = 1;
    if (children.size() == 0) return num_nodes;
    for (auto ch: children){
        num_nodes += ch->getNumNodes();
    }
    return num_nodes;
}

void stringSplit (std::string const& s, char delim, std::vector<std::string>& words) {
    size_t start_pos = 0, end_pos = 0, temp_pos = 0;
    while ((end_pos = s.find(delim, start_pos)) != std::string::npos) {
        if (end_pos >= s.length()) {
            break;
        }
        std::string sub;
        if (temp_pos == 0) {
            sub = s.substr(start_pos, end_pos-start_pos);
            if (std::count(sub.begin(), sub.end(), '\'') % 2 == 1) {
                temp_pos = start_pos;
            }
            else {
                words.emplace_back(sub);
            }
        }
        else {
            sub = s.substr(temp_pos, end_pos-temp_pos);
            if (std::count(sub.begin(), sub.end(), '\'') % 2 == 0) {
                temp_pos = 0;
                words.emplace_back(sub);
            }
        }
        // words.emplace_back(s.substr(start_pos, end_pos-start_pos));
        start_pos = end_pos+1;
    }
    auto last = s.substr(start_pos, s.size()-start_pos);
    if (last != "") {
        words.push_back(std::move(last));
    }
}

std::string stripString(std::string s){
    while(s.length() && s[s.length() - 1] == ' '){
        s.pop_back();
    }
    for(size_t i = 0; i < s.length(); i++){
        if(s[i] != ' '){
            return s.substr(i);
        }
    }
    return s;
}

void Tree::dfsExpansion(Node* node,
                                     std::vector< Node* >& vec) {
    vec.push_back(node);
    for(auto child: node->children) {
        dfsExpansion(child, vec);
    }
}


int dfsExpansionSize(Node* node) {
    int c = 0;
    if (node->children.size() == 0) {
        c++;
        return c;
    }
    for (auto &n: node->children){
        c += dfsExpansionSize(n);
    }
    std::cout << node->name << "\t" << c << std::endl;
    return c;
}

std::string Tree::getNewickString(Node* node) {

    // traversal to print each node subtree size
    // int s = dfsExpansionSize(node);
    // exit(0);
    std::vector< Node* > traversal;
    dfsExpansion(node, traversal);
    std::string newick;

    if (traversal.size() == 1) {
        newick += node->name;
        return newick;
    }

    size_t level_offset = node->level-1;
    size_t curr_level = 0;
    bool prev_open = true;

    std::stack<std::string> node_stack;
    std::stack<float> branch_length_stack;

    for (auto n: traversal) {
        size_t level = n->level-level_offset;
        float branch_length = n->bl;

        if(curr_level < level) {
            if (!prev_open) {
                newick += ',';
            }
            size_t l = level - 1;
            if (curr_level > 1) {
                l = level - curr_level;
            }
            for (size_t i=0; i < l; i++) {
                newick += '(';
                prev_open = true;
            }
            if (n->children.size() == 0) {

                newick += n->name;

                if (branch_length >= 0) {
                    newick += ':';
                    newick += std::to_string(branch_length);
                }
                prev_open = false;
            } else {
                node_stack.push(n->name);
                branch_length_stack.push(branch_length);
            }
        } else if (curr_level > level) {
            prev_open = false;
            for (size_t i = level; i < curr_level; i++) {
                newick += ')';

                newick += node_stack.top();

                if (branch_length_stack.top() >= 0) {
                    newick += ':';
                    newick += std::to_string(branch_length_stack.top());
                }
                node_stack.pop();
                branch_length_stack.pop();
            }
            if (n->children.size() == 0) {
                newick += ',';
                newick += n->name;

                if (branch_length >= 0) {
                    newick += ':';
                    newick += std::to_string(branch_length);
                }
            } else {
                node_stack.push(n->name);
                branch_length_stack.push(branch_length);
            }
        } else {
            prev_open = false;
            if (n->children.size() == 0) {

                newick += ',';
                newick += n->name;

                if (branch_length >= 0) {
                    newick += ':';
                    newick += std::to_string(branch_length);
                }
            } else {
                node_stack.push(n->name);
                branch_length_stack.push(branch_length);
            }
        }
        curr_level = level;
    }
    size_t remaining = node_stack.size();
    for (size_t i = 0; i < remaining; i++) {
        newick += ')';
        newick += node_stack.top();

        if (branch_length_stack.top() >= 0) {
            newick += ':';
            newick += std::to_string(branch_length_stack.top());
        }
        node_stack.pop();
        branch_length_stack.pop();
    }

    newick += ';';
    std::cerr << newick << std::endl;
    return newick;
}

Tree::Tree(std::string newickString, size_t totalLeaves) {
    newickString = stripString(newickString);

    Node* treeRoot = nullptr;

    std::vector<std::string> leaves;
    std::vector<size_t> numOpen;
    std::vector<size_t> numClose;
    std::vector<std::queue<float>> branchLen (128);  // will be resized later if needed
    size_t level = 0;

    std::vector<std::string> s1;
    stringSplit(newickString, ',', s1);

    numOpen.reserve(s1.size());
    numClose.reserve(s1.size());
    

    for (auto s: s1) {
        size_t no = 0;
        size_t nc = 0;
        size_t leafDepth = 0;

        bool stop = false;
        bool branchStart = false;
        bool nameZone = false;
        bool hasApo = false;
        std::string leaf = "";
        std::string branch = "";

        for (auto c: s) {
            if (nameZone) {
                leaf += c;
                if (c == '\'') nameZone = false;
            } else if (c == '\'' && !nameZone) {
                nameZone = true;
                hasApo = true;
                leaf += c;
            } else if (c == ':') {
                stop = true;
                branch = "";
                branchStart = true;
            } else if (c == '(') {
                no++;
                level++;
                if (branchLen.size() <= level) {
                    branchLen.resize(level*2);
                }
            } else if (c == ')') {
                stop = true;
                nc++;
                // float len = (branch.size() > 0) ? std::stof(branch) : -1.0;
                float len = (branch.size() >= 0) ? std::stof(branch) : 1.0;
                // if (len == 0) len = 1.0;
                branchLen[level].push(len);
                level--;
                branchStart = false;
            } else if (!stop) {
                leaf += c;
                branchStart = false;
                leafDepth = level;

            } else if (branchStart) {
                if (isdigit(c)  || c == '.' || c == 'e' || c == 'E' || c == '-' || c == '+') {
                    branch += c;
                }
            }
        }
        if (hasApo && leaf[0] == '\'' && leaf[leaf.length()-1] == '\'') leaf = leaf.substr(1, leaf.length()-2);
        leaves.push_back(std::move(leaf));
        numOpen.push_back(no);
        numClose.push_back(nc);
        // float len = (branch.size() > 0) ? std::stof(branch) : -1.0;
        float len = (branch.size() >= 0) ? std::stof(branch) : 1.0;
        // if (len == 0) len = 1.0;
        branchLen[level].push(len);

        // Adjusting max and mean depths
        m_maxDepth = std::max(m_maxDepth, leafDepth);
        m_meanDepth += leafDepth;

    }

    m_meanDepth /= leaves.size();

    if (level != 0) {
        fprintf(stderr, "ERROR: incorrect Newick format!\n");
        exit(1);
    }

    m_numLeaves = leaves.size();

    m_currInternalNode = totalLeaves - 1; // Internal nodes start from the last leaf index
    std::stack<Node*> parentStack;

    for (size_t i=0; i<leaves.size(); i++) {
        auto leaf = leaves[i];
        auto no = numOpen[i];
        auto nc = numClose[i];
        for (size_t j=0; j<no; j++) {
            // int idx = m_currInternalNode+1;
            int idx = m_currInternalNode+1;
            std::string nid = newInternalNodeId();
            Node* newNode = nullptr;
            if (parentStack.size() == 0) {
                // nid = "node_"  + std::to_string(2*m_numLeaves+totalLeaves-2);
                // int idx = m_numLeaves+totalLeaves-2;
                newNode = new Node(nid, branchLen[level].front());
                newNode->idx = idx;
                treeRoot = newNode;
            } else {
                // int idx = m_currInternalNode+1;
                // nid = newInternalNodeId();
                newNode = new Node(nid, parentStack.top(), branchLen[level].front());
                newNode->idx = idx;
                // if (branchLen[level].front() < 0.001) std::cout << nid << '\t' << branchLen[level].front() << '\n';
        
            }
            branchLen[level].pop();
            level++;

            allNodes[nid] = newNode;
            parentStack.push(newNode);
        }
        Node* leafNode = new Node(leaf, parentStack.top(), branchLen[level].front());
        leafNode->idx = m_numLeafID++;
        if (leafNode->name=="T22") std::cerr << "Leaf T22 found with len " << branchLen[level].front() << '\n';
        // if (branchLen[level].front() < 0.001) std::cout << leaf << '\t' << branchLen[level].front() << '\n';
        /* Group Id */
        allNodes[leaf] = leafNode;

        branchLen[level].pop();
        for (size_t j=0; j<nc; j++) {
            parentStack.pop();
            level--;
        }
    }

    if (treeRoot == nullptr) {
        fprintf(stderr, "WARNING: Tree found empty!\n");
    }

    treeRoot->bl = 0;
    root = treeRoot;

}


Tree::~Tree() {
    for (auto n: this->allNodes) {
        delete n.second;
    }
    this->allNodes.clear();
    this->root = nullptr;
}