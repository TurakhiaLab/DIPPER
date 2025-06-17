#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <random>
#include <algorithm>
#include <regex>
#include <zlib.h>

using namespace std;

struct FastaSeq {
    string header;
    string sequence;
};

vector<FastaSeq> readFastaGz(const string& filename) {
    vector<FastaSeq> sequences;
    gzFile file = gzopen(filename.c_str(), "rb");
    if (!file) {
        cerr << "Error: cannot open gzipped FASTA file " << filename << endl;
        exit(1);
    }

    const int BUF_SIZE = 8192;
    char buffer[BUF_SIZE];
    string currentHeader, currentSeq;

    while (gzgets(file, buffer, BUF_SIZE)) {
        string line(buffer);
        if (!line.empty() && line.back() == '\n') line.pop_back();
        if (!line.empty() && line.back() == '\r') line.pop_back();

        if (line.empty()) continue;

        if (line[0] == '>') {
            if (!currentHeader.empty()) {
                sequences.push_back({currentHeader, currentSeq});
                currentSeq.clear();
            }
            currentHeader = line.substr(1);
        } else {
            currentSeq += line;
        }
    }

    if (!currentHeader.empty()) {
        sequences.push_back({currentHeader, currentSeq});
    }

    gzclose(file);
    return sequences;
}

vector<FastaSeq> readFastaPlain(const string& filename) {
    ifstream in(filename);
    if (!in) {
        cerr << "Error: cannot open FASTA file " << filename << endl;
        exit(1);
    }

    vector<FastaSeq> sequences;
    string line, currentHeader, currentSeq;

    while (getline(in, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (!currentHeader.empty()) {
                sequences.push_back({currentHeader, currentSeq});
                currentSeq.clear();
            }
            currentHeader = line.substr(1); // remove '>'
        } else {
            currentSeq += line;
        }
    }

    if (!currentHeader.empty()) {
        sequences.push_back({currentHeader, currentSeq});
    }

    return sequences;
}

unordered_set<string> parseNewickLeaves(const string& filename) {
    ifstream in(filename);
    if (!in) {
        cerr << "Error: cannot open Newick file " << filename << endl;
        exit(1);
    }

    stringstream buffer;
    buffer << in.rdbuf();
    string tree = buffer.str();

    // Remove whitespace
    tree.erase(remove_if(tree.begin(), tree.end(), ::isspace), tree.end());

    unordered_set<string> leaves;
    regex leafRegex("([A-Za-z0-9_\\.-]+)[):,]");  // captures names before ):,
    smatch match;

    string::const_iterator searchStart(tree.cbegin());
    while (regex_search(searchStart, tree.cend(), match, leafRegex)) {
        leaves.insert(match[1]);
        searchStart = match.suffix().first;
    }

    return leaves;
}

vector<FastaSeq> filterByLeaves(const vector<FastaSeq>& sequences, const unordered_set<string>& leafNames) {
    vector<FastaSeq> filtered;
    for (const auto& seq : sequences) {
        string headerBase = seq.header;
        size_t spacePos = headerBase.find_first_of(" \t");
        if (spacePos != string::npos)
            headerBase = headerBase.substr(0, spacePos);

        if (leafNames.count(headerBase)) {
            filtered.push_back(seq);
        }
    }
    return filtered;
}

bool isGzipFile(const string& filename) {
    return filename.size() > 3 && filename.substr(filename.size() - 3) == ".gz";
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " input.fasta[.gz] tree.nwk" << endl;
        return 1;
    }

    string fastaFile = argv[1];
    string newickFile = argv[2];

    vector<FastaSeq> sequences;
    if (isGzipFile(fastaFile)) {
        sequences = readFastaGz(fastaFile);
    } else {
        sequences = readFastaPlain(fastaFile);
    }

    unordered_set<string> leafNames = parseNewickLeaves(newickFile);
    vector<FastaSeq> filtered = filterByLeaves(sequences, leafNames);

    for (const auto& seq : filtered) {
        cout << ">" << seq.header << "\n" << seq.sequence << "\n";
    }

    return 0;
}
