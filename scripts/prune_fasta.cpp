#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <random>
#include <algorithm>
#include <zlib.h>

using namespace std;

struct FastaSeq {
    string header;
    string sequence;
};

// Read from gzipped file
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

// Read from plain text file
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

// Random sampling
void printRandomSequences(const vector<FastaSeq>& sequences, size_t n) {
    if (n > sequences.size()) {
        cerr << "Error: Requested more sequences than available (" << sequences.size() << ")." << endl;
        exit(1);
    }

    vector<size_t> indices(sequences.size());
    iota(indices.begin(), indices.end(), 0);
    shuffle(indices.begin(), indices.end(), mt19937(random_device{}()));

    for (size_t i = 0; i < n; ++i) {
        const auto& seq = sequences[indices[i]];
        cout << ">" << seq.header << "\n" << seq.sequence << "\n";
    }
}

bool isGzipFile(const string& filename) {
    return filename.size() > 3 &&
           filename.substr(filename.size() - 3) == ".gz";
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " input.fasta[.gz] N" << endl;
        return 1;
    }

    string fastaFile = argv[1];
    size_t n = stoi(argv[2]);

    vector<FastaSeq> sequences;

    if (isGzipFile(fastaFile)) {
        sequences = readFastaGz(fastaFile);
    } else {
        sequences = readFastaPlain(fastaFile);
    }

    printRandomSequences(sequences, n);
    return 0;
}

