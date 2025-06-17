#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <sstream>

using namespace std;

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <input_file>\n";
        return 1;
    }

    ifstream infile(argv[1]);

    if (!infile ) {
        cerr << "Error opening file.\n";
        return 1;
    }

    int num_taxa;
    infile >> num_taxa;
    infile.ignore(); // Skip to next line

    vector<string> taxa_names(num_taxa);
    vector<vector<double>> matrix(num_taxa, vector<double>(num_taxa, 0.0));

    int i = 0;
    string line;
    while (getline(infile, line) && i < num_taxa) {
        istringstream iss(line);
        iss >> taxa_names[i];

        for (int j = 0; j < i; ++j) {
            double dist;
            if (!(iss >> dist)) {
                cerr << "Error reading distance at line " << i + 1 << endl;
                return 1;
            }
            matrix[i][j] = dist;
            matrix[j][i] = dist; // fill upper triangle
        }
        i++;
    }
    std::cerr << "Converting to full matrix format..." << std::endl;

    // Write full matrix to output
    std::cout << num_taxa << "\n";
    for (int i = 0; i < num_taxa; ++i) {
        std::cout << left << setw(10) << taxa_names[i];
        for (int j = 0; j < num_taxa; ++j) {
            std::cout << fixed << setprecision(6) << setw(10) << matrix[i][j];
        }
        std::cout << "\n";
    }

    infile.close();

    cout << "Conversion complete. Output written to " << argv[2] << endl;
    return 0;
}
