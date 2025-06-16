#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>

int main(int argc, char* argv[]) {
    std::ifstream infile(argv[1]); 
    if (!infile.is_open()) {
        std::cerr << "Error opening file.\n";
        return 1;
    }

    int n;
    infile >> n;
    infile.ignore();  

    std::vector<std::string> labels;
    std::vector<std::vector<double>> matrix;

    std::string line;
    for (int i = 0; i < n; ++i) {
        std::getline(infile, line);
        std::istringstream iss(line);

        std::string label;
        iss >> label;
        labels.push_back(label);

        std::vector<double> row;
        double value;
        while (iss >> value) {
            row.push_back(value);
        }
        matrix.push_back(row);
    }

    // Output as tab-delimited
    std::cout << "\t";
    for (const auto& label : labels) {
        std::cout << label << "\t";
    }
    std::cout << "\n";

    for (int i = 0; i < n; ++i) {
        std::cout << labels[i] << "\t";
        for (double value : matrix[i]) {
            std::cout << std::fixed << std::setprecision(10) << value << "\t";
        }
        std::cout << "\n";
    }

    return 0;
}
