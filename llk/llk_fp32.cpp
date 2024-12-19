#ifndef CORE_LIKELIHOOD_FP32
#include "core_likelihood_fp32.hpp"
#endif
#include "core_likelihood_fp32.cpp"
#include "tree_fp32.cpp"
#include <chrono>

int main(int argc, char **argv)
{
    if (argc != 3)
        std::cerr << "Usage: " << argv[0] << " [alignment file] [newick tree file]" << std::endl;

    utility::msa_seq(argv[1]);      // Store MSA data into data-structure
    utility::subs_param.resize(10); // Set Subs parameter
    for (size_t i = 0; i < 10; i++)
        utility::subs_param[i] = 1;
    utility::rate_matrix_calc(); // Find the rate matrix

    if (1) // Print matrix_exp Matrix
    {
        for (size_t i = 0; i < utility::rate_matrix.size(); i++)
        {
            for (size_t j = 0; j < utility::rate_matrix[0].size(); j++)
                std::cout << utility::rate_matrix[i][j] << "\t";
            std::cout << "\n";
        }
        for (auto &p : utility::pi)
            std::cout << p << "\n";
    }

    std::ifstream newickTree(argv[2]);
    std::string newickString;
    newickTree >> newickString;
    utility::Tree tree(newickString);

    auto start = std::chrono::high_resolution_clock::now();

    multiplication_values_file.open("addition_values_SUM1.txt");
    multiplication_values_file << std::setprecision(15);

    utility::felsenstein_pruning(tree);

    multiplication_values_file.close();

    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "Felsensteins pruning execution took " << duration.count() << " microseconds." << std::endl;

    std::ofstream outFile("bottom_values_SUM1.txt");

    for (auto &p : tree.allNodes)
    {
        utility::Node *node = p.second;
        for (size_t i = 0; i < node->bottom.size(); ++i)
        {
            for (size_t j = 0; j < node->bottom[i].size(); ++j)
            {
                outFile << node->bottom[i][j] << "\n";
            }
        }
    }
    outFile.close();

    std::ofstream outFile1("bottom_values_by_layer_SUM1.csv");
    outFile1 << "layer,bottom_value\n";

    for (auto &p : tree.allNodes)
    {
        utility::Node *node = p.second;
        size_t layer = node->level;
        for (size_t i = 0; i < node->bottom.size(); ++i)
        {
            for (size_t j = 0; j < node->bottom[i].size(); ++j)
            {
                float bottom_value = node->bottom[i][j];
                outFile1 << layer << "," << bottom_value << "\n";
            }
        }
    }
    outFile.close();

    // // print leaf nodes
    // for (auto &p: tree.allNodes) {
    //     if (p.second->children.size() == 0) {
    //         printf("Leaf Node\n");
    //         printf("branchLength: %f\n", p.second->branchLength);
    //         printf("identifier: %s\n", p.second->identifier.c_str());
    //         printf("bottom matrix:\n");
    //         for (int i = 0; i < 10; ++i) {
    //             for (int j = 0; j < 5; ++j) {
    //                 printf("%f ", p.second->bottom[j][i]);
    //             }
    //             printf("\n");
    //         }
    //     }
    // }

    // printf("Root Node\n");
    // printf("branchLength: %f\n", tree.root->branchLength);

    // printf("Bottom Matrix:\n");
    // for (int i = 0; i < 10; ++i) {
    //     for (int j = 0; j < 5; ++j) {
    //         printf("%f ", tree.root->bottom[i][j]);
    //     }
    //     printf("\n");
    // }

    // utility::bottom_up(tree);
    // std::cout << "\n";
    // utility::top_down(tree);
    // std::cout << "\n";
    // utility::marginal(tree);

    // utility::fitch(tree);

    // if (PRINT_INFERENCE)
    // {
    //     std::string tip="Mon";
    //     auto search = tree.allNodes.find(tip);
    //     if (search != tree.allNodes.end()) utility::printLikelihoodInference(tree, search->second);
    //     else {std::cerr<<"Tip not found";}
    // }

    return 0;
}
