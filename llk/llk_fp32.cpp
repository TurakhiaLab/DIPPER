// #include "core_likelihood_fp32.hpp"
// #include "core_likelihood_fp32.cpp"
// #include "tree_fp32.cpp"

// int main(int argc, char** argv)
// {
//     if (argc != 3)
//         std::cerr << "Usage: "<< argv[0] <<" [alignment file] [newick tree file]" << std::endl;

//     utility::msa_seq(argv[1]);    // Store MSA data into data-structure
//     utility::subs_param.resize(10); // Set Subs parameter
//     for (size_t i = 0; i < 10; i++) utility::subs_param[i] = 1;
//     utility::rate_matrix_calc(); // Find the rate matrix

    
//     if (1) // Print matrix_exp Matrix
//     {
//         for (size_t i = 0; i < utility::rate_matrix.size(); i++)
//         {
//             for (size_t j = 0; j < utility::rate_matrix[0].size(); j++) {
//                 printf("%.20f ", utility::rate_matrix[i][j]);
//                 // std::cout << utility::rate_matrix[i][j] << "\t";
//             }
//             std::cout << "\n";
//         }
//         for (auto &p: utility::pi) std::cout << p << "\n";
//     }


//     std::ifstream newickTree(argv[2]);
//     std::string newickString;
//     newickTree >> newickString;
//     utility::Tree tree(newickString);

//     utility::felsenstein_pruning(tree);

//     // utility::bottom_up(tree);
//     // std::cout << "\n";
//     // utility::top_down(tree);
//     // std::cout << "\n";
//     // utility::marginal(tree);


//     // utility::fitch(tree);

//     // if (PRINT_INFERENCE)
//     // {
//     //     std::string tip="Mon";
//     //     auto search = tree.allNodes.find(tip);
//     //     if (search != tree.allNodes.end()) utility::printLikelihoodInference(tree, search->second);
//     //     else {std::cerr<<"Tip not found";}
//     // }

//     return 0;
// }




#ifndef CORE_LIKELIHOOD_FP32
#include "core_likelihood_fp32.hpp"
#endif
#include "core_likelihood_fp32.cpp"
#include "tree_fp32.cpp"
#include <chrono>

int main(int argc, char** argv)
{
    if (argc != 3)
        std::cerr << "Usage: "<< argv[0] <<" [alignment file] [newick tree file]" << std::endl;

    utility::msa_seq(argv[1]);    // Store MSA data into data-structure
    utility::subs_param.resize(10); // Set Subs parameter
    for (size_t i = 0; i < 10; i++) utility::subs_param[i] = 1;
    utility::rate_matrix_calc(); // Find the rate matrix

    
    if (0) // Print matrix_exp Matrix
    {
        for (size_t i = 0; i < utility::rate_matrix.size(); i++)
        {
            for (size_t j = 0; j < utility::rate_matrix[0].size(); j++) std::cout << utility::rate_matrix[i][j] << "\t";
            std::cout << "\n";
        }
        for (auto &p: utility::pi) {
            printf("%.50f\n", p);
        }
    }


    std::ifstream newickTree(argv[2]);
    std::string newickString;
    newickTree >> newickString;
    utility::Tree tree(newickString);

    auto start = std::chrono::high_resolution_clock::now();

    utility::felsenstein_pruning(tree);

    auto end = std::chrono::high_resolution_clock::now();  

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "Felsensteins pruning execution took " << duration.count() << " microseconds." << std::endl;
   
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

