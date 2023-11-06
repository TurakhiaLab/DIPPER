/*
Naive matrix multiplication algorithm 
CPU-Only, for benchmarking purposes.

Compute e^R using Maclaurin series approximation for 100 iterations. 
*/

#include <iostream>
#include <vector>
#include <cstdlib> // for random number generation
#include <ctime>   // for seeding random number generator

using namespace std;

// generate a random sparse 2D matrix
vector<vector<int>> generateSparseMatrix(int rows, int cols, double sparsity) {
    vector<vector<int>> matrix(rows, vector<int>(cols, 0));

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double randomValue = (double)rand() / RAND_MAX; // Generate a random value between 0 and 1

            if (randomValue > sparsity) {
                // Set a non-zero value in the matrix
                matrix[i][j] = rand() % 100; // You can adjust the range of non-zero values as needed
            }
        }
    }

    return matrix;
}

// compute the factorial of a number
long factorial(const int n) {
    long f = 1;
    for (int i=1; i<=n; i++)
        f *= i;
    return f;
}

// multiply two square matrices together
vector<vector<int>> mMultiply(const vector<vector<int>>& A, const vector<vector<int>>& B) {
    int n = A.size();
    vector<vector<int>> result(n, vector<int>(n, 0));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return result;
}

// add two square matrices together
vector<vector<int>> mAdd(const vector<vector<int>>& A, const vector<vector<int>>& B) {
    int n = A.size();
    vector<vector<int>> result(n, vector<int>(n, 0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i][j] = A[i][j] + B[i][j];
        }
    }

    return result;      
}

// divide a vector by a scalar.
vector<vector<int>> sDivide(float f, const vector<vector<int>>& matrix) {
    int n = matrix.size();
    vector<vector<int>> result(n, vector<int>(n, 0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i][j] = matrix[i][j]/f;
        }
    }

    return result;

}

// compute the Maclaurin series to approximate e^R
vector<vector<int>> computeSeries(int iterations, const vector<vector<int>>& matrix) {
    vector<vector<int>> matrixProduct = matrix;

    std::vector<std::vector<int>> ones(matrix.size(), std::vector<int>(matrix.size(), 1));

    vector<vector<int>> output = mAdd(ones, matrix); // terms 0 and 1 of the series.

    for (int i = 2; i <= iterations; i++) {
        matrixProduct = mMultiply(matrix, matrixProduct);
        output = mAdd(output, sDivide(factorial(i), matrixProduct));
    }

    return output;
}



int main() {
    srand(55); // Seed the random number generator

    int rows = 10000; // Set the number of rows
    int cols = 10000; // Set the number of columns
    double sparsity = 0.7; // Set the sparsity factor (higher value = more sparse)

    vector<vector<int>> sparseMatrix = generateSparseMatrix(rows, cols, sparsity);

    vector<vector<int>> computed = computeSeries(100, sparseMatrix); // compute the series (e^R)

    // Print the generated matrix
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            cout << computed[i][j] << " ";
        }
        cout << endl;
    }

    return 0;
}
