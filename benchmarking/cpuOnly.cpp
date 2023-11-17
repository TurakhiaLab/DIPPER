#include "generateMatrix.hpp"

__half hA[N*N];
__half hB[N*N];
float result[N * N];

int main() {
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            float sum  = 0.0f;
            for (int k1 = 0; k1 < N; k1++) {
                auto posA = i * N + k1;
                auto posB = k1 * N + j;
                sum      += static_cast<float>(hA[posA]) *  // [i][k]
                            static_cast<float>(hB[posB]);   // [k][j]
            }
            auto posC       = i * N + j;
            result[posC] = sum;  // [i][j]
        }
    }

    return 0;

}