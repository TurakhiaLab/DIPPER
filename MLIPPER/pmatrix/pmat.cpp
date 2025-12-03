#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
// #include <lapacke.h>
#include "pmat.h"

extern "C" {
    void dsyevd_(char* jobz, char* uplo, int* n,
                 double* a, int* lda,
                 double* w,
                 double* work, int* lwork,
                 int* iwork, int* liwork,
                 int* info);
}
static void eye(double* I, int n){
    std::fill(I, I+n*n, 0.0);
    for(int i=0;i<n;++i) I[i*n+i] = 1.0;
}

static void matmul_nn(const double* A, const double* B, double* C, int n){
    std::fill(C, C+n*n, 0.0);
    for(int i=0;i<n;++i){
        for(int k=0;k<n;++k){
            double aik = A[i*n + k];
            for(int j=0;j<n;++j){
                C[i*n + j] += aik * B[k*n + j];
            }
        }
    }
}

static void matadd_inplace(double* A, const double* B, int n){
    for(int i=0;i<n*n;++i) A[i] += B[i];
}

static double inf_norm(const std::vector<double>& M, int n){
    double m = 0.0;
    for(int i=0;i<n*n;++i) m = std::max(m, std::abs(M[i]));
    return m;
}

static double check_rowsum_max_dev(const double* P, int n){
    double worst = 0.0;
    for(int i=0;i<n;++i){
        double s = 0.0;
        for(int j=0;j<n;++j) s += P[i*n + j];
        worst = std::max(worst, std::abs(s - 1.0));
    }
    return worst;
}

static void clamp_neg_and_row_norm(double* P, int n){
    for(int i=0;i<n;++i){
        double s = 0.0;
        for(int j=0;j<n;++j){
            double &x = P[i*n + j];
            if (x < 0.0 && x > -1e-14) x = 0.0;
            s += x;
        }
        if (s != 0.0){
            for(int j=0;j<n;++j) P[i*n + j] /= s;
        }
    }
}

void pmatrix_from_triple(const double* Vinv, const double* V,
                                const double* lamb, double r, double t, double p,
                                double* P, int n)
{
    // I
    std::vector<double> I(n*n, 0.0); 
    for (int i=0;i<n;++i) I[i*n+i] = 1.0;

    // D = diag( expm1(λ r t) )
    std::vector<double> D(n*n, 0.0);
    for (int j=0;j<n;++j) D[j*n + j] = std::expm1(lamb[j] * r * t);

    // T = V * D
    std::vector<double> T(n*n, 0.0);
    for (int i=0;i<n;++i){
        for (int k=0;k<n;++k){
            double vik = V[i*n + k];
            for (int j=0;j<n;++j){
                T[i*n + j] += vik * D[k*n + j];
            }
        }
    }

    // P = I + T * Vinv  (= I + V * [e^{Λt}-I] * V^{-1})
    std::fill(P, P + n*n, 0.0);
    for (int i=0;i<n;++i){
        for (int k=0;k<n;++k){
            double tik = T[i*n + k];
            for (int j=0;j<n;++j){
                P[i*n + j] += tik * Vinv[k*n + j];
            }
        }
    }
    for (int i=0;i<n;++i) P[i*n + i] += 1.0;

    // invariants: P = (1-p) P + p I   （只做線性混合，不要動 exponent）
    if (p > 0){
        for (int i=0;i<n*n;++i) P[i] = (1.0 - p)*P[i] + p*I[i];
    }
}



void pmatrix_direct(const double* U, const double* pi,
                           const double* lamb, double r, double t, double p,
                           double* P, int n)
{
    std::vector<double> Dsqrt(n), Dsqrt_inv(n);
    for(int i=0;i<n;++i){ Dsqrt[i]=std::sqrt(pi[i]); Dsqrt_inv[i]=1.0/Dsqrt[i]; }
    // E = diag(exp(lamb*r*t))
    std::vector<double> E(n*n,0.0), temp(n*n,0.0), M(n*n,0.0), I(n*n,0.0);
    for(int j=0;j<n;++j) E[j*n + j] = std::exp(lamb[j]*r*t);
    // temp = U * E
    matmul_nn(U, E.data(), temp.data(), n);
    // M = temp * U^T  （用索引展開 U^T）
    std::fill(M.begin(), M.end(), 0.0);
    for(int i=0;i<n;++i){
        for(int k=0;k<n;++k){
            double tik = temp[i*n + k];
            for(int j=0;j<n;++j){
                M[i*n + j] += tik * U[j*n + k]; // U^T[k,j] = U[j,k]
            }
        }
    }
    // P0 = D^{-1} * M * D
    std::vector<double> P0(n*n, 0.0);
    for(int i=0;i<n;++i){
        for(int j=0;j<n;++j){
            P0[i*n + j] = Dsqrt_inv[i] * M[i*n + j] * Dsqrt[j];
        }
    }
    eye(I.data(), n);
    for(int i=0;i<n*n;++i) P[i] = (1.0 - p)*P0[i] + p*I[i];
}


EigResult gtr_eigendecomp_cpu(
    const double* Q_rowmajor,   // Q（row-major，大小 n*n）
    const double* pi,           // 平衡頻率（長度 n，sum=1）
    int n)
{
    if (n <= 0) throw std::invalid_argument("n must be > 0");

    // 1) D = diag(sqrt(pi))
    std::vector<double> sqrtpi(n);
    for (int i = 0; i < n; ++i) {
        if (pi[i] <= 0.0) throw std::invalid_argument("pi must be > 0");
        sqrtpi[i] = std::sqrt(pi[i]);
    }

    // 2) S = D * Q * D^{-1}  （S 對稱）
    //   S[i,j] = sqrt(pi[i]) * Q[i,j] / sqrt(pi[j])
    std::vector<double> S(n * n);
    for (int i = 0; i < n; ++i) {
        const double si = sqrtpi[i];
        for (int j = 0; j < n; ++j) {
            S[i*n + j] = si * Q_rowmajor[i*n + j] / sqrtpi[j];
        }
    }

    // 3) 對稱特徵分解：S = U * Λ * U^T
    //    LAPACKE_dsyevd：Divide & Conquer，速度/穩定性俱佳
    //    Row-major + 上三角（'U'）慣例
    for (int i=0;i<n;++i){
        for (int j=i+1;j<n;++j){
            double sym = 0.5 * (S[i*n + j] + S[j*n + i]);
            S[i*n + j] = sym;
            S[j*n + i] = sym;
        }
    }
    std::vector<double> lambdas(n);
    // int info = LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'V', 'U', n, S.data(), n, lambdas.data());

    std::vector<double> A(n * n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            // row-major S(i,j) -> column-major A(i,j)
            A[j * n + i] = S[i * n + j];
        }
    }

    char jobz  = 'V';  // 計算 eigenvalues + eigenvectors
    char uplo  = 'U';  // 用上三角（跟原本 LAPACKE 呼叫一致）
    int  N     = n;
    int  lda   = n;
    int  info  = 0;

    // workspace query
    int    lwork  = -1;
    int    liwork = -1;
    double work_query;
    int    iwork_query;

    dsyevd_(&jobz, &uplo, &N,
        A.data(), &lda,
        lambdas.data(),
        &work_query, &lwork,
        &iwork_query, &liwork,
        &info);
    lwork  = static_cast<int>(work_query);   // ✔ 賦值，不是重新宣告
    liwork = iwork_query;      
    std::vector<double> work(lwork);
    std::vector<int>    iwork(liwork);

    dsyevd_(&jobz, &uplo, &N,
        A.data(), &lda,
        lambdas.data(),
        work.data(), &lwork,
        iwork.data(), &liwork,
        &info);
    if (info != 0) throw std::runtime_error("dsyevd failed, info=" + std::to_string(info));
    // 現在 S 裡面存的是 U（row-major；第 j 欄是第 j 個特徵向量）
    std::vector<double> U(n*n);
    for (int i=0;i<n;++i)
        for (int j=0;j<n;++j)
            U[i*n + j] = A[j*n + i];   // U(i,j) = S^T(i,j)
    // 4) 構造 V = D^{-1} U、Vinv = U^T D
    EigResult out;
    out.lambdas = std::move(lambdas);
    out.V.resize(n * n);
    out.Vinv.resize(n * n);
    out.U = U;  // <= 這行把 U 留存起來，給 pmatrix_direct 使用

    // V[i,j] = U[i,j] / sqrt(pi[i])
    for (int i = 0; i < n; ++i) {
        const double invsqrt = 1.0 / sqrtpi[i];
        for (int j = 0; j < n; ++j) {
            out.V[i*n + j] = U[i*n + j] * invsqrt;
        }
    }
    // Vinv[i,j] = U^T D = U(j,i) * sqrt(pi[j])
    for (int i = 0; i < n; ++i) {
        const double sp = sqrtpi[i];
        for (int j = 0; j < n; ++j) {
            out.Vinv[j*n + i] = U[i*n + j] * sp; // 這行等價於 Vinv[j,i] = U(i,j)*sqrtpi[i]
        }
    }

    // 5) 數值修正（可選）：把極小正值夾到 0，避免後續 exp 產生意外
    //    Rate matrix 的特徵值應該 ≤ 0 且有一個 0（對應穩態）
    for (int j = 0; j < n; ++j) {
        if (out.lambdas[j] > -1e-15 && out.lambdas[j] < 1e-15) {
            out.lambdas[j] = 0.0;
        }
    }

    return out;
}
