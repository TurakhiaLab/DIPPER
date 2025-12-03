// test_gtr_pmatrix.cpp
#include <vector>
#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
#include <algorithm>
#include <cassert>

#include "pmat.h"   // 你的檔案；需提供 EigResult { lambdas,V,Vinv,U }

// ---------- 小工具 ----------

// ---------- 隨機產生合法的 GTR Q ----------
static std::vector<double> make_Q_GTR(const std::vector<double>& pi, std::mt19937_64& rng){
    const int n = (int)pi.size();
    std::uniform_real_distribution<double> U(0.05, 2.0);
    // symmetric exchange rates a_ij
    std::vector<double> a(n*n, 0.0);
    for(int i=0;i<n;++i){
        for(int j=i+1;j<n;++j){
            double aij = U(rng);
            a[i*n + j] = aij;
            a[j*n + i] = aij;
        }
    }
    // Q_ij = a_ij * pi_j ; Q_ii = -sum_{j!=i} Q_ij
    std::vector<double> Q(n*n, 0.0);
    for(int i=0;i<n;++i){
        double rowsum = 0.0;
        for(int j=0;j<n;++j){
            if (i==j) continue;
            double qij = a[i*n + j] * pi[j];
            Q[i*n + j] = qij;
            rowsum += qij;
        }
        Q[i*n + i] = -rowsum;
    }
    // scale to mean rate = 1: sum_i pi_i * (-Q_ii) = 1
    double rate = 0.0;
    for(int i=0;i<n;++i) rate += pi[i] * (-Q[i*n + i]);
    for(double& x : Q) x /= rate;

    // sanity: row sums ~ 0
    for(int i=0;i<n;++i){
        double s=0.0;
        for(int j=0;j<n;++j) s += Q[i*n + j];
        assert(std::abs(s) < 1e-12);
    }
    return Q;
}


int main(int argc, char** argv){
    std::cout.setf(std::ios::scientific);
    std::cout << std::setprecision(3);

    const int n = 4;           // 4 or 20
    const int rate_cats = 4;
    const int NUM_TEST = 500;  // 可調大到 1000

    // 產生隨機 pi
    std::random_device rd;
    std::mt19937_64 rng(rd());
    std::uniform_real_distribution<double> U01(0.0, 1.0);
    std::uniform_real_distribution<double> Urate(0.2, 2.0);
    std::uniform_real_distribution<double> Utlen(0.0, 2.0);
    std::uniform_real_distribution<double> Upinv(0.0, 0.2);

    std::vector<double> pi(n);
    {
        double s=0.0;
        for(int i=0;i<n;++i){ pi[i] = U01(rng); s += pi[i]; }
        for(int i=0;i<n;++i) pi[i] /= s;
    }

    // 產生 Q，並 sanity check 列和
    auto Q = make_Q_GTR(pi, rng);

    // 分解
    EigResult ER = gtr_eigendecomp_cpu(Q.data(), pi.data(), n);

    // 檢查 Vinv*V ≈ I
    {
        std::vector<double> I(n*n), VV(n*n);
        matmul_nn(ER.Vinv.data(), ER.V.data(), VV.data(), n);
        eye(I.data(), n);
        for(int i=0;i<n*n;++i) VV[i] -= I[i];
        double vinv_v_err = inf_norm(VV, n);
        std::cout << "[CHK] ||Vinv*V - I||_inf = " << vinv_v_err << "\n";
    }

    // rates
    std::vector<double> rates(rate_cats);
    for (int i=0;i<rate_cats;++i) rates[i] = Urate(rng);

    double max_diff_all = 0.0, mean_max_diff = 0.0;
    double worst_rowsum_before = 0.0, worst_rowsum_after = 0.0;

    for(int it=0; it<NUM_TEST; ++it){
        double r = rates[it % rate_cats];
        double t = Utlen(rng);
        double p = Upinv(rng);

        std::vector<double> P1(n*n), P2(n*n);
        pmatrix_from_triple(ER.Vinv.data(), ER.V.data(), ER.lambdas.data(),
                            r, t, p, P1.data(), n);
        pmatrix_direct(ER.U.data(), pi.data(), ER.lambdas.data(),
                       r, t, p, P2.data(), n);

        // diff
        double md = 0.0;
        for(int i=0;i<n*n;++i) md = std::max(md, std::abs(P1[i]-P2[i]));
        max_diff_all = std::max(max_diff_all, md);
        mean_max_diff += md;

        // rowsum before fix
        double dev1 = check_rowsum_max_dev(P1.data(), n);
        double dev2 = check_rowsum_max_dev(P2.data(), n);
        worst_rowsum_before = std::max(worst_rowsum_before, std::max(dev1, dev2));

        // optional: clamp small negatives + renorm (工程保險絲)
        clamp_neg_and_row_norm(P1.data(), n);
        clamp_neg_and_row_norm(P2.data(), n);

        // rowsum after fix
        double dev1a = check_rowsum_max_dev(P1.data(), n);
        double dev2a = check_rowsum_max_dev(P2.data(), n);
        worst_rowsum_after = std::max(worst_rowsum_after, std::max(dev1a, dev2a));
    }

    mean_max_diff /= NUM_TEST;

    std::cout << "[RES] tests=" << NUM_TEST
              << "  max_abs_diff=" << max_diff_all
              << "  mean_max_abs_diff=" << mean_max_diff << "\n";
    std::cout << "[ROW] max dev before fix = " << worst_rowsum_before
              << " ; after fix = " << worst_rowsum_after << "\n";

    // 期望：double 精度下 max_abs_diff 通常在 1e-12 ~ 1e-10
    // 若你的 BLAS/LAPACK 實作或 pi 很極端，可能更鬆一些
    if (max_diff_all > 1e-9) {
        std::cerr << "[FAIL] difference too large.\n";
        return 1;
    }
    return 0;
}
