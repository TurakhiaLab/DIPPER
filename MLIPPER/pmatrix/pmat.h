#pragma once
#include <vector>
#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
#include <algorithm>
#include <cassert>

struct EigResult {
    std::vector<double> lambdas;
    std::vector<double> V;
    std::vector<double> Vinv;
    std::vector<double> U;   // <= 新增：保存 U，供直接公式使用
};

EigResult gtr_eigendecomp_cpu(
    const double* Q_rowmajor,   // Q（row-major，大小 n*n）
    const double* pi,           // 平衡頻率（長度 n，sum=1）
    int n);

void pmatrix_from_triple(const double* Vinv, const double* V,
                                const double* lamb, double r, double t, double p,
                                double* P, int n);

void pmatrix_direct(const double* U, const double* pi,
                           const double* lamb, double r, double t, double p,
                           double* P, int n);