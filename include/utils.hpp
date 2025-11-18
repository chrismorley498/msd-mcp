#pragma once

#include<eigen3/Eigen/Dense>

Eigen::VectorXd polyfit(const std::vector<double>& x,
                        const std::vector<double>& y,
                        int degree)
{
    const int N = x.size();
    const int M = degree + 1;

    Eigen::MatrixXd A(N, M);
    Eigen::VectorXd b(N);

    for (int i = 0; i < N; ++i) {
        double xi = 1.0;
        for (int j = 0; j < M; ++j) {
            A(i, j) = xi;
            xi *= x[i];
        }
        b(i) = y[i];
    }
    return A.colPivHouseholderQr().solve(b);
}

double polyeval(const Eigen::VectorXd coeffs, double x){
    double val=0;
    for(int i=0;i<coeffs.size();++i)
        val += std::pow(x,i) * coeffs(i);
    return val;
}