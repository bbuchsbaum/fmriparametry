#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

// helper to compute convolution of a signal with a kernel where the kernel
// vector is already reversed (as expected by stats::convolve with rev())
static inline void convolve_open_rev(const Eigen::VectorXd &signal,
                                     const Eigen::VectorXd &kernel_rev,
                                     Eigen::VectorXd &out) {
    int n_time = out.size();
    int k_len = kernel_rev.size();
    out.setZero();
    for (int i = 0; i < n_time; ++i) {
        double sum = 0.0;
        for (int k = 0; k < k_len; ++k) {
            int idx = i - k;
            if (idx >= 0 && idx < signal.size()) {
                sum += signal(idx) * kernel_rev(k);
            }
        }
        out(i) += sum;
    }
}

// [[Rcpp::export]]
Rcpp::List compute_standard_errors_bulk_cpp(
        const Rcpp::List &derivs_list,
        const Rcpp::List &hrf_list,
        const Eigen::Map<Eigen::MatrixXd> &Y,
        const Eigen::Map<Eigen::MatrixXd> &S,
        const Eigen::VectorXd &beta0) {
    int n_time = Y.rows();
    int n_vox = Y.cols();
    int n_params = Rcpp::as<Rcpp::NumericMatrix>(derivs_list[0]).ncol();
    int n_hrf = Rcpp::as<Rcpp::NumericMatrix>(derivs_list[0]).nrow();
    int n_reg = S.cols();

    Rcpp::NumericMatrix se_theta_hat(n_vox, n_params);
    Rcpp::NumericVector se_beta0(n_vox);

#pragma omp parallel for
    for (int v = 0; v < n_vox; ++v) {
        Rcpp::NumericMatrix derivs_rcpp = derivs_list[v];
        Rcpp::NumericVector hrf_rcpp = hrf_list[v];
        Eigen::Map<Eigen::MatrixXd> derivs(derivs_rcpp.begin(), n_hrf, n_params);
        Eigen::Map<Eigen::VectorXd> hrf(hrf_rcpp.begin(), n_hrf);

        Eigen::VectorXd x_hrf(n_time);
        x_hrf.setZero();
        Eigen::VectorXd tmp(n_time);
        Eigen::VectorXd hrf_rev = hrf.reverse();
        for (int r = 0; r < n_reg; ++r) {
            Eigen::VectorXd signal = S.col(r);
            convolve_open_rev(signal, hrf_rev, tmp);
            x_hrf += tmp;
        }

        Eigen::MatrixXd X_derivs(n_time, n_params);
        for (int p = 0; p < n_params; ++p) {
            Eigen::VectorXd col(n_time);
            col.setZero();
            Eigen::VectorXd deriv_rev = derivs.col(p).reverse();
            for (int r = 0; r < n_reg; ++r) {
                Eigen::VectorXd signal = S.col(r);
                convolve_open_rev(signal, deriv_rev, tmp);
                col += tmp;
            }
            col *= beta0(v);
            X_derivs.col(p) = col;
        }

        Eigen::VectorXd residuals = Y.col(v) - beta0(v) * x_hrf;
        double sigma2 = residuals.squaredNorm() / std::max(1, n_time - 1);

        Eigen::MatrixXd fisher = (X_derivs.transpose() * X_derivs) / sigma2;
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(fisher, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::VectorXd s = svd.singularValues();
        Eigen::VectorXd inv_s = s;
        for (int i = 0; i < s.size(); ++i) {
            inv_s(i) = (s(i) > 1e-12) ? 1.0 / s(i) : 0.0;
        }
        Eigen::MatrixXd fisher_inv = svd.matrixV() * inv_s.asDiagonal() * svd.matrixU().transpose();

        for (int p = 0; p < n_params; ++p) {
            double val = fisher_inv(p, p);
            if (val < 0.0) val = 0.0;
            se_theta_hat(v, p) = std::sqrt(val);
        }

        double denom = x_hrf.squaredNorm();
        se_beta0(v) = std::sqrt(sigma2 / std::max(1e-12, denom));
    }

    return Rcpp::List::create(Rcpp::Named("se_theta_hat") = se_theta_hat,
                              Rcpp::Named("se_beta0") = se_beta0);
}

