#include <RcppEigen.h>
#include <limits>
#include <vector>
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
    const int n_time = Y.rows();
    const int n_vox = Y.cols();

    if (n_time <= 0 || n_vox <= 0) {
        Rcpp::stop("Y must be a non-empty matrix");
    }
    if (S.rows() != n_time) {
        Rcpp::stop("S and Y must have matching row counts");
    }
    if (derivs_list.size() != n_vox || hrf_list.size() != n_vox) {
        Rcpp::stop("derivs_list/hrf_list lengths must match number of voxels in Y");
    }
    if (beta0.size() != n_vox) {
        Rcpp::stop("beta0 length must match number of voxels in Y");
    }

    int n_params = -1;
    int n_hrf = -1;

    std::vector<Eigen::MatrixXd> derivs_cache;
    std::vector<Eigen::VectorXd> hrf_cache;
    derivs_cache.reserve(n_vox);
    hrf_cache.reserve(n_vox);

    for (int v = 0; v < n_vox; ++v) {
        Rcpp::NumericMatrix derivs_rcpp = Rcpp::as<Rcpp::NumericMatrix>(derivs_list[v]);
        Rcpp::NumericVector hrf_rcpp = Rcpp::as<Rcpp::NumericVector>(hrf_list[v]);
        const int this_n_hrf = derivs_rcpp.nrow();
        const int this_n_params = derivs_rcpp.ncol();

        if (this_n_hrf <= 0 || this_n_params <= 0) {
            Rcpp::stop("Each derivatives matrix must have positive dimensions");
        }
        if (hrf_rcpp.size() != this_n_hrf) {
            Rcpp::stop("Each HRF vector length must match derivative rows");
        }

        if (v == 0) {
            n_hrf = this_n_hrf;
            n_params = this_n_params;
        } else {
            if (this_n_hrf != n_hrf || this_n_params != n_params) {
                Rcpp::stop("All derivative/HRF entries must share common dimensions");
            }
        }

        Eigen::Map<Eigen::MatrixXd> derivs_map(derivs_rcpp.begin(), this_n_hrf, this_n_params);
        Eigen::Map<Eigen::VectorXd> hrf_map(hrf_rcpp.begin(), this_n_hrf);
        derivs_cache.emplace_back(derivs_map);
        hrf_cache.emplace_back(hrf_map);
    }

    const int n_reg = S.cols();
    std::vector<Eigen::VectorXd> signal_cache;
    signal_cache.reserve(n_reg);
    for (int r = 0; r < n_reg; ++r) {
        signal_cache.emplace_back(S.col(r));
    }

    std::vector<double> se_theta_vals(static_cast<size_t>(n_vox) * static_cast<size_t>(n_params), NA_REAL);
    std::vector<double> se_beta_vals(static_cast<size_t>(n_vox), NA_REAL);

#pragma omp parallel for
    for (int v = 0; v < n_vox; ++v) {
        const Eigen::MatrixXd &derivs = derivs_cache[v];
        const Eigen::VectorXd &hrf = hrf_cache[v];
        const double beta = beta0(v);
        if (!std::isfinite(beta)) {
            continue;
        }

        Eigen::VectorXd x_hrf(n_time);
        x_hrf.setZero();
        Eigen::VectorXd tmp(n_time);
        const Eigen::VectorXd hrf_rev = hrf.reverse();
        for (int r = 0; r < n_reg; ++r) {
            convolve_open_rev(signal_cache[r], hrf_rev, tmp);
            x_hrf += tmp;
        }

        Eigen::MatrixXd X_derivs(n_time, n_params);
        for (int p = 0; p < n_params; ++p) {
            Eigen::VectorXd col(n_time);
            col.setZero();
            const Eigen::VectorXd deriv_rev = derivs.col(p).reverse();
            for (int r = 0; r < n_reg; ++r) {
                convolve_open_rev(signal_cache[r], deriv_rev, tmp);
                col += tmp;
            }
            col *= beta;
            X_derivs.col(p) = col;
        }

        const Eigen::VectorXd residuals = Y.col(v) - beta * x_hrf;
        const int dof = std::max(1, n_time - (n_params + 1));
        const double sigma2 = residuals.squaredNorm() / static_cast<double>(dof);
        if (!std::isfinite(sigma2) || sigma2 <= std::numeric_limits<double>::epsilon()) {
            continue;
        }

        const Eigen::MatrixXd fisher = (X_derivs.transpose() * X_derivs) / sigma2;
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(fisher, Eigen::ComputeFullU | Eigen::ComputeFullV);
        const Eigen::VectorXd s = svd.singularValues();
        if (s.size() == 0) {
            continue;
        }
        const double s_max = s.maxCoeff();
        const double svd_tol = std::numeric_limits<double>::epsilon() *
            static_cast<double>(std::max(fisher.rows(), fisher.cols())) * s_max;

        Eigen::VectorXd inv_s(s.size());
        for (int i = 0; i < s.size(); ++i) {
            inv_s(i) = (s(i) > svd_tol) ? (1.0 / s(i)) : 0.0;
        }
        const Eigen::MatrixXd fisher_inv =
            svd.matrixV() * inv_s.asDiagonal() * svd.matrixU().transpose();

        for (int p = 0; p < n_params; ++p) {
            double val = fisher_inv(p, p);
            if (!std::isfinite(val) || val < 0.0) {
                val = 0.0;
            }
            se_theta_vals[static_cast<size_t>(v) * static_cast<size_t>(n_params) + static_cast<size_t>(p)] =
                std::sqrt(val);
        }

        const double denom = x_hrf.squaredNorm();
        if (std::isfinite(denom) && denom > 1e-12) {
            se_beta_vals[static_cast<size_t>(v)] = std::sqrt(sigma2 / denom);
        }
    }

    Rcpp::NumericMatrix se_theta_hat(n_vox, n_params);
    Rcpp::NumericVector se_beta0(n_vox);
    for (int v = 0; v < n_vox; ++v) {
        se_beta0(v) = se_beta_vals[static_cast<size_t>(v)];
        for (int p = 0; p < n_params; ++p) {
            se_theta_hat(v, p) =
                se_theta_vals[static_cast<size_t>(v) * static_cast<size_t>(n_params) + static_cast<size_t>(p)];
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("se_theta_hat") = se_theta_hat,
        Rcpp::Named("se_beta0") = se_beta0
    );
}
