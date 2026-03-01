#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
Eigen::MatrixXd ridge_linear_solve_cpp(const Eigen::Map<Eigen::MatrixXd> &X,
                                       const Eigen::Map<Eigen::MatrixXd> &Y,
                                       double lambda) {
  if (!std::isfinite(lambda) || lambda < 0.0) {
    Rcpp::stop("lambda must be finite and non-negative");
  }
  if (X.rows() != Y.rows()) {
    Rcpp::stop("X and Y must have the same number of rows");
  }

  const int n = X.rows();
  const int p = X.cols();
  const int v = Y.cols();
  const bool use_ridge = (lambda > 0.0);
  const int n_aug = use_ridge ? (n + p) : n;

  Eigen::MatrixXd X_aug(n_aug, p);
  X_aug.topRows(n) = X;
  if (use_ridge) {
    X_aug.bottomRows(p).setZero();
    X_aug.bottomRows(p).diagonal().array() = std::sqrt(lambda);
  }

  Eigen::MatrixXd Y_aug(n_aug, v);
  Y_aug.topRows(n) = Y;
  if (use_ridge) {
    Y_aug.bottomRows(p).setZero();
  }

  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(X_aug);
  if (qr.rank() < p) {
    Rcpp::stop("ridge_linear_solve_cpp: design matrix is rank deficient");
  }

  Eigen::MatrixXd coeffs = qr.solve(Y_aug);
  if (!coeffs.allFinite()) {
    Rcpp::stop("ridge_linear_solve_cpp: failed to compute finite coefficients");
  }
  return coeffs;
}
