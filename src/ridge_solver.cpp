#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
Eigen::MatrixXd ridge_linear_solve_cpp(const Eigen::Map<Eigen::MatrixXd> &X,
                                       const Eigen::Map<Eigen::MatrixXd> &Y,
                                       double lambda) {
  const int p = X.cols();
  const int v = Y.cols();

  Eigen::MatrixXd XtX = X.transpose() * X;
  if (lambda > 0) {
    XtX.diagonal().array() += lambda;
  }
  Eigen::LDLT<Eigen::MatrixXd> solver(XtX);

  Eigen::MatrixXd XtY = X.transpose() * Y;
  Eigen::MatrixXd coeffs(p, v);

  #pragma omp parallel for if(v > 1)
  for (int i = 0; i < v; ++i) {
    coeffs.col(i) = solver.solve(XtY.col(i));
  }

  return coeffs;
}
