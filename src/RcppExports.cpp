#include <Rcpp.h>
#include <RcppParallel.h>
#include <RcppEigen.h>

using namespace Rcpp;

// Declarations
NumericMatrix fast_batch_convolution_cpp(NumericVector signal, NumericMatrix kernels, int output_length);
Eigen::MatrixXd ridge_linear_solve_cpp(const Eigen::Map<Eigen::MatrixXd> &X,
                                       const Eigen::Map<Eigen::MatrixXd> &Y,
                                       double lambda);

// Wrappers
extern "C" SEXP _fmriparametric_fast_batch_convolution_cpp(SEXP signalSEXP, SEXP kernelsSEXP, SEXP output_lengthSEXP) {
  BEGIN_RCPP
  NumericVector signal(signalSEXP);
  NumericMatrix kernels(kernelsSEXP);
  int output_length = as<int>(output_lengthSEXP);
  NumericMatrix res = fast_batch_convolution_cpp(signal, kernels, output_length);
  return res;
  END_RCPP
}

extern "C" SEXP _fmriparametric_ridge_linear_solve_cpp(SEXP XSEXP, SEXP YSEXP, SEXP lambdaSEXP) {
  BEGIN_RCPP
  Eigen::Map<Eigen::MatrixXd> X = as<Eigen::Map<Eigen::MatrixXd> >(XSEXP);
  Eigen::Map<Eigen::MatrixXd> Y = as<Eigen::Map<Eigen::MatrixXd> >(YSEXP);
  double lambda = as<double>(lambdaSEXP);
  Eigen::MatrixXd res = ridge_linear_solve_cpp(X, Y, lambda);
  return Rcpp::wrap(res);
  END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
  {"_fmriparametric_fast_batch_convolution_cpp", (DL_FUNC) &_fmriparametric_fast_batch_convolution_cpp, 3},
  {"_fmriparametric_ridge_linear_solve_cpp", (DL_FUNC) &_fmriparametric_ridge_linear_solve_cpp, 3},
  {NULL, NULL, 0}
};

extern "C" void R_init_fmriparametric(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
