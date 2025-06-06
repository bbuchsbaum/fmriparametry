// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// fast_batch_convolution_cpp
NumericMatrix fast_batch_convolution_cpp(NumericVector signal, NumericMatrix kernels, int output_length);
RcppExport SEXP _fmriparametric_fast_batch_convolution_cpp(SEXP signalSEXP, SEXP kernelsSEXP, SEXP output_lengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type signal(signalSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type kernels(kernelsSEXP);
    Rcpp::traits::input_parameter< int >::type output_length(output_lengthSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_batch_convolution_cpp(signal, kernels, output_length));
    return rcpp_result_gen;
END_RCPP
}
// ridge_linear_solve_cpp
Eigen::MatrixXd ridge_linear_solve_cpp(const Eigen::Map<Eigen::MatrixXd>& X, const Eigen::Map<Eigen::MatrixXd>& Y, double lambda);
RcppExport SEXP _fmriparametric_ridge_linear_solve_cpp(SEXP XSEXP, SEXP YSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(ridge_linear_solve_cpp(X, Y, lambda));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fmriparametric_fast_batch_convolution_cpp", (DL_FUNC) &_fmriparametric_fast_batch_convolution_cpp, 3},
    {"_fmriparametric_ridge_linear_solve_cpp", (DL_FUNC) &_fmriparametric_ridge_linear_solve_cpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_fmriparametric(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
