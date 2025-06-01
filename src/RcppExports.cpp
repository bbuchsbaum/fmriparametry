#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;

// Declaration of C++ function
NumericMatrix fast_batch_convolution_cpp(NumericVector signal, NumericMatrix kernels, int output_length);

// Wrapper to call from R
extern "C" SEXP _fmriparametric_fast_batch_convolution_cpp(SEXP signalSEXP, SEXP kernelsSEXP, SEXP output_lengthSEXP) {
  BEGIN_RCPP
  NumericVector signal(signalSEXP);
  NumericMatrix kernels(kernelsSEXP);
  int output_length = as<int>(output_lengthSEXP);
  NumericMatrix res = fast_batch_convolution_cpp(signal, kernels, output_length);
  return res;
  END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
  {"_fmriparametric_fast_batch_convolution_cpp", (DL_FUNC) &_fmriparametric_fast_batch_convolution_cpp, 3},
  {NULL, NULL, 0}
};

extern "C" void R_init_fmriparametric(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

