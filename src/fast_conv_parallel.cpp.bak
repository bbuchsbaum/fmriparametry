#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

struct ConvWorker : public Worker {
  RVector<double> signal;
  RMatrix<double> kernels;
  RMatrix<double> result;
  int output_length;
  int kernel_length;

  ConvWorker(NumericVector signal, NumericMatrix kernels, NumericMatrix result, int output_length)
    : signal(signal), kernels(kernels), result(result), output_length(output_length), kernel_length(kernels.nrow()) {}

  void operator()(size_t begin, size_t end) {
    for(size_t j = begin; j < end; ++j) {
      for(int i = 0; i < output_length; ++i) {
        double sum = 0.0;
        for(int k = 0; k < kernel_length; ++k) {
          int idx = i + k;
          if(idx < signal.length()) {
            sum += signal[idx] * kernels(kernel_length - 1 - k, j);
          }
        }
        result(i, j) = sum;
      }
    }
  }
};

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::export]]
NumericMatrix fast_batch_convolution_cpp(NumericVector signal, NumericMatrix kernels, int output_length) {
  int n_kernels = kernels.ncol();
  NumericMatrix result(output_length, n_kernels);
  ConvWorker worker(signal, kernels, result, output_length);
  parallelFor(0, n_kernels, worker);
  return result;
}
