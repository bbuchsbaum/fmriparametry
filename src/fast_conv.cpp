#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix fast_batch_convolution_cpp(NumericVector signal, NumericMatrix kernels, int output_length) {
  int signal_length = signal.size();
  int kernel_length = kernels.nrow();
  int n_kernels = kernels.ncol();
  
  NumericMatrix result(output_length, n_kernels);
  
  // Simple nested loop implementation
  // This implements the convolution y[i] = sum_k x[i-k] * h[k]
  // where h is already reversed (as expected by stats::convolve with rev(kernel))
  for(int j = 0; j < n_kernels; ++j) {
    for(int i = 0; i < output_length; ++i) {
      double sum = 0.0;
      for(int k = 0; k < kernel_length; ++k) {
        int idx = i - k;
        if(idx >= 0 && idx < signal_length) {
          sum += signal[idx] * kernels(k, j);
        }
      }
      result(i, j) = sum;
    }
  }
  
  return result;
}