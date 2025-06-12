#' High-performance batch convolution
#'
#' Performs convolution of a signal with multiple kernels using an
#' FFT-based method for large problems and a parallel C++ fallback for
#' smaller cases. This helper is exported but considered internal.
#'
#' @param signal Numeric vector signal to convolve.
#' @param kernels Numeric matrix with one kernel per column.
#' @param output_length Integer length of output time series.
#' @return Matrix of convolved signals (time x kernels).
#' @keywords internal
#' @export
.fast_batch_convolution <- function(signal, kernels, output_length) {
  n_kernels <- ncol(kernels)
  kernel_length <- nrow(kernels)

  use_fft <- (output_length > 200 && kernel_length > 20)

  if (use_fft) {
    n_fft <- nextn(output_length + kernel_length - 1, factors = 2)
    signal_fft <- fft(c(signal, rep(0, n_fft - length(signal))))
    kernels_padded <- rbind(kernels, matrix(0, n_fft - kernel_length, n_kernels))
    kernels_fft <- mvfft(kernels_padded)
    conv_fft <- signal_fft * kernels_fft
    conv_time <- mvfft(conv_fft, inverse = TRUE) / n_fft
    return(Re(conv_time[seq_len(output_length), , drop = FALSE]))
  } else {
    # Ensure inputs are double-precision for C++
    storage.mode(signal) <- "double"
    storage.mode(kernels) <- "double"
    return(fast_batch_convolution_cpp(signal, kernels, output_length))
  }
}

