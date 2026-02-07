#' High-performance batch convolution
#'
#' Performs convolution of a signal with multiple kernels using an
#' FFT-based method for large problems and a parallel C++ fallback for
#' smaller cases. This helper is exported but considered internal.
#'
#' @param signal Numeric vector signal to convolve.
#' @param kernels Numeric matrix with one kernel per column.
#' @param output_length Integer length of output time series.
#' @param conv_context Optional precomputed convolution context from
#'   `.prepare_fft_convolution_context()`.
#' @return Matrix of convolved signals (time x kernels).
#' @keywords internal
#' @export
.prepare_fft_convolution_context <- function(signal, output_length, kernel_length) {
  if (is.null(dim(signal))) {
    signal <- as.numeric(signal)
  } else if (is.matrix(signal) && ncol(signal) > 1L) {
    signal <- rowSums(signal)
  } else {
    signal <- as.numeric(signal[, 1])
  }

  use_fft <- (output_length > 200 && kernel_length > 20)
  if (!use_fft) {
    return(list(use_fft = FALSE))
  }

  n_fft <- nextn(output_length + kernel_length - 1L, factors = 2)
  signal_fft <- fft(c(signal, rep(0, n_fft - length(signal))))

  list(
    use_fft = TRUE,
    n_fft = n_fft,
    signal_fft = signal_fft,
    signal_length = length(signal),
    output_length = output_length,
    kernel_length = kernel_length
  )
}

.fast_batch_convolution <- function(signal, kernels, output_length, conv_context = NULL) {
  n_kernels <- ncol(kernels)
  kernel_length <- nrow(kernels)

  use_fft <- (output_length > 200 && kernel_length > 20)

  if (use_fft) {
    use_cached_fft <- !is.null(conv_context) &&
      isTRUE(conv_context$use_fft) &&
      identical(conv_context$output_length, output_length) &&
      identical(conv_context$kernel_length, kernel_length) &&
      identical(conv_context$signal_length, length(signal)) &&
      !is.null(conv_context$signal_fft) &&
      !is.null(conv_context$n_fft)

    if (use_cached_fft) {
      n_fft <- conv_context$n_fft
      signal_fft <- conv_context$signal_fft
    } else {
      n_fft <- nextn(output_length + kernel_length - 1L, factors = 2)
      signal_fft <- fft(c(signal, rep(0, n_fft - length(signal))))
    }

    kernels_padded <- matrix(0, nrow = n_fft, ncol = n_kernels)
    kernels_padded[seq_len(kernel_length), ] <- kernels
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
