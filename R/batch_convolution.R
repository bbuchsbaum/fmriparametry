#' Batch convolution helper
#'
#' Performs convolution of multiple stimulus regressors with a set of
#' kernel columns and sums the results. Used by the parametric engine
#' to build design matrices efficiently.
#'
#' @param signals Matrix of stimulus regressors (time x regressors)
#' @param kernels Matrix of kernels/basis functions (time x kernels)
#' @param output_length Length of the desired output time series
#' @keywords internal
.batch_convolution <- function(signals, kernels, output_length) {
  # Ensure signals is a matrix
  if (is.null(dim(signals))) {
    signals <- matrix(signals, ncol = 1)
  }
  
  # Validate inputs
  if (!is.matrix(kernels)) {
    stop(".batch_convolution: kernels must be a matrix")
  }
  if (ncol(kernels) == 0 || nrow(kernels) == 0) {
    stop(".batch_convolution: kernels matrix is empty")
  }
  
  # Initialize result matrix
  result <- matrix(0, nrow = output_length, ncol = ncol(kernels))
  
  # Convolve each signal with kernels and sum
  for (s in seq_len(ncol(signals))) {
    conv_mat <- .fast_batch_convolution(signals[, s], kernels, output_length)
    if (is.null(conv_mat)) {
      stop(sprintf(".batch_convolution: .fast_batch_convolution returned NULL for signal %d", s))
    }
    if (!is.matrix(conv_mat)) {
      stop(sprintf(".batch_convolution: .fast_batch_convolution returned non-matrix for signal %d", s))
    }
    result <- result + conv_mat
  }
  
  # Ensure we return a matrix
  if (!is.matrix(result)) {
    stop(".batch_convolution: result is not a matrix after processing")
  }
  
  return(result)
}
