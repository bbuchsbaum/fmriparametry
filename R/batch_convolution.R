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
  if (is.null(dim(signals))) {
    signals <- matrix(signals, ncol = 1)
  }
  result <- matrix(0, nrow = output_length, ncol = ncol(kernels))
  for (s in seq_len(ncol(signals))) {
    conv_mat <- .fast_batch_convolution(signals[, s], kernels, output_length)
    result <- result + conv_mat
  }
  result
}
