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
  
  # DEBUG: Check inputs (commented out for production)
  debug_mode <- getOption("fmriparametric.debug", FALSE)
  # if (debug_mode || all(abs(signals) < 1e-10)) {
  #   cat("\n=== BATCH CONVOLUTION DEBUG ===\n")
  #   cat("signals dim:", dim(signals), "\n")
  #   cat("signals range:", range(signals), "\n")
  #   cat("signals sum:", sum(abs(signals)), "\n")
  #   cat("kernels dim:", dim(kernels), "\n")
  #   cat("kernels[,1] (HRF) range:", range(kernels[,1]), "\n")
  #   cat("output_length:", output_length, "\n")
  # }
  
  # Optimize: If summing signals, convolution is linear so we can sum first
  # This is much faster than looping
  if (ncol(signals) > 1) {
    # Sum all signal columns first
    summed_signal <- rowSums(signals)
    
    if (debug_mode) {
      cat("Summed signal range:", range(summed_signal), "\n")
    }
    
    # Single convolution with summed signal
    result <- .fast_batch_convolution(summed_signal, kernels, output_length)
  } else {
    # Single signal column - direct convolution
    result <- .fast_batch_convolution(signals[, 1], kernels, output_length)
  }
  
  if (is.null(result)) {
    stop(".batch_convolution: .fast_batch_convolution returned NULL")
  }
  if (!is.matrix(result)) {
    stop(".batch_convolution: .fast_batch_convolution returned non-matrix")
  }
  
  # if (debug_mode || all(abs(signals) < 1e-10)) {
  #   cat("Final result range:", range(result), "\n")
  #   cat("================================\n\n")
  # }
  
  # Ensure we return a matrix
  if (!is.matrix(result)) {
    stop(".batch_convolution: result is not a matrix after processing")
  }
  
  return(result)
}
