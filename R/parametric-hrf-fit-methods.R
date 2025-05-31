#' Print a parametric_hrf_fit object
#'
#' Displays a concise summary of the fit including the HRF model,
#' number of voxels and basic parameter statistics.
#'
#' @param x An object of class `parametric_hrf_fit`.
#' @param ... Additional arguments passed to methods. Ignored.
#' @return The input object `x` invisibly.
#' @export
print.parametric_hrf_fit <- function(x, ...) {
  cat("Parametric HRF Fit\n")
  cat("Model:", x$hrf_model, "\n")
  cat("Voxels:", nrow(x$estimated_parameters), "\n")
  cat("\nParameter Summary:\n")
  print(summary(x$estimated_parameters))
  invisible(x)
}

#' Extract coefficients from a parametric_hrf_fit object
#'
#' Returns the matrix of estimated HRF parameters for each voxel.
#'
#' @param object An object of class `parametric_hrf_fit`.
#' @param ... Additional arguments passed to methods. Ignored.
#' @return Numeric matrix of parameter estimates.
#' @export
coef.parametric_hrf_fit <- function(object, ...) {
  object$estimated_parameters
}

#' Summarize a parametric_hrf_fit object
#'
#' Produces summary statistics of the estimated parameters and amplitudes.
#'
#' @param object An object of class `parametric_hrf_fit`.
#' @param ... Additional arguments passed to methods. Ignored.
#' @return A list with parameter and amplitude summaries.
#' @export
summary.parametric_hrf_fit <- function(object, ...) {
  param_sum <- apply(object$estimated_parameters, 2, summary)
  amp_sum <- summary(object$amplitudes)
  list(
    parameter_summary = param_sum,
    amplitude_summary = amp_sum,
    hrf_model = object$hrf_model,
    n_voxels = n_voxels(object)
  )
}
