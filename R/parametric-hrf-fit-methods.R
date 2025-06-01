#' Print a parametric_hrf_fit object
#'
#' Displays a concise summary of the parametric HRF fit including the HRF model 
#' used, number of voxels analyzed, and basic summary statistics for each 
#' estimated parameter.
#'
#' @param x An object of class \code{parametric_hrf_fit}, typically from 
#'   \code{\link{estimate_parametric_hrf}}.
#' @param ... Additional arguments (currently ignored).
#' 
#' @return The input object \code{x} invisibly.
#' 
#' @examples
#' \dontrun{
#' # Fit a model
#' fit <- estimate_parametric_hrf(fmri_data, event_model)
#' 
#' # Print summary
#' print(fit)
#' # or simply:
#' fit
#' }
#' 
#' @seealso \code{\link{estimate_parametric_hrf}}, \code{\link{summary.parametric_hrf_fit}}
#' 
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
#' Returns the matrix of estimated HRF parameters for each voxel. Each row 
#' corresponds to a voxel and each column to a parameter.
#'
#' @param object An object of class \code{parametric_hrf_fit}, typically from 
#'   \code{\link{estimate_parametric_hrf}}.
#' @param ... Additional arguments (currently ignored).
#' 
#' @return A numeric matrix with dimensions (n_voxels × n_parameters). For the 
#'   LWU model, columns are named \code{"tau"}, \code{"sigma"}, and \code{"rho"}.
#' 
#' @examples
#' \dontrun{
#' # Fit a model
#' fit <- estimate_parametric_hrf(fmri_data, event_model)
#' 
#' # Extract parameters
#' params <- coef(fit)
#' 
#' # Get parameters for specific voxel
#' voxel_100_params <- params[100, ]
#' 
#' # Summary statistics for each parameter
#' apply(params, 2, summary)
#' }
#' 
#' @seealso \code{\link{estimate_parametric_hrf}}, \code{\link{summary.parametric_hrf_fit}}
#' 
#' @export
coef.parametric_hrf_fit <- function(object, ...) {
  object$estimated_parameters
}

#' Summarize a parametric_hrf_fit object
#'
#' Produces comprehensive summary statistics of the estimated HRF parameters and 
#' response amplitudes across all voxels.
#'
#' @param object An object of class \code{parametric_hrf_fit}, typically from 
#'   \code{\link{estimate_parametric_hrf}}.
#' @param ... Additional arguments (currently ignored).
#' 
#' @return A list containing:
#'   \item{parameter_summary}{Matrix with summary statistics (min, 1st quartile, 
#'     median, mean, 3rd quartile, max) for each parameter}
#'   \item{amplitude_summary}{Vector of summary statistics for response amplitudes}
#'   \item{hrf_model}{Character string identifying the HRF model used}
#'   \item{n_voxels}{Integer number of voxels analyzed}
#' 
#' @examples
#' \dontrun{
#' # Fit a model
#' fit <- estimate_parametric_hrf(fmri_data, event_model)
#' 
#' # Get summary
#' summ <- summary(fit)
#' 
#' # View parameter summaries
#' summ$parameter_summary
#' 
#' # Check amplitude distribution
#' summ$amplitude_summary
#' }
#' 
#' @seealso \code{\link{estimate_parametric_hrf}}, \code{\link{coef.parametric_hrf_fit}}
#' 
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
