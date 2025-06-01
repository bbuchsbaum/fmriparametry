#' Construct a parametric_hrf_fit object
#'
#' Creates a new S3 object storing results from parametric HRF estimation.
#' This constructor validates that all required fields are present and
#' returns an object of class `"parametric_hrf_fit"`.
#'
#' @param estimated_parameters numeric matrix of parameter estimates
#' @param amplitudes numeric vector of fitted amplitudes
#' @param parameter_names character vector naming the parameters
#' @param hrf_model character string identifying the HRF model
#' @param r_squared optional numeric vector of R-squared values
#' @param residuals optional numeric matrix of residuals
#' @param parameter_ses optional numeric matrix of parameter standard errors
#' @param convergence_info list of convergence diagnostics
#' @param metadata list containing additional metadata such as the call,
#'   number of voxels and time points, the parameter seed and bounds
#'
#' @return An object of class `parametric_hrf_fit`
#' @keywords internal
new_parametric_hrf_fit <- function(
  estimated_parameters,
  amplitudes,
  parameter_names,
  hrf_model = "lwu",
  r_squared = NULL,
  residuals = NULL,
  parameter_ses = NULL,
  convergence_info = list(),
  metadata = list(),
  ...  # Ignore any extra arguments for backward compatibility
) {
  assertthat::assert_that(is.matrix(estimated_parameters))
  assertthat::assert_that(is.numeric(amplitudes))
  assertthat::assert_that(nrow(estimated_parameters) == length(amplitudes))
  assertthat::assert_that(is.character(parameter_names))
  assertthat::assert_that(ncol(estimated_parameters) == length(parameter_names))
  assertthat::assert_that(is.character(hrf_model), length(hrf_model) == 1)
  if (!is.null(r_squared)) {
    assertthat::assert_that(is.numeric(r_squared),
                            length(r_squared) == nrow(estimated_parameters))
  }
  if (!is.null(residuals)) {
    assertthat::assert_that(is.matrix(residuals),
                            ncol(residuals) == nrow(estimated_parameters))
  }
  if (!is.null(parameter_ses)) {
    assertthat::assert_that(is.matrix(parameter_ses),
                            nrow(parameter_ses) == nrow(estimated_parameters),
                            ncol(parameter_ses) == length(parameter_names))
  }
  assertthat::assert_that(is.list(convergence_info))
  assertthat::assert_that(is.list(metadata))

  meta_defaults <- list(
    call = NULL,
    n_voxels = nrow(estimated_parameters),
    n_timepoints = NA_integer_,
    theta_seed = rep(NA_real_, length(parameter_names)),
    theta_bounds = list(lower = rep(NA_real_, length(parameter_names)),
                        upper = rep(NA_real_, length(parameter_names)))
  )
  metadata <- utils::modifyList(meta_defaults, metadata)
  if (is.null(metadata$call)) {
    metadata$call <- sys.call(-1)
  }

  obj <- list(
    estimated_parameters = estimated_parameters,
    amplitudes = as.numeric(amplitudes),
    parameter_names = parameter_names,
    hrf_model = hrf_model,
    r_squared = r_squared,
    residuals = residuals,
    parameter_ses = parameter_ses,
    convergence_info = convergence_info,
    metadata = metadata
  )
  # Backward compatibility aliases
  obj$parameters <- obj$estimated_parameters
  obj$convergence <- obj$convergence_info
  class(obj) <- "parametric_hrf_fit"
  obj
}

#' Number of voxels in a fit object
#' @param x parametric_hrf_fit
#' @keywords internal
n_voxels <- function(x) {
  x$metadata$n_voxels
}

#' Number of time points used during fitting
#' @param x parametric_hrf_fit
#' @keywords internal
n_timepoints <- function(x) {
  x$metadata$n_timepoints
}
