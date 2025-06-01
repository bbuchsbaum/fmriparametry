#' Construct a parametric_hrf_fit object (Sprint 2 enhanced version)
#'
#' Creates a new S3 object storing results from parametric HRF estimation.
#' This enhanced constructor includes additional fields for diagnostics and
#' standard errors introduced in Sprint 2.
#'
#' @param estimated_parameters numeric matrix of parameter estimates (voxels x parameters)
#' @param amplitudes numeric vector of fitted amplitudes
#' @param parameter_names character vector naming the parameters
#' @param hrf_model character string identifying the HRF model
#' @param r_squared numeric vector of R-squared values for each voxel
#' @param residuals numeric matrix of residuals (timepoints x voxels) or NULL
#' @param parameter_ses numeric matrix of standard errors (voxels x parameters) or NULL
#' @param convergence_info list of convergence diagnostics including trajectory
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
  metadata = list()
) {
  # Basic validation (as before)
  assertthat::assert_that(is.matrix(estimated_parameters))
  assertthat::assert_that(is.numeric(amplitudes))
  assertthat::assert_that(nrow(estimated_parameters) == length(amplitudes))
  assertthat::assert_that(is.character(parameter_names))
  assertthat::assert_that(ncol(estimated_parameters) == length(parameter_names))
  assertthat::assert_that(is.character(hrf_model), length(hrf_model) == 1)
  
  # New field validation
  n_vox <- nrow(estimated_parameters)
  if (!is.null(r_squared)) {
    assertthat::assert_that(is.numeric(r_squared), length(r_squared) == n_vox)
  }
  
  if (!is.null(residuals)) {
    assertthat::assert_that(is.matrix(residuals), ncol(residuals) == n_vox)
  }
  
  if (!is.null(parameter_ses)) {
    assertthat::assert_that(is.matrix(parameter_ses),
                            nrow(parameter_ses) == n_vox,
                            ncol(parameter_ses) == length(parameter_names))
  }
  
  assertthat::assert_that(is.list(convergence_info))
  assertthat::assert_that(is.list(metadata))
  
  # Set column names
  colnames(estimated_parameters) <- parameter_names
  if (!is.null(parameter_ses)) {
    colnames(parameter_ses) <- parameter_names
  }
  
  # Handle metadata defaults
  meta_defaults <- list(
    call = NULL,
    n_voxels = n_vox,
    n_timepoints = if(!is.null(residuals)) nrow(residuals) else NA_integer_,
    theta_seed = rep(NA_real_, length(parameter_names)),
    theta_bounds = list(lower = rep(NA_real_, length(parameter_names)),
                        upper = rep(NA_real_, length(parameter_names))),
    recenter_global_passes = 0,
    coeffs = NULL,
    theta_expansion = NULL
  )
  metadata <- utils::modifyList(meta_defaults, metadata)
  if (is.null(metadata$call)) {
    metadata$call <- sys.call(-1)
  }
  
  # Construct object with enhanced fields
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
  
  # For backward compatibility, also include old fields
  obj$convergence <- convergence_info  # Alias for backward compatibility
  
  class(obj) <- "parametric_hrf_fit"
  obj
}

#' Check if a parametric_hrf_fit has Sprint 2 enhancements
#' @param x parametric_hrf_fit object
#' @return Logical indicating if object has Sprint 2 fields
#' @keywords internal
is_v2_fit <- function(x) {
  !is.null(x$r_squared) || !is.null(x$residuals) || !is.null(x$parameter_ses)
}

#' Get fitted values from parametric_hrf_fit
#' @param x parametric_hrf_fit object  
#' @param Y_proj Original projected Y data (required if residuals not stored)
#' @return Matrix of fitted values
#' @keywords internal
get_fitted_values <- function(x, Y_proj = NULL) {
  if (!is.null(x$residuals)) {
    # If we have residuals, fitted = Y - residuals
    if (is.null(Y_proj)) {
      stop("Y_proj required to compute fitted values from residuals")
    }
    return(Y_proj - x$residuals)
  } else {
    stop("Cannot compute fitted values without residuals or design matrix")
  }
}