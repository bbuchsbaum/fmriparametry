#' Construct a parametric_hrf_fit object (Sprint 2 enhanced version)
#'
#' Creates a new S3 object storing results from parametric HRF estimation.
#' This enhanced constructor includes additional fields for diagnostics and
#' standard errors introduced in Sprint 2.
#'
#' @param estimated_parameters numeric matrix of parameter estimates (voxels x parameters)
#' @param amplitudes numeric vector of fitted amplitudes
#' @param parameter_names character vector naming the parameters
#' @param parametric_model character string identifying the HRF model
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
  parametric_model = "lwu",
  r_squared = NULL,
  residuals = NULL,
  parameter_ses = NULL,
  convergence_info = list(),
  metadata = list(),
  hrf_shape = NULL
) {
  # Basic validation (as before)
  assertthat::assert_that(is.matrix(estimated_parameters))
  assertthat::assert_that(is.numeric(amplitudes))
  assertthat::assert_that(nrow(estimated_parameters) == length(amplitudes))
  assertthat::assert_that(is.character(parameter_names))
  assertthat::assert_that(ncol(estimated_parameters) == length(parameter_names))
  assertthat::assert_that(is.character(parametric_model), length(parametric_model) == 1)
  
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
    method_used = "parametric_taylor",
    design_info = list(
      n_time = if(!is.null(residuals)) nrow(residuals) else NA_integer_,
      n_vox = n_vox,
      n_cond = NA_integer_,
      basis_dim = ncol(estimated_parameters),
      projected = NA
    ),
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
  
  # Construct object with enhanced fields using new structure
  obj <- list(
    # New standardized structure
    model_specific = list(
      parameters = estimated_parameters,
      parameter_names = parameter_names,
      model = parametric_model,
      standard_errors = parameter_ses
    ),
    
    # Pre-computed HRF shapes
    hrf_shape = hrf_shape,
    
    # Top-level fields
    amplitudes = as.numeric(amplitudes),
    r_squared = r_squared,
    residuals = residuals,
    convergence_info = convergence_info,
    metadata = metadata
  )

  # Set convenient top-level access
  obj$parameter_names <- obj$model_specific$parameter_names
  obj$parametric_model <- obj$model_specific$model
  obj$design_info <- metadata$design_info
  obj$gof_per_voxel <- r_squared
  obj$method_used <- metadata$method_used
  
  # Add backward compatibility field
  obj$estimated_parameters <- obj$model_specific$parameters
  
  # Add fields expected by tests but not included in constructor
  obj$standard_errors <- parameter_ses  # Use the passed value or NULL
  obj$se_amplitudes <- NULL
  obj$fit_quality <- NULL
  obj$refinement_info <- list(applied = FALSE)  # Always present, non-NULL
  
  class(obj) <- "parametric_hrf_fit"
  obj
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

#' Generate HRF shapes from parameters
#' 
#' Pre-computes HRF curves for all voxels for efficient plotting
#' 
#' @param parameters Matrix of HRF parameters (voxels x parameters)
#' @param amplitudes Vector of amplitude values
#' @param times Time points for HRF evaluation
#' @param model Character string specifying HRF model
#' @param hrf_interface Optional pre-created HRF interface
#' @return Matrix of HRF curves (time points x voxels)
#' @keywords internal
.generate_hrf_shape <- function(parameters, amplitudes, times, model, 
                               hrf_interface = NULL) {
  if (is.null(hrf_interface)) {
    hrf_interface <- .create_hrf_interface(model)
  }
  
  n_vox <- nrow(parameters)
  n_time <- length(times)
  hrf_shape <- matrix(NA_real_, n_time, n_vox)
  
  for (v in seq_len(n_vox)) {
    hrf_vals <- hrf_interface$hrf_function(times, parameters[v, ])
    hrf_shape[, v] <- hrf_vals * amplitudes[v]
  }
  
  hrf_shape
}
