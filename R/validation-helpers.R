#' Validation helpers for fmriparametric
#'
#' This module provides legacy validation wrappers that delegate to the
#' unified validation core. These functions maintain backward compatibility
#' while using the consolidated validation logic.

#' Validate fMRI data
#' 
#' @param fmri_data fMRI data (matrix or fmri_dataset object)
#' @param caller Name of calling function
#' @return List with validated data properties
#' @keywords internal
.validate_fmri_data <- function(fmri_data, caller = "function") {
  if (is.null(fmri_data)) {
    stop(sprintf("In %s: fmri_data cannot be NULL", caller), call. = FALSE)
  }
  
  # Check if it's an fmri_dataset object
  if (inherits(fmri_data, "fmri_dataset")) {
    # Get dimensions - handle mock objects that might not have proper dim method
    dims <- tryCatch(dim(fmri_data), error = function(e) NULL)
    
    if (is.null(dims)) {
      # For mock objects, try to get dimensions from data component
      if (!is.null(fmri_data$data) && is.matrix(fmri_data$data)) {
        n_time <- nrow(fmri_data$data)
        n_voxels <- ncol(fmri_data$data)
      } else {
        stop(sprintf("In %s: Cannot determine dimensions of fmri_dataset object", caller), 
             call. = FALSE)
      }
    } else {
      n_dims <- length(dims)
      n_time <- dims[n_dims]  # Last dimension is time
      n_voxels <- prod(dims[-n_dims])  # All other dimensions are spatial
    }
    
    return(list(
      data = fmri_data,
      is_fmri_dataset = TRUE,
      n_time = n_time,
      n_voxels = n_voxels
    ))
  }
  
  # Otherwise must be a matrix
  if (!is.matrix(fmri_data)) {
    stop(sprintf("In %s: fmri_data must be a matrix or fmri_dataset object. Got: %s",
                 caller, class(fmri_data)[1]), call. = FALSE)
  }
  
  # Check for too few time points first (this handles empty matrices too)
  if (nrow(fmri_data) < 10) {
    stop(sprintf("In %s: Insufficient time points (%d). Need at least 10 time points for meaningful HRF estimation.",
                 caller, nrow(fmri_data)), call. = FALSE)
  }
  
  # Check for no voxels
  if (ncol(fmri_data) == 0) {
    stop(sprintf("In %s: No voxels found in fMRI data", caller), call. = FALSE)
  }
  
  # Check if matrix is numeric (after we know it's not empty)
  if (!is.numeric(fmri_data)) {
    stop(sprintf("In %s: fmri_data must be a matrix or fmri_dataset object. Got: %s",
                 caller, class(fmri_data)[1]), call. = FALSE)
  }
  
  # Validate matrix properties
  spec <- list(
    name = "fmri_data",
    type = "matrix",
    finite = TRUE
  )
  
  .validate_data(fmri_data, spec, safety_mode = "balanced", caller = caller)
  
  # Check for constant voxels
  voxel_sds <- apply(fmri_data, 2, sd, na.rm = TRUE)
  n_constant <- sum(voxel_sds < .Machine$double.eps, na.rm = TRUE)
  if (n_constant > 0) {
    warning(sprintf("In %s: %d voxels have constant values and will produce invalid results",
                    caller, n_constant), call. = FALSE)
  }
  
  list(
    data = fmri_data,
    is_fmri_dataset = FALSE,
    n_time = nrow(fmri_data),
    n_voxels = ncol(fmri_data)
  )
}

#' Validate event model
#'
#' @param event_model Event model (matrix or event_model object)
#' @param n_time Number of time points
#' @param caller Name of calling function
#' @keywords internal
.validate_event_model <- function(event_model, n_time, caller = "function") {
  if (is.null(event_model)) {
    stop(sprintf("In %s: event_model cannot be NULL", caller), call. = FALSE)
  }
  
  # Handle event_model objects
  if (inherits(event_model, "event_model")) {
    # Extract design matrix from event_model
    if (requireNamespace("fmrireg", quietly = TRUE)) {
      design <- tryCatch({
        as.matrix(fmrireg::design_matrix(event_model))
      }, error = function(e) {
        stop(sprintf("In %s: Cannot extract design matrix from event_model: %s",
                     caller, e$message), call. = FALSE)
      })
      
      return(list(
        design = design,
        type = "event_model",
        n_events = ncol(design)
      ))
    } else {
      stop(sprintf("In %s: fmrireg package required to handle event_model objects",
                   caller), call. = FALSE)
    }
  }
  
  # Handle numeric vectors by converting to matrix
  if (is.numeric(event_model) && is.vector(event_model)) {
    event_model <- matrix(event_model, ncol = 1)
  }
  
  # Otherwise must be a matrix
  if (!is.matrix(event_model)) {
    stop(sprintf("In %s: event_model must be a matrix or event_model object. Got: %s",
                 caller, class(event_model)[1]), call. = FALSE)
  }
  
  # Check dimensions
  if (nrow(event_model) != n_time) {
    stop(sprintf("In %s: event_model rows (%d) must match fmri_data time points (%d)",
                 caller, nrow(event_model), n_time), call. = FALSE)
  }
  
  if (ncol(event_model) == 0) {
    stop(sprintf("In %s: event_model has no events", caller), call. = FALSE)
  }
  
  # Check for empty events
  event_density <- mean(event_model != 0)
  if (event_density < 0.01) {
    perc <- round(event_density * 100, 1)
    warning(sprintf("In %s: Very low event density (%g%% of time points). Consider if you have enough events for reliable estimation.",
                    caller, perc), call. = FALSE)
  }
  
  # Return structured result
  list(
    design = event_model,
    type = "matrix",
    n_events = ncol(event_model),
    event_density = event_density
  )
}

#' Validate theta bounds (delegates to unified validation)
#'
#' @param theta_bounds List with lower and upper bounds
#' @param n_params Number of parameters
#' @param param_names Optional parameter names
#' @param caller Calling function name
#' @return Validated bounds or NULL
#' @keywords internal
.validate_theta_bounds <- function(theta_bounds, n_params, param_names = NULL, 
                                   caller = "estimate_parametric_hrf") {
  .validate_bounds(theta_bounds, n_params, param_names, caller)
}

#' Validate numeric parameter (legacy wrapper)
#'
#' @param x Numeric value to validate
#' @param name Parameter name
#' @param min_val Minimum value
#' @param max_val Maximum value
#' @param allow_null Whether NULL is allowed
#' @param default Default value if NULL
#' @param caller Calling function name
#' @return Validated value
#' @keywords internal
.validate_numeric_param <- function(x, name, min_val = -Inf, max_val = Inf, 
                                    allow_null = TRUE, default = NULL,
                                    caller = "estimate_parametric_hrf") {
  # Handle NULL with default early
  if (is.null(x) && allow_null) {
    return(default)
  }
  
  # Special handling for NA
  if (length(x) == 1 && is.na(x)) {
    stop(sprintf("In %s: %s must be a single numeric value. Got: %s", 
                 caller, name, class(x)[1]), call. = FALSE)
  }
  
  # Special handling for infinite values
  if (length(x) == 1 && is.infinite(x)) {
    stop(sprintf("In %s: %s cannot be infinite", caller, name), call. = FALSE)
  }
  
  constraints <- list(
    null_ok = allow_null,
    finite = TRUE
  )
  
  if (is.finite(min_val) || is.finite(max_val)) {
    constraints$range <- c(min_val, max_val)
  }
  
  .validate_numeric(x, name, constraints, caller)
}

#' Master validation function for comprehensive mode
#'
#' @param fmri_data fMRI data
#' @param event_model Event model
#' @param parametric_model Model name
#' @param theta_seed Seed parameters
#' @param theta_bounds Parameter bounds
#' @param hrf_span HRF span
#' @param lambda_ridge Ridge parameter
#' @param recenter_global_passes Global passes
#' @param recenter_epsilon Convergence epsilon
#' @param r2_threshold R-squared threshold
#' @param mask Optional mask
#' @param verbose Verbose output
#' @param caller Calling function
#' @keywords internal
.rock_solid_validate_inputs <- function(
  fmri_data,
  event_model,
  parametric_model,
  theta_seed,
  theta_bounds,
  hrf_span,
  lambda_ridge,
  recenter_global_passes,
  recenter_epsilon,
  r2_threshold,
  mask = NULL,
  verbose = TRUE,
  caller = "estimate_parametric_hrf"
) {
  
  if (verbose) {
    cat("  [Validation] Starting comprehensive input validation...\n")
  }
  
  # fMRI data validation
  fmri_info <- .validate_fmri_data(fmri_data, caller)
  n_time <- fmri_info$n_time
  n_voxels <- fmri_info$n_voxels
  
  # Event model validation
  event_info <- .validate_event_model(event_model, n_time, caller)
  
  # Model validation
  if (!identical(tolower(parametric_model), "lwu")) {
    stop(sprintf("In %s: Only 'lwu' model currently supported. Got: '%s'", 
                 caller, parametric_model), call. = FALSE)
  }
  
  # Numeric parameters
  hrf_span <- .validate_numeric_param(hrf_span, "hrf_span", 
                                      min_val = 5, max_val = 60,
                                      allow_null = TRUE, default = 30, caller = caller)
  
  lambda_ridge <- .validate_numeric_param(lambda_ridge, "lambda_ridge",
                                          min_val = 0, max_val = 10,
                                          allow_null = TRUE, default = 0.01, caller = caller)
  
  recenter_global_passes <- .validate_numeric_param(recenter_global_passes, 
                                                    "recenter_global_passes",
                                                    min_val = 0, max_val = 20,
                                                    allow_null = TRUE, default = 3, caller = caller)
  
  recenter_epsilon <- .validate_numeric_param(recenter_epsilon, "recenter_epsilon",
                                              min_val = 1e-6, max_val = 1,
                                              allow_null = TRUE, default = 0.01, caller = caller)
  
  r2_threshold <- .validate_numeric_param(r2_threshold, "r2_threshold",
                                          min_val = 0, max_val = 1,
                                          allow_null = TRUE, default = 0.1, caller = caller)
  
  # Bounds validation
  if (!is.null(theta_bounds)) {
    theta_bounds <- .validate_theta_bounds(theta_bounds, n_params = 3, 
                                          param_names = c("tau", "sigma", "rho"),
                                          caller = caller)
  }
  
  # Seed validation
  if (!is.null(theta_seed) && !identical(theta_seed, "data_driven")) {
    if (!is.numeric(theta_seed) || length(theta_seed) != 3) {
      stop(caller, ": theta_seed must be numeric vector of length 3 or 'data_driven'. ",
           "Got: ", class(theta_seed)[1], " of length ", length(theta_seed), 
           call. = FALSE)
    }
    
    if (any(!is.finite(theta_seed))) {
      stop(caller, ": theta_seed must contain finite values.", call. = FALSE)
    }
  }
  
  # Mask validation
  if (!is.null(mask)) {
    if (!is.logical(mask) && !is.numeric(mask)) {
      stop(caller, ": mask must be logical or numeric. Got: ", 
           class(mask)[1], call. = FALSE)
    }
    
    if (length(mask) != n_voxels) {
      stop(caller, ": mask length (", length(mask), 
           ") must match number of voxels (", n_voxels, ")", 
           call. = FALSE)
    }
    
    # Convert to logical
    if (is.numeric(mask)) {
      mask <- mask != 0
    }
    
    if (sum(mask) == 0) {
      stop(caller, ": mask excludes all voxels", call. = FALSE)
    }
  }
  
  if (verbose) {
    cat("  [Validation] All inputs validated successfully\n")
    cat(sprintf("    - Data: %d timepoints x %d voxels\n", n_time, n_voxels))
    if (!is.null(mask)) {
      cat(sprintf("    - Mask: %d active voxels\n", sum(mask)))
    }
  }
  
  # Return validated data structure
  return(list(
    fmri_data = list(
      data = fmri_info$data,
      is_fmri_dataset = fmri_info$is_fmri_dataset,
      n_time = n_time,
      n_vox = n_voxels
    ),
    event_model = event_info$design,
    parametric_model = tolower(parametric_model),
    theta_seed = theta_seed,
    theta_bounds = theta_bounds,
    hrf_span = hrf_span,
    lambda_ridge = lambda_ridge,
    recenter_global_passes = recenter_global_passes,
    recenter_epsilon = recenter_epsilon,
    r2_threshold = r2_threshold,
    mask = mask,
    verbose = verbose
  ))
}