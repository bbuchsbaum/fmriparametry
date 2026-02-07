#' Bounds utilities for fmriparametric
#'
#' This module provides centralized functions for managing parameter bounds
#' throughout the estimation pipeline, ensuring consistency and proper
#' enforcement at all stages.

#' Apply parameter bounds
#'
#' Enforces parameter bounds using vectorized pmin/pmax operations.
#' Works with both vectors and matrices of parameters.
#'
#' @param params Numeric vector or matrix of parameters
#' @param bounds List with 'lower' and 'upper' numeric vectors
#' @param epsilon Optional small value to avoid exact boundary values
#' @return Bounded parameters of same shape as input
#' @keywords internal
.apply_bounds <- function(params, bounds, epsilon = NULL) {
  if (is.null(bounds)) return(params)
  
  # Handle both vector and matrix inputs
  if (is.matrix(params)) {
    n_vox <- nrow(params)
    n_params <- ncol(params)
    
    # Expand bounds to match matrix dimensions
    lower <- matrix(bounds$lower, nrow = n_vox, ncol = n_params, byrow = TRUE)
    upper <- matrix(bounds$upper, nrow = n_vox, ncol = n_params, byrow = TRUE)
    
    # Apply epsilon if requested (for numerical derivatives)
    if (!is.null(epsilon)) {
      if (length(epsilon) == 1) epsilon <- rep(epsilon, n_params)
      eps_mat <- matrix(epsilon, nrow = n_vox, ncol = n_params, byrow = TRUE)
      lower <- lower + eps_mat
      upper <- upper - eps_mat
    }
    
    # Apply bounds
    pmin(upper, pmax(lower, params))
  } else {
    # Vector case
    lower <- bounds$lower
    upper <- bounds$upper
    
    if (!is.null(epsilon)) {
      if (length(epsilon) == 1) epsilon <- rep(epsilon, length(params))
      lower <- lower + epsilon
      upper <- upper - epsilon
    }
    
    pmax(lower, pmin(params, upper))
  }
}

#' Check if parameters are within bounds
#'
#' @param params Numeric vector or matrix of parameters
#' @param bounds List with 'lower' and 'upper' numeric vectors
#' @param tolerance Small value for numerical tolerance
#' @return Logical vector/matrix indicating which parameters are in bounds
#' @keywords internal
.check_bounds <- function(params, bounds, tolerance = .Machine$double.eps^0.5) {
  if (is.null(bounds)) return(rep(TRUE, length(params)))
  
  if (is.matrix(params)) {
    n_vox <- nrow(params)
    n_params <- ncol(params)
    
    # Expand bounds to match matrix dimensions
    lower <- matrix(bounds$lower, nrow = n_vox, ncol = n_params, byrow = TRUE)
    upper <- matrix(bounds$upper, nrow = n_vox, ncol = n_params, byrow = TRUE)
    
    # Check bounds with tolerance
    (params >= lower - tolerance) & (params <= upper + tolerance)
  } else {
    (params >= bounds$lower - tolerance) & (params <= bounds$upper + tolerance)
  }
}

#' Get active bounds from HRF interface or config
#'
#' Retrieves the authoritative bounds to use, preferring active_bounds
#' from the HRF interface if available.
#'
#' @param hrf_interface HRF interface list
#' @param config Configuration list
#' @return List with 'lower' and 'upper' bounds
#' @keywords internal
.get_active_bounds <- function(hrf_interface, config = NULL) {
  # Prefer active_bounds from interface
  if (!is.null(hrf_interface$active_bounds)) {
    return(hrf_interface$active_bounds)
  }
  
  # Fall back to config bounds
  if (!is.null(config) && !is.null(config$theta_bounds)) {
    return(config$theta_bounds)
  }
  
  # Last resort: get default bounds from interface
  if (!is.null(hrf_interface$default_bounds)) {
    return(hrf_interface$default_bounds())
  }
  
  stop("No parameter bounds available", call. = FALSE)
}

#' Create bounds-aware parameter update function
#'
#' Returns a function that updates parameters while respecting bounds.
#' Useful for optimization routines.
#'
#' @param bounds List with 'lower' and 'upper' bounds
#' @return Function that takes params and delta, returns bounded update
#' @keywords internal
.create_bounded_update <- function(bounds) {
  force(bounds)  # Capture bounds in closure
  
  function(params, delta, step_size = 1.0) {
    updated <- params + step_size * delta
    .apply_bounds(updated, bounds)
  }
}

#' Validate and adjust bounds for numerical stability
#'
#' Ensures bounds are suitable for the specific HRF model and
#' numerical methods being used.
#'
#' @param bounds List with 'lower' and 'upper' bounds
#' @param model Character string specifying the HRF model
#' @param warn Logical whether to issue warnings
#' @return Adjusted bounds list
#' @keywords internal
.adjust_bounds_for_stability <- function(bounds, model = "lwu", warn = TRUE) {
  adjusted_bounds <- bounds
  
  if (model == "lwu") {
    # Ensure sigma > 0.05 for numerical stability in fmrihrf
    if (bounds$lower[2] < 0.051) {
      if (warn) {
        warning("Adjusting sigma lower bound from ", bounds$lower[2], 
                " to 0.051 for numerical stability", call. = FALSE)
      }
      adjusted_bounds$lower[2] <- 0.051
    }
    
    # Ensure reasonable upper bounds
    if (bounds$upper[1] > 50) {
      if (warn) {
        warning("tau upper bound > 50 is extreme, consider reducing", 
                call. = FALSE)
      }
    }
    
    if (bounds$upper[2] > 30) {
      if (warn) {
        warning("sigma upper bound > 30 is extreme, consider reducing", 
                call. = FALSE)
      }
    }
  }
  
  adjusted_bounds
}

#' Report bounds violations
#'
#' Diagnostic function to report which parameters are outside bounds
#' and by how much.
#'
#' @param params Parameters to check
#' @param bounds Bounds list
#' @param param_names Optional parameter names
#' @param voxel_idx Optional voxel indices for reporting
#' @return Invisible NULL, prints diagnostic information
#' @keywords internal
.report_bounds_violations <- function(params, bounds, param_names = NULL, 
                                     voxel_idx = NULL) {
  if (is.null(bounds)) return(invisible(NULL))
  
  violations <- !.check_bounds(params, bounds)
  
  if (any(violations)) {
    if (is.matrix(params)) {
      # Matrix case - report by voxel
      vox_with_violations <- which(rowSums(violations) > 0)
      
      cat("Bounds violations detected in", length(vox_with_violations), "voxels:\n")
      
      for (v in head(vox_with_violations, 5)) {  # Show first 5
        v_idx <- if (!is.null(voxel_idx)) voxel_idx[v] else v
        param_violations <- which(violations[v, ])
        
        for (p in param_violations) {
          p_name <- if (!is.null(param_names)) param_names[p] else paste0("param", p)
          cat(sprintf("  Voxel %d: %s = %.3f (bounds: [%.3f, %.3f])\n",
                      v_idx, p_name, params[v, p], 
                      bounds$lower[p], bounds$upper[p]))
        }
      }
      
      if (length(vox_with_violations) > 5) {
        cat("  ... and", length(vox_with_violations) - 5, "more voxels\n")
      }
    } else {
      # Vector case
      param_violations <- which(violations)
      for (p in param_violations) {
        p_name <- if (!is.null(param_names)) param_names[p] else paste0("param", p)
        cat(sprintf("  %s = %.3f (bounds: [%.3f, %.3f])\n",
                    p_name, params[p], bounds$lower[p], bounds$upper[p]))
      }
    }
  }
  
  invisible(NULL)
}