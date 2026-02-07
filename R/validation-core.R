#' Core validation functions for fmriparametric
#'
#' This module provides unified validation functions that consolidate
#' redundant validation code across the package.

#' Unified data validation with performance levels
#'
#' @param x Data to validate
#' @param spec Validation specification list
#' @param safety_mode Safety mode: "maximum", "balanced", or "performance"
#' @param caller Name of calling function for error messages
#' @return Validated data (may be modified)
#' @keywords internal
.validate_data <- function(x, spec, safety_mode = "balanced", caller = NULL) {
  # Quick return for performance mode
  if (safety_mode == "performance" && !is.null(x)) {
    return(x)
  }
  
  # Check NULL
  if (is.null(x) && !isTRUE(spec$null_ok)) {
    .diag_abort(
      sprintf("%s cannot be NULL", spec$name %||% "Input"),
      caller = caller
    )
  }
  
  if (is.null(x)) return(x)
  
  # Type checking
  if (!is.null(spec$type)) {
    type_check <- switch(spec$type,
      numeric = is.numeric(x),
      matrix = is.matrix(x),
      list = is.list(x),
      character = is.character(x),
      logical = is.logical(x),
      TRUE
    )
    
    if (!type_check) {
      .diag_abort(
        sprintf("%s must be of type %s", spec$name %||% "Input", spec$type),
        caller = caller
      )
    }
  }
  
  # Dimension checking for matrices
  if (is.matrix(x) && !is.null(spec$dims)) {
    if (!is.null(spec$dims$nrow) && nrow(x) != spec$dims$nrow) {
      .diag_abort(
        sprintf("%s must have %d rows, got %d", 
                spec$name %||% "Matrix", spec$dims$nrow, nrow(x)),
        caller = caller
      )
    }
    if (!is.null(spec$dims$ncol) && ncol(x) != spec$dims$ncol) {
      .diag_abort(
        sprintf("%s must have %d columns, got %d", 
                spec$name %||% "Matrix", spec$dims$ncol, ncol(x)),
        caller = caller
      )
    }
  }
  
  # Numeric constraints (only in balanced/maximum mode)
  if (safety_mode != "performance" && is.numeric(x)) {
    # Single-pass validation for efficiency - vectorized operations
    if (is.matrix(x)) {
      # For matrices, use vectorized column operations
      if (!is.null(spec$finite) && spec$finite) {
        finite_check <- is.finite(x)
        if (!all(finite_check)) {
          bad_cols <- which(colSums(!finite_check) > 0)
          .diag_abort(
            sprintf("%s contains non-finite values in columns: %s", 
                    spec$name %||% "Input", 
                    paste(head(bad_cols, 5), collapse = ", ")),
            caller = caller
          )
        }
      }
      
      if (!is.null(spec$positive) && spec$positive) {
        pos_check <- x > 0 | is.na(x)
        if (!all(pos_check)) {
          .diag_abort(
            sprintf("%s must contain only positive values", spec$name %||% "Input"),
            caller = caller
          )
        }
      }
      
      if (!is.null(spec$range)) {
        # Vectorized range check
        range_check <- x >= spec$range[1] & x <= spec$range[2]
        if (!all(range_check, na.rm = TRUE)) {
          .diag_abort(
            sprintf("%s values must be in range [%g, %g]", 
                    spec$name %||% "Input", spec$range[1], spec$range[2]),
            caller = caller
          )
        }
      }
    } else {
      # Vector case - original logic but with early exit
      if (!is.null(spec$finite) && spec$finite && any(!is.finite(x))) {
        .diag_abort(
          sprintf("%s contains non-finite values", spec$name %||% "Input"),
          caller = caller
        )
      }
      
      if (!is.null(spec$positive) && spec$positive && any(x <= 0, na.rm = TRUE)) {
        .diag_abort(
          sprintf("%s must contain only positive values", spec$name %||% "Input"),
          caller = caller
        )
      }
      
      if (!is.null(spec$range)) {
        if (any(x < spec$range[1] | x > spec$range[2], na.rm = TRUE)) {
          .diag_abort(
            sprintf("%s values must be in range [%g, %g]", 
                    spec$name %||% "Input", spec$range[1], spec$range[2]),
            caller = caller
          )
        }
      }
    }
  }
  
  # Length constraints
  if (!is.null(spec$length)) {
    actual_length <- if (is.matrix(x)) nrow(x) * ncol(x) else length(x)
    if (actual_length != spec$length) {
      .diag_abort(
        sprintf("%s must have length %d, got %d", 
                spec$name %||% "Input", spec$length, actual_length),
        caller = caller
      )
    }
  }
  
  x
}

#' Validate fMRI data with tiered validation
#'
#' @param fmri_data fMRI data (matrix or fmri_dataset)
#' @param caller Calling function name
#' @param safety_mode Validation level
#' @return List with validated data and metadata
#' @keywords internal
.validate_fmri_data_tiered <- function(fmri_data, caller = NULL, 
                                      safety_mode = "balanced") {
  
  # Minimal validation - just check it exists and has right structure
  if (safety_mode == "performance") {
    if (is.null(fmri_data)) {
      stop("fmri_data cannot be NULL")
    }
    
    if (inherits(fmri_data, "fmri_dataset")) {
      return(list(
        data = fmri_data,
        is_fmri_dataset = TRUE,
        n_time = dim(fmri_data)[length(dim(fmri_data))],
        n_voxels = prod(dim(fmri_data)[-length(dim(fmri_data))])
      ))
    } else if (is.matrix(fmri_data)) {
      return(list(
        data = fmri_data,
        is_fmri_dataset = FALSE,
        n_time = nrow(fmri_data),
        n_voxels = ncol(fmri_data)
      ))
    } else {
      stop("fmri_data must be a matrix or fmri_dataset object")
    }
  }
  
  # Use existing validation for balanced/maximum modes
  .validate_fmri_data(fmri_data, caller)
}

#' Unified numeric parameter validation
#'
#' Replaces both .validate_numeric_param and .validate_constraints
#'
#' @param x Numeric value to validate
#' @param name Parameter name for error messages
#' @param constraints List of constraints (range, finite, positive, integer)
#' @param caller Calling function name
#' @return Validated numeric value
#' @keywords internal
.validate_numeric <- function(x, name, constraints = list(), caller = NULL) {
  # NULL handling
  if (is.null(x) && !isTRUE(constraints$null_ok)) {
    .diag_abort(sprintf("%s cannot be NULL", name), caller = caller)
  }
  
  if (is.null(x)) return(x)
  
  # Type and length check
  if (!is.numeric(x)) {
    .diag_abort(sprintf("%s must be a single numeric value. Got: %s", name, class(x)[1]), caller = caller)
  }
  
  if (length(x) != 1) {
    .diag_abort(sprintf("%s must be a single numeric value (length 1). Got length %d", name, length(x)), caller = caller)
  }
  
  # Apply constraints efficiently
  if (isTRUE(constraints$integer) && x != round(x)) {
    .diag_abort(sprintf("%s must be an integer", name), caller = caller)
  }
  
  if (isTRUE(constraints$finite) && !is.finite(x)) {
    .diag_abort(sprintf("%s must be finite", name), caller = caller)
  }
  
  if (isTRUE(constraints$positive) && x <= 0) {
    .diag_abort(sprintf("%s must be positive", name), caller = caller)
  }
  
  if (!is.null(constraints$range)) {
    if (x < constraints$range[1] || x > constraints$range[2]) {
      .diag_abort(
        sprintf("%s: %s is outside valid range [%g, %g]. Got: %g", 
                caller %||% "validate_numeric", name, 
                constraints$range[1], constraints$range[2], x),
        caller = NULL
      )
    }
  }
  
  x
}

#' Validate parameter bounds with enhanced checking
#'
#' @param bounds List with lower and upper bounds
#' @param n_params Expected number of parameters
#' @param param_names Optional parameter names for better errors
#' @param caller Calling function name
#' @return Validated bounds
#' @keywords internal
.validate_bounds <- function(bounds, n_params, param_names = NULL, caller = NULL) {
  if (is.null(bounds)) return(NULL)
  
  # Check if it's a list first
  if (!is.list(bounds)) {
    .diag_abort(
      "theta_bounds must be a list with 'lower' and 'upper' elements",
      caller = caller
    )
  }
  
  # Structure validation
  if (!all(c("lower", "upper") %in% names(bounds))) {
    missing <- setdiff(c("lower", "upper"), names(bounds))
    .diag_abort(
      sprintf("theta_bounds missing required elements: %s", paste(missing, collapse = ", ")),
      caller = caller
    )
  }
  
  # Length validation
  if (length(bounds$lower) != n_params || length(bounds$upper) != n_params) {
    .diag_abort(
      sprintf("theta_bounds dimensions incorrect: lower (%d), upper (%d). Expected %d", 
              length(bounds$lower), length(bounds$upper), n_params),
      caller = caller
    )
  }
  
  # Numeric validation
  if (!is.numeric(bounds$lower) || !is.numeric(bounds$upper)) {
    .diag_abort("theta_bounds must contain numeric values", caller = caller)
  }
  
  # Order validation
  if (any(bounds$lower >= bounds$upper)) {
    bad_idx <- which(bounds$lower >= bounds$upper)
    param_info <- if (!is.null(param_names)) {
      paste(param_names[bad_idx], collapse = ", ")
    } else {
      paste(bad_idx, collapse = ", ")
    }
    .diag_abort(
      sprintf("theta_bounds: lower must be less than upper for %s", 
              param_info),
      caller = caller
    )
  }
  
  # Physiological plausibility checks (if param names available)
  if (!is.null(param_names)) {
    # Check tau bounds
    if ("tau" %in% param_names) {
      tau_idx <- which(param_names == "tau")
      if (bounds$lower[tau_idx] < 0) {
        warning(sprintf("%s: tau lower bound < 0 is non-physiological. Consider using >= 0.", 
                        caller %||% "validate_bounds"), call. = FALSE)
      }
      if (bounds$upper[tau_idx] > 30) {
        warning(sprintf("%s: tau upper bound > 30s is unusually high for HRF peak time.", 
                        caller %||% "validate_bounds"), call. = FALSE)
      }
    }
    
    # Check sigma bounds
    if ("sigma" %in% param_names) {
      sigma_idx <- which(param_names == "sigma")
      if (bounds$lower[sigma_idx] < 0.1) {
        warning(sprintf("%s: sigma lower bound < 0.1 may cause numerical instability.", 
                        caller %||% "validate_bounds"), call. = FALSE)
      }
      if (bounds$upper[sigma_idx] > 20) {
        warning(sprintf("%s: sigma upper bound > 20s is unusually wide for HRF.", 
                        caller %||% "validate_bounds"), call. = FALSE)
      }
    }
    
    # Check rho bounds
    if ("rho" %in% param_names) {
      rho_idx <- which(param_names == "rho")
      if (bounds$lower[rho_idx] < 0) {
        warning(sprintf("%s: rho lower bound < 0 prevents undershoot modeling.", 
                        caller %||% "validate_bounds"), call. = FALSE)
      }
      if (bounds$upper[rho_idx] > 2) {
        warning(sprintf("%s: rho upper bound > 2 is non-physiological for undershoot ratio.", 
                        caller %||% "validate_bounds"), call. = FALSE)
      }
    }
  }
  
  bounds
}

#' Create validation cache for repeated validations
#'
#' @return Environment for caching validation results
#' @keywords internal
.create_validation_cache <- function() {
  cache <- new.env(parent = emptyenv())
  cache$results <- list()
  cache$hits <- 0
  cache$misses <- 0
  
  cache
}

#' Get or compute validation result with caching
#'
#' @param cache Validation cache environment
#' @param key Cache key (usually data digest)
#' @param validate_fn Function to call if not cached
#' @param ... Arguments to pass to validate_fn
#' @return Validation result
#' @keywords internal
.cached_validate <- function(cache, key, validate_fn, ...) {
  if (key %in% names(cache$results)) {
    cache$hits <- cache$hits + 1
    return(cache$results[[key]])
  }
  
  cache$misses <- cache$misses + 1
  result <- validate_fn(...)
  cache$results[[key]] <- result
  
  # Limit cache size to prevent memory bloat
  if (length(cache$results) > 100) {
    # Remove oldest entries
    cache$results <- tail(cache$results, 50)
  }
  
  result
}

#' Create a digest key for validation caching
#'
#' @param x Data to create key for
#' @param spec Validation spec
#' @return Character digest key
#' @keywords internal
.validation_digest <- function(x, spec) {
  # For large matrices, sample to create digest
  if (is.matrix(x) && length(x) > 10000) {
    # Sample deterministically
    idx <- seq(1, length(x), length.out = 1000)
    sample_data <- x[idx]
  } else {
    sample_data <- x
  }
  
  digest::digest(list(
    data = sample_data,
    spec = spec,
    dims = dim(x),
    class = class(x)
  ), algo = "xxhash32")
}