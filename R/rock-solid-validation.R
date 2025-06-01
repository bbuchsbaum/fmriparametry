#' Rock Solid Input Validation Functions
#'
#' Comprehensive input validation that NEVER lets bad data through.
#' Every validation function provides actionable error messages.

#' Validate fMRI data input with extreme prejudice
#' @keywords internal
.validate_fmri_data <- function(fmri_data, caller = "estimate_parametric_hrf") {
  # Type validation
  if (is.null(fmri_data)) {
    stop(caller, ": fmri_data cannot be NULL. ",
         "Provide a matrix, fmri_dataset, or matrix_dataset object.", 
         call. = FALSE)
  }
  
  # Handle different input types
  if (inherits(fmri_data, c("fmri_dataset", "matrix_dataset"))) {
    # Extract data matrix
    if ("data" %in% names(fmri_data)) {
      data_matrix <- fmri_data$data
    } else {
      stop(caller, ": fmri_data object missing 'data' field. ",
           "Ensure your dataset object is properly formatted.", 
           call. = FALSE)
    }
  } else if (is.matrix(fmri_data)) {
    data_matrix <- fmri_data
  } else if (is.data.frame(fmri_data)) {
    warning(caller, ": Converting data.frame to matrix. ",
            "Consider providing a matrix directly for better performance.")
    data_matrix <- as.matrix(fmri_data)
  } else {
    stop(caller, ": fmri_data must be a matrix, fmri_dataset, or matrix_dataset. ",
         "Got: ", class(fmri_data)[1], 
         call. = FALSE)
  }
  
  # Dimension checks
  if (!is.matrix(data_matrix)) {
    stop(caller, ": Failed to extract data matrix from fmri_data.", 
         call. = FALSE)
  }
  
  dims <- dim(data_matrix)
  if (length(dims) != 2) {
    stop(caller, ": fmri_data must be 2D (time x voxels). ",
         "Got ", length(dims), "D data.", 
         call. = FALSE)
  }
  
  n_time <- dims[1]
  n_vox <- dims[2]
  
  if (n_time < 10) {
    stop(caller, ": Insufficient time points (", n_time, "). ",
         "Need at least 10 time points for meaningful HRF estimation.", 
         call. = FALSE)
  }
  
  if (n_vox < 1) {
    stop(caller, ": No voxels found in data (columns = ", n_vox, ").", 
         call. = FALSE)
  }
  
  # Data quality checks
  if (!is.numeric(data_matrix)) {
    stop(caller, ": fmri_data must contain numeric values. ",
         "Found: ", typeof(data_matrix), 
         call. = FALSE)
  }
  
  # Check for NA/NaN/Inf
  n_na <- sum(is.na(data_matrix))
  n_inf <- sum(is.infinite(data_matrix))
  
  if (n_na > 0) {
    prop_na <- n_na / length(data_matrix)
    if (prop_na > 0.5) {
      stop(caller, ": More than 50% of data is NA (", 
           round(prop_na * 100, 1), "%). ",
           "Check your data preprocessing.", 
           call. = FALSE)
    } else if (prop_na > 0.1) {
      warning(caller, ": ", round(prop_na * 100, 1), 
              "% of data contains NA values. ",
              "These will be handled but may affect results.")
    }
  }
  
  if (n_inf > 0) {
    stop(caller, ": Data contains ", n_inf, " infinite values. ",
         "Check for numerical overflow in preprocessing.", 
         call. = FALSE)
  }
  
  # Check for constant/zero voxels
  voxel_sds <- apply(data_matrix, 2, sd, na.rm = TRUE)
  n_constant <- sum(voxel_sds < .Machine$double.eps)
  n_zero <- sum(colSums(abs(data_matrix), na.rm = TRUE) < .Machine$double.eps)
  
  if (n_constant > 0) {
    warning(caller, ": Found ", n_constant, " constant voxels (no variation). ",
            "These will produce undefined results.")
  }
  
  if (n_zero > 0) {
    warning(caller, ": Found ", n_zero, " all-zero voxels. ",
            "These will produce zero parameter estimates.")
  }
  
  # Return validated data
  list(
    data = data_matrix,
    n_time = n_time,
    n_vox = n_vox,
    n_na = n_na,
    n_constant = n_constant,
    n_zero = n_zero,
    type = class(fmri_data)[1]
  )
}

#' Validate event model with comprehensive checks
#' @keywords internal
.validate_event_model <- function(event_model, n_time, caller = "estimate_parametric_hrf") {
  if (is.null(event_model)) {
    stop(caller, ": event_model cannot be NULL. ",
         "Provide an event_model object or stimulus matrix.", 
         call. = FALSE)
  }
  
  # Handle different input types
  if (inherits(event_model, "event_model")) {
    # Extract design matrix
    if ("terms" %in% names(event_model)) {
      if (length(event_model$terms) == 0) {
        stop(caller, ": event_model contains no terms. ",
             "Check your event model specification.", 
             call. = FALSE)
      }
      design_matrix <- event_model$terms[[1]]
    } else {
      stop(caller, ": event_model missing 'terms' field.", 
           call. = FALSE)
    }
  } else if (is.matrix(event_model) || is.numeric(event_model)) {
    design_matrix <- as.matrix(event_model)
  } else {
    stop(caller, ": event_model must be an event_model object or numeric matrix. ",
         "Got: ", class(event_model)[1], 
         call. = FALSE)
  }
  
  # Dimension checks
  if (nrow(design_matrix) != n_time) {
    stop(caller, ": event_model time points (", nrow(design_matrix), 
         ") don't match fmri_data time points (", n_time, ").", 
         call. = FALSE)
  }
  
  # Content checks
  if (sum(abs(design_matrix), na.rm = TRUE) < .Machine$double.eps) {
    stop(caller, ": event_model contains no events (all zeros). ",
         "Check your event timing specification.", 
         call. = FALSE)
  }
  
  event_density <- mean(design_matrix > 0, na.rm = TRUE)
  if (event_density > 0.9) {
    warning(caller, ": Very high event density (", 
            round(event_density * 100), "% of time points). ",
            "This may lead to poor HRF estimation.")
  } else if (event_density < 0.01) {
    warning(caller, ": Very low event density (", 
            round(event_density * 100, 2), "% of time points). ",
            "Consider if you have enough events for reliable estimation.")
  }
  
  list(
    design = design_matrix,
    n_events = sum(design_matrix > 0),
    event_density = event_density,
    type = class(event_model)[1]
  )
}

#' Validate parameter bounds with physiological constraints
#' @keywords internal
.validate_theta_bounds <- function(theta_bounds, n_params, param_names = NULL, 
                                   caller = "estimate_parametric_hrf") {
  if (is.null(theta_bounds)) {
    return(NULL)  # Use defaults
  }
  
  if (!is.list(theta_bounds)) {
    stop(caller, ": theta_bounds must be a list with 'lower' and 'upper' elements. ",
         "Got: ", class(theta_bounds)[1], 
         call. = FALSE)
  }
  
  if (!all(c("lower", "upper") %in% names(theta_bounds))) {
    missing <- setdiff(c("lower", "upper"), names(theta_bounds))
    stop(caller, ": theta_bounds missing required elements: ",
         paste(missing, collapse = ", "), 
         call. = FALSE)
  }
  
  lower <- theta_bounds$lower
  upper <- theta_bounds$upper
  
  # Length checks
  if (length(lower) != n_params || length(upper) != n_params) {
    stop(caller, ": theta_bounds dimensions incorrect. Expected ", n_params,
         " parameters, got lower = ", length(lower), ", upper = ", length(upper), 
         call. = FALSE)
  }
  
  # Numeric checks
  if (!is.numeric(lower) || !is.numeric(upper)) {
    stop(caller, ": theta_bounds must contain numeric values.", 
         call. = FALSE)
  }
  
  # Order checks
  violations <- which(lower >= upper)
  if (length(violations) > 0) {
    param_info <- if (!is.null(param_names)) {
      paste0(param_names[violations], " (", violations, ")")
    } else {
      as.character(violations)
    }
    stop(caller, ": theta_bounds: lower must be less than upper for all parameters. ",
         "Violations at: ", paste(param_info, collapse = ", "), 
         call. = FALSE)
  }
  
  # Physiological plausibility checks for LWU
  if (!is.null(param_names) && length(param_names) == 3 && 
      all(param_names == c("tau", "sigma", "rho"))) {
    
    # Tau (lag) checks
    if (lower[1] < 0) {
      warning(caller, ": tau lower bound < 0 is non-physiological. Consider using >= 0.")
    }
    if (upper[1] > 30) {
      warning(caller, ": tau upper bound > 30s is unusually high for HRF peak time.")
    }
    
    # Sigma (width) checks  
    if (lower[2] < 0.1) {
      warning(caller, ": sigma lower bound < 0.1 may cause numerical instability.")
    }
    if (upper[2] > 20) {
      warning(caller, ": sigma upper bound > 20s is unusually wide for HRF.")
    }
    
    # Rho (undershoot) checks
    if (lower[3] < 0) {
      warning(caller, ": rho lower bound < 0 prevents undershoot modeling.")
    }
    if (upper[3] > 2) {
      warning(caller, ": rho upper bound > 2 is non-physiological for undershoot ratio.")
    }
  }
  
  theta_bounds
}

#' Validate numeric parameters with range and sanity checks
#' @keywords internal
.validate_numeric_param <- function(x, name, min_val = -Inf, max_val = Inf, 
                                    allow_null = TRUE, default = NULL,
                                    caller = "estimate_parametric_hrf") {
  if (is.null(x)) {
    if (allow_null) {
      return(default)
    } else {
      stop(caller, ": ", name, " cannot be NULL.", call. = FALSE)
    }
  }
  
  if (!is.numeric(x) || length(x) != 1) {
    stop(caller, ": ", name, " must be a single numeric value. ",
         "Got: ", class(x)[1], " of length ", length(x), 
         call. = FALSE)
  }
  
  if (is.na(x)) {
    stop(caller, ": ", name, " cannot be NA.", call. = FALSE)
  }
  
  if (is.infinite(x)) {
    stop(caller, ": ", name, " cannot be infinite.", call. = FALSE)
  }
  
  if (x < min_val || x > max_val) {
    stop(caller, ": ", name, " = ", x, " is outside valid range [",
         min_val, ", ", max_val, "].", 
         call. = FALSE)
  }
  
  x
}

#' Master validation function that orchestrates all checks
#' @keywords internal
.rock_solid_validate_inputs <- function(
  fmri_data,
  event_model,
  parametric_hrf,
  theta_seed,
  theta_bounds,
  hrf_span,
  lambda_ridge,
  recenter_global_passes,
  recenter_epsilon,
  r2_threshold,
  mask,
  verbose,
  caller = "estimate_parametric_hrf"
) {
  
  if (verbose) cat("Performing rock-solid input validation...\n")
  
  # Core data validation
  fmri_valid <- .validate_fmri_data(fmri_data, caller)
  event_valid <- .validate_event_model(event_model, fmri_valid$n_time, caller)
  
  # Model validation
  if (!is.character(parametric_hrf) || length(parametric_hrf) != 1) {
    stop(caller, ": parametric_hrf must be a single character string.", 
         call. = FALSE)
  }
  
  if (tolower(parametric_hrf) != "lwu") {
    stop(caller, ": Only 'lwu' model currently supported. Got: '", parametric_hrf, "'", 
         call. = FALSE)
  }
  
  # Numeric parameter validation
  hrf_span <- .validate_numeric_param(hrf_span, "hrf_span", 
                                      min_val = 5, max_val = 60, 
                                      default = 30, caller = caller)
  
  lambda_ridge <- .validate_numeric_param(lambda_ridge, "lambda_ridge",
                                          min_val = 0, max_val = 10,
                                          default = 0.01, caller = caller)
  
  recenter_global_passes <- .validate_numeric_param(recenter_global_passes, 
                                                    "recenter_global_passes",
                                                    min_val = 0, max_val = 20,
                                                    default = 3, caller = caller)
  
  recenter_epsilon <- .validate_numeric_param(recenter_epsilon, "recenter_epsilon",
                                              min_val = 1e-10, max_val = 1,
                                              default = 0.01, caller = caller)
  
  r2_threshold <- .validate_numeric_param(r2_threshold, "r2_threshold",
                                          min_val = -1, max_val = 1,
                                          default = 0.1, caller = caller)
  
  # Return validated inputs
  list(
    fmri_data = fmri_valid,
    event_model = event_valid,
    parametric_hrf = tolower(parametric_hrf),
    hrf_span = hrf_span,
    lambda_ridge = lambda_ridge,
    recenter_global_passes = as.integer(recenter_global_passes),
    recenter_epsilon = recenter_epsilon,
    r2_threshold = r2_threshold,
    verbose = isTRUE(verbose)
  )
}