# Helper to check if baseline_model is "intercept" string
.is_intercept_baseline <- function(baseline_model) {
  !is.null(baseline_model) && 
  is.character(baseline_model) && 
  length(baseline_model) == 1 && 
  baseline_model == "intercept"
}

# Identify voxels with negligible temporal variance.
.detect_constant_voxels <- function(Y, tol = .Machine$double.eps^0.5) {
  if (is.null(dim(Y)) || ncol(Y) == 0) {
    return(integer(0))
  }
  voxel_sd <- apply(Y, 2, stats::sd, na.rm = TRUE)
  which(!is.finite(voxel_sd) | voxel_sd < tol)
}

# Force stable outputs for known degenerate voxels.
.apply_constant_voxel_guards <- function(r_squared, amplitudes, constant_idx) {
  if (length(constant_idx) == 0) {
    return(list(r_squared = r_squared, amplitudes = amplitudes))
  }
  r_squared[constant_idx] <- 0
  amplitudes[constant_idx] <- 0
  list(r_squared = r_squared, amplitudes = amplitudes)
}

# Build stage-0 output from precomputed projected inputs.
.stage0_use_prepared_inputs <- function(prepared_inputs, baseline_model = NULL, verbose = FALSE) {
  required_fields <- c("Y_proj", "S_target_proj", "hrf_eval_times")
  missing_fields <- setdiff(required_fields, names(prepared_inputs))
  if (length(missing_fields) > 0) {
    stop(
      "prepared_inputs is missing required fields: ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.matrix(prepared_inputs$Y_proj) || !is.numeric(prepared_inputs$Y_proj)) {
    stop("prepared_inputs$Y_proj must be a numeric matrix", call. = FALSE)
  }
  if (!is.matrix(prepared_inputs$S_target_proj) || !is.numeric(prepared_inputs$S_target_proj)) {
    stop("prepared_inputs$S_target_proj must be a numeric matrix", call. = FALSE)
  }
  if (nrow(prepared_inputs$Y_proj) != nrow(prepared_inputs$S_target_proj)) {
    stop("prepared_inputs$Y_proj and prepared_inputs$S_target_proj must have same number of rows", call. = FALSE)
  }
  if (!is.numeric(prepared_inputs$hrf_eval_times) || length(prepared_inputs$hrf_eval_times) < 2) {
    stop("prepared_inputs$hrf_eval_times must be numeric with length >= 2", call. = FALSE)
  }

  n_time <- nrow(prepared_inputs$Y_proj)
  n_vox <- ncol(prepared_inputs$Y_proj)

  if (is.null(prepared_inputs$Y_raw)) {
    prepared_inputs$Y_raw <- prepared_inputs$Y_proj
  }
  if (is.null(prepared_inputs$S_target)) {
    prepared_inputs$S_target <- prepared_inputs$S_target_proj
  }
  if (is.null(prepared_inputs$scan_times)) {
    prepared_inputs$scan_times <- seq_len(n_time)
  }
  if (is.null(prepared_inputs$Z)) {
    prepared_inputs$Z <- NULL
  }
  if (is.null(prepared_inputs$baseline_model)) {
    prepared_inputs$baseline_model <- baseline_model
  }

  storage.mode(prepared_inputs$Y_proj) <- "double"
  storage.mode(prepared_inputs$S_target_proj) <- "double"
  storage.mode(prepared_inputs$S_target) <- "double"

  constant_voxel_idx <- .detect_constant_voxels(prepared_inputs$Y_proj)
  if (verbose && length(constant_voxel_idx) > 0) {
    cat(sprintf("  Detected %d constant/near-constant voxels\n", length(constant_voxel_idx)))
  }

  list(
    inputs = prepared_inputs,
    n_vox = n_vox,
    n_time = n_time,
    constant_voxel_idx = constant_voxel_idx,
    validation_level = "prepared",
    profiler = NULL
  )
}

#' Stage 0: Validate inputs and prepare data
#'
#' @param fmri_data Raw fMRI data
#' @param event_model Event model
#' @param hrf_interface HRF interface object
#' @param config List containing all configuration parameters
#' @return List with validated and prepared data
#' @noRd
.stage0_validate_and_prepare <- function(fmri_data, event_model, hrf_interface, config) {
  # Extract config parameters
  confound_formula <- config$confound_formula
  baseline_model <- config$baseline_model
  hrf_eval_times <- config$hrf_eval_times
  hrf_span <- config$hrf_span
  mask <- config$mask
  safety_mode <- config$safety_mode
  verbose <- config$verbose
  
  if (verbose) cat("-> Stage 0: Input validation and safety checks...\n")
  
  # Create stage profiler if verbose
  profiler <- if (verbose) .diag_stage_profiler() else NULL
  if (!is.null(profiler)) profiler$start_stage("validation")
  
  # Validate inputs based on safety mode
  validation_level <- switch(safety_mode,
    maximum = "comprehensive",
    balanced = "standard",
    performance = "minimal"
  )
  
  # Run validation
  if (validation_level == "comprehensive") {
    .rock_solid_validate_inputs(
      fmri_data = fmri_data,
      event_model = event_model,
      parametric_model = config$parametric_model,
      theta_seed = config$theta_seed,
      theta_bounds = config$theta_bounds,
      hrf_span = hrf_span,
      lambda_ridge = config$lambda_ridge,
      recenter_global_passes = config$global_passes,
      recenter_epsilon = config$convergence_epsilon,
      r2_threshold = 0.1,
      mask = mask,
      verbose = verbose,
      caller = "estimate_parametric_hrf"
    )
  } else if (validation_level == "standard") {
    fmri_ok <- .validate_fmri_data(fmri_data, "estimate_parametric_hrf")
    .validate_event_model(event_model, fmri_ok$n_time, "estimate_parametric_hrf")
  }
  
  # Prepare data
  if (verbose) cat("-> Preparing data matrices...\n")
  inputs <- .prepare_parametric_inputs(
    fmri_data = fmri_data,
    event_model = event_model,
    confound_formula = confound_formula,
    baseline_model = baseline_model,
    hrf_eval_times = hrf_eval_times,
    hrf_span = hrf_span,
    mask = mask
  )
  
  n_vox <- ncol(inputs$Y_proj)
  n_time <- nrow(inputs$Y_proj)
  constant_voxel_idx <- .detect_constant_voxels(inputs$Y_proj)
  if (verbose && length(constant_voxel_idx) > 0) {
    cat(sprintf("  Detected %d constant/near-constant voxels\n", length(constant_voxel_idx)))
  }
  
  # End profiling
  if (!is.null(profiler)) {
    profiler$end_stage("validation", list(n_vox = n_vox, n_time = n_time))
  }
  
  list(
    inputs = inputs,
    n_vox = n_vox,
    n_time = n_time,
    constant_voxel_idx = constant_voxel_idx,
    validation_level = validation_level,
    profiler = profiler
  )
}

#' Stage 1: Initialize parameters
#'
#' @param prepared_data Output from stage 0
#' @param hrf_interface HRF interface object
#' @param config Configuration list
#' @return List with initialized parameters
#' @noRd
.stage1_initialize_parameters <- function(prepared_data, hrf_interface, config) {
  verbose <- config$verbose
  if (verbose) cat("\n-> Stage 1: Parameter initialization...\n")
  
  n_vox <- prepared_data$n_vox
  n_params <- length(hrf_interface$parameter_names)
  
  # Initialize theta matrix
  theta_current <- matrix(NA_real_, n_vox, n_params)
  colnames(theta_current) <- hrf_interface$parameter_names
  
  # Handle theta_seed
  theta_seed <- config$theta_seed
  if (is.null(theta_seed)) {
    theta_seed <- hrf_interface$default_seed()
  } else if (identical(theta_seed, "data_driven")) {
    if (verbose) cat("  Computing data-driven initialization...\n")
    theta_seed <- .compute_data_driven_seed(
      Y = prepared_data$inputs$Y_proj,
      S = prepared_data$inputs$S_target_proj,
      hrf_interface = hrf_interface,
      theta_bounds = config$theta_bounds
    )
  } else {
    # Validate theta_seed
    if (!is.numeric(theta_seed)) {
      stop("theta_seed must be numeric", call. = FALSE)
    }
    if (length(theta_seed) != n_params) {
      stop(sprintf("theta_seed must have length %d", n_params), call. = FALSE)
    }
  }
  
  # Handle theta_bounds - use active_bounds from interface if available
  if (!is.null(hrf_interface$active_bounds)) {
    theta_bounds <- hrf_interface$active_bounds
  } else if (!is.null(config$theta_bounds)) {
    theta_bounds <- config$theta_bounds
  } else {
    theta_bounds <- hrf_interface$default_bounds()
  }
  
  # Validate bounds if they came from config
  if (!is.null(config$theta_bounds) && is.null(hrf_interface$active_bounds)) {
    # Validate theta_bounds
    if (!is.list(theta_bounds) || 
        !all(c("lower", "upper") %in% names(theta_bounds))) {
      stop("theta_bounds missing required elements 'lower' and 'upper'", call. = FALSE)
    }
    if (length(theta_bounds$lower) != n_params || 
        length(theta_bounds$upper) != n_params) {
      stop(sprintf("theta_bounds components must have length %d", n_params), call. = FALSE)
    }
    if (!is.numeric(theta_bounds$lower) || !is.numeric(theta_bounds$upper)) {
      stop("theta_bounds values must be numeric", call. = FALSE)
    }
  }
  
  # Ensure theta_seed is within safe bounds
  eps <- c(0.01, 0.01, 0.01)
  theta_seed <- as.numeric(.clamp_theta_to_bounds(
    theta = theta_seed,
    theta_bounds = theta_bounds,
    epsilon = eps
  ))
  
  # Set initial values
  for (j in seq_len(ncol(theta_current))) {
    theta_current[, j] <- theta_seed[j]
  }
  
  # K-means initialization if requested
  kmeans_iterations <- 0
  if (config$kmeans_refinement && config$kmeans_k > 1) {
    if (verbose) cat("  Performing K-means clustering for initialization...\n")
    kmeans_result <- .perform_kmeans_initialization(
      Y = prepared_data$inputs$Y_proj,
      S = prepared_data$inputs$S_target_proj,
      k = config$kmeans_k,
      hrf_interface = hrf_interface,
      theta_bounds = theta_bounds
    )
    kmeans_iterations <- kmeans_result$iterations
    
    # Update seeds based on clusters
    for (k in seq_len(config$kmeans_k)) {
      cluster_voxels <- which(kmeans_result$cluster == k)
      if (length(cluster_voxels) > 0) {
        theta_current[cluster_voxels, ] <- matrix(
          kmeans_result$centers[k, ],
          nrow = length(cluster_voxels),
          ncol = ncol(theta_current),
          byrow = TRUE
        )
      }
    }
  }
  
  list(
    theta_current = theta_current,
    theta_seed = theta_seed,
    theta_bounds = theta_bounds,
    kmeans_iterations = kmeans_iterations
  )
}

#' Stage 2: Run core estimation
#'
#' @param prepared_data Output from stage 0
#' @param init_params Output from stage 1
#' @param hrf_interface HRF interface object
#' @param config Configuration list
#' @return List with core estimation results
#' @noRd
.stage2_core_estimation <- function(prepared_data, init_params, hrf_interface, config) {
  verbose <- config$verbose
  if (verbose) cat("\n-> Stage 2: Core parametric estimation...\n")
  
  # Parallel processing not yet implemented
  parallel_config <- NULL
  process_function <- .parametric_engine
  
  # Run core engine
  args <- list(
    Y_proj = prepared_data$inputs$Y_proj,
    S_target_proj = prepared_data$inputs$S_target_proj,
    hrf_eval_times = prepared_data$inputs$hrf_eval_times,
    hrf_interface = hrf_interface,
    theta_seed = init_params$theta_seed,
    theta_bounds = init_params$theta_bounds,
    lambda_ridge = config$lambda_ridge,
    baseline_model = config$baseline_model
  )
  
  # Parallel config would be added here when implemented
  
  core_result <- do.call(process_function, args)
  
  # Update current estimates
  theta_current <- core_result$theta_hat
  amplitudes <- core_result$beta0
  
  # Ensure theta_current is always a matrix
  if (!is.matrix(theta_current) || length(dim(theta_current)) != 2) {
    theta_current <- matrix(theta_current, 
                            nrow = prepared_data$n_vox, 
                            ncol = length(hrf_interface$parameter_names))
    colnames(theta_current) <- hrf_interface$parameter_names
  }
  
  r_squared <- core_result$r_squared
  guarded <- .apply_constant_voxel_guards(
    r_squared = r_squared,
    amplitudes = amplitudes,
    constant_idx = prepared_data$constant_voxel_idx
  )
  r_squared <- guarded$r_squared
  amplitudes <- guarded$amplitudes
  
  if (verbose) {
    cat(sprintf("  Initial fit: Mean R^2 = %.3f (range: %.3f - %.3f)\n",
                mean(r_squared), min(r_squared), max(r_squared)))
  }
  
  list(
    theta_current = theta_current,
    amplitudes = amplitudes,
    r_squared = r_squared,
    residuals = core_result$residuals,  # Pass through residuals
    intercepts = if (!is.null(core_result$coeffs) && 
                     .is_intercept_baseline(config$baseline_model)) {
                   core_result$coeffs[1, ]  # First row contains intercepts
                 } else NULL,
    core_result = core_result,
    parallel_config = parallel_config
  )
}
