# Helper to check if baseline_model is "intercept" string
.is_intercept_baseline <- function(baseline_model) {
  !is.null(baseline_model) && 
  is.character(baseline_model) && 
  length(baseline_model) == 1 && 
  baseline_model == "intercept"
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
  
  # End profiling
  if (!is.null(profiler)) {
    profiler$end_stage("validation", list(n_vox = n_vox, n_time = n_time))
  }
  
  list(
    inputs = inputs,
    n_vox = n_vox,
    n_time = n_time,
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
  theta_seed <- pmax(theta_bounds$lower + eps,
                     pmin(theta_seed, theta_bounds$upper - eps))
  
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

#' Stage 3: Apply global refinement
#'
#' @param prepared_data Output from stage 0
#' @param core_results Output from stage 2
#' @param hrf_interface HRF interface object
#' @param config Configuration list
#' @return List with refined parameters
#' @noRd
.stage3_global_refinement <- function(prepared_data, core_results, hrf_interface, config) {
  # Check if refinement is enabled
  use_global_refinement <- isTRUE(getOption("fmriparametric.refine_global", TRUE))
  if (!use_global_refinement || !config$global_refinement || config$global_passes <= 0) {
    return(list(
      theta_current = core_results$theta_current,
      amplitudes = core_results$amplitudes,
      r_squared = core_results$r_squared,
      residuals = core_results$residuals,
      intercepts = core_results$intercepts,
      convergence_info = list()
    ))
  }
  
  verbose <- config$verbose
  if (verbose) cat("\n-> Stage 3: Global iterative refinement...\n")
  
  # Initialize tracking variables
  theta_current <- core_results$theta_current
  amplitudes <- core_results$amplitudes
  r_squared <- core_results$r_squared
  
  best_theta <- theta_current
  best_amplitudes <- amplitudes
  best_r2 <- mean(r_squared)
  
  convergence_info <- list(
    global_iterations = 0,
    converged = FALSE
  )
  
  # Use non-parallel processing for now
  process_function <- .parametric_engine
  
  # Refinement loop
  for (iter in seq_len(config$global_passes)) {
    if (verbose) cat(sprintf("  Iteration %d/%d: ", iter, config$global_passes))
    
    # Store previous parameters
    theta_prev <- theta_current
    r_squared_prev <- r_squared
    
    # Re-center globally
    theta_center <- apply(theta_current, 2, median)
    
    # Ensure within bounds - use active bounds from interface
    theta_bounds <- if (!is.null(hrf_interface$active_bounds)) {
      hrf_interface$active_bounds
    } else if (!is.null(core_results$core_result$theta_bounds)) {
      core_results$core_result$theta_bounds
    } else {
      NULL
    }
    
    if (!is.null(theta_bounds)) {
      eps <- 1e-6
      theta_center <- pmax(theta_bounds$lower + eps,
                           pmin(theta_center, theta_bounds$upper - eps))
    }
    
    # Re-run engine with new center
    args <- list(
      Y_proj = prepared_data$inputs$Y_proj,
      S_target_proj = prepared_data$inputs$S_target_proj,
      hrf_eval_times = prepared_data$inputs$hrf_eval_times,
      hrf_interface = hrf_interface,
      theta_seed = theta_center,
      theta_bounds = theta_bounds,
      lambda_ridge = config$lambda_ridge,
      baseline_model = config$baseline_model
    )
    
    # Parallel config would be added here when implemented
    
    iter_result <- do.call(process_function, args)
    
    # Update estimates
    theta_current <- iter_result$theta_hat
    
    # Ensure theta_current remains a matrix (fix for single voxel case)
    if (!is.matrix(theta_current) || length(dim(theta_current)) != 2) {
      theta_current <- matrix(theta_current, 
                              nrow = prepared_data$n_vox, 
                              ncol = length(hrf_interface$parameter_names))
    }
    
    # Enforce bounds on the refined parameters
    if (!is.null(theta_bounds)) {
      lower_mat <- matrix(theta_bounds$lower, 
                          nrow = nrow(theta_current), 
                          ncol = ncol(theta_current), 
                          byrow = TRUE)
      upper_mat <- matrix(theta_bounds$upper, 
                          nrow = nrow(theta_current), 
                          ncol = ncol(theta_current), 
                          byrow = TRUE)
      theta_current <- pmax(lower_mat, pmin(theta_current, upper_mat))
    }
    
    amplitudes <- iter_result$beta0
    r_squared_new <- iter_result$r_squared
    
    # Check convergence using unified criteria
    convergence_config <- .create_convergence_config(
      param_tol = config$convergence_epsilon,
      r2_tol = 1e-5
    )
    
    conv_check <- .check_convergence(
      current = theta_current,
      previous = theta_prev,
      config = convergence_config
    )
    
    mean_r2_change <- mean(r_squared_new) - mean(r_squared_prev)
    
    if (mean_r2_change < -convergence_config$r2_tol) {
      if (verbose) cat("No improvement, rolling back\n")
      theta_current <- theta_prev
      amplitudes <- core_results$amplitudes
      r_squared <- r_squared_prev
      break
    } else {
      r_squared <- r_squared_new
      if (mean(r_squared) > best_r2) {
        best_r2 <- mean(r_squared)
        best_theta <- theta_current
        best_amplitudes <- amplitudes
      }
      
      if (verbose) {
        max_change <- max(abs(theta_current - theta_prev))
        cat(sprintf("Max Deltatheta = %.4f, Mean R^2 = %.3f (Delta = %+.4f)\n",
                    max_change, mean(r_squared), mean_r2_change))
      }
      
      if (conv_check$converged) {
        if (verbose) cat("  [OK] Converged!\n")
        convergence_info$converged <- TRUE
        convergence_info$convergence_reason <- conv_check$reason
        break
      }
    }
  }
  
  convergence_info$global_iterations <- iter
  
  list(
    theta_current = best_theta,
    amplitudes = best_amplitudes,
    r_squared = r_squared,
    residuals = core_results$residuals,  # Pass through from core
    intercepts = core_results$intercepts,  # Pass through from core
    se_amplitudes = core_results$se_amplitudes,  # Pass through if exists
    convergence_info = convergence_info
  )
}

#' Stage 4: Apply tiered refinement
#'
#' @param prepared_data Output from stage 0
#' @param refined_results Output from stage 3
#' @param hrf_interface HRF interface object
#' @param config Configuration list
#' @return List with tiered refinement results
#' @noRd
.stage4_tiered_refinement <- function(prepared_data, refined_results, hrf_interface, config) {
  if (config$tiered_refinement == "none") {
    return(list(
      theta_current = refined_results$theta_current,
      amplitudes = refined_results$amplitudes,
      se_amplitudes = refined_results$se_amplitudes,  # Pass through
      r_squared = refined_results$r_squared,
      residuals = refined_results$residuals,
      intercepts = refined_results$intercepts,
      refinement_info = list(),
      se_theta = NULL,
      convergence_info = refined_results$convergence_info
    ))
  }
  
  verbose <- config$verbose
  if (verbose) cat("\n-> Stage 4: Tiered voxel refinement...\n")
  
  theta_current <- refined_results$theta_current
  amplitudes <- refined_results$amplitudes
  r_squared <- refined_results$r_squared
  n_vox <- prepared_data$n_vox
  
  # Classify voxels using R^2-only mode (Stage 5 will compute real SEs)
  voxel_classes <- .classify_voxels_for_refinement(
    r_squared = r_squared,
    se_theta = NULL,  # NULL triggers R^2-only classification
    thresholds = config$refinement_thresholds
  )
  
  if (verbose) {
    cat(sprintf("  Voxel classification:\n"))
    cat(sprintf("    Easy (high R^2): %d voxels\n", sum(voxel_classes == "easy")))
    cat(sprintf("    Moderate: %d voxels\n", sum(voxel_classes == "moderate")))
    cat(sprintf("    Hard (low R^2): %d voxels\n", sum(voxel_classes == "hard")))
  }
  
  refinement_info <- list(
    classification = voxel_classes,
    n_easy = sum(voxel_classes == "easy"),
    n_moderate = sum(voxel_classes == "moderate"),
    n_hard = sum(voxel_classes == "hard")
  )
  
  # Apply refinements based on classification
  
  # Moderate voxels: Local re-centering
  moderate_idx <- which(voxel_classes == "moderate")
  if (length(moderate_idx) > 0 && config$tiered_refinement %in% c("moderate", "aggressive")) {
    if (verbose) cat("  Refining moderate voxels with local re-centering...\n")
    
    moderate_result <- .refine_moderate_voxels_grouped(
      voxel_idx = moderate_idx,
      Y_proj = prepared_data$inputs$Y_proj,
      S_target_proj = prepared_data$inputs$S_target_proj,
      theta_current = theta_current,
      r_squared = r_squared,
      hrf_interface = hrf_interface,
      hrf_eval_times = prepared_data$inputs$hrf_eval_times,
      theta_bounds = if (!is.null(hrf_interface$active_bounds)) hrf_interface$active_bounds else config$theta_bounds,
      lambda_ridge = config$lambda_ridge
    )
    
    # Update results
    theta_current[moderate_idx, ] <- moderate_result$theta_refined
    amplitudes[moderate_idx] <- moderate_result$amplitudes
    refinement_info$moderate_refined <- length(moderate_idx)
    refinement_info$moderate_improved <- moderate_result$n_improved
  }
  
  # Hard voxels: Gauss-Newton optimization
  hard_idx <- which(voxel_classes == "hard")
  if (length(hard_idx) > 0 && config$tiered_refinement == "aggressive") {
    if (verbose) cat("  Refining hard voxels with Gauss-Newton optimization...\n")
    
    # Create queue labels for Gauss-Newton refinement
    queue_labels <- rep("easy", ncol(prepared_data$inputs$Y_proj))
    queue_labels[hard_idx] <- "hard_GN"
    
    # Get the projected stimulus matrix
    S_target_proj <- prepared_data$inputs$S_target_proj
    
    # Ensure it's a matrix
    if (is.null(dim(S_target_proj))) {
      S_target_proj <- matrix(S_target_proj, ncol = 1)
    }
    
    # Ensure theta_current is a proper matrix
    if (!is.matrix(theta_current)) {
      theta_current <- matrix(theta_current, ncol = length(hrf_interface$parameter_names))
    }
    
    hard_result <- .gauss_newton_refinement(
      theta_hat_voxel = theta_current,
      r2_voxel = r_squared,
      Y_proj = prepared_data$inputs$Y_proj,
      S_target_proj = S_target_proj,
      scan_times = seq_len(nrow(prepared_data$inputs$Y_proj)),
      hrf_eval_times = prepared_data$inputs$hrf_eval_times,
      hrf_interface = hrf_interface,
      theta_bounds = if (!is.null(hrf_interface$active_bounds)) hrf_interface$active_bounds else config$theta_bounds,
      queue_labels = queue_labels,
      max_iter_gn = config$refinement_thresholds$gauss_newton_maxiter,
      verbose = FALSE
    )
    
    # Update results - only update the hard voxels
    theta_current[hard_idx, ] <- hard_result$theta_hat[hard_idx, ]
    # Recompute amplitudes for hard voxels
    hard_amplitudes <- .compute_amplitudes_for_voxels(
      voxel_idx = hard_idx,
      theta_hat = theta_current,
      Y_proj = prepared_data$inputs$Y_proj,
      S_target_proj = prepared_data$inputs$S_target_proj,
      hrf_interface = hrf_interface,
      hrf_eval_times = prepared_data$inputs$hrf_eval_times
    )
    amplitudes[hard_idx] <- hard_amplitudes
    
    refinement_info$hard_refined <- length(hard_idx)
    refinement_info$n_converged <- hard_result$n_converged
    refinement_info$n_improved <- hard_result$n_improved
  }
  
  # Recompute R-squared after refinement
  if (refinement_info$n_moderate > 0 || refinement_info$n_hard > 0) {
    if (verbose) cat("  Recomputing fit quality metrics...\n")
    
    # Generate HRF basis using refined parameters
    n_vox <- ncol(prepared_data$inputs$Y_proj)
    n_time <- nrow(prepared_data$inputs$Y_proj)
    
    # Vectorized HRF computation
    S_target_refined <- vapply(seq_len(n_vox), function(v) {
      hrf_interface$hrf_function(prepared_data$inputs$hrf_eval_times, theta_current[v, ])
    }, numeric(length(prepared_data$inputs$hrf_eval_times)))
    
    # Vectorized convolution
    if (ncol(prepared_data$inputs$S_target_proj) == 1) {
      # Single stimulus column
      Y_conv <- .fast_batch_convolution(
        signal = prepared_data$inputs$S_target_proj[, 1],
        kernels = S_target_refined,
        output_length = n_time
      )
    } else {
      # Multiple stimulus columns - sum convolutions
      Y_conv <- matrix(0, n_time, n_vox)
      for (j in seq_len(ncol(prepared_data$inputs$S_target_proj))) {
        Y_conv_j <- .fast_batch_convolution(
          signal = prepared_data$inputs$S_target_proj[, j],
          kernels = S_target_refined,
          output_length = n_time
        )
        Y_conv <- Y_conv + Y_conv_j
      }
    }
    
    # Scale by amplitudes
    Y_fitted <- sweep(Y_conv, 2, amplitudes, "*")
    
    # Add intercept if present
    if (.is_intercept_baseline(config$baseline_model)) {
      # Intercepts should be available from the parametric engine results
      if (!is.null(refined_results$intercepts)) {
        Y_fitted <- sweep(Y_fitted, 2, refined_results$intercepts, "+")
      }
    }
    
    # Calculate R-squared
    r_squared <- .compute_r_squared(
      prepared_data$inputs$Y_proj, 
      Y_fitted, 
      has_intercept = .is_intercept_baseline(config$baseline_model)
    )
    
    if (verbose) {
      cat(sprintf("  Final fit after refinement: Mean R^2 = %.3f\n", mean(r_squared)))
    }
  }
  
  # CRITICAL: Enforce bounds on final parameters after all refinement
  # This ensures no parameters violate user-specified bounds
  # Use active_bounds from interface if available, otherwise use config
  final_bounds <- if (!is.null(hrf_interface$active_bounds)) {
    hrf_interface$active_bounds
  } else if (!is.null(config$theta_bounds)) {
    config$theta_bounds
  } else {
    NULL
  }
  
  if (!is.null(final_bounds)) {
    lower_mat <- matrix(final_bounds$lower, 
                        nrow = nrow(theta_current), 
                        ncol = ncol(theta_current), 
                        byrow = TRUE)
    upper_mat <- matrix(final_bounds$upper, 
                        nrow = nrow(theta_current), 
                        ncol = ncol(theta_current), 
                        byrow = TRUE)
    theta_current <- pmax(lower_mat, pmin(theta_current, upper_mat))
  }
  
  # Return fitted values if we computed them
  result <- list(
    theta_current = theta_current,
    amplitudes = amplitudes,
    r_squared = r_squared,
    residuals = refined_results$residuals,  # Pass through
    intercepts = refined_results$intercepts,  # Pass through
    refinement_info = refinement_info,
    se_theta = NULL,  # Let Stage 5 compute real SEs
    se_amplitudes = refined_results$se_amplitudes,  # Pass through if exists
    convergence_info = refined_results$convergence_info
  )
  
  # Add fitted values if computed (overrides residuals from core)
  if (exists("Y_fitted") && !is.null(Y_fitted)) {
    result$fitted_values <- Y_fitted
    # Compute updated residuals
    result$residuals <- prepared_data$inputs$Y_proj - Y_fitted
  }
  
  result
}

#' Stage 5: Compute final standard errors
#'
#' @param tiered_results Output from stage 4
#' @param prepared_data Output from stage 0
#' @param hrf_interface HRF interface object
#' @param config Configuration list
#' @return List with standard errors
#' @noRd
.stage5_statistical_inference <- function(tiered_results, prepared_data, hrf_interface, config) {
  if (!config$compute_se || !is.null(tiered_results$se_theta)) {
    return(tiered_results)
  }
  
  verbose <- config$verbose
  if (verbose) cat("\n-> Stage 5: Computing standard errors...\n")
  
  se_result <- .compute_standard_errors_delta(
    theta_hat = tiered_results$theta_current,
    beta0 = tiered_results$amplitudes,
    Y_proj = prepared_data$inputs$Y_proj,
    S_target_proj = prepared_data$inputs$S_target_proj,
    hrf_interface = hrf_interface,
    hrf_eval_times = prepared_data$inputs$hrf_eval_times
  )
  
  tiered_results$se_theta <- se_result$se_theta_hat
  tiered_results$se_amplitudes <- se_result$se_beta0
  
  tiered_results
}

#' Package results into final output object
#'
#' @param final_results Combined results from all stages
#' @param prepared_data Output from stage 0
#' @param hrf_interface HRF interface object
#' @param config Configuration list
#' @param total_time Total computation time
#' @return parametric_hrf_fit object
#' @noRd
.package_final_results <- function(final_results, prepared_data, hrf_interface, config, total_time) {
  n_vox <- prepared_data$n_vox
  n_time <- prepared_data$n_time
  
  # Compute residuals if not already available
  residuals <- NULL
  rmse <- NULL
  
  if (!is.null(final_results$residuals)) {
    # Use residuals from the engine if available
    residuals <- final_results$residuals
  } else if (!is.null(final_results$fitted_values)) {
    # Compute residuals from fitted values
    residuals <- prepared_data$inputs$Y_proj - final_results$fitted_values
  } else {
    # Need to compute fitted values and residuals
    # This replicates the computation from Stage 4
    # Use the hrf_interface that was passed as a parameter
    
    # Generate HRF basis using final parameters
    S_target_refined <- vapply(seq_len(n_vox), function(v) {
      hrf_interface$hrf_function(prepared_data$inputs$hrf_eval_times, 
                                 final_results$theta_current[v, ])
    }, numeric(length(prepared_data$inputs$hrf_eval_times)))
    
    # Vectorized convolution
    if (ncol(prepared_data$inputs$S_target_proj) == 1) {
      Y_conv <- .fast_batch_convolution(
        signal = prepared_data$inputs$S_target_proj[, 1],
        kernels = S_target_refined,
        output_length = n_time
      )
    } else {
      Y_conv <- matrix(0, n_time, n_vox)
      for (j in seq_len(ncol(prepared_data$inputs$S_target_proj))) {
        Y_conv_j <- .fast_batch_convolution(
          signal = prepared_data$inputs$S_target_proj[, j],
          kernels = S_target_refined,
          output_length = n_time
        )
        Y_conv <- Y_conv + Y_conv_j
      }
    }
    
    # Scale by amplitudes
    fitted_values <- sweep(Y_conv, 2, final_results$amplitudes, "*")
    
    # Add intercept if present
    if (.is_intercept_baseline(config$baseline_model) &&
        !is.null(final_results$intercepts)) {
      fitted_values <- sweep(fitted_values, 2, final_results$intercepts, "+")
    }
    
    # Compute residuals
    residuals <- prepared_data$inputs$Y_proj - fitted_values
  }
  
  # Compute RMSE if we have residuals
  if (!is.null(residuals)) {
    # RMSE per voxel
    rmse <- sqrt(colMeans(residuals^2))
  }
  
  # Create fit quality object
  fit_quality <- list(
    r_squared = final_results$r_squared,
    mean_r2 = mean(final_results$r_squared),
    min_r2 = min(final_results$r_squared),
    max_r2 = max(final_results$r_squared),
    rmse = rmse,
    mean_rmse = if (!is.null(rmse)) mean(rmse) else NA
  )
  
  # Create metadata
  metadata <- list(
    call = config$call,
    n_voxels = n_vox,
    n_timepoints = prepared_data$n_time,
    parametric_model = config$parametric_model,
    theta_seed = config$theta_seed,
    theta_bounds = final_results$theta_bounds,
    settings = list(
      global_refinement = config$global_refinement,
      global_passes = config$global_passes,
      kmeans_refinement = config$kmeans_refinement,
      kmeans_k = config$kmeans_k,
      tiered_refinement = config$tiered_refinement,
      parallel = config$parallel,
      n_cores = config$n_cores,
      safety_mode = config$safety_mode
    ),
    timing = list(
      total_seconds = total_time,
      voxels_per_second = n_vox / total_time
    ),
    version = as.character(utils::packageVersion("fmriparametric"))
  )
  
  # Ensure convergence_info is a list
  convergence_info <- final_results$convergence_info
  if (is.null(convergence_info)) {
    convergence_info <- list()
  }
  
  # Compute HRF shapes for each voxel
  hrf_shape <- NULL
  if (!is.null(hrf_interface) && !is.null(hrf_interface$hrf_function)) {
    # Standard time grid for HRF evaluation (0 to 24 seconds, 0.1s resolution)
    time_grid <- seq(0, 24, by = 0.1)
    n_voxels <- nrow(final_results$theta_current)
    hrf_curves <- matrix(NA_real_, nrow = length(time_grid), ncol = n_voxels)
    
    # Compute HRF curve for each voxel
    for (i in seq_len(n_voxels)) {
      tryCatch({
        hrf_curves[, i] <- hrf_interface$hrf_function(time_grid, final_results$theta_current[i, ])
      }, error = function(e) {
        # If HRF computation fails for a voxel, leave as NA
        warning("Failed to compute HRF for voxel ", i, ": ", e$message, call. = FALSE)
      })
    }
    
    hrf_shape <- list(
      time_grid = time_grid,
      curves = hrf_curves
    )
  }

  # Create and return parametric_hrf_fit object
  fit <- new_parametric_hrf_fit(
    estimated_parameters = final_results$theta_current,
    amplitudes = as.numeric(final_results$amplitudes),
    parameter_names = if (!is.null(colnames(final_results$theta_current))) {
      colnames(final_results$theta_current)
    } else if (!is.null(hrf_interface) && !is.null(hrf_interface$parameter_names)) {
      hrf_interface$parameter_names
    } else {
      paste0("param", seq_len(ncol(final_results$theta_current)))
    },
    parametric_model = config$parametric_model,
    r_squared = final_results$r_squared,
    residuals = residuals,
    parameter_ses = final_results$se_theta,
    convergence_info = convergence_info,
    metadata = metadata,
    hrf_shape = hrf_shape
  )
  
  fit$standard_errors <- final_results$se_theta
  fit$se_amplitudes <- final_results$se_amplitudes
  fit$fit_quality <- fit_quality
  fit$refinement_info <- final_results$refinement_info
  
  fit
}

# ========== Stage-specific helper functions ==========

#' Refine moderate voxels using grouped local recentering
#'
#' Groups voxels with identical parameters and refines them together
#' using local Taylor expansion. This is more efficient than individual
#' voxel refinement.
#'
#' @noRd
.refine_moderate_voxels_grouped <- function(
  voxel_idx, Y_proj, S_target_proj, theta_current, r_squared,
  hrf_interface, hrf_eval_times, theta_bounds = NULL,
  lambda_ridge = 0.01
) {
  
  if (length(voxel_idx) == 0) {
    return(list(
      theta_refined = theta_current[voxel_idx, , drop = FALSE],
      amplitudes = numeric(0),
      n_improved = 0
    ))
  }
  
  # Use default bounds if not provided
  if (is.null(theta_bounds)) {
    bounds <- hrf_interface$default_bounds()
  } else {
    bounds <- theta_bounds
  }
  
  n_time <- nrow(Y_proj)
  n_params <- length(hrf_interface$parameter_names)
  
  # Extract parameters for moderate voxels
  theta_block <- theta_current[voxel_idx, , drop = FALSE]
  
  # Group voxels with identical parameters
  id <- apply(theta_block, 1, paste, collapse = ",")
  groups <- split(seq_along(id), id)
  
  # Initialize output
  theta_out <- theta_block
  amps_out <- numeric(length(voxel_idx))
  n_improved <- 0
  
  # Process each group
  for (g in groups) {
    # Get shared parameters for this group
    theta_v <- theta_block[g[1], ]
    
    # Construct Taylor basis
    basis <- hrf_interface$taylor_basis(theta_v, hrf_eval_times)
    if (!is.matrix(basis)) {
      basis <- matrix(basis, ncol = n_params + 1)
    }
    
    # Create design matrix via convolution
    X <- sapply(seq_len(ncol(basis)), function(j) {
      conv_full <- stats::convolve(S_target_proj[, 1], rev(basis[, j]), type = "open")
      conv_full[seq_len(n_time)]
    })
    
    # QR decomposition for stability
    qr_decomp <- qr(X)
    Q <- qr.Q(qr_decomp)
    R <- qr.R(qr_decomp)
    
    # Extract data for this group
    Y_block <- Y_proj[, voxel_idx[g], drop = FALSE]
    
    # Solve with ridge regularization
    coeffs <- solve(R + lambda_ridge * diag(ncol(R)), t(Q) %*% Y_block)
    
    # Extract amplitudes
    beta0_new <- as.numeric(coeffs[1, ])
    
    # Compute parameter updates
    delta <- matrix(0, length(beta0_new), n_params)
    valid <- abs(beta0_new) >= 1e-6
    
    if (any(valid)) {
      # Normalize by amplitude to get parameter deltas
      delta[valid, ] <- t(coeffs[2:(n_params + 1), valid, drop = FALSE]) / beta0_new[valid]
    }
    
    # Update parameters
    theta_new <- sweep(delta, 2, theta_v, FUN = "+")
    
    # Apply bounds
    theta_new <- pmax(bounds$lower, pmin(bounds$upper, theta_new))
    
    # Ensure theta_new is always a matrix
    if (!is.matrix(theta_new)) {
      theta_new <- matrix(theta_new, nrow = length(beta0_new), ncol = n_params, byrow = TRUE)
    }
    
    # Compute new R-squared
    fitted <- X %*% coeffs
    r2_new <- .compute_r_squared(Y_block, fitted, has_intercept = FALSE)
    
    # Only keep improvements
    keep <- !is.na(r2_new) & r2_new > r_squared[voxel_idx[g]] & valid
    if (any(keep)) {
      theta_out[g[keep], ] <- theta_new[keep, , drop = FALSE]
      n_improved <- n_improved + sum(keep)
    }
    
    # Store amplitudes
    amps_out[g] <- beta0_new
  }
  
  list(
    theta_refined = theta_out,
    amplitudes = amps_out,
    n_improved = n_improved
  )
}

#' Compute amplitudes for specific voxels
#'
#' Computes optimal amplitudes given HRF parameters
#'
#' @noRd
.compute_amplitudes_for_voxels <- function(
  voxel_idx, theta_hat, Y_proj, S_target_proj,
  hrf_interface, hrf_eval_times
) {
  n_time <- nrow(Y_proj)
  
  amp_fun <- function(v) {
    hrf_vals <- hrf_interface$hrf_function(hrf_eval_times, theta_hat[v, ])
    conv_full <- stats::convolve(S_target_proj[, 1], rev(hrf_vals), type = "open")
    x_pred <- conv_full[seq_len(n_time)]
    as.numeric(crossprod(x_pred, Y_proj[, v])) / sum(x_pred^2)
  }
  
  vapply(voxel_idx, amp_fun, numeric(1))
}