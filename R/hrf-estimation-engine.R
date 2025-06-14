#' HRF estimation engine - modular implementation
#'
#' This is the modular implementation of HRF parameter estimation that uses
#' stage functions instead of a monolithic implementation for better
#' maintainability and testability.
#'
#' @inherit estimate_parametric_hrf params return
#' @keywords internal
.run_hrf_estimation_engine <- function(
  fmri_data,
  event_model,
  parametric_hrf = "lwu",
  # Basic parameters
  theta_seed = NULL,
  theta_bounds = NULL,
  confound_formula = NULL,
  baseline_model = "intercept",
  hrf_eval_times = NULL,
  hrf_span = 30,
  lambda_ridge = 0.01,
  mask = NULL,
  # Global refinement
  global_refinement = TRUE,
  global_passes = 3,
  convergence_epsilon = 0.01,
  # K-means refinement
  kmeans_refinement = FALSE,
  kmeans_k = 5,
  kmeans_passes = 2,
  # Tiered refinement
  tiered_refinement = c("none", "moderate", "aggressive"),
  refinement_thresholds = list(
    r2_easy = 0.7,
    r2_hard = 0.3,
    se_low = 0.3,
    se_high = 0.7,
    gauss_newton_maxiter = 10
  ),
  # Parallel processing
  parallel = FALSE,
  n_cores = NULL,
  # Output options
  compute_se = TRUE,
  # Safety and diagnostics
  safety_mode = c("balanced", "maximum", "performance"),
  progress = TRUE,
  verbose = TRUE
) {
  
  # Start timing
  total_start <- Sys.time()
  
  # Match arguments
  tiered_refinement <- match.arg(tiered_refinement)
  safety_mode <- match.arg(safety_mode)
  
  # Validate basic inputs
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("verbose must be logical TRUE/FALSE", call. = FALSE)
  }
  if (!is.logical(progress) || length(progress) != 1) {
    stop("progress must be logical TRUE/FALSE", call. = FALSE)
  }
  
  # Initialize progress tracking
  if (progress && verbose) {
    cat("╔══════════════════════════════════════════════════════════════╗\n")
    cat("║          Parametric HRF estimation in progress (v2)          ║\n")
    cat("╚══════════════════════════════════════════════════════════════╝\n\n")
  }
  
  # Get base HRF interface from registry
  base_interface <- .get_hrf_interface(parametric_hrf)
  
  # Create HRF interface with user-specified bounds
  if (parametric_hrf == "lwu") {
    # For LWU, we need to wrap the functions to use the specified bounds
    final_bounds <- if (!is.null(theta_bounds)) {
      # Ensure sigma lower bound is at least 0.051
      if (theta_bounds$lower[2] < 0.051) {
        warning("Adjusting sigma lower bound to 0.051 for numerical stability", call. = FALSE)
        theta_bounds$lower[2] <- 0.051
      }
      theta_bounds
    } else {
      base_interface$default_bounds()
    }
    
    hrf_interface <- list(
      hrf_function = function(t, params_vector, ...) {
        base_interface$hrf_function(t, params_vector, bounds = final_bounds, ...)
      },
      taylor_basis = function(params_vector0, t_hrf_eval, ...) {
        base_interface$taylor_basis(params_vector0, t_hrf_eval, bounds = final_bounds, ...)
      },
      parameter_names = base_interface$parameter_names,
      default_seed = base_interface$default_seed,
      default_bounds = base_interface$default_bounds,
      active_bounds = final_bounds
    )
  } else {
    # For other models, use the base interface as-is
    hrf_interface <- base_interface
    hrf_interface$active_bounds <- theta_bounds
  }
  
  # Package all configuration into a single list
  config <- list(
    # Original parameters
    parametric_hrf = parametric_hrf,
    theta_seed = theta_seed,
    theta_bounds = theta_bounds,
    confound_formula = confound_formula,
    baseline_model = baseline_model,
    hrf_eval_times = hrf_eval_times,
    hrf_span = hrf_span,
    lambda_ridge = lambda_ridge,
    mask = mask,
    global_refinement = global_refinement,
    global_passes = global_passes,
    convergence_epsilon = convergence_epsilon,
    kmeans_refinement = kmeans_refinement,
    kmeans_k = kmeans_k,
    kmeans_passes = kmeans_passes,
    tiered_refinement = tiered_refinement,
    refinement_thresholds = refinement_thresholds,
    parallel = parallel,
    n_cores = n_cores,
    compute_se = compute_se,
    safety_mode = safety_mode,
    progress = progress,
    verbose = verbose,
    # Additional metadata
    call = match.call()
  )
  
  # Execute stages in sequence
  tryCatch({
    # Stage 0: Validation and data preparation
    stage0_results <- .stage0_validate_and_prepare(
      fmri_data = fmri_data,
      event_model = event_model,
      hrf_interface = hrf_interface,
      config = config
    )
    
    # Stage 1: Parameter initialization
    stage1_results <- .stage1_initialize_parameters(
      prepared_data = stage0_results,
      hrf_interface = hrf_interface,
      config = config
    )
    
    # Stage 2: Core estimation
    stage2_results <- .stage2_core_estimation(
      prepared_data = stage0_results,
      init_params = stage1_results,
      hrf_interface = hrf_interface,
      config = config
    )
    
    # Stage 3: Global refinement
    stage3_results <- .stage3_global_refinement(
      prepared_data = stage0_results,
      core_results = stage2_results,
      hrf_interface = hrf_interface,
      config = config
    )
    
    # Stage 4: Tiered refinement
    # Update config with the bounds from stage 1 (in case they were defaulted)
    config_with_bounds <- config
    config_with_bounds$theta_bounds <- stage1_results$theta_bounds
    
    stage4_results <- .stage4_tiered_refinement(
      prepared_data = stage0_results,
      refined_results = stage3_results,
      hrf_interface = hrf_interface,
      config = config_with_bounds
    )
    
    # Stage 5: Statistical inference
    final_results <- .stage5_statistical_inference(
      tiered_results = stage4_results,
      prepared_data = stage0_results,
      hrf_interface = hrf_interface,
      config = config
    )
    
    # Clean up parallel backend if used
    if (!is.null(stage2_results$parallel_config)) {
      stage2_results$parallel_config$cleanup()
    }
    
  }, error = function(e) {
    # Clean up on error
    if (exists("stage2_results") && !is.null(stage2_results$parallel_config)) {
      stage2_results$parallel_config$cleanup()
    }
    stop(e)
  })
  
  # Calculate total time
  total_time <- as.numeric(difftime(Sys.time(), total_start, units = "secs"))
  
  if (verbose) {
    cat("\n╔══════════════════════════════════════════════════════════════╗\n")
    cat("║                    ESTIMATION COMPLETE                       ║\n")
    cat("╚══════════════════════════════════════════════════════════════╝\n")
    cat(sprintf("Total time: %.2f seconds (%.0f voxels/second)\n",
                total_time, stage0_results$n_vox / total_time))
    cat(sprintf("Final mean R²: %.3f\n", mean(final_results$r_squared)))
  }
  
  # Package and return results
  .package_final_results(
    final_results = final_results,
    prepared_data = stage0_results,
    hrf_interface = hrf_interface,
    config = config,
    total_time = total_time
  )
}

#' Create HRF interface (helper function)
#' @keywords internal
.create_hrf_interface <- function(model = "lwu") {
  .get_hrf_interface(model)
}