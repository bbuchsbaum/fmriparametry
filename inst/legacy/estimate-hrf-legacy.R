#' Legacy implementation of parametric HRF estimation
#'
#' This is the original monolithic implementation of HRF parameter estimation.
#' It is being preserved during the transition to the refactored modular version.
#' 
#' @inherit estimate_parametric_hrf params return
#' @keywords internal
.estimate_hrf_legacy <- function(
  fmri_data,
  event_model,
  parametric_hrf = "lwu",
  # Basic parameters (Sprint 1)
  theta_seed = NULL,
  theta_bounds = NULL,
  confound_formula = NULL,
  baseline_model = "intercept",
  hrf_eval_times = NULL,
  hrf_span = 30,
  lambda_ridge = 0.01,
  mask = NULL,
  # Global refinement (Sprint 2)
  global_refinement = TRUE,
  global_passes = 3,
  convergence_epsilon = 0.01,
  # K-means refinement (Sprint 3)
  kmeans_refinement = FALSE,
  kmeans_k = 5,
  kmeans_passes = 2,
  # Tiered refinement (Sprint 3)
  tiered_refinement = c("none", "moderate", "aggressive"),
  refinement_thresholds = list(
    r2_easy = 0.7,
    r2_hard = 0.3,
    se_low = 0.3,
    se_high = 0.7,
    gauss_newton_maxiter = 10
  ),
  # Parallel processing (Sprint 3)
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

  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("verbose must be logical TRUE/FALSE", call. = FALSE)
  }
  if (!is.logical(progress) || length(progress) != 1) {
    stop("progress must be logical TRUE/FALSE", call. = FALSE)
  }
  
  # Initialize progress tracking if requested
  if (progress && verbose) {
    cat("╔══════════════════════════════════════════════════════════════╗\n")
    cat("║             Parametric HRF estimation in progress            ║\n")
    cat("╚══════════════════════════════════════════════════════════════╝\n\n")
  }
  
  # ========== STAGE 0: VALIDATION (Rock-solid safety) ==========
  if (verbose) cat("→ Stage 0: Input validation and safety checks...\n")
  
  # Validate inputs based on safety mode
  validation_level <- switch(safety_mode,
    maximum = "comprehensive",
    balanced = "standard",
    performance = "minimal"
  )

  # Run input validation helpers depending on level

  caller <- "estimate_parametric_hrf"

  if (validation_level == "comprehensive") {
    .rock_solid_validate_inputs(
      fmri_data = fmri_data,
      event_model = event_model,
      parametric_hrf = parametric_hrf,
      theta_seed = theta_seed,
      theta_bounds = theta_bounds,
      hrf_span = hrf_span,
      lambda_ridge = lambda_ridge,
      recenter_global_passes = global_passes,
      recenter_epsilon = convergence_epsilon,
      r2_threshold = 0.1,
      mask = mask,
      verbose = verbose,
      caller = caller
    )
  } else if (validation_level == "standard") {
    fmri_ok <- .validate_fmri_data(fmri_data, caller)
    .validate_event_model(event_model, fmri_ok$n_time, caller)
    theta_bounds <- .validate_theta_bounds(
      theta_bounds,
      length(.create_hrf_interface(parametric_hrf)$parameter_names),
      .create_hrf_interface(parametric_hrf)$parameter_names,
      caller
    )
  } else {
    if (is.null(fmri_data)) {
      stop(caller, ": fmri_data cannot be NULL.", call. = FALSE)
    }
    if (is.null(event_model)) {
      stop(caller, ": event_model cannot be NULL.", call. = FALSE)
    }
  }

  # ========== STAGE 0A: Data preparation ==========
  
  # Create HRF interface
  hrf_interface <- .create_hrf_interface(parametric_hrf)
  
  # Prepare data
  if (verbose) cat("→ Preparing data matrices...\n")
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
  
  # Initialize results storage
  theta_current <- matrix(NA_real_, n_vox, length(hrf_interface$parameter_names))
  colnames(theta_current) <- hrf_interface$parameter_names
  
  # ========== STAGE 1: INITIALIZATION ==========
  if (verbose) cat("\n→ Stage 1: Parameter initialization...\n")
  
  # Handle theta_seed
  if (is.null(theta_seed)) {
    theta_seed <- hrf_interface$default_seed()
  } else if (identical(theta_seed, "data_driven")) {
    if (verbose) cat("  Computing data-driven initialization...\n")
    theta_seed <- .compute_data_driven_seed(
      Y = inputs$Y_proj,
      S = inputs$S_target_proj,
      hrf_interface = hrf_interface,
      theta_bounds = theta_bounds
    )
  } else {
    if (!is.numeric(theta_seed)) {
      stop("theta_seed must be numeric", call. = FALSE)
    }
    expected_len <- length(hrf_interface$parameter_names)
    if (length(theta_seed) != expected_len) {
      stop(sprintf("theta_seed must have length %d", expected_len), call. = FALSE)
    }
  }
  
  # Handle theta_bounds
  if (is.null(theta_bounds)) {
    theta_bounds <- hrf_interface$default_bounds()
  } else {
    if (!is.list(theta_bounds) ||
        !all(c("lower", "upper") %in% names(theta_bounds))) {
      stop("theta_bounds must have both 'lower' and 'upper' components",
           call. = FALSE)
    }
    expected_len <- length(hrf_interface$parameter_names)
    if (length(theta_bounds$lower) != expected_len ||
        length(theta_bounds$upper) != expected_len) {
      stop(sprintf("theta_bounds components must have length %d",
                   expected_len), call. = FALSE)
    }
    if (!is.numeric(theta_bounds$lower) || !is.numeric(theta_bounds$upper)) {
      stop("theta_bounds values must be numeric", call. = FALSE)
    }
  }
  
  # Ensure theta_seed is within safe bounds for numerical derivatives
  eps <- c(0.01, 0.01, 0.01)
  theta_seed <- pmax(theta_bounds$lower + eps,
                     pmin(theta_seed, theta_bounds$upper - eps))
  
  # Set initial values
  for (j in seq_len(ncol(theta_current))) {
    theta_current[, j] <- theta_seed[j]
  }
  
  # K-means initialization if requested
  kmeans_iterations <- 0
  if (kmeans_refinement && kmeans_k > 1) {
    if (verbose) cat("  Performing K-means clustering for initialization...\n")
    kmeans_result <- .perform_kmeans_initialization(
      Y = inputs$Y_proj,
      S = inputs$S_target_proj,
      k = kmeans_k,
      hrf_interface = hrf_interface,
      theta_bounds = theta_bounds
    )
    kmeans_iterations <- kmeans_result$iterations
    # Update seeds based on clusters
    for (k in seq_len(kmeans_k)) {
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
  
  # ========== STAGE 2: CORE ESTIMATION ==========
  if (verbose) cat("\n→ Stage 2: Core parametric estimation...\n")
  
  # Setup parallel backend if requested
  parallel_config <- NULL
  if (parallel) {
    parallel_config <- .setup_parallel_backend(n_cores = n_cores, verbose = verbose)
    n_cores <- parallel_config$n_cores
    process_function <- .parametric_engine_parallel
  } else {
    process_function <- .parametric_engine
  }
  
  # Run core engine (with proper parameter handling)
  if (parallel) {
    core_result <- process_function(
      Y_proj = inputs$Y_proj,
      S_target_proj = inputs$S_target_proj,
      hrf_eval_times = inputs$hrf_eval_times,
      hrf_interface = hrf_interface,
      theta_seed = theta_seed,
      theta_bounds = theta_bounds,
      lambda_ridge = lambda_ridge,
      parallel_config = parallel_config
    )
  } else {
    core_result <- process_function(
      Y_proj = inputs$Y_proj,
      S_target_proj = inputs$S_target_proj,
      hrf_eval_times = inputs$hrf_eval_times,
      hrf_interface = hrf_interface,
      theta_seed = theta_seed,
      theta_bounds = theta_bounds,
      lambda_ridge = lambda_ridge
    )
  }
  
  # Update current estimates
  theta_current <- core_result$theta_hat
  amplitudes <- core_result$beta0
  
  # Ensure theta_current is always a matrix
  if (!is.matrix(theta_current) || length(dim(theta_current)) != 2) {
    theta_current <- matrix(theta_current, nrow = n_vox, ncol = length(hrf_interface$parameter_names))
    colnames(theta_current) <- hrf_interface$parameter_names
  }
  
  # Use R-squared from the parametric engine
  r_squared <- core_result$r_squared
  
  if (verbose) {
    cat(sprintf("  Initial fit: Mean R² = %.3f (range: %.3f - %.3f)\n",
                mean(r_squared), min(r_squared), max(r_squared)))
  }
  
  # ========== STAGE 3: ITERATIVE REFINEMENT ==========
  convergence_info <- list()

  use_global_refinement <- isTRUE(getOption("fmriparametric.refine_global", TRUE))
  
  
  if (use_global_refinement && global_refinement && global_passes > 0) {
    # Only initialize these fields if we're actually doing global refinement
    convergence_info$global_iterations <- 0
    convergence_info$converged <- TRUE
    if (verbose) cat("\n→ Stage 3: Global iterative refinement...\n")

    best_theta <- theta_current
    best_amplitudes <- amplitudes
    best_r2 <- mean(r_squared)
    global_iter_count <- 0
    global_converged <- FALSE

    for (iter in seq_len(global_passes)) {
      global_iter_count <- iter
      if (verbose) cat(sprintf("  Iteration %d/%d: ", iter, global_passes))

      # Store previous parameters
      theta_prev <- theta_current
      amplitudes_prev <- amplitudes
      r_squared_prev <- r_squared
      
      # Re-center globally with bounds enforcement
      if (!is.matrix(theta_current) || length(dim(theta_current)) != 2) {
        theta_current <- matrix(
          theta_current,
          nrow = n_vox,
          ncol = length(hrf_interface$parameter_names)
        )
        colnames(theta_current) <- hrf_interface$parameter_names
      }

      theta_center <- apply(theta_current, 2, median)
      
      # Ensure theta_center is within bounds before using as seed
      if (!is.null(theta_bounds)) {
        theta_center <- pmax(theta_bounds$lower, pmin(theta_center, theta_bounds$upper))
        # Stay slightly inside bounds to avoid issues in derivative calculations
        eps <- 1e-6
        theta_center <- pmax(theta_bounds$lower + eps,
                             pmin(theta_center, theta_bounds$upper - eps))
      }
      
      # Re-run engine with new center
      if (parallel) {
        iter_result <- process_function(
          Y_proj = inputs$Y_proj,
          S_target_proj = inputs$S_target_proj,
          hrf_eval_times = inputs$hrf_eval_times,
          hrf_interface = hrf_interface,
          theta_seed = theta_center,
          theta_bounds = theta_bounds,
          lambda_ridge = lambda_ridge,
          parallel_config = parallel_config
        )
      } else {
        iter_result <- process_function(
          Y_proj = inputs$Y_proj,
          S_target_proj = inputs$S_target_proj,
          hrf_eval_times = inputs$hrf_eval_times,
          hrf_interface = hrf_interface,
          theta_seed = theta_center,
          theta_bounds = theta_bounds,
          lambda_ridge = lambda_ridge
        )
      }
      
      # Update estimates
      theta_current <- iter_result$theta_hat
      amplitudes <- iter_result$beta0
      
      # Ensure theta_current is always a matrix (critical for apply() calls)
      if (!is.matrix(theta_current) || length(dim(theta_current)) != 2) {
        theta_current <- matrix(theta_current, nrow = n_vox, ncol = length(hrf_interface$parameter_names))
        colnames(theta_current) <- hrf_interface$parameter_names
      }
      
      # Additional safety check - verify dimensions
      if (nrow(theta_current) != n_vox || ncol(theta_current) != length(hrf_interface$parameter_names)) {
        warning("theta_current dimensions incorrect, reshaping...")
        theta_current <- matrix(as.vector(theta_current), nrow = n_vox, ncol = length(hrf_interface$parameter_names))
        colnames(theta_current) <- hrf_interface$parameter_names
      }
      
      # Apply bounds to prevent parameter drift
      if (!is.null(theta_bounds)) {
        theta_current <- pmax(theta_bounds$lower, pmin(theta_current, theta_bounds$upper))
      }

      # Check convergence and fit quality
      max_change <- max(abs(theta_current - theta_prev))
      r_squared_new <- iter_result$r_squared

      mean_r2_change <- mean(r_squared_new) - mean(r_squared_prev)

      if (mean_r2_change < 0) {
        if (verbose) cat("No improvement, rolling back\n")
        theta_current <- theta_prev
        amplitudes <- amplitudes_prev
        r_squared <- r_squared_prev
        convergence_info[[paste0("global_iter_", iter)]] <- list(
          max_param_change = max_change,
          mean_r2 = mean(r_squared_prev),
          r2_improvement = mean_r2_change,
          rolled_back = TRUE
        )
        break
      } else {
        amplitudes <- iter_result$beta0
        r_squared <- r_squared_new

        if (mean(r_squared) > best_r2) {
          best_r2 <- mean(r_squared)
          best_theta <- theta_current
          best_amplitudes <- amplitudes
        }

        if (verbose) {
          cat(sprintf("Max Δθ = %.4f, Mean R² = %.3f (Δ = %+.4f)\n",
                      max_change, mean(r_squared), mean_r2_change))
        }

        convergence_info[[paste0("global_iter_", iter)]] <- list(
          max_param_change = max_change,
          mean_r2 = mean(r_squared),
          r2_improvement = mean_r2_change,
          rolled_back = FALSE
        )

        if (max_change < convergence_epsilon) {
          if (verbose) cat("  ✓ Converged!\n")
          global_converged <- TRUE
          break
        }
      }
    }

    theta_current <- best_theta
    amplitudes <- best_amplitudes
    convergence_info$global_iterations <- global_iter_count
    convergence_info$converged <- global_converged
  }
  
  # ========== STAGE 4: TIERED REFINEMENT ==========
  refinement_info <- list()
  se_theta <- NULL
  
  if (tiered_refinement != "none") {
    if (verbose) cat("\n→ Stage 4: Tiered voxel refinement...\n")
    
    # First compute standard errors if needed for classification
    if (compute_se) {
      se_result <- .compute_standard_errors_delta(
        theta_hat = theta_current,
        beta0 = amplitudes,
        Y_proj = inputs$Y_proj,
        S_target_proj = inputs$S_target_proj,
        hrf_interface = hrf_interface,
        hrf_eval_times = inputs$hrf_eval_times
      )
      se_theta <- se_result$se_theta_hat
    } else {
      se_theta <- matrix(0, n_vox, ncol(theta_current))
    }
    
    # Classify voxels
    voxel_classes <- .classify_voxels_for_refinement(
      r_squared = r_squared,
      se_theta = se_theta,
      thresholds = refinement_thresholds
    )
    
    if (verbose) {
      cat(sprintf("  Voxel classification:\n"))
      cat(sprintf("    Easy (high R²): %d voxels\n", sum(voxel_classes == "easy")))
      cat(sprintf("    Moderate: %d voxels\n", sum(voxel_classes == "moderate")))
      cat(sprintf("    Hard (low R²): %d voxels\n", sum(voxel_classes == "hard")))
    }
    
    # Apply refinement based on classification
    refinement_info <- list(
      classification = voxel_classes,
      n_easy = sum(voxel_classes == "easy"),
      n_moderate = sum(voxel_classes == "moderate"),
      n_hard = sum(voxel_classes == "hard")
    )
    
    # Moderate voxels: Local re-centering
    moderate_idx <- which(voxel_classes == "moderate")
    if (length(moderate_idx) > 0 && tiered_refinement %in% c("moderate", "aggressive")) {
      if (verbose) cat("  Refining moderate voxels with local re-centering...\n")
      
      moderate_result <- .refine_moderate_voxels(
        voxel_idx = moderate_idx,
        Y_proj = inputs$Y_proj,
        S_target_proj = inputs$S_target_proj,
        theta_current = theta_current,
        r_squared = r_squared,
        hrf_interface = hrf_interface,
        hrf_eval_times = inputs$hrf_eval_times,
        theta_bounds = theta_bounds,
        parallel = parallel,
        n_cores = if (parallel) parallel_config$n_cores else 1
      )
      
      # Update results
      theta_current[moderate_idx, ] <- moderate_result$theta_refined
      amplitudes[moderate_idx] <- moderate_result$amplitudes
      refinement_info$moderate_refined <- length(moderate_idx)
    }
    
    # Hard voxels: Gauss-Newton
    hard_idx <- which(voxel_classes == "hard")
    if (length(hard_idx) > 0 && tiered_refinement == "aggressive") {
      if (verbose) cat("  Refining hard voxels with Gauss-Newton optimization...\n")
      
      hard_result <- .refine_hard_voxels(
        voxel_idx = hard_idx,
        Y_proj = inputs$Y_proj,
        S_target_proj = inputs$S_target_proj,
        theta_current = theta_current,
        r_squared = r_squared,
        hrf_interface = hrf_interface,
        hrf_eval_times = inputs$hrf_eval_times,
        theta_bounds = theta_bounds,
        max_iter = refinement_thresholds$gauss_newton_maxiter,
        parallel = parallel,
        n_cores = if (parallel) parallel_config$n_cores else 1
      )
      
      # Update results
      theta_current[hard_idx, ] <- hard_result$theta_refined
      amplitudes[hard_idx] <- hard_result$amplitudes
      refinement_info$hard_refined <- length(hard_idx)
      refinement_info$n_converged <- hard_result$gn_info$n_converged
      refinement_info$n_improved <- hard_result$gn_info$n_improved
      convergence_info$gauss_newton <- list(
        n_hard = length(hard_idx),
        n_converged = hard_result$gn_info$n_converged,
        n_improved = hard_result$gn_info$n_improved,
        mean_iterations = mean(hard_result$gn_info$iteration_counts)
      )
    }
    
    # Recompute final R-squared after refinement
    # Use vectorized approach for efficiency
    
    # Generate HRF basis using refined parameters
    # Use vapply for safe vectorization if hrf_function isn't vectorized
    S_target_refined <- vapply(1:n_vox, function(v) {
      hrf_interface$hrf_function(inputs$hrf_eval_times, theta_current[v, ])
    }, numeric(length(inputs$hrf_eval_times)))
    
    # Vectorized convolution approach
    if (ncol(inputs$S_target) == 1) {
      # Single stimulus column - direct batch convolution
      Y_conv <- .fast_batch_convolution(
        signal = inputs$S_target[, 1],
        kernels = S_target_refined,
        output_length = n_time
      )
    } else {
      # Multiple stimulus columns - sum convolutions
      # First convolve each stimulus with all HRFs
      Y_conv <- matrix(0, n_time, n_vox)
      for (j in seq_len(ncol(inputs$S_target))) {
        Y_conv_j <- .fast_batch_convolution(
          signal = inputs$S_target[, j],
          kernels = S_target_refined,
          output_length = n_time
        )
        Y_conv <- Y_conv + Y_conv_j
      }
    }
    
    # Scale by amplitudes using efficient broadcasting
    Y_fitted <- sweep(Y_conv, 2, amplitudes, "*")
    
    # Add intercept/baseline if the model includes one
    if (!is.null(baseline_model) && baseline_model == "intercept") {
      # The fitted values need to include the intercept term
      # Get the intercept from the last engine result (stored in coeffs)
      # Note: With intercept, coeffs[1,] contains intercepts
      if (exists("core_result") && !is.null(core_result$coeffs)) {
        if (nrow(core_result$coeffs) > ncol(inputs$S_target_proj)) {
          # First row is intercept
          intercepts <- core_result$coeffs[1, ]
          Y_fitted <- sweep(Y_fitted, 2, intercepts, "+")
        }
      }
    }
    
    # Calculate R-squared using proper baseline handling
    r_squared <- .compute_r_squared(inputs$Y_proj, Y_fitted, 
                                    has_intercept = (!is.null(baseline_model) && baseline_model == "intercept"))
    
    if (verbose) {
      cat(sprintf("  Final fit after refinement: Mean R² = %.3f\n", mean(r_squared)))
    }
  }
  
  # ========== STAGE 5: STATISTICAL INFERENCE ==========
  if (verbose && compute_se) cat("\n→ Stage 5: Computing standard errors...\n")
  
  # Compute standard errors if not already done
  if (compute_se && is.null(se_theta)) {
    se_result <- .compute_standard_errors_delta(
      theta_hat = theta_current,
      beta0 = amplitudes,
      Y_proj = inputs$Y_proj,
      S_target_proj = inputs$S_target_proj,
      hrf_interface = hrf_interface,
      hrf_eval_times = inputs$hrf_eval_times
    )
    se_theta <- se_result$se_theta_hat
    se_amplitudes <- se_result$se_beta0
  } else {
    se_theta <- NULL
    se_amplitudes <- NULL
  }
  
  # ========== FINALIZE RESULTS ==========
  total_time <- as.numeric(difftime(Sys.time(), total_start, units = "secs"))
  
  if (verbose) {
    cat("\n╔══════════════════════════════════════════════════════════════╗\n")
    cat("║                    ESTIMATION COMPLETE                       ║\n")
    cat("╚══════════════════════════════════════════════════════════════╝\n")
    cat(sprintf("Total time: %.2f seconds (%.0f voxels/second)\n",
                total_time, n_vox / total_time))
    cat(sprintf("Final mean R²: %.3f\n", mean(r_squared)))
  }
  
  # Create fit quality object
  fit_quality <- list(
    r_squared = r_squared,
    mean_r2 = mean(r_squared),
    min_r2 = min(r_squared),
    max_r2 = max(r_squared),
    rmse = NA  # TODO: Compute RMSE properly using fitted values
  )
  
  # Create comprehensive metadata
  metadata <- list(
    call = match.call(),
    n_voxels = n_vox,
    n_timepoints = n_time,
    hrf_model = parametric_hrf,
    theta_seed = if(is.character(theta_seed)) theta_seed else as.numeric(theta_seed),
    theta_bounds = theta_bounds,
    settings = list(
      global_refinement = global_refinement,
      global_passes = global_passes,
      kmeans_refinement = kmeans_refinement,
      kmeans_k = kmeans_k,
      tiered_refinement = tiered_refinement,
      parallel = parallel,
      n_cores = n_cores,
      safety_mode = safety_mode
    ),
    parallel_info = if (!is.null(parallel_config)) list(
      backend = parallel_config$backend,
      n_cores = parallel_config$n_cores
    ) else list(
      backend = "sequential",
      n_cores = 1
    ),
    timing = list(
      total_seconds = total_time,
      voxels_per_second = n_vox / total_time
    ),
    version = "ultimate_impeccable_v1.0",
    kmeans_info = list(
      applied = kmeans_refinement && kmeans_k > 1,
      n_clusters = kmeans_k,
      total_iterations = kmeans_iterations
    )
  )

  # Clean up parallel backend
  if (!is.null(parallel_config)) {
    parallel_config$cleanup()
  }

  # Ensure theta_current is a proper matrix
  if (!is.matrix(theta_current) || length(dim(theta_current)) != 2) {
    theta_current <- matrix(theta_current, nrow = n_vox, ncol = length(hrf_interface$parameter_names))
  }
  
  # Ensure proper column names
  if (is.null(colnames(theta_current))) {
    colnames(theta_current) <- hrf_interface$parameter_names
  }
  
  # Ensure theta_current is still a matrix before creating fit object
  if (!is.matrix(theta_current) || length(dim(theta_current)) != 2) {
    theta_current <- matrix(theta_current, nrow = n_vox, ncol = length(hrf_interface$parameter_names))
    colnames(theta_current) <- hrf_interface$parameter_names
  }
  
  # Create and return parametric_hrf_fit object
  fit <- new_parametric_hrf_fit(
    estimated_parameters = theta_current,
    amplitudes = as.numeric(amplitudes),
    parameter_names = hrf_interface$parameter_names,
    hrf_model = parametric_hrf,
    r_squared = r_squared,
    residuals = core_result$residuals,
    parameter_ses = se_theta,
    convergence_info = convergence_info,
    metadata = metadata
  )
  
  # Add additional components
  fit$standard_errors <- se_theta
  fit$se_amplitudes <- se_amplitudes
  fit$fit_quality <- fit_quality
  # Only overwrite refinement_info if we have a non-NULL value
  if (!is.null(refinement_info)) {
    fit$refinement_info <- refinement_info
  }
  
  return(fit)
}


# ========== HELPER FUNCTION IMPLEMENTATIONS ==========

# MISSING FUNCTION 1: Parallel engine

.parametric_engine_parallel <- function(Y_proj, S_target_proj, hrf_eval_times,
                                       hrf_interface, theta_seed, theta_bounds,
                                       lambda_ridge = 0.01, parallel_config,
                                       epsilon_beta = 1e-6) {

  n_vox <- ncol(Y_proj)
  n_params <- length(hrf_interface$parameter_names)

  # Processing function for a chunk of voxels
  process_chunk <- function(voxel_idx) {
    res <- .parametric_engine(
      Y_proj = Y_proj[, voxel_idx, drop = FALSE],
      S_target_proj = S_target_proj,
      hrf_eval_times = hrf_eval_times,
      hrf_interface = hrf_interface,
      theta_seed = theta_seed,
      theta_bounds = theta_bounds,
      lambda_ridge = lambda_ridge,
      epsilon_beta = epsilon_beta,
      verbose = FALSE,
      baseline_model = baseline_model
    )
    list(list(indices = voxel_idx, theta_hat = res$theta_hat, beta0 = res$beta0))
  }

  # Run using the generic parallel backend
  chunk_results <- .parallel_voxel_processing(
    voxel_indices = seq_len(n_vox),
    process_function = process_chunk,
    parallel_config = parallel_config,
    chunk_size = "auto",
    progress = FALSE
  )

  # Combine results
  theta_hat <- matrix(NA_real_, n_vox, n_params)
  beta0 <- numeric(n_vox)

  for (res in chunk_results) {
    theta_hat[res$indices, ] <- res$theta_hat
    beta0[res$indices] <- res$beta0
  }

  list(theta_hat = theta_hat, beta0 = beta0)
}

# MISSING FUNCTION 2: Standard errors via Delta method
.compute_standard_errors_delta <- function(theta_hat, beta0, Y_proj, S_target_proj,
                                          hrf_interface, hrf_eval_times) {
  n_vox <- ncol(Y_proj)
  n_params <- length(hrf_interface$parameter_names)
  n_time <- nrow(Y_proj)

  # Ensure theta_hat is a matrix
  if (!is.matrix(theta_hat) || length(dim(theta_hat)) != 2) {
    theta_hat <- matrix(theta_hat, nrow = n_vox, ncol = n_params,
                        byrow = FALSE)
    colnames(theta_hat) <- hrf_interface$parameter_names
  }
  
  basis_list <- vector("list", n_vox)
  hrf_list <- vector("list", n_vox)
  for (v in seq_len(n_vox)) {
    theta_v <- theta_hat[v, ]
    taylor_basis <- hrf_interface$taylor_basis(theta_v, hrf_eval_times)
    if (!is.matrix(taylor_basis)) {
      taylor_basis <- matrix(taylor_basis, ncol = n_params + 1)
    }
    basis_list[[v]] <- taylor_basis[, -1, drop = FALSE]
    hrf_list[[v]] <- hrf_interface$hrf_function(hrf_eval_times, theta_v)
  }

  res <- compute_standard_errors_bulk_cpp(
    basis_list, hrf_list, Y_proj, S_target_proj, as.numeric(beta0)
  )

  return(res)
}

# Removed duplicate functions that are now in utils-hrf-estimation.R

# Removed - now using unified .compute_r_squared from fit-metrics.R

.classify_voxels_for_refinement <- function(r_squared, se_theta, thresholds) {
  classes <- rep("moderate", length(r_squared))
  
  # Easy: high R² and low SE
  easy_mask <- r_squared > thresholds$r2_easy & 
               rowMeans(se_theta) < thresholds$se_low
  classes[easy_mask] <- "easy"
  
  # Hard: low R² or high SE
  hard_mask <- r_squared < thresholds$r2_hard | 
               rowMeans(se_theta) > thresholds$se_high
  classes[hard_mask] <- "hard"
  
  return(classes)
}

# Refinement functions
.refine_moderate_voxels <- function(voxel_idx, Y_proj, S_target_proj,
                                    theta_current, r_squared,
                                    hrf_interface, hrf_eval_times,
                                    theta_bounds = NULL,
                                    parallel = FALSE,
                                    n_cores = 1) {

  if (length(voxel_idx) == 0) {
    return(list(theta_refined = theta_current[voxel_idx, , drop = FALSE],
                amplitudes = numeric(0)))
  }

  if (is.null(theta_bounds)) {
    bounds <- hrf_interface$default_bounds()
  } else {
    bounds <- theta_bounds
  }
  n_time <- nrow(Y_proj)
  n_params <- length(hrf_interface$parameter_names)

  theta_block <- theta_current[voxel_idx, , drop = FALSE]
  id <- apply(theta_block, 1, paste, collapse = ",")
  groups <- split(seq_along(id), id)

  theta_out <- theta_block
  amps_out <- numeric(length(voxel_idx))

  for (g in groups) {
    theta_v <- theta_block[g[1], ]
    basis <- hrf_interface$taylor_basis(theta_v, hrf_eval_times)
    if (!is.matrix(basis)) {
      basis <- matrix(basis, ncol = n_params + 1)
    }
    X <- sapply(seq_len(ncol(basis)), function(j) {
      conv_full <- stats::convolve(S_target_proj[, 1], rev(basis[, j]), type = "open")
      conv_full[seq_len(n_time)]
    })
    qr_decomp <- qr(X)
    Q <- qr.Q(qr_decomp)
    R <- qr.R(qr_decomp)
    Y_block <- Y_proj[, voxel_idx[g], drop = FALSE]
    coeffs <- solve(R + 0.01 * diag(ncol(R)), t(Q) %*% Y_block)
    beta0_new <- as.numeric(coeffs[1, ])

    delta <- matrix(0, length(beta0_new), n_params)
    valid <- abs(beta0_new) >= 1e-6
    if (any(valid)) {
      delta[valid, ] <- t(coeffs[2:(n_params + 1), valid, drop = FALSE]) / beta0_new[valid]
    }
    theta_new <- sweep(delta, 2, theta_v, FUN = "+")
    theta_new <- pmax(bounds$lower, pmin(bounds$upper, theta_new))

    fitted <- X %*% coeffs
    r2_new <- .compute_r_squared(Y_block, fitted, has_intercept = TRUE)

    keep <- !is.na(r2_new) & r2_new > r_squared[voxel_idx[g]] & valid
    theta_out[g[keep], ] <- theta_new[keep, , drop = FALSE]
    amps_out[g] <- beta0_new
  }

  list(theta_refined = theta_out, amplitudes = amps_out)
}

.refine_hard_voxels <- function(voxel_idx, Y_proj, S_target_proj,
                                theta_current, r_squared,
                                hrf_interface, hrf_eval_times,
                                theta_bounds = NULL,
                                max_iter = 5, parallel = FALSE,
                                n_cores = 1) {

  if (length(voxel_idx) == 0) {
    return(list(theta_refined = theta_current[voxel_idx, , drop = FALSE],
                amplitudes = numeric(0)))
  }

  n_vox <- ncol(Y_proj)
  if (is.null(theta_bounds)) {
    bounds <- hrf_interface$default_bounds()
  } else {
    bounds <- theta_bounds
  }
  queue_labels <- rep("easy", n_vox)
  queue_labels[voxel_idx] <- "hard_GN"

  gn <- .gauss_newton_refinement(
    theta_hat_voxel = theta_current,
    r2_voxel = r_squared,
    Y_proj = Y_proj,
    S_target_proj = S_target_proj,
    scan_times = seq_len(nrow(Y_proj)),
    hrf_eval_times = hrf_eval_times,
    hrf_interface = hrf_interface,
    theta_bounds = bounds,
    queue_labels = queue_labels,
    max_iter_gn = max_iter,
    verbose = FALSE
  )

  theta_updated <- gn$theta_hat

  n_time <- nrow(Y_proj)
  amp_fun <- function(v) {
    hrf_vals <- hrf_interface$hrf_function(hrf_eval_times, theta_updated[v, ])
    conv_full <- stats::convolve(S_target_proj[, 1], rev(hrf_vals), type = "open")
    x_pred <- conv_full[seq_len(n_time)]
    as.numeric(crossprod(x_pred, Y_proj[, v])) / sum(x_pred^2)
  }

  amps <- vapply(voxel_idx, amp_fun, numeric(1))

  list(
    theta_refined = theta_updated[voxel_idx, , drop = FALSE],
    amplitudes = amps,
    gn_info = list(
      n_refined = gn$n_refined,
      n_converged = gn$n_converged,
      n_improved = gn$n_improved,
      convergence_status = gn$convergence_status,
      iteration_counts = gn$iteration_counts
    )
  )
}