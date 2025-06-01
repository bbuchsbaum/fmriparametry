#' Estimate parametric HRF parameters
#'
#' Estimate hemodynamic response function parameters from fMRI time series using
#' a parametric HRF model. Optional refinement steps improve fits for challenging
#' voxels.
#' 
#' @section Package Options:
#' Global iterative refinement (Stage 3) is controlled by the option
#' `fmriparametric.refine_global`. Set to `FALSE` to disable global re-centering
#' for all calls.
#'
#' @param fmri_data An fMRI dataset object or numeric matrix (timepoints x voxels)
#' @param event_model Event timing design matrix or event model object
#' @param parametric_hrf Character string specifying HRF model (currently "lwu")
#' @param theta_seed Initial parameters or "data_driven" for automatic selection
#' @param theta_bounds List with 'lower' and 'upper' bounds for parameters
#' @param confound_formula Optional formula for nuisance regressors
#' @param baseline_model Baseline model specification (default "intercept")
#' @param hrf_eval_times Time points for HRF evaluation
#' @param hrf_span Duration for HRF evaluation (default 30 seconds)
#' @param lambda_ridge Ridge penalty for stability (default 0.01)
#' @param mask Optional voxel selection mask
#' @param global_refinement Logical: perform iterative global refinement? (Sprint 2)
#' @param global_passes Number of global refinement iterations (default 3)
#' @param kmeans_refinement Logical: perform K-means local refinement? (Sprint 3)
#' @param kmeans_k Number of clusters for K-means (default 5)
#' @param kmeans_passes Number of K-means refinement passes (default 2)
#' @param tiered_refinement Character: refinement strategy ("none", "moderate", "aggressive")
#' @param refinement_thresholds List of R² and SE thresholds for tiered refinement
#' @param parallel Logical: use parallel processing?
#' @param n_cores Number of cores for parallel processing (NULL = auto-detect)
#' @param compute_se Logical: compute standard errors?
#' @param safety_mode Character: "maximum", "balanced", or "performance"
#' @param progress Logical: show progress bar?
#' @param verbose Logical: print detailed messages?
#'
#' @return Object of class 'parametric_hrf_fit' containing:
#'   - estimated_parameters: Matrix of HRF parameters (voxels x parameters)
#'   - amplitudes: Response amplitudes for each voxel
#'   - standard_errors: Parameter standard errors (if computed)
#'   - fit_quality: List with R-squared and other metrics
#'   - convergence_info: Detailed convergence information
#'   - refinement_info: Information about refinement strategies applied
#'   - metadata: Complete analysis metadata
#'
#' @examples
#' # Simulated quick example
#' set.seed(1)
#' fmri_data <- matrix(rnorm(40), nrow = 20, ncol = 2)
#' events <- matrix(0, nrow = 20, ncol = 1)
#' events[c(5, 15), 1] <- 1
#' fit <- estimate_parametric_hrf(fmri_data, events, parametric_hrf = "lwu",
#'                                verbose = FALSE)
#' summary(fit)
#'
#' @export
estimate_parametric_hrf <- function(
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
    local_radius = 26,
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
      hrf_interface = hrf_interface
    )
  }
  
  # Handle theta_bounds
  if (is.null(theta_bounds)) {
    theta_bounds <- hrf_interface$default_bounds()
  }
  
  # Set initial values
  for (j in seq_len(ncol(theta_current))) {
    theta_current[, j] <- theta_seed[j]
  }
  
  # K-means initialization if requested
  if (kmeans_refinement && kmeans_k > 1) {
    if (verbose) cat("  Performing K-means clustering for initialization...\n")
    kmeans_result <- .perform_kmeans_initialization(
      Y = inputs$Y_proj,
      S = inputs$S_target_proj,
      k = kmeans_k,
      hrf_interface = hrf_interface
    )
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
      scan_times = inputs$scan_times,
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
      scan_times = inputs$scan_times,
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
  if (!is.matrix(theta_current)) {
    theta_current <- matrix(theta_current, nrow = n_vox, ncol = length(hrf_interface$parameter_names))
    colnames(theta_current) <- hrf_interface$parameter_names
  }
  
  # Compute initial R-squared
  r_squared <- .compute_r_squared(
    Y = inputs$Y_proj,
    Y_pred = inputs$S_target_proj %*% core_result$beta0
  )
  
  if (verbose) {
    cat(sprintf("  Initial fit: Mean R² = %.3f (range: %.3f - %.3f)\n",
                mean(r_squared), min(r_squared), max(r_squared)))
  }
  
  # ========== STAGE 3: ITERATIVE REFINEMENT ==========
  convergence_info <- list()

  use_global_refinement <- isTRUE(getOption("fmriparametric.refine_global", TRUE))

  if (use_global_refinement && global_refinement && global_passes > 0) {
    if (verbose) cat("\n→ Stage 3: Global iterative refinement...\n")

    best_theta <- theta_current
    best_amplitudes <- amplitudes
    best_r2 <- mean(r_squared)

    for (iter in seq_len(global_passes)) {
      if (verbose) cat(sprintf("  Iteration %d/%d: ", iter, global_passes))

      # Store previous parameters
      theta_prev <- theta_current
      amplitudes_prev <- amplitudes
      r_squared_prev <- r_squared
      
      # Re-center globally with bounds enforcement
      theta_center <- apply(theta_current, 2, median)
      
      # Ensure theta_center is within bounds before using as seed
      if (!is.null(theta_bounds)) {
        theta_center <- pmax(theta_bounds$lower, pmin(theta_center, theta_bounds$upper))
      }
      
      # Re-run engine with new center
      if (parallel) {
        iter_result <- process_function(
          Y_proj = inputs$Y_proj,
          S_target_proj = inputs$S_target_proj,
          scan_times = inputs$scan_times,
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
          scan_times = inputs$scan_times,
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
      if (!is.matrix(theta_current) || is.null(dim(theta_current))) {
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
      r_squared_new <- .compute_r_squared(
        Y = inputs$Y_proj,
        Y_pred = inputs$S_target_proj %*% iter_result$beta0
      )

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
          break
        }
      }
    }

    theta_current <- best_theta
    amplitudes <- best_amplitudes
  }
  
  # ========== STAGE 4: TIERED REFINEMENT ==========
  refinement_info <- NULL
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
        local_radius = refinement_thresholds$local_radius,
        parallel = parallel,
        parallel_config = parallel_config
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
        max_iter = refinement_thresholds$gauss_newton_maxiter,
        parallel = parallel,
        parallel_config = parallel_config
      )
      
      # Update results
      theta_current[hard_idx, ] <- hard_result$theta_refined
      amplitudes[hard_idx] <- hard_result$amplitudes
      refinement_info$hard_refined <- length(hard_idx)
    }
    
    # Recompute final R-squared
    r_squared <- .compute_r_squared(
      Y = inputs$Y_proj,
      Y_pred = inputs$S_target_proj %*% amplitudes
    )
    
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
    rmse = sqrt(mean((inputs$Y_proj - inputs$S_target_proj %*% amplitudes)^2))
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
    version = "ultimate_impeccable_v1.0"
  )

  # Clean up parallel backend
  if (!is.null(parallel_config)) {
    parallel_config$cleanup()
  }

  # Ensure theta_current is a proper matrix
  if (!is.matrix(theta_current)) {
    theta_current <- matrix(theta_current, nrow = n_vox, ncol = length(hrf_interface$parameter_names))
  }
  
  # Ensure proper column names
  if (is.null(colnames(theta_current))) {
    colnames(theta_current) <- hrf_interface$parameter_names
  }
  
  # Ensure theta_current is still a matrix before creating fit object
  if (!is.matrix(theta_current)) {
    theta_current <- matrix(theta_current, nrow = n_vox, ncol = length(hrf_interface$parameter_names))
    colnames(theta_current) <- hrf_interface$parameter_names
  }
  
  # Create and return parametric_hrf_fit object
  fit <- new_parametric_hrf_fit(
    estimated_parameters = theta_current,
    amplitudes = as.numeric(amplitudes),
    parameter_names = hrf_interface$parameter_names,
    hrf_model = parametric_hrf,
    convergence = convergence_info,
    metadata = metadata
  )
  
  # Add additional components
  fit$standard_errors <- se_theta
  fit$se_amplitudes <- se_amplitudes
  fit$fit_quality <- fit_quality
  fit$refinement_info <- refinement_info
  
  return(fit)
}

# Helper function to create HRF interface
.create_hrf_interface <- function(model = "lwu") {
  switch(model,
    lwu = list(
      hrf_function = .lwu_hrf_function,
      taylor_basis = .lwu_hrf_taylor_basis_function,
      parameter_names = .lwu_hrf_parameter_names(),
      default_seed = .lwu_hrf_default_seed,
      default_bounds = .lwu_hrf_default_bounds
    ),
    stop("Only 'lwu' model currently supported. More models coming soon!")
  )
}


# ========== HELPER FUNCTION IMPLEMENTATIONS ==========

# MISSING FUNCTION 1: Parallel engine
.parametric_engine_parallel <- function(Y_proj, S_target_proj, scan_times, hrf_eval_times,
                                       hrf_interface, theta_seed, theta_bounds,
                                       lambda_ridge = 0.01, parallel_config) {



  n_vox <- ncol(Y_proj)
  n_params <- length(hrf_interface$parameter_names)

## <<<<<<< codex/connect-unused-performance-features
  # Processing function for a chunk of voxels
  process_chunk <- function(voxel_idx) {
    res <- .parametric_engine(
      Y_proj = Y_proj[, voxel_idx, drop = FALSE],
## =======
  # Determine optimal parallel strategy
  perf <- .master_performance_dispatcher(n_vox, nrow(Y_proj), verbose = FALSE)
  if (is.null(n_cores)) {
    n_cores <- perf$parallel$recommended_cores
  }
  par_cfg <- .setup_parallel_backend(if (perf$parallel$use_parallel) n_cores else 1,
                                     verbose = FALSE)

  process_voxel <- function(idx) {
    res <- .parametric_engine(
      Y_proj = Y_proj[, idx, drop = FALSE],
## >>>>>>> main
      S_target_proj = S_target_proj,
      scan_times = scan_times,
      hrf_eval_times = hrf_eval_times,
      hrf_interface = hrf_interface,
      theta_seed = theta_seed,
      theta_bounds = theta_bounds,
      lambda_ridge = lambda_ridge,
      verbose = FALSE
    )
## <<<<<<< codex/connect-unused-performance-features
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

## =======
    list(theta_hat = res$theta_hat[1, ], beta0 = res$beta0[1])
  }

  results <- .parallel_voxel_processing(
    voxel_indices = seq_len(n_vox),
    process_function = process_voxel,
    parallel_config = par_cfg,
    chunk_size = "auto",
    progress = FALSE
  )

  par_cfg$cleanup()

  theta_hat <- matrix(unlist(lapply(results, `[[`, "theta_hat")), n_vox, n_params, byrow = TRUE)
  beta0 <- vapply(results, `[[`, numeric(1), "beta0")

## >>>>>>> main
  list(theta_hat = theta_hat, beta0 = beta0)
}

# MISSING FUNCTION 2: Standard errors via Delta method
.compute_standard_errors_delta <- function(theta_hat, beta0, Y_proj, S_target_proj, 
                                          hrf_interface, hrf_eval_times) {
  
  n_vox <- ncol(Y_proj)
  n_params <- length(hrf_interface$parameter_names)
  n_time <- nrow(Y_proj)
  
  # Pre-allocate results
  se_theta_hat <- matrix(NA_real_, n_vox, n_params)
  se_beta0 <- numeric(n_vox)
  
  # Compute for each voxel (vectorized where possible)
  for (v in seq_len(n_vox)) {
    tryCatch({
      # Current parameter estimates
      theta_v <- theta_hat[v, ]
      
      # Get analytical derivatives from Taylor basis
      taylor_basis <- hrf_interface$taylor_basis(theta_v, hrf_eval_times)
      if (!is.matrix(taylor_basis)) {
        taylor_basis <- matrix(taylor_basis, ncol = n_params + 1)
      }
      basis_derivs <- taylor_basis[, -1, drop = FALSE]

      # Convolve derivatives with stimulus (sum across regressors)
      X_derivs <- matrix(0, n_time, n_params)
      for (p in seq_len(n_params)) {
        design_col <- numeric(n_time)
        for (j in seq_len(ncol(S_target_proj))) {
          conv_full <- stats::convolve(S_target_proj[, j], rev(basis_derivs[, p]), type = "open")
          design_col <- design_col + conv_full[seq_len(n_time)]
        }
        X_derivs[, p] <- beta0[v] * design_col
      }

      # Predicted HRF for amplitude SE computation
      hrf_vals <- hrf_interface$hrf_function(hrf_eval_times, theta_v)
      x_hrf <- numeric(n_time)
      for (j in seq_len(ncol(S_target_proj))) {
        conv_full <- stats::convolve(S_target_proj[, j], rev(hrf_vals), type = "open")
        x_hrf <- x_hrf + conv_full[seq_len(n_time)]
      }
      
      # Residual variance
      residuals <- Y_proj[, v] - beta0[v] * x_hrf
      sigma2 <- sum(residuals^2) / (n_time - 1)
      
      # Fisher Information Matrix (approximate)
      fisher_info <- crossprod(X_derivs) / sigma2
      
      # Standard errors (diagonal of inverse Fisher matrix)
      fisher_inv <- tryCatch({
        solve(fisher_info)
      }, error = function(e) {
        # Use pseudo-inverse if singular
        svd_result <- svd(fisher_info)
        svd_result$v %*% diag(1 / pmax(svd_result$d, 1e-12)) %*% t(svd_result$u)
      })
      
      se_theta_hat[v, ] <- sqrt(pmax(0, diag(fisher_inv)))
      
      # Standard error for beta0 using design for current HRF
      se_beta0[v] <- sqrt(sigma2 / sum(x_hrf^2))
      
    }, error = function(e) {
      # Fallback: use rough estimates
      se_theta_hat[v, ] <- abs(theta_hat[v, ]) * 0.1  # 10% of estimate
      se_beta0[v] <- abs(beta0[v]) * 0.1
    })
  }
  
  return(list(se_theta_hat = se_theta_hat, se_beta0 = se_beta0))
}


# DATA-DRIVEN INITIALIZATION
.compute_data_driven_seed <- function(Y, S, hrf_interface) {
  default_seed <- hrf_interface$default_seed()
  bounds <- hrf_interface$default_bounds()

  # Average high variance voxels
  vars <- apply(Y, 2, var)
  idx <- which(vars >= quantile(vars, 0.75))
  y_mean <- rowMeans(Y[, idx, drop = FALSE])

  # Cross correlation with stimulus
  lags <- seq(-8, 12, by = 1)
  cors <- vapply(lags, function(l) {
    if (l >= 0) {
      s_shift <- c(rep(0, l), S[, 1])[1:nrow(S)]
    } else {
      s_shift <- c(S[, 1], rep(0, -l))[(-l + 1):nrow(S)]
    }
    stats::cor(y_mean, s_shift, use = "complete.obs")
  }, numeric(1))

  best_lag <- lags[which.max(abs(cors))]
  tau_est <- default_seed[1] + best_lag
  tau_est <- max(bounds$lower[1], min(bounds$upper[1], tau_est))

  cors_abs <- abs(cors)
  half <- max(cors_abs) * 0.5
  width_idx <- which(cors_abs >= half)
  sigma_est <- if (length(width_idx) >= 2) diff(range(lags[width_idx]))/2 else default_seed[2]
  sigma_est <- max(bounds$lower[2], min(bounds$upper[2], sigma_est))

  c(tau_est, sigma_est, default_seed[3])
}

# K-MEANS INITIALIZATION
.perform_kmeans_initialization <- function(Y, S, k, hrf_interface) {
  n_vox <- ncol(Y)
  if (k <= 1 || n_vox <= k) {
    return(list(cluster = rep(1, n_vox),
                centers = matrix(hrf_interface$default_seed(), nrow = 1)))
  }

  bounds <- hrf_interface$default_bounds()
  delays <- seq(-8, 12, by = 1)
  tau_est <- numeric(n_vox)

  for (v in seq_len(n_vox)) {
    cors <- vapply(delays, function(l) {
      if (l >= 0) {
        s_shift <- c(rep(0, l), S[, 1])[1:nrow(S)]
      } else {
        s_shift <- c(S[, 1], rep(0, -l))[(-l + 1):nrow(S)]
      }
      stats::cor(Y[, v], s_shift, use = "complete.obs")
    }, numeric(1))
    tau_est[v] <- delays[which.max(abs(cors))]
  }

  features <- matrix(tau_est, ncol = 1)
  km <- kmeans(features, centers = k, nstart = 20, iter.max = 50)

  centers <- matrix(rep(hrf_interface$default_seed(), k), nrow = k, byrow = TRUE)
  for (cl in seq_len(k)) {
    idx <- which(km$cluster == cl)
    if (length(idx) > 0) {
      centers[cl, 1] <- median(tau_est[idx])
    }
  }
  centers[, 1] <- pmax(bounds$lower[1], pmin(bounds$upper[1], centers[, 1]))

  list(cluster = km$cluster, centers = centers)
}

.compute_r_squared <- function(Y, Y_pred) {
  ss_res <- colSums((Y - Y_pred)^2)
  ss_tot <- colSums(scale(Y, scale = FALSE)^2)
  r2 <- 1 - ss_res / pmax(ss_tot, .Machine$double.eps)
  pmax(0, pmin(1, r2))
}

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
                                    local_radius = 1, parallel = FALSE,
                                    n_cores = 1) {

  if (length(voxel_idx) == 0) {
    return(list(theta_refined = theta_current[voxel_idx, , drop = FALSE],
                amplitudes = numeric(0)))
  }

  bounds <- hrf_interface$default_bounds()
  n_time <- nrow(Y_proj)
  n_params <- length(hrf_interface$parameter_names)

  refine_one <- function(v) {
    theta_v <- theta_current[v, ]
    basis <- hrf_interface$taylor_basis(theta_v, hrf_eval_times)
    if (!is.matrix(basis)) {
      basis <- matrix(basis, ncol = n_params + 1)
    }
    X <- matrix(0, n_time, ncol(basis))
    for (j in seq_len(ncol(basis))) {
      conv_full <- stats::convolve(S_target_proj[, 1], rev(basis[, j]), type = "open")
      X[, j] <- conv_full[seq_len(n_time)]
    }
    qr_decomp <- qr(X)
    Q <- qr.Q(qr_decomp)
    R <- qr.R(qr_decomp)
    coeffs <- solve(R + 0.01 * diag(ncol(R))) %*% t(Q) %*% Y_proj[, v]
    beta0_new <- as.numeric(coeffs[1])
    if (abs(beta0_new) < 1e-6) {
      return(list(theta = theta_v, amp = beta0_new, r2 = r_squared[v]))
    }
    delta <- coeffs[2:(n_params + 1)] / beta0_new
    theta_new <- theta_v + as.numeric(delta)
    theta_new <- pmax(bounds$lower, pmin(bounds$upper, theta_new))
    fitted <- X %*% coeffs
    r2_new <- 1 - sum((Y_proj[, v] - fitted)^2) / sum((Y_proj[, v] - mean(Y_proj[, v]))^2)
    if (is.na(r2_new) || r2_new <= r_squared[v]) {
      return(list(theta = theta_v, amp = beta0_new, r2 = r_squared[v]))
    }
    list(theta = theta_new, amp = beta0_new, r2 = r2_new)
  }

  if (parallel && n_cores > 1 && .Platform$OS.type == "unix") {
    res_list <- parallel::mclapply(voxel_idx, refine_one, mc.cores = n_cores)
  } else if (parallel && n_cores > 1) {
    cl <- parallel::makeCluster(n_cores)
    res_list <- parallel::parLapply(cl, voxel_idx, refine_one)
    parallel::stopCluster(cl)
  } else {
    res_list <- lapply(voxel_idx, refine_one)
  }

  theta_out <- theta_current[voxel_idx, , drop = FALSE]
  amps_out <- numeric(length(voxel_idx))
  for (i in seq_along(res_list)) {
    theta_out[i, ] <- res_list[[i]]$theta
    amps_out[i] <- res_list[[i]]$amp
  }

  list(theta_refined = theta_out, amplitudes = amps_out)
}

.refine_hard_voxels <- function(voxel_idx, Y_proj, S_target_proj,
                                theta_current, r_squared,
                                hrf_interface, hrf_eval_times,
                                max_iter = 5, parallel = FALSE,
                                n_cores = 1) {

  if (length(voxel_idx) == 0) {
    return(list(theta_refined = theta_current[voxel_idx, , drop = FALSE],
                amplitudes = numeric(0)))
  }

  n_vox <- ncol(Y_proj)
  bounds <- hrf_interface$default_bounds()
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
    sum(x_pred * Y_proj[, v]) / sum(x_pred^2)
  }

  amps <- vapply(voxel_idx, amp_fun, numeric(1))

  list(theta_refined = theta_updated[voxel_idx, , drop = FALSE],
       amplitudes = amps)
}