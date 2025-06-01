#' Estimate parametric HRF parameters (Sprint 3 version)
#'
#' Main function for estimating voxel-wise parameters of parametric Hemodynamic 
#' Response Function (HRF) models from fMRI data. This version includes all
#' Sprint 3 enhancements: K-means re-centering, tiered refinement, and parallel
#' processing capabilities.
#'
#' @param fmri_data An fMRI dataset object or numeric matrix
#' @param event_model An event model object or numeric matrix
#' @param parametric_hrf Character string specifying the parametric HRF model
#' @param theta_seed Initial parameter values, "data_median", or NULL
#' @param theta_bounds List with elements `lower` and `upper`
#' @param confound_formula Formula specifying confound regressors
#' @param baseline_model Character string specifying baseline model
#' @param hrf_eval_times Numeric vector of HRF evaluation time points
#' @param hrf_span Duration in seconds for HRF evaluation
#' @param lambda_ridge Ridge regularization penalty
#' @param recenter_global_passes Number of global re-centering iterations
#' @param recenter_epsilon Convergence tolerance for re-centering
#' @param r2_threshold R-squared threshold for selecting good voxels
#' @param recenter_kmeans_passes Number of K-means re-centering passes
#' @param kmeans_k Number of clusters for K-means
#' @param r2_threshold_kmeans R² threshold for K-means clustering
#' @param refinement_opts List of refinement options
#' @param compute_se Logical whether to compute standard errors
#' @param n_cores Number of cores for parallel processing (NULL = auto)
#' @param progress Logical whether to show progress
#' @param mask Optional mask for spatial subsetting
#' @param verbose Logical for progress messages
#'
#' @return An object of class `parametric_hrf_fit` with Sprint 3 enhancements
#' @export
estimate_parametric_hrf_v3 <- function(
  fmri_data,
  event_model,
  parametric_hrf = "lwu",
  theta_seed = NULL,
  theta_bounds = NULL,
  confound_formula = NULL,
  baseline_model = "intercept",
  hrf_eval_times = NULL,
  hrf_span = 30,
  lambda_ridge = 0.01,
  recenter_global_passes = 3,
  recenter_epsilon = 0.01,
  r2_threshold = 0.1,
  recenter_kmeans_passes = 2,
  kmeans_k = 5,
  r2_threshold_kmeans = 0.7,
  refinement_opts = list(
    apply_refinement = TRUE,
    r2_threshold_hard = 0.3,
    r2_threshold_moderate = 0.7,
    se_threshold_hard = 0.5,
    se_threshold_moderate = 0.3,
    max_iter_gn = 5
  ),
  compute_se = TRUE,
  n_cores = NULL,
  progress = TRUE,
  mask = NULL,
  verbose = FALSE
) {
  # Input validation (same as v2)
  assertthat::assert_that(!missing(fmri_data), !missing(event_model))
  assertthat::assert_that(is.character(parametric_hrf), length(parametric_hrf) == 1)
  parametric_hrf <- tolower(parametric_hrf)
  
  if (identical(parametric_hrf, "lwu")) {
    hrf_interface <- list(
      hrf_function = .lwu_hrf_function,
      taylor_basis = .lwu_hrf_taylor_basis_function,
      parameter_names = .lwu_hrf_parameter_names(),
      default_seed = .lwu_hrf_default_seed(),
      default_bounds = .lwu_hrf_default_bounds()
    )
  } else {
    stop("Unsupported parametric_hrf: ", parametric_hrf, call. = FALSE)
  }
  
  # Handle theta_bounds
  if (is.null(theta_bounds)) {
    theta_bounds <- hrf_interface$default_bounds
  } else {
    assertthat::assert_that(
      is.list(theta_bounds),
      all(c("lower", "upper") %in% names(theta_bounds)),
      is.numeric(theta_bounds$lower),
      is.numeric(theta_bounds$upper),
      length(theta_bounds$lower) == length(hrf_interface$parameter_names),
      length(theta_bounds$upper) == length(hrf_interface$parameter_names),
      msg = "theta_bounds must be a list with numeric 'lower' and 'upper' vectors"
    )
  }
  
  # Validate K-means parameters
  assertthat::assert_that(is.numeric(recenter_kmeans_passes), 
                          length(recenter_kmeans_passes) == 1,
                          recenter_kmeans_passes >= 0)
  assertthat::assert_that(is.numeric(kmeans_k), 
                          length(kmeans_k) == 1,
                          kmeans_k >= 2)
  assertthat::assert_that(is.numeric(r2_threshold_kmeans), 
                          length(r2_threshold_kmeans) == 1,
                          r2_threshold_kmeans >= 0,
                          r2_threshold_kmeans <= 1)
  
  # Validate refinement options
  assertthat::assert_that(is.list(refinement_opts))
  refinement_defaults <- list(
    apply_refinement = TRUE,
    r2_threshold_hard = 0.3,
    r2_threshold_moderate = 0.7,
    se_threshold_hard = 0.5,
    se_threshold_moderate = 0.3,
    max_iter_gn = 5
  )
  refinement_opts <- utils::modifyList(refinement_defaults, refinement_opts)
  
  # Prepare inputs
  if (verbose) message("Preparing inputs...")
  prepared <- .prepare_parametric_inputs(
    fmri_data = fmri_data,
    event_model = event_model,
    confound_formula = confound_formula,
    baseline_model = baseline_model,
    hrf_eval_times = hrf_eval_times,
    hrf_span = hrf_span,
    mask = mask
  )
  
  # Handle theta_seed including "data_median" option (same as v2)
  if (is.null(theta_seed)) {
    theta_seed_final <- hrf_interface$default_seed
  } else if (identical(theta_seed, "data_median")) {
    if (verbose) message("Computing data-driven seed...")
    
    # Load iterative engine
    source(file.path(dirname(getwd()), "R", "parametric-engine-iterative.R"), local = TRUE)
    
    # Preliminary pass
    prelim_res <- .parametric_engine_iterative(
      Y_proj = prepared$Y_proj,
      S_target_proj = prepared$S_target_proj,
      scan_times = prepared$scan_times,
      hrf_eval_times = prepared$hrf_eval_times,
      hrf_interface = hrf_interface,
      theta_seed = hrf_interface$default_seed,
      theta_bounds = theta_bounds,
      lambda_ridge = lambda_ridge,
      recenter_global_passes = 1,
      compute_residuals = FALSE,
      compute_se = FALSE,
      verbose = FALSE
    )
    
    # Select good voxels
    r2_threshold_prelim <- quantile(prelim_res$r_squared, 0.75, na.rm = TRUE)
    idx_good <- which(prelim_res$r_squared >= r2_threshold_prelim)
    
    if (length(idx_good) >= 10) {
      theta_seed_final <- apply(
        prelim_res$theta_hat[idx_good, , drop = FALSE],
        2,
        median,
        na.rm = TRUE
      )
      if (verbose) {
        message("Data-driven seed from ", length(idx_good), " good voxels: ", 
                paste(round(theta_seed_final, 3), collapse = ", "))
      }
    } else {
      warning("Too few good voxels for data-driven seed; using default")
      theta_seed_final <- hrf_interface$default_seed
    }
  } else {
    assertthat::assert_that(
      is.numeric(theta_seed),
      length(theta_seed) == length(hrf_interface$parameter_names),
      msg = "theta_seed must be numeric and match number of parameters"
    )
    theta_seed_final <- theta_seed
  }
  
  # Validate seed is within bounds
  assertthat::assert_that(
    all(theta_seed_final >= theta_bounds$lower),
    all(theta_seed_final <= theta_bounds$upper),
    msg = "theta_seed must fall within theta_bounds"
  )
  
  # Set up parallel processing if requested
  parallel_info <- NULL
  parallel_config <- NULL
  if (!is.null(n_cores) && n_cores > 1) {
    if (verbose) message("Setting up parallel processing...")
    source(file.path(dirname(getwd()), "R", "parallel-processing.R"), local = TRUE)
    parallel_config <- .setup_parallel_backend(n_cores = n_cores, verbose = verbose)
    parallel_info <- list(
      n_cores = parallel_config$n_cores, 
      backend = parallel_config$backend
    )
  }
  
  # Run main estimation with K-means
  if (verbose) message("Running parametric engine with iterative refinement and K-means...")
  
  # Load required functions
  source(file.path(dirname(getwd()), "R", "parametric-engine-iterative.R"), local = TRUE)
  source(file.path(dirname(getwd()), "R", "kmeans-recentering.R"), local = TRUE)
  
  fit_res <- .parametric_engine_iterative(
    Y_proj = prepared$Y_proj,
    S_target_proj = prepared$S_target_proj,
    scan_times = prepared$scan_times,
    hrf_eval_times = prepared$hrf_eval_times,
    hrf_interface = hrf_interface,
    theta_seed = theta_seed_final,
    theta_bounds = theta_bounds,
    lambda_ridge = lambda_ridge,
    recenter_global_passes = recenter_global_passes,
    recenter_epsilon = recenter_epsilon,
    r2_threshold = r2_threshold,
    compute_residuals = TRUE,
    compute_se = compute_se,
    recenter_kmeans_passes = recenter_kmeans_passes,
    kmeans_k = kmeans_k,
    r2_threshold_kmeans = r2_threshold_kmeans,
    verbose = verbose
  )
  
  # Apply tiered refinement if requested
  refinement_info <- NULL
  if (refinement_opts$apply_refinement) {
    if (verbose) message("Applying tiered refinement...")
    
    # Load refinement functions
    source(file.path(dirname(getwd()), "R", "refinement-queue.R"), local = TRUE)
    source(file.path(dirname(getwd()), "R", "gauss-newton-refinement.R"), local = TRUE)
    source(file.path(dirname(getwd()), "R", "parametric-engine.R"), local = TRUE)
    
    # Classify voxels into refinement queues
    queue_result <- .classify_refinement_queue(
      r2_voxel = fit_res$r_squared,
      se_theta_hat_voxel = fit_res$se_theta_hat,
      refinement_opts = refinement_opts
    )
    
    # Print summary
    .print_refinement_summary(queue_result, verbose = verbose)
    
    if (queue_result$refinement_needed) {
      # Apply local re-centering for moderate voxels
      idx_moderate <- which(queue_result$queue_labels == "moderate_local_recenter")
      if (length(idx_moderate) > 0) {
        if (verbose) message("Applying local re-centering for ", length(idx_moderate), " moderate voxels...")
        
        if (!is.null(parallel_config) && parallel_config$n_cores > 1) {
          # Parallel processing
          recenter_results <- .parallel_local_recentering(
            moderate_indices = idx_moderate,
            fit_data = fit_res,
            prepared_data = prepared,
            hrf_interface = hrf_interface,
            parallel_config = parallel_config,
            theta_bounds = theta_bounds,
            lambda_ridge = lambda_ridge,
            verbose = verbose
          )
          
          # Apply updates
          n_improved <- 0
          for (res in recenter_results) {
            if (res$improved) {
              v <- res$voxel_idx
              fit_res$theta_hat[v, ] <- res$theta_new
              fit_res$beta0[v] <- res$beta_new
              fit_res$r_squared[v] <- res$r2_new
              if (!is.null(fit_res$residuals)) {
                hrf_vals <- hrf_interface$hrf_function(prepared$hrf_eval_times, res$theta_new)
                conv_full <- stats::convolve(prepared$S_target_proj[, 1], rev(hrf_vals), type = "open")
                x_pred <- conv_full[seq_len(nrow(prepared$Y_proj))]
                fit_res$residuals[, v] <- prepared$Y_proj[, v] - res$beta_new * x_pred
              }
              n_improved <- n_improved + 1
            }
          }
          
          if (verbose) message("Local re-centering complete: ", n_improved, " voxels improved")
        } else {
          # Sequential processing
          n_improved <- 0
          for (v in idx_moderate) {
            # Use voxel's current estimate as expansion point
            theta_local <- fit_res$theta_hat[v, ]
            
            # Single Taylor pass for this voxel
            local_res <- .parametric_engine(
              Y_proj = prepared$Y_proj[, v, drop = FALSE],
              S_target_proj = prepared$S_target_proj,
              scan_times = prepared$scan_times,
              hrf_eval_times = prepared$hrf_eval_times,
              hrf_interface = hrf_interface,
              theta_seed = theta_local,
              theta_bounds = theta_bounds,
              lambda_ridge = lambda_ridge,
              verbose = FALSE
            )
            
            # Update only if R² improves
            if (local_res$r_squared[1] > fit_res$r_squared[v]) {
              fit_res$theta_hat[v, ] <- local_res$theta_hat[1, ]
              fit_res$beta0[v] <- local_res$beta0[1]
              fit_res$r_squared[v] <- local_res$r_squared[1]
              if (!is.null(fit_res$residuals)) {
                fit_res$residuals[, v] <- prepared$Y_proj[, v] - local_res$beta0[1] * local_res$X_conv[, 1]
              }
              n_improved <- n_improved + 1
            }
          }
          
          if (verbose) message("Local re-centering complete: ", n_improved, " voxels improved")
        }
      }
      
      # Apply Gauss-Newton for hard voxels
      if (any(queue_result$queue_labels == "hard_GN")) {
        if (verbose) message("Applying Gauss-Newton refinement...")
        
        idx_hard <- which(queue_result$queue_labels == "hard_GN")
        
        if (!is.null(parallel_config) && parallel_config$n_cores > 1 && length(idx_hard) > 1) {
          # Parallel Gauss-Newton
          gn_results <- .parallel_gauss_newton(
            hard_indices = idx_hard,
            fit_data = fit_res,
            prepared_data = prepared,
            hrf_interface = hrf_interface,
            parallel_config = parallel_config,
            theta_bounds = theta_bounds,
            max_iter_gn = refinement_opts$max_iter_gn,
            tol_gn = 1e-4,
            lambda_ridge = lambda_ridge,
            verbose = verbose
          )
          
          # Apply updates
          n_converged <- 0
          n_improved <- 0
          for (res in gn_results) {
            v <- res$voxel_idx
            if (res$converged) n_converged <- n_converged + 1
            if (res$improved) {
              fit_res$theta_hat[v, ] <- res$theta_new
              fit_res$beta0[v] <- res$beta_new
              fit_res$r_squared[v] <- res$r2_new
              queue_result$queue_labels[v] <- "easy"  # Mark as resolved
              n_improved <- n_improved + 1
              
              if (!is.null(fit_res$residuals)) {
                hrf_vals <- hrf_interface$hrf_function(prepared$hrf_eval_times, res$theta_new)
                conv_full <- stats::convolve(prepared$S_target_proj[, 1], rev(hrf_vals), type = "open")
                x_pred <- conv_full[seq_len(nrow(prepared$Y_proj))]
                fit_res$residuals[, v] <- prepared$Y_proj[, v] - res$beta_new * x_pred
              }
            }
          }
          
          if (verbose) {
            message("Gauss-Newton complete: ", n_converged, " converged, ", 
                    n_improved, " improved out of ", length(idx_hard))
          }
          
          # Store GN result info
          gn_result <- list(
            n_refined = length(idx_hard),
            n_converged = n_converged,
            n_improved = n_improved
          )
        } else {
          # Sequential Gauss-Newton
          gn_result <- .gauss_newton_refinement(
            theta_hat_voxel = fit_res$theta_hat,
            r2_voxel = fit_res$r_squared,
            Y_proj = prepared$Y_proj,
            S_target_proj = prepared$S_target_proj,
            scan_times = prepared$scan_times,
            hrf_eval_times = prepared$hrf_eval_times,
            hrf_interface = hrf_interface,
            theta_bounds = theta_bounds,
            queue_labels = queue_result$queue_labels,
            max_iter_gn = refinement_opts$max_iter_gn,
            lambda_ridge = lambda_ridge,
            verbose = verbose
          )
          
          # Update results
          fit_res$theta_hat <- gn_result$theta_hat
          fit_res$r_squared <- gn_result$r2
          queue_result$queue_labels <- gn_result$queue_labels
          
          # Re-estimate amplitudes for refined voxels
          idx_refined <- which(gn_result$r2 != fit_res$r_squared)
          if (length(idx_refined) > 0) {
            for (v in idx_refined) {
              hrf_vals <- hrf_interface$hrf_function(prepared$hrf_eval_times, fit_res$theta_hat[v, ])
              conv_full <- stats::convolve(prepared$S_target_proj[, 1], rev(hrf_vals), type = "open")
              x_pred <- conv_full[seq_len(nrow(prepared$Y_proj))]
              fit_res$beta0[v] <- sum(x_pred * prepared$Y_proj[, v]) / sum(x_pred^2)
            }
          }
        }
      }
      
      # Recalculate standard errors for refined voxels if needed
      if (compute_se && !is.null(fit_res$se_theta_hat)) {
        idx_refined_any <- which(queue_result$queue_labels != "easy")
        if (length(idx_refined_any) > 0) {
          if (verbose) message("Recalculating standard errors for refined voxels...")
          
          # Recalculate SEs using final parameter estimates
          for (v in idx_refined_any) {
            taylor_basis <- hrf_interface$taylor_basis(fit_res$theta_hat[v, ], prepared$hrf_eval_times)
            if (!is.matrix(taylor_basis)) {
              taylor_basis <- matrix(taylor_basis, ncol = length(fit_res$theta_hat[v, ]) + 1)
            }
            
            # Convolve derivatives
            X_deriv <- matrix(0, nrow = nrow(prepared$Y_proj), ncol = ncol(taylor_basis) - 1)
            for (j in seq_len(ncol(X_deriv))) {
              conv_full <- stats::convolve(prepared$S_target_proj[, 1], rev(taylor_basis[, j + 1]), type = "open")
              X_deriv[, j] <- conv_full[seq_len(nrow(prepared$Y_proj))]
            }
            
            # Calculate SE via Delta method
            var_beta <- if (!is.null(fit_res$residuals)) {
              sum(fit_res$residuals[, v]^2) / (nrow(prepared$Y_proj) - ncol(X_deriv) - 1)
            } else {
              1
            }
            
            cov_theta <- tryCatch({
              var_beta * solve(crossprod(X_deriv) + lambda_ridge * diag(ncol(X_deriv)))
            }, error = function(e) {
              matrix(NA, ncol(X_deriv), ncol(X_deriv))
            })
            
            fit_res$se_theta_hat[v, ] <- sqrt(diag(cov_theta))
          }
        }
      }
    }
    
    # Create refinement info summary
    refinement_info <- list(
      applied = TRUE,
      queue_result = queue_result,
      n_moderate_refined = length(idx_moderate),
      n_hard_refined = if (exists("gn_result")) gn_result$n_refined else 0,
      n_converged = if (exists("gn_result")) gn_result$n_converged else 0,
      n_improved = if (exists("gn_result")) gn_result$n_improved else 0,
      final_queue_summary = table(queue_result$queue_labels)
    )
  }
  
  # Clean up parallel backend if used
  if (!is.null(parallel_config)) {
    parallel_config$cleanup()
  }
  
  # Load enhanced class constructor
  source(file.path(dirname(getwd()), "R", "parametric-hrf-fit-class-v2.R"), local = TRUE)
  
  # Construct enhanced output
  if (verbose) message("Constructing output...")
  new_parametric_hrf_fit(
    estimated_parameters = fit_res$theta_hat,
    amplitudes = fit_res$beta0,
    parameter_names = hrf_interface$parameter_names,
    hrf_model = parametric_hrf,
    r_squared = fit_res$r_squared,
    residuals = fit_res$residuals,
    parameter_ses = fit_res$se_theta_hat,
    convergence_info = fit_res$convergence_info,
    metadata = list(
      call = match.call(),
      n_voxels = ncol(prepared$Y_raw),
      n_timepoints = nrow(prepared$Y_raw),
      theta_seed = if(identical(theta_seed, "data_median")) "data_median" else theta_seed_final,
      theta_bounds = theta_bounds,
      recenter_global_passes = recenter_global_passes,
      recenter_kmeans_passes = recenter_kmeans_passes,
      kmeans_k = kmeans_k,
      kmeans_info = fit_res$kmeans_info,
      refinement_info = refinement_info,
      parallel_info = parallel_info,
      coeffs = fit_res$coeffs,
      theta_expansion = fit_res$theta_expansion
    )
  )
}