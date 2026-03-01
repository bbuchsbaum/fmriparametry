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
      r_squared_raw = core_results$r_squared_raw %||% core_results$r_squared,
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
  r_squared_raw <- core_results$r_squared_raw %||% core_results$r_squared
  guarded <- .apply_constant_voxel_guards(
    r_squared = r_squared,
    amplitudes = amplitudes,
    constant_idx = prepared_data$constant_voxel_idx
  )
  r_squared <- guarded$r_squared
  amplitudes <- guarded$amplitudes
  if (length(prepared_data$constant_voxel_idx) > 0) {
    r_squared_raw[prepared_data$constant_voxel_idx] <- 0
  }
  
  convergence_info <- list(
    global_iterations = 0,
    converged = FALSE
  )
  
  # Use non-parallel processing for now
  process_function <- .parametric_engine
  convergence_config <- .create_convergence_config(
    param_tol = config$convergence_epsilon,
    r2_tol = 1e-5
  )
  theta_bounds <- .resolve_optional_theta_bounds(hrf_interface, config)
  
  iter <- 0
  default_seed <- hrf_interface$default_seed()

  # Refinement loop
  for (iter in seq_len(config$global_passes)) {
    if (verbose) cat(sprintf("  Iteration %d/%d: ", iter, config$global_passes))

    # Store previous parameters
    theta_prev <- theta_current
    amplitudes_prev <- amplitudes
    r_squared_prev <- r_squared

    # Initialize candidate arrays
    n_vox <- prepared_data$n_vox
    n_params <- length(hrf_interface$parameter_names)
    theta_candidate <- theta_current
    amplitudes_candidate <- amplitudes
    r_squared_candidate <- r_squared
    r_squared_candidate_raw <- r_squared_raw

    # 1. Always run global center pass (baseline, same as original algorithm)
    theta_center <- .compute_global_refinement_center(
      theta_current = theta_current,
      default_seed = default_seed,
      theta_bounds = theta_bounds
    )

    global_result <- process_function(
      Y_proj = prepared_data$inputs$Y_proj,
      S_target_proj = prepared_data$inputs$S_target_proj,
      hrf_eval_times = prepared_data$inputs$hrf_eval_times,
      hrf_interface = hrf_interface,
      theta_seed = theta_center,
      theta_bounds = theta_bounds,
      lambda_ridge = config$lambda_ridge,
      baseline_model = config$baseline_model
    )

    g_theta <- global_result$theta_hat
    if (!is.matrix(g_theta) || length(dim(g_theta)) != 2) {
      g_theta <- matrix(g_theta, nrow = n_vox, ncol = n_params)
    }
    g_theta <- .clamp_theta_to_bounds(g_theta, theta_bounds)
    theta_candidate <- g_theta
    amplitudes_candidate <- as.numeric(global_result$beta0)
    r_squared_candidate <- as.numeric(global_result$r_squared)
    r_squared_candidate_raw <- as.numeric(global_result$r_squared_raw %||% global_result$r_squared)

    # 2. Group voxels by coarse theta grid and run per-group recentering
    groups <- .compute_theta_grid_groups(
      theta_current = theta_current,
      theta_bounds = theta_bounds
    )
    n_groups <- length(groups)
    if (verbose) {
      cat(sprintf("%d groups ", n_groups))
    }

    # Only run grouped recentering if it adds value (more than 1 unique group
    # and groups aren't all singletons with the same centroid as the global center)
    if (n_groups > 1 && n_groups < n_vox) {
      for (g in groups) {
        idx <- g$idx
        Y_group <- prepared_data$inputs$Y_proj[, idx, drop = FALSE]

        group_result <- process_function(
          Y_proj = Y_group,
          S_target_proj = prepared_data$inputs$S_target_proj,
          hrf_eval_times = prepared_data$inputs$hrf_eval_times,
          hrf_interface = hrf_interface,
          theta_seed = g$centroid,
          theta_bounds = theta_bounds,
          lambda_ridge = config$lambda_ridge,
          baseline_model = config$baseline_model
        )

        gr_theta <- group_result$theta_hat
        if (!is.matrix(gr_theta) || length(dim(gr_theta)) != 2) {
          gr_theta <- matrix(gr_theta, nrow = length(idx), ncol = n_params)
        }
        gr_theta <- .clamp_theta_to_bounds(gr_theta, theta_bounds)
        gr_r2 <- as.numeric(group_result$r_squared)
        gr_r2_raw <- as.numeric(group_result$r_squared_raw %||% group_result$r_squared)
        gr_amp <- as.numeric(group_result$beta0)

        # For this group's voxels, keep whichever is better: global or grouped
        better <- is.finite(gr_r2) & (gr_r2 > r_squared_candidate[idx])
        if (any(better)) {
          update_idx <- idx[better]
          theta_candidate[update_idx, ] <- gr_theta[better, , drop = FALSE]
          amplitudes_candidate[update_idx] <- gr_amp[better]
          r_squared_candidate[update_idx] <- gr_r2[better]
          r_squared_candidate_raw[update_idx] <- gr_r2_raw[better]
        }
      }
    }

    guarded <- .apply_constant_voxel_guards(
      r_squared = r_squared_candidate,
      amplitudes = amplitudes_candidate,
      constant_idx = prepared_data$constant_voxel_idx
    )
    r_squared_candidate <- guarded$r_squared
    amplitudes_candidate <- guarded$amplitudes
    if (length(prepared_data$constant_voxel_idx) > 0) {
      r_squared_candidate_raw[prepared_data$constant_voxel_idx] <- 0
    }

    # Keep per-voxel improvements only
    improved_idx <- .select_global_refinement_updates(
      theta_prev = theta_prev,
      theta_candidate = theta_candidate,
      r_squared_prev = r_squared_prev,
      r_squared_candidate = r_squared_candidate,
      default_seed = default_seed,
      theta_bounds = theta_bounds,
      r2_tol = convergence_config$r2_tol
    )
    if (length(improved_idx) > 0) {
      theta_current[improved_idx, ] <- theta_candidate[improved_idx, , drop = FALSE]
      amplitudes[improved_idx] <- amplitudes_candidate[improved_idx]
      r_squared[improved_idx] <- r_squared_candidate[improved_idx]
      r_squared_raw[improved_idx] <- r_squared_candidate_raw[improved_idx]
    } else {
      theta_current <- theta_prev
      amplitudes <- amplitudes_prev
      r_squared <- r_squared_prev
    }

    conv_check <- .check_convergence(
      current = theta_current,
      previous = theta_prev,
      config = convergence_config
    )

    mean_r2_change <- mean(r_squared) - mean(r_squared_prev)
    if (verbose) {
      max_change <- max(abs(theta_current - theta_prev))
      cat(sprintf(
        "Max Deltatheta = %.4f, Mean R^2 = %.3f (Delta = %+.4f, improved=%d)\n",
        max_change, mean(r_squared), mean_r2_change, length(improved_idx)
      ))
    }

    if (length(improved_idx) == 0) {
      convergence_info$converged <- TRUE
      convergence_info$convergence_reason <- "no_improvement"
      break
    }

    if (conv_check$converged) {
      if (verbose) cat("  [OK] Converged!\n")
      convergence_info$converged <- TRUE
      convergence_info$convergence_reason <- conv_check$reason
      break
    }
  }
  
  convergence_info$global_iterations <- iter
  convergence_info$iterations <- iter  # Add for backward compatibility with print methods

  # Recompute linear terms to keep residuals/intercepts consistent with final theta.
  linear_terms <- .recompute_linear_terms(
    theta_hat = theta_current,
    Y_proj = prepared_data$inputs$Y_proj,
    S_target_proj = prepared_data$inputs$S_target_proj,
    hrf_interface = hrf_interface,
    hrf_eval_times = prepared_data$inputs$hrf_eval_times,
    baseline_model = config$baseline_model,
    amplitudes = amplitudes
  )
  
  list(
    theta_current = theta_current,
    amplitudes = amplitudes,
    r_squared = r_squared,
    r_squared_raw = r_squared_raw,
    residuals = linear_terms$residuals,
    intercepts = linear_terms$intercepts,
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
      r_squared_raw = refined_results$r_squared_raw %||% refined_results$r_squared,
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
  r_squared_raw <- refined_results$r_squared_raw %||% refined_results$r_squared
  constant_idx <- prepared_data$constant_voxel_idx
  guarded <- .apply_constant_voxel_guards(
    r_squared = r_squared,
    amplitudes = amplitudes,
    constant_idx = constant_idx
  )
  r_squared <- guarded$r_squared
  amplitudes <- guarded$amplitudes
  if (length(constant_idx) > 0) {
    r_squared_raw[constant_idx] <- 0
  }

  n_vox <- prepared_data$n_vox
  
  # Classify voxels using R^2-only mode (Stage 5 will compute real SEs)
  voxel_classes <- .classify_voxels_for_refinement(
    r_squared = r_squared,
    se_theta = NULL,  # NULL triggers R^2-only classification
    thresholds = config$refinement_thresholds
  )
  if (length(constant_idx) > 0) {
    voxel_classes[constant_idx] <- "easy"
  }
  
  if (verbose) {
    cat(sprintf("  Voxel classification:\n"))
    cat(sprintf("    Easy (high R^2): %d voxels\n", sum(voxel_classes == "easy")))
    cat(sprintf("    Moderate: %d voxels\n", sum(voxel_classes == "moderate")))
    cat(sprintf("    Hard (low R^2): %d voxels\n", sum(voxel_classes == "hard")))
    if (length(constant_idx) > 0) {
      cat(sprintf("    Constant/degenerate: %d voxels\n", length(constant_idx)))
    }
  }
  
  refinement_info <- list(
    classification = voxel_classes,
    n_easy = sum(voxel_classes == "easy"),
    n_moderate = sum(voxel_classes == "moderate"),
    n_hard = sum(voxel_classes == "hard"),
    n_constant = length(constant_idx)
  )
  did_update <- FALSE
  
  # Apply refinements based on classification

  # Easy voxels: one pass of local re-centering to tighten estimates
  easy_idx <- which(voxel_classes == "easy")
  # Exclude constant voxels from easy refinement
  easy_idx <- setdiff(easy_idx, constant_idx)
  if (length(easy_idx) > 0 && config$tiered_refinement %in% c("moderate", "aggressive")) {
    if (verbose) cat("  Refining easy voxels with local re-centering...\n")

    easy_result <- .refine_moderate_voxels_grouped(
      voxel_idx = easy_idx,
      Y_proj = prepared_data$inputs$Y_proj,
      S_target_proj = prepared_data$inputs$S_target_proj,
      theta_current = theta_current,
      amplitudes_current = amplitudes,
      r_squared = r_squared,
      hrf_interface = hrf_interface,
      hrf_eval_times = prepared_data$inputs$hrf_eval_times,
      theta_bounds = if (!is.null(hrf_interface$active_bounds)) hrf_interface$active_bounds else config$theta_bounds,
      lambda_ridge = config$lambda_ridge,
      baseline_model = config$baseline_model
    )

    theta_current[easy_idx, ] <- easy_result$theta_refined
    amplitudes[easy_idx] <- easy_result$amplitudes
    r_squared[easy_idx] <- easy_result$r_squared
    r_squared_raw[easy_idx] <- easy_result$r_squared
    refinement_info$easy_refined <- length(easy_idx)
    refinement_info$easy_improved <- easy_result$n_improved
    if (easy_result$n_improved > 0) {
      did_update <- TRUE
    }
    if (verbose) {
      cat(sprintf("    Easy voxels improved: %d/%d\n",
                  easy_result$n_improved, length(easy_idx)))
    }
  }

  # Moderate voxels: Local re-centering
  moderate_idx <- which(voxel_classes == "moderate")
  if (length(moderate_idx) > 0 && config$tiered_refinement %in% c("moderate", "aggressive")) {
    if (verbose) cat("  Refining moderate voxels with local re-centering...\n")
    
    moderate_result <- .refine_moderate_voxels_grouped(
      voxel_idx = moderate_idx,
      Y_proj = prepared_data$inputs$Y_proj,
      S_target_proj = prepared_data$inputs$S_target_proj,
      theta_current = theta_current,
      amplitudes_current = amplitudes,
      r_squared = r_squared,
      hrf_interface = hrf_interface,
      hrf_eval_times = prepared_data$inputs$hrf_eval_times,
      theta_bounds = if (!is.null(hrf_interface$active_bounds)) hrf_interface$active_bounds else config$theta_bounds,
      lambda_ridge = config$lambda_ridge,
      baseline_model = config$baseline_model
    )
    
    # Update results
    theta_current[moderate_idx, ] <- moderate_result$theta_refined
    amplitudes[moderate_idx] <- moderate_result$amplitudes
    r_squared[moderate_idx] <- moderate_result$r_squared
    r_squared_raw[moderate_idx] <- moderate_result$r_squared
    refinement_info$moderate_refined <- length(moderate_idx)
    refinement_info$moderate_improved <- moderate_result$n_improved
    if (moderate_result$n_improved > 0) {
      did_update <- TRUE
    }
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
    
    gn_parallel_min_voxels <- config$refinement_thresholds$gauss_newton_parallel_min_voxels
    if (is.null(gn_parallel_min_voxels)) {
      gn_parallel_min_voxels <- 200L
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
      parallel = isTRUE(config$parallel),
      n_cores = config$n_cores,
      parallel_min_voxels = gn_parallel_min_voxels,
      baseline_model = config$baseline_model,
      verbose = FALSE
    )
    
    improved_hard_idx <- hard_idx[hard_result$r2[hard_idx] > (r_squared[hard_idx] + 1e-5)]
    if (length(improved_hard_idx) > 0) {
      theta_current[improved_hard_idx, ] <- hard_result$theta_hat[improved_hard_idx, , drop = FALSE]
      r_squared[improved_hard_idx] <- hard_result$r2[improved_hard_idx]
      r_squared_raw[improved_hard_idx] <- hard_result$r2[improved_hard_idx]
      
      # Recompute amplitudes only for improved hard voxels.
      hard_amplitudes <- .compute_amplitudes_for_voxels(
        voxel_idx = improved_hard_idx,
        theta_hat = theta_current,
        Y_proj = prepared_data$inputs$Y_proj,
        S_target_proj = prepared_data$inputs$S_target_proj,
        hrf_interface = hrf_interface,
        hrf_eval_times = prepared_data$inputs$hrf_eval_times
      )
      amplitudes[improved_hard_idx] <- hard_amplitudes
      did_update <- TRUE
    }
    
    refinement_info$hard_refined <- length(hard_idx)
    refinement_info$n_converged <- hard_result$n_converged
    refinement_info$hard_improved <- length(improved_hard_idx)
  }
  
  if (did_update) {
    guarded <- .apply_constant_voxel_guards(
      r_squared = r_squared,
      amplitudes = amplitudes,
      constant_idx = constant_idx
    )
    r_squared <- guarded$r_squared
    amplitudes <- guarded$amplitudes
    if (length(constant_idx) > 0) {
      r_squared_raw[constant_idx] <- 0
    }
    if (verbose) {
      cat(sprintf("  Final fit after refinement: Mean R^2 = %.3f\n", mean(r_squared)))
    }
  }
  
  # CRITICAL: Enforce bounds on final parameters after all refinement
  # This ensures no parameters violate user-specified bounds
  # Use active_bounds from interface if available, otherwise use config
  final_bounds <- .resolve_optional_theta_bounds(hrf_interface, config)
  theta_current <- .clamp_theta_to_bounds(theta_current, final_bounds)

  # Recompute amplitudes/intercepts/residuals/R^2 for final theta estimates.
  linear_terms <- .recompute_linear_terms(
    theta_hat = theta_current,
    Y_proj = prepared_data$inputs$Y_proj,
    S_target_proj = prepared_data$inputs$S_target_proj,
    hrf_interface = hrf_interface,
    hrf_eval_times = prepared_data$inputs$hrf_eval_times,
    baseline_model = config$baseline_model,
    amplitudes = amplitudes
  )
  r_squared_raw <- linear_terms$r_squared_raw %||% r_squared
  
  # Return fitted values if we computed them
  result <- list(
    theta_current = theta_current,
    amplitudes = amplitudes,
    r_squared = r_squared,
    r_squared_raw = r_squared_raw,
    residuals = linear_terms$residuals,
    intercepts = linear_terms$intercepts,
    refinement_info = refinement_info,
    se_theta = NULL,  # Let Stage 5 compute real SEs
    se_amplitudes = refined_results$se_amplitudes,  # Pass through if exists
    convergence_info = refined_results$convergence_info
  )
  
  result
}

# ========== Stage-specific helper functions ==========
