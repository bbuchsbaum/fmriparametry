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
  guarded <- .apply_constant_voxel_guards(
    r_squared = r_squared,
    amplitudes = amplitudes,
    constant_idx = prepared_data$constant_voxel_idx
  )
  r_squared <- guarded$r_squared
  amplitudes <- guarded$amplitudes
  
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
  
  # Refinement loop
  for (iter in seq_len(config$global_passes)) {
    if (verbose) cat(sprintf("  Iteration %d/%d: ", iter, config$global_passes))
    
    # Store previous parameters
    theta_prev <- theta_current
    amplitudes_prev <- amplitudes
    r_squared_prev <- r_squared
    
    # Re-center globally at a stable interior anchor.
    default_seed <- hrf_interface$default_seed()
    theta_center <- .compute_global_refinement_center(
      theta_current = theta_current,
      default_seed = default_seed,
      theta_bounds = theta_bounds
    )
    
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
    
    theta_candidate <- iter_result$theta_hat
    
    # Ensure theta_candidate remains a matrix (fix for single voxel case)
    if (!is.matrix(theta_candidate) || length(dim(theta_candidate)) != 2) {
      theta_candidate <- matrix(theta_candidate, 
                                nrow = prepared_data$n_vox, 
                                ncol = length(hrf_interface$parameter_names))
    }
    
    # Enforce bounds on the refined parameters
    theta_candidate <- .clamp_theta_to_bounds(theta_candidate, theta_bounds)
    
    amplitudes_candidate <- as.numeric(iter_result$beta0)
    r_squared_candidate <- as.numeric(iter_result$r_squared)
    guarded <- .apply_constant_voxel_guards(
      r_squared = r_squared_candidate,
      amplitudes = amplitudes_candidate,
      constant_idx = prepared_data$constant_voxel_idx
    )
    r_squared_candidate <- guarded$r_squared
    amplitudes_candidate <- guarded$amplitudes
    
    # Keep per-voxel improvements only; never replace a clearly better fit.
    # If R^2 is effectively tied, prefer candidates that are closer to the
    # model's default interior seed and not pinned to bounds.
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
  
  list(
    theta_current = theta_current,
    amplitudes = amplitudes,
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
  constant_idx <- prepared_data$constant_voxel_idx
  guarded <- .apply_constant_voxel_guards(
    r_squared = r_squared,
    amplitudes = amplitudes,
    constant_idx = constant_idx
  )
  r_squared <- guarded$r_squared
  amplitudes <- guarded$amplitudes
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
      verbose = FALSE
    )
    
    improved_hard_idx <- hard_idx[hard_result$r2[hard_idx] > (r_squared[hard_idx] + 1e-5)]
    if (length(improved_hard_idx) > 0) {
      theta_current[improved_hard_idx, ] <- hard_result$theta_hat[improved_hard_idx, , drop = FALSE]
      r_squared[improved_hard_idx] <- hard_result$r2[improved_hard_idx]
      
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
    if (verbose) {
      cat(sprintf("  Final fit after refinement: Mean R^2 = %.3f\n", mean(r_squared)))
    }
  }
  
  # CRITICAL: Enforce bounds on final parameters after all refinement
  # This ensures no parameters violate user-specified bounds
  # Use active_bounds from interface if available, otherwise use config
  final_bounds <- .resolve_optional_theta_bounds(hrf_interface, config)
  theta_current <- .clamp_theta_to_bounds(theta_current, final_bounds)
  
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
  
  result
}

# ========== Stage-specific helper functions ==========
