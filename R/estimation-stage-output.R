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
  
  # Get bounds from hrf_interface or config
  theta_bounds <- if (!is.null(hrf_interface$active_bounds)) {
    hrf_interface$active_bounds
  } else if (!is.null(config$theta_bounds)) {
    config$theta_bounds
  } else {
    hrf_interface$default_bounds()
  }
  
  # Create metadata
  metadata <- list(
    call = config$call,
    n_voxels = n_vox,
    n_timepoints = prepared_data$n_time,
    parametric_model = config$parametric_model,
    method_used = "parametric_taylor",
    theta_seed = config$theta_seed,
    theta_bounds = theta_bounds,
    S_target_proj = prepared_data$inputs$S_target_proj,
    hrf_eval_times = prepared_data$inputs$hrf_eval_times,
    design_info = list(
      n_time = prepared_data$n_time,
      n_vox = n_vox,
      n_cond = ncol(prepared_data$inputs$S_target_proj),
      basis_dim = ncol(final_results$theta_current),
      projected = !is.null(prepared_data$inputs$Z)
    ),
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
