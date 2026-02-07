#' S3 methods for parametric_hrf_fit objects
#'
#' These methods provide printing, summarizing and plotting functionality
#' for objects returned by `estimate_parametric_hrf()` and related
#' constructors.  They correspond to the final Sprint 3 versions which
#' include refinement diagnostics and advanced visualisations.
#'
#' @name parametric_hrf_fit-methods
NULL


#' Print a parametric_hrf_fit object
#'
#' Displays a comprehensive summary including refinement information,
#' convergence diagnostics and parallel processing details.
#'
#' @param x A `parametric_hrf_fit` object
#' @param ... Additional arguments (ignored)
#' @return The input object invisibly
#' @export
#' @examples
#' \dontrun{
#' # Assuming 'fit' is a parametric_hrf_fit object from estimate_parametric_hrf()
#' print(fit)
#' 
#' # Shows summary information about the fit including:
#' # - Number of voxels and parameters
#' # - Mean R-squared
#' # - Parameter ranges
#' }
print.parametric_hrf_fit <- function(x, ...) {
  cat("Parametric HRF Fit\n")
  cat("==================\n")
  cat("Model:", x$parametric_model, "\n")
  cat("Voxels:", nrow(get_parameters(x)), "\n")

  # Basic parameter summary
  cat("\nParameter Summary:\n")
  params <- get_parameters(x)
  param_means <- colMeans(params, na.rm = TRUE)
  param_sds <- apply(params, 2, sd, na.rm = TRUE)
  param_summary <- data.frame(
    Parameter = x$parameter_names,
    Mean = round(param_means, 3),
    SD = round(param_sds, 3)
  )
  print(param_summary, row.names = FALSE)

  if (!is.null(x$r_squared)) {
    cat("\nModel Fit (R^2):\n")
    r2_summary <- summary(x$r_squared)
    cat("  Min:", round(r2_summary[1], 3), "\n")
    cat("  Median:", round(r2_summary[3], 3), "\n")
    cat("  Mean:", round(mean(x$r_squared, na.rm = TRUE), 3), "\n")
    cat("  Max:", round(r2_summary[6], 3), "\n")
  }
  
  if (!is.null(x$convergence_info)) {
    cat("\nConvergence:\n")
    cat("  Iterations:", x$convergence_info$iterations, "\n")
    if (!is.null(x$convergence_info$converged)) {
      cat("  Converged:", x$convergence_info$converged, "\n")
    }
    if (!is.null(x$convergence_info$final_change)) {
      cat("  Final parameter change:", 
          format(x$convergence_info$final_change, scientific = TRUE, digits = 3), 
          "\n")
    }
  }
  
  # Sprint 3: Refinement information
  if (!is.null(x$refinement_info) && !is.null(x$refinement_info$applied) && x$refinement_info$applied) {
    cat("\nRefinement Applied:\n")
    if (!is.null(x$refinement_info$voxels_refined)) {
      cat("  Voxels refined:", x$refinement_info$voxels_refined, "\n")
      
      if (!is.null(x$refinement_info$r2_improvement)) {
        cat("  Mean R^2 improvement:", 
            round(x$refinement_info$r2_improvement, 4), "\n")
      }
    }
    if (!is.null(x$refinement_info$strategy)) {
      cat("  Strategy:", x$refinement_info$strategy, "\n")
    }
    if (!is.null(x$refinement_info$duration)) {
      cat("  Refinement time:", 
          round(x$refinement_info$duration, 2), "seconds\n")
    }
  }
  
  invisible(x)
}

#' Summarize a parametric_hrf_fit object
#'
#' Provides detailed statistical summaries of the fitted HRF parameters and
#' model quality metrics. The Sprint 3 implementation includes refinement
#' statistics and stage-wise timing information.
#'
#' @param object A `parametric_hrf_fit` object
#' @param ... Additional arguments (ignored)
#' @return An object of class `summary_parametric_hrf_fit` containing:
#'   - `parameter_summary`: Summary statistics for each HRF parameter
#'   - `amplitude_summary`: Summary of response amplitudes
#'   - `r_squared_summary`: Distribution of R-squared values
#'   - `convergence_summary`: Convergence diagnostics
#'   - `refinement_summary`: Refinement statistics (if applied)
#'   - `performance_summary`: Computational performance metrics
#' @export
#' @examples
#' \dontrun{
#' s <- summary(fit)
#' print(s)
#' }
summary.parametric_hrf_fit <- function(object, ...) {
  # Get parameters using helper function
  params <- get_parameters(object)
  
  # Parameter summary
  parameter_summary <- data.frame(
    Parameter = object$parameter_names,
    Mean = colMeans(params, na.rm = TRUE),
    SD = apply(params, 2, sd, na.rm = TRUE),
    Min = apply(params, 2, min, na.rm = TRUE),
    Max = apply(params, 2, max, na.rm = TRUE)
  )
  
  # Amplitude summary
  amplitude_summary <- NULL
  if (!is.null(object$amplitudes)) {
    amplitude_summary <- summary(object$amplitudes)
  }
  
  # R-squared summary
  r_squared_summary <- NULL
  if (!is.null(object$r_squared)) {
    r_squared_summary <- summary(object$r_squared)
  }
  
  # Convergence summary
  convergence_summary <- object$convergence_info
  
  # Sprint 3: Refinement summary
  refinement_summary <- NULL
  if (!is.null(object$refinement_info) && !is.null(object$refinement_info$applied) && object$refinement_info$applied) {
    refinement_summary <- object$refinement_info
    
    # Add classification breakdown if available
    if (!is.null(object$metadata$voxel_classification)) {
      refinement_summary$classification <- table(object$metadata$voxel_classification)
    }
  }
  
  # Performance summary
  performance_summary <- NULL
  if (!is.null(object$metadata$computational_details)) {
    comp_details <- object$metadata$computational_details
    performance_summary <- list(
      total_time = comp_details$total_duration,
      parallel = comp_details$parallel_used,
      n_cores = comp_details$n_cores
    )
    
    # Add stage timings if available
    if (!is.null(comp_details$stage_timings)) {
      performance_summary$stage_timings <- comp_details$stage_timings
    }
  }
  
  result <- list(
    parameter_summary = parameter_summary,
    amplitude_summary = amplitude_summary,
    r_squared_summary = r_squared_summary,
    convergence_summary = convergence_summary,
    refinement_summary = refinement_summary,
    performance_summary = performance_summary,
    n_voxels = nrow(params),
    model = object$parametric_model
  )
  
  class(result) <- "summary_parametric_hrf_fit"
  result
}

#' @export
print.summary_parametric_hrf_fit <- function(x, ...) {
  cat("Summary of Parametric HRF Fit\n")
  cat("=============================\n")
  cat("Model:", x$model, "\n")
  cat("Voxels:", x$n_voxels, "\n\n")
  
  cat("Parameter Estimates:\n")
  print(x$parameter_summary, row.names = FALSE, digits = 3)
  
  if (!is.null(x$amplitude_summary)) {
    cat("\nAmplitude Summary:\n")
    print(x$amplitude_summary)
  }
  
  if (!is.null(x$r_squared_summary)) {
    cat("\nModel Fit (R^2) Summary:\n")
    print(x$r_squared_summary)
  }
  
  if (!is.null(x$convergence_summary)) {
    cat("\nConvergence Information:\n")
    cat("  Iterations:", x$convergence_summary$iterations, "\n")
    if (!is.null(x$convergence_summary$converged)) {
      cat("  Converged:", x$convergence_summary$converged, "\n")
    }
  }
  
  if (!is.null(x$refinement_summary) && !is.null(x$refinement_summary$applied) && x$refinement_summary$applied) {
    cat("\nRefinement Information:\n")
    cat("  Strategy:", x$refinement_summary$strategy, "\n")
    cat("  Voxels refined:", x$refinement_summary$voxels_refined, "\n")
    
    if (!is.null(x$refinement_summary$classification)) {
      cat("\nVoxel Classification:\n")
      print(x$refinement_summary$classification)
    }
  }
  
  if (!is.null(x$performance_summary)) {
    cat("\nPerformance:\n")
    cat("  Total time:", round(x$performance_summary$total_time, 2), "seconds\n")
    if (!is.null(x$performance_summary$parallel)) {
      cat("  Parallel:", x$performance_summary$parallel, "\n")
      if (x$performance_summary$parallel && !is.null(x$performance_summary$n_cores)) {
        cat("  Cores used:", x$performance_summary$n_cores, "\n")
      }
    }
  }
  
  invisible(x)
}

#' Extract coefficients from a parametric_hrf_fit object
#'
#' Returns HRF parameters, amplitudes, or standard errors depending on
#' the requested type. Supports both v1 and v2 object structures.
#'
#' @param object A `parametric_hrf_fit` object
#' @param type Character string specifying what to extract:
#'   - `"parameters"`: HRF parameters (default)
#'   - `"amplitudes"`: Response amplitudes (supports alias "amplitude")
#'   - `"se"`: Parameter standard errors
#' @param ... Additional arguments (ignored)
#' @return A matrix or vector of the requested coefficients
#' @export
#' @examples
#' \dontrun{
#' # Get HRF parameters (tau, sigma, rho for LWU model)
#' params <- coef(fit)
#' head(params)
#' 
#' # Get amplitudes
#' amps <- coef(fit, type = "amplitudes")
#' 
#' # Get standard errors
#' ses <- coef(fit, type = "se")
#' }
coef.parametric_hrf_fit <- function(object,
                                   type = c("parameters", "amplitudes", "amplitude", "se"),
                                   ...) {
  type <- match.arg(type)
  
  # Handle amplitude/amplitudes alias
  if (type == "amplitude") type <- "amplitudes"
  
  result <- switch(type,
    parameters = object$model_specific$parameters,
    amplitudes = object$amplitudes,
    se = {
      if (is.null(object$model_specific$standard_errors)) {
        warning("Standard errors not available")
        NULL
      } else {
        object$model_specific$standard_errors
      }
    })
  
  # Add row names for parameters
  if (type == "parameters" && !is.null(result)) {
    if (is.null(rownames(result))) {
      rownames(result) <- paste0("Voxel_", seq_len(nrow(result)))
    }
  }
  
  result
}

#' Compute residuals from a parametric_hrf_fit object
#'
#' Returns residuals if available, otherwise computes them from fitted values.
#' Requires the original data if residuals are not stored in the object.
#'
#' @param object A `parametric_hrf_fit` object
#' @param Y_proj Optional projected Y data. Required if residuals not stored.
#' @param ... Additional arguments (ignored)
#' @return Matrix of residuals (timepoints x voxels)
#' @export
residuals.parametric_hrf_fit <- function(object, Y_proj = NULL, ...) {
  # Check if residuals are already computed
  if (!is.null(object$residuals)) {
    return(object$residuals)
  }
  
  # Need Y_proj to compute residuals
  if (is.null(Y_proj)) {
    stop("Y_proj must be provided when residuals are not stored in the fit object")
  }
  
  # Compute fitted values
  fitted_vals <- fitted(object, Y_proj)
  
  # Return residuals
  Y_proj - fitted_vals
}

#' Plot a parametric_hrf_fit object
#'
#' Create various visualizations of the HRF fit results including individual
#' HRF curves, parameter distributions, and model quality metrics.
#'
#' @param x A `parametric_hrf_fit` object
#' @param type Character string specifying plot type:
#'   - `"hrf"`: Plot HRF curves (default)
#'   - `"parameters"`: Parameter distribution plots
#'   - `"quality"`: Model quality metrics (R-squared distribution)
#'   - `"refinement"`: Refinement diagnostics (if applicable)
#' @param voxel_indices Integer vector of voxel indices to plot (for type="hrf")
#' @param n_voxels Number of random voxels to plot if voxel_indices not specified
#' @param add_mean Logical, whether to add mean HRF curve
#' @param ... Additional graphical parameters
#' @return The input object invisibly
#' @export
plot.parametric_hrf_fit <- function(x,
                                    type = c("hrf", "parameters", "quality", "refinement"),
                                    voxel_indices = NULL,
                                    n_voxels = 10,
                                    add_mean = TRUE,
                                    ...) {
  type <- match.arg(type)
  
  switch(type,
         hrf = .plot_hrf(x, voxel_indices, n_voxels, add_mean, ...),
         parameters = .plot_parameters(x, ...),
         quality = .plot_quality(x, ...),
         refinement = .plot_refinement(x, ...)
  )
  
  invisible(x)
}

# Helper function to plot HRF curves
.plot_hrf <- function(x, voxel_indices = NULL, n_voxels = 10, add_mean = TRUE, ...) {
  # Use pre-computed HRF shapes if available
  if (!is.null(x$hrf_shape)) {
    .plot_hrf_shape(x$hrf_shape, x$metadata$hrf_eval_times, 
                    voxel_indices, n_voxels, add_mean, ...)
  } else {
    # Compute HRF curves on the fly if not pre-computed
    .plot_hrf_compute(x, voxel_indices, n_voxels, add_mean, ...)
  }
}

# Plot pre-computed HRF shapes
.plot_hrf_shape <- function(hrf_shape, times, voxel_indices = NULL, 
                           n_voxels = 10, add_mean = TRUE, ...) {
  # Check if hrf_shape is valid
  if (is.null(hrf_shape) || !is.matrix(hrf_shape) || ncol(hrf_shape) == 0) {
    warning("No HRF shapes available to plot")
    return(invisible(NULL))
  }
  
  # Check if times is valid
  if (is.null(times) || length(times) == 0) {
    warning("No time points available to plot")
    return(invisible(NULL))
  }
  
  n_total <- ncol(hrf_shape)
  
  if (is.null(voxel_indices)) {
    if (n_voxels >= n_total) {
      voxel_indices <- seq_len(n_total)
    } else {
      voxel_indices <- sample(n_total, n_voxels)
    }
  }
  
  # Setup plot
  y_range <- range(hrf_shape[, voxel_indices], na.rm = TRUE)
  plot(times, hrf_shape[, voxel_indices[1]], type = "l",
       ylim = y_range, col = "gray70",
       xlab = "Time (s)", ylab = "HRF Amplitude",
       main = "Hemodynamic Response Functions", ...)
  
  # Add other voxels
  if (length(voxel_indices) > 1) {
    for (i in 2:length(voxel_indices)) {
      lines(times, hrf_shape[, voxel_indices[i]], col = "gray70")
    }
  }
  
  # Add mean curve
  if (add_mean) {
    mean_hrf <- rowMeans(hrf_shape, na.rm = TRUE)
    lines(times, mean_hrf, col = "red", lwd = 2)
    legend("topright", c("Individual", "Mean"), 
           col = c("gray70", "red"), lwd = c(1, 2))
  }
}

# Compute HRF curves on the fly
.plot_hrf_compute <- function(x, voxel_indices = NULL, n_voxels = 10, 
                            add_mean = TRUE, ...) {
  params <- get_parameters(x)
  n_total <- nrow(params)
  
  if (is.null(voxel_indices)) {
    if (n_voxels >= n_total) {
      voxel_indices <- seq_len(n_total)
    } else {
      voxel_indices <- sample(n_total, n_voxels)
    }
  }
  
  # Get HRF interface
  hrf_interface <- .create_hrf_interface(x$parametric_model)
  times <- x$metadata$hrf_eval_times
  if (is.null(times)) {
    times <- seq(0, 30, length.out = 61)
  }
  
  # Compute HRF for selected voxels
  hrf_curves <- matrix(NA, length(times), length(voxel_indices))
  for (i in seq_along(voxel_indices)) {
    v <- voxel_indices[i]
    theta_v <- params[v, ]
    hrf_curves[, i] <- hrf_interface$hrf_function(times, theta_v)
  }
  
  # Plot
  y_range <- range(hrf_curves, na.rm = TRUE)
  plot(times, hrf_curves[, 1], type = "l", ylim = y_range, col = "gray70",
       xlab = "Time (s)", ylab = "HRF Amplitude",
       main = "Hemodynamic Response Functions", ...)
  
  if (ncol(hrf_curves) > 1) {
    for (i in 2:ncol(hrf_curves)) {
      lines(times, hrf_curves[, i], col = "gray70")
    }
  }
  
  if (add_mean) {
    mean_params <- colMeans(params, na.rm = TRUE)
    mean_hrf <- hrf_interface$hrf_function(times, mean_params)
    lines(times, mean_hrf, col = "red", lwd = 2)
    legend("topright", c("Individual", "Mean"), 
           col = c("gray70", "red"), lwd = c(1, 2))
  }
}

# Helper function to plot parameter distributions
.plot_parameters <- function(x, ...) {
  params <- get_parameters(x)
  n_params <- ncol(params)
  
  # Setup multi-panel plot
  old_par <- par(mfrow = c(1, n_params), mar = c(4, 4, 2, 1))
  on.exit(par(old_par))
  
  for (i in seq_len(n_params)) {
    hist(params[, i], 
         main = x$parameter_names[i],
         xlab = x$parameter_names[i],
         col = "lightblue", ...)
    
    # Add mean line
    abline(v = mean(params[, i], na.rm = TRUE), 
           col = "red", lwd = 2, lty = 2)
  }
}

# Helper function to plot model quality metrics
.plot_quality <- function(x, ...) {
  if (is.null(x$r_squared)) {
    warning("No R-squared values available")
    return(invisible(NULL))
  }
  
  # R-squared histogram
  hist(x$r_squared, 
       main = "Model Fit Quality (R^2)",
       xlab = "R-squared",
       col = "lightgreen",
       breaks = 30, ...)
  
  # Add summary lines
  abline(v = mean(x$r_squared, na.rm = TRUE), col = "red", lwd = 2, lty = 1)
  abline(v = median(x$r_squared, na.rm = TRUE), col = "blue", lwd = 2, lty = 2)
  
  legend("topleft", 
         c(paste("Mean:", round(mean(x$r_squared, na.rm = TRUE), 3)),
           paste("Median:", round(median(x$r_squared, na.rm = TRUE), 3))),
         col = c("red", "blue"), lty = c(1, 2), lwd = 2)
}

# Helper function to plot refinement diagnostics
.plot_refinement <- function(x, ...) {
  if (is.null(x$refinement_info) || is.null(x$refinement_info$applied) || !x$refinement_info$applied) {
    warning("No refinement information available")
    return(invisible(NULL))
  }
  
  # Check for before/after R-squared
  if (!is.null(x$refinement_info$r2_before) && 
      !is.null(x$refinement_info$r2_after)) {
    
    # Create before/after comparison
    r2_change <- x$refinement_info$r2_after - x$refinement_info$r2_before
    
    # Plot improvement
    hist(r2_change[r2_change != 0], 
         main = "R^2 Improvement from Refinement",
         xlab = "Change in R^2",
         col = "lightcoral", ...)
    
    abline(v = 0, col = "gray", lwd = 2, lty = 2)
    abline(v = mean(r2_change[r2_change != 0], na.rm = TRUE), 
           col = "red", lwd = 2)
    
    # Add summary
    improved <- sum(r2_change > 0, na.rm = TRUE)
    total_refined <- sum(r2_change != 0, na.rm = TRUE)
    
    legend("topright", 
           c(paste("Improved:", improved, "/", total_refined),
             paste("Mean improvement:", 
                   round(mean(r2_change[r2_change > 0], na.rm = TRUE), 4))),
           bty = "n")
  }
}

#' Get fitted values from parametric_hrf_fit
#' @param object parametric_hrf_fit object  
#' @param Y_proj Original projected Y data (required if residuals not stored)
#' @param ... Additional arguments (ignored)
#' @return Matrix of fitted values
#' @export
fitted.parametric_hrf_fit <- function(object, Y_proj = NULL, ...) {
  # Implementation depends on whether residuals are stored
  if (!is.null(object$residuals) && !is.null(Y_proj)) {
    return(Y_proj - object$residuals)
  }
  
  # Otherwise need Y_proj to compute fitted values
  if (is.null(Y_proj)) {
    stop("Y_proj required")
  }
  
  # Need design matrix to recompute
  if (is.null(object$metadata$S_target_proj)) {
    stop("Cannot compute fitted values without design matrix")
  }
  
  params <- get_parameters(object)
  amplitudes <- object$amplitudes
  S <- object$metadata$S_target_proj
  times <- object$metadata$hrf_eval_times
  
  # Get HRF interface
  hrf_interface <- .create_hrf_interface(object$parametric_model)
  
  # Compute fitted values for each voxel
  n_time <- nrow(S)
  n_vox <- nrow(params)
  fitted_vals <- matrix(0, n_time, n_vox)
  
  for (v in seq_len(n_vox)) {
    hrf_v <- hrf_interface$hrf_function(times, params[v,])
    conv_v <- .batch_convolution(S, matrix(hrf_v, ncol = 1), n_time)
    fitted_vals[, v] <- amplitudes[v] * conv_v[, 1]
  }
  
  fitted_vals
}

#' Predict method for parametric_hrf_fit
#'
#' Generate predictions for new stimulus designs using fitted HRF parameters.
#'
#' @param object A `parametric_hrf_fit` object
#' @param newdata Optional stimulus design matrix in projected space
#'   (`S_target_proj` format). If `NULL`, uses the design stored in
#'   `object$metadata$S_target_proj` when available.
#' @param voxel_indices Integer vector of voxel indices to predict. Defaults to
#'   all voxels.
#' @param ... Additional arguments (ignored)
#' @return Numeric matrix of predicted time series (timepoints x voxels)
#' @export
predict.parametric_hrf_fit <- function(object,
                                       newdata = NULL,
                                       voxel_indices = NULL,
                                       ...) {
  params <- get_parameters(object)
  
  if (is.null(voxel_indices)) {
    voxel_indices <- seq_len(nrow(params))
  }

  if (is.null(newdata)) {
    newdata <- object$metadata$S_target_proj
    if (is.null(newdata)) {
      stop("newdata must be supplied when original design is unavailable")
    }
  }

  if (is.null(dim(newdata))) {
    newdata <- matrix(newdata, ncol = 1)
  }

  hrf_eval_times <- object$metadata$hrf_eval_times
  if (is.null(hrf_eval_times)) {
    hrf_eval_times <- seq(0, 30, length.out = 61)
  }

  # Use factory instead of registry
  hrf_interface <- .create_hrf_interface(object$parametric_model, 
                                        user_bounds = object$metadata$theta_bounds)
  
  n_time <- nrow(newdata)
  result <- matrix(NA_real_, nrow = n_time, ncol = length(voxel_indices))

  for (i in seq_along(voxel_indices)) {
    v <- voxel_indices[i]
    theta_v <- params[v, ]
    beta0_v <- object$amplitudes[v]
    hrf_vals <- hrf_interface$hrf_function(hrf_eval_times, theta_v)
    conv <- .batch_convolution(newdata, matrix(hrf_vals, ncol = 1), n_time)
    result[, i] <- beta0_v * conv[, 1]
  }

  result
}

#' Get parameters from parametric_hrf_fit
#' 
#' @param fit A parametric_hrf_fit object
#' @return Matrix of estimated parameters
#' @export
get_parameters <- function(fit) {
  fit$model_specific$parameters
}

#' Get parameter standard errors from parametric_hrf_fit
#' 
#' @param fit A parametric_hrf_fit object
#' @return Matrix of parameter standard errors
#' @export
get_parameter_ses <- function(fit) {
  fit$model_specific$standard_errors
}

#' Get model name from parametric_hrf_fit
#' 
#' @param fit A parametric_hrf_fit object
#' @return Character string of model name
#' @export
get_model_name <- function(fit) {
  fit$parametric_model
}

#' Get design metadata from parametric_hrf_fit
#'
#' Returns normalized design information aligned with other fit-object based
#' modeling packages.
#'
#' @param fit A parametric_hrf_fit object.
#' @return A named list with `n_time`, `n_vox`, `n_cond`, `basis_dim`, and
#'   `projected`.
#' @export
get_design_info <- function(fit) {
  info <- fit$metadata$design_info
  if (is.null(info)) {
    info <- list(
      n_time = fit$metadata$n_timepoints,
      n_vox = fit$metadata$n_voxels,
      n_cond = NA_integer_,
      basis_dim = ncol(get_parameters(fit)),
      projected = NA
    )
  }
  info
}

#' Get per-voxel goodness-of-fit from parametric_hrf_fit
#'
#' @param fit A parametric_hrf_fit object.
#' @return Numeric vector of per-voxel goodness-of-fit (`R^2`).
#' @export
get_gof_per_voxel <- function(fit) {
  fit$r_squared
}

#' Get estimation method label from parametric_hrf_fit
#'
#' @param fit A parametric_hrf_fit object.
#' @return Character scalar with method label.
#' @export
get_method_used <- function(fit) {
  method <- fit$metadata$method_used
  if (is.null(method) || !is.character(method) || length(method) != 1) {
    return("parametric_taylor")
  }
  method
}
