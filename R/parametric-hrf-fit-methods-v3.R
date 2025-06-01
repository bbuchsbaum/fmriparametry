#' Enhanced S3 methods for parametric_hrf_fit with refinement information
#'
#' These methods provide enhanced functionality for Sprint 3 features including
#' refinement diagnostics, parallel processing info, and advanced visualizations.

#' Print a parametric_hrf_fit object (Sprint 3 version)
#'
#' Displays a comprehensive summary including refinement information, convergence
#' diagnostics, and parallel processing details.
#'
#' @param x An object of class \code{parametric_hrf_fit}
#' @param ... Additional arguments (currently ignored)
#' 
#' @return The input object invisibly
#' @export
print.parametric_hrf_fit <- function(x, ...) {
  cat("Parametric HRF Fit\n")
  cat("==================\n")
  cat("Model:", x$hrf_model, "\n")
  cat("Voxels:", nrow(x$estimated_parameters), "\n")
  
  # Basic parameter summary
  cat("\nParameter Summary:\n")
  param_means <- colMeans(x$estimated_parameters, na.rm = TRUE)
  param_sds <- apply(x$estimated_parameters, 2, sd, na.rm = TRUE)
  param_summary <- data.frame(
    Parameter = x$parameter_names,
    Mean = round(param_means, 3),
    SD = round(param_sds, 3)
  )
  print(param_summary, row.names = FALSE)
  
  # R-squared summary
  if (!is.null(x$r_squared)) {
    cat("\nModel Fit (R²):\n")
    r2_summary <- summary(x$r_squared)
    cat("  Min:", round(r2_summary[1], 3), "\n")
    cat("  Median:", round(r2_summary[3], 3), "\n")
    cat("  Mean:", round(mean(x$r_squared, na.rm = TRUE), 3), "\n")
    cat("  Max:", round(r2_summary[6], 3), "\n")
    cat("  % R² > 0.5:", round(100 * mean(x$r_squared > 0.5, na.rm = TRUE), 1), "%\n")
  }
  
  # Refinement information
  if (!is.null(x$metadata$refinement_info) && x$metadata$refinement_info$applied) {
    cat("\nRefinement Applied:\n")
    ref_info <- x$metadata$refinement_info
    cat("  Moderate voxels refined:", ref_info$n_moderate_refined, "\n")
    cat("  Hard voxels refined:", ref_info$n_hard_refined, "\n")
    if (ref_info$n_hard_refined > 0) {
      cat("  Gauss-Newton converged:", ref_info$n_converged,
          "(", round(100 * ref_info$n_converged / ref_info$n_hard_refined, 1), "%)\n")
      cat("  Improved after refinement:", ref_info$n_improved,
          "(", round(100 * ref_info$n_improved / ref_info$n_hard_refined, 1), "%)\n")
    }
  }
  
  # Convergence information
  if (!is.null(x$convergence_info)) {
    cat("\nConvergence:\n")
    cat("  Global iterations:", x$convergence_info$global_iterations, "\n")
    if (!is.null(x$convergence_info$converged)) {
      cat("  Converged:", ifelse(x$convergence_info$converged, "Yes", "No"), "\n")
    }
    if (!is.null(x$metadata$kmeans_info) && x$metadata$kmeans_info$applied) {
      cat("  K-means clusters:", x$metadata$kmeans_info$n_clusters, "\n")
      cat("  K-means iterations:", x$metadata$kmeans_info$total_iterations, "\n")
    }
  }
  
  # Parallel processing info
  if (!is.null(x$metadata$parallel_info)) {
    cat("\nParallel Processing:\n")
    cat("  Backend:", x$metadata$parallel_info$backend, "\n")
    cat("  Cores used:", x$metadata$parallel_info$n_cores, "\n")
  }
  
  invisible(x)
}

#' Summary method with refinement information
#'
#' Produces comprehensive summary statistics including refinement diagnostics
#' and queue information.
#'
#' @param object An object of class \code{parametric_hrf_fit}
#' @param ... Additional arguments
#' 
#' @return A list with comprehensive summary information
#' @export
summary.parametric_hrf_fit <- function(object, ...) {
  # Basic summaries
  param_sum <- apply(object$estimated_parameters, 2, summary)
  amp_sum <- summary(object$amplitudes)
  
  # R-squared summary
  r2_sum <- if (!is.null(object$r_squared)) {
    list(
      summary = summary(object$r_squared),
      prop_good = mean(object$r_squared > 0.5, na.rm = TRUE),
      prop_excellent = mean(object$r_squared > 0.7, na.rm = TRUE),
      n_failed = sum(object$r_squared < 0.1, na.rm = TRUE)
    )
  } else NULL
  
  # Standard error summary
  se_sum <- if (!is.null(object$parameter_ses)) {
    apply(object$parameter_ses, 2, function(x) summary(x[is.finite(x)]))
  } else NULL
  
  # Refinement summary
  refinement_sum <- if (!is.null(object$metadata$refinement_info) && 
                        object$metadata$refinement_info$applied) {
    ref_info <- object$metadata$refinement_info
    queue_sum <- ref_info$final_queue_summary
    
    list(
      applied = TRUE,
      queue_summary = queue_sum,
      queue_proportions = prop.table(queue_sum),
      n_moderate_refined = ref_info$n_moderate_refined,
      n_hard_refined = ref_info$n_hard_refined,
      n_converged = ref_info$n_converged,
      n_improved = ref_info$n_improved,
      improvement_rate = if (ref_info$n_hard_refined > 0) {
        ref_info$n_improved / ref_info$n_hard_refined
      } else NA
    )
  } else {
    list(applied = FALSE)
  }
  
  # Construct output
  structure(
    list(
      parameter_summary = param_sum,
      amplitude_summary = amp_sum,
      r_squared_summary = r2_sum,
      se_summary = se_sum,
      refinement_summary = refinement_sum,
      hrf_model = object$hrf_model,
      n_voxels = nrow(object$estimated_parameters),
      n_timepoints = object$metadata$n_timepoints,
      computation_time = object$metadata$computation_time,
      parallel_info = object$metadata$parallel_info
    ),
    class = "summary.parametric_hrf_fit"
  )
}

#' Print summary.parametric_hrf_fit
#' 
#' @param x A summary.parametric_hrf_fit object
#' @param ... Additional arguments
#' @export
print.summary.parametric_hrf_fit <- function(x, ...) {
  cat("Summary of Parametric HRF Fit\n")
  cat("=============================\n")
  cat("\nModel:", x$hrf_model, "\n")
  cat("Voxels:", x$n_voxels, "\n")
  cat("Timepoints:", x$n_timepoints, "\n")
  
  cat("\nParameter Estimates:\n")
  print(round(x$parameter_summary, 3))
  
  cat("\nAmplitude Summary:\n")
  print(round(x$amplitude_summary, 3))
  
  if (!is.null(x$r_squared_summary)) {
    cat("\nModel Fit (R²):\n")
    print(round(x$r_squared_summary$summary, 3))
    cat("Proportion R² > 0.5:", round(100 * x$r_squared_summary$prop_good, 1), "%\n")
    cat("Proportion R² > 0.7:", round(100 * x$r_squared_summary$prop_excellent, 1), "%\n")
    cat("Failed voxels (R² < 0.1):", x$r_squared_summary$n_failed, "\n")
  }
  
  if (!is.null(x$se_summary)) {
    cat("\nStandard Errors:\n")
    print(round(x$se_summary, 3))
  }
  
  if (x$refinement_summary$applied) {
    cat("\nRefinement Summary:\n")
    cat("Queue distribution:\n")
    queue_df <- data.frame(
      Queue = names(x$refinement_summary$queue_summary),
      Count = as.numeric(x$refinement_summary$queue_summary),
      Proportion = round(100 * as.numeric(x$refinement_summary$queue_proportions), 1)
    )
    print(queue_df, row.names = FALSE)
    
    cat("\nRefinement results:\n")
    cat("  Moderate voxels refined:", x$refinement_summary$n_moderate_refined, "\n")
    cat("  Hard voxels refined:", x$refinement_summary$n_hard_refined, "\n")
    if (x$refinement_summary$n_hard_refined > 0) {
      cat("  Gauss-Newton converged:", x$refinement_summary$n_converged, "\n")
      cat("  Improved after refinement:", x$refinement_summary$n_improved,
          "(", round(100 * x$refinement_summary$improvement_rate, 1), "%)\n")
    }
  }
  
  if (!is.null(x$parallel_info)) {
    cat("\nComputation:\n")
    cat("  Parallel backend:", x$parallel_info$backend, "\n")
    cat("  Cores used:", x$parallel_info$n_cores, "\n")
  }
  
  if (!is.null(x$computation_time)) {
    cat("  Total time:", round(x$computation_time, 1), "seconds\n")
  }
  
  invisible(x)
}

#' Plot method with refinement diagnostics
#'
#' Enhanced plotting method that can visualize HRFs, parameter distributions,
#' diagnostic plots, and refinement information.
#'
#' @param x An object of class \code{parametric_hrf_fit}
#' @param type Plot type: "hrf", "parameters", "diagnostic", "refinement"
#' @param voxels Which voxels to plot (for "hrf" type)
#' @param ... Additional arguments passed to plotting functions
#' 
#' @export
plot.parametric_hrf_fit <- function(x, 
                                    type = c("hrf", "parameters", "diagnostic", "refinement"),
                                    voxels = NULL,
                                    ...) {
  type <- match.arg(type)
  
  # Ensure we have required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting")
  }
  
  switch(type,
    hrf = .plot_hrf(x, voxels, ...),
    parameters = .plot_parameters(x, ...),
    diagnostic = .plot_diagnostic(x, ...),
    refinement = .plot_refinement(x, ...)
  )
}

#' Plot HRF curves
#' @keywords internal
.plot_hrf <- function(x, voxels = NULL, n_curves = 10, ...) {
  # Select voxels to plot
  if (is.null(voxels)) {
    # Sample across R² range
    r2_quantiles <- quantile(x$r_squared, probs = seq(0.1, 0.9, length.out = n_curves))
    voxels <- sapply(r2_quantiles, function(q) {
      which.min(abs(x$r_squared - q))[1]
    })
  }
  
  # Generate HRF curves
  t_hrf <- seq(0, 30, by = 0.1)
  hrf_curves <- matrix(NA, nrow = length(t_hrf), ncol = length(voxels))
  
  # Load HRF function
  if (x$hrf_model == "lwu") {
    source(system.file("R", "hrf-interface-lwu.R", package = "fmriparametric"), local = TRUE)
    hrf_fn <- .lwu_hrf_function
  }
  
  for (i in seq_along(voxels)) {
    v <- voxels[i]
    theta <- x$estimated_parameters[v, ]
    hrf_curves[, i] <- x$amplitudes[v] * hrf_fn(t_hrf, theta)
  }
  
  # Create plot data
  plot_data <- data.frame(
    time = rep(t_hrf, length(voxels)),
    hrf = as.vector(hrf_curves),
    voxel = factor(rep(voxels, each = length(t_hrf))),
    r2 = rep(x$r_squared[voxels], each = length(t_hrf))
  )
  
  ggplot2::ggplot(plot_data, ggplot2::aes(x = time, y = hrf, 
                                           color = r2, group = voxel)) +
    ggplot2::geom_line(alpha = 0.7) +
    ggplot2::scale_color_viridis_c(name = expression(R^2)) +
    ggplot2::labs(
      title = "Estimated HRF Curves",
      x = "Time (seconds)",
      y = "Response"
    ) +
    ggplot2::theme_minimal()
}

#' Plot parameter distributions
#' @keywords internal
.plot_parameters <- function(x, ...) {
  # Prepare data
  param_data <- as.data.frame(x$estimated_parameters)
  colnames(param_data) <- x$parameter_names
  param_data$r2 <- x$r_squared
  
  # Reshape for faceting
  param_long <- reshape(param_data, 
                        direction = "long",
                        varying = x$parameter_names,
                        v.names = "value",
                        timevar = "parameter",
                        times = x$parameter_names)
  
  ggplot2::ggplot(param_long, ggplot2::aes(x = value, fill = parameter)) +
    ggplot2::geom_histogram(bins = 30, alpha = 0.7) +
    ggplot2::facet_wrap(~ parameter, scales = "free") +
    ggplot2::labs(
      title = "Parameter Distributions",
      x = "Parameter Value",
      y = "Count"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
}

#' Plot diagnostics
#' @keywords internal
.plot_diagnostic <- function(x, ...) {
  if (is.null(x$r_squared)) {
    stop("No R-squared values available for diagnostic plot")
  }
  
  # Create diagnostic data
  diag_data <- data.frame(
    voxel = seq_len(nrow(x$estimated_parameters)),
    r2 = x$r_squared,
    amplitude = x$amplitudes
  )
  
  # Add parameter values
  for (i in seq_along(x$parameter_names)) {
    diag_data[[x$parameter_names[i]]] <- x$estimated_parameters[, i]
  }
  
  # R² histogram
  p1 <- ggplot2::ggplot(diag_data, ggplot2::aes(x = r2)) +
    ggplot2::geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    ggplot2::geom_vline(xintercept = c(0.3, 0.5, 0.7), 
                        linetype = "dashed", alpha = 0.5) +
    ggplot2::labs(title = "R² Distribution", x = expression(R^2), y = "Count") +
    ggplot2::theme_minimal()
  
  # Parameter vs R² scatter
  param_plots <- list()
  for (param in x$parameter_names) {
    param_plots[[param]] <- ggplot2::ggplot(diag_data, 
                                             ggplot2::aes_string(x = param, y = "r2")) +
      ggplot2::geom_point(alpha = 0.3, size = 0.5) +
      ggplot2::geom_smooth(method = "loess", se = FALSE, color = "red") +
      ggplot2::labs(title = paste(param, "vs R²"), x = param, y = expression(R^2)) +
      ggplot2::theme_minimal()
  }
  
  # Combine plots
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(p1, grobs = param_plots, ncol = 2)
  } else {
    p1
  }
}

#' Plot refinement information
#' @keywords internal
.plot_refinement <- function(x, ...) {
  if (is.null(x$metadata$refinement_info) || !x$metadata$refinement_info$applied) {
    stop("No refinement information available")
  }
  
  ref_info <- x$metadata$refinement_info
  queue_result <- ref_info$queue_result
  
  # Create refinement data
  ref_data <- data.frame(
    voxel = seq_len(nrow(x$estimated_parameters)),
    r2_initial = x$r_squared,  # This would be post-refinement
    queue = queue_result$queue_labels
  )
  
  # Queue distribution pie chart
  queue_summary <- as.data.frame(table(ref_data$queue))
  colnames(queue_summary) <- c("Queue", "Count")
  
  p1 <- ggplot2::ggplot(queue_summary, ggplot2::aes(x = "", y = Count, fill = Queue)) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::coord_polar("y", start = 0) +
    ggplot2::labs(title = "Refinement Queue Distribution") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.title = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank())
  
  # R² by queue
  p2 <- ggplot2::ggplot(ref_data, ggplot2::aes(x = queue, y = r2_initial, fill = queue)) +
    ggplot2::geom_boxplot(alpha = 0.7) +
    ggplot2::labs(
      title = "R² by Refinement Queue",
      x = "Queue",
      y = expression(R^2)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
  
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(p1, p2, ncol = 2)
  } else {
    p1
  }
}

#' Extract residuals
#' 
#' @param object A parametric_hrf_fit object
#' @param ... Additional arguments
#' @export
residuals.parametric_hrf_fit <- function(object, ...) {
  if (is.null(object$residuals)) {
    warning("No residuals stored in the fit object")
    return(NULL)
  }
  object$residuals
}

#' Extract fitted values
#' 
#' @param object A parametric_hrf_fit object
#' @param ... Additional arguments
#' @export
fitted.parametric_hrf_fit <- function(object, ...) {
  if (is.null(object$residuals)) {
    warning("Cannot compute fitted values without residuals")
    return(NULL)
  }
  
  # Fitted = Original - Residuals
  # This assumes Y_proj is available, which it isn't in the fit object
  # So we'd need to reconstruct or store fitted values
  warning("Fitted values not directly available. Use residuals and original data.")
  NULL
}