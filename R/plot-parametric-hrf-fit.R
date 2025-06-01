#' Plot method for parametric_hrf_fit objects
#'
#' Visualizes the estimated HRF shape for selected voxels or summary statistics.
#' Can display individual voxel HRFs or aggregate HRFs.
#'
#' @param x An object of class \code{parametric_hrf_fit}
#' @param voxel_indices Integer vector of voxel indices to plot. If NULL,
#'   plots the median parameter HRF
#' @param type Character string: "hrf" for HRF shape, "parameters" for 
#'   parameter distributions
#' @param hrf_time_max Maximum time in seconds for HRF plot
#' @param add_amplitude Logical whether to scale HRF by amplitude
#' @param ... Additional arguments passed to plotting functions
#' 
#' @return Invisible NULL. Called for side effect of plotting.
#' 
#' @examples
#' \dontrun{
#' # Plot median HRF
#' plot(fit)
#' 
#' # Plot specific voxels
#' plot(fit, voxel_indices = c(100, 200, 300))
#' 
#' # Plot parameter distributions
#' plot(fit, type = "parameters")
#' }
#' 
#' @export
plot.parametric_hrf_fit <- function(x, 
                                    voxel_indices = NULL, 
                                    type = c("hrf", "parameters"),
                                    hrf_time_max = 30,
                                    add_amplitude = TRUE,
                                    ...) {
  type <- match.arg(type)
  
  if (type == "parameters") {
    # Parameter distribution plots
    plot_parameter_distributions(x, ...)
  } else {
    # HRF shape plots
    plot_hrf_shapes(x, voxel_indices, hrf_time_max, add_amplitude, ...)
  }
  
  invisible(NULL)
}

#' Plot HRF shapes for selected voxels
#' @keywords internal
plot_hrf_shapes <- function(x, voxel_indices, hrf_time_max, add_amplitude, ...) {
  # Get HRF interface based on model
  if (x$hrf_model == "lwu") {
    hrf_function <- .lwu_hrf_function
  } else {
    stop("Unknown HRF model: ", x$hrf_model)
  }
  
  # Time points for evaluation
  t_plot <- seq(0, hrf_time_max, by = 0.1)
  
  # Determine which voxels to plot
  if (is.null(voxel_indices)) {
    # Use median parameters
    median_params <- apply(x$estimated_parameters, 2, median, na.rm = TRUE)
    median_amp <- median(x$amplitudes, na.rm = TRUE)
    
    # Generate HRF
    hrf_shape <- hrf_function(t_plot, median_params)
    if (add_amplitude) {
      hrf_shape <- median_amp * hrf_shape
    }
    
    # Plot
    plot(t_plot, hrf_shape, type = "l", lwd = 2,
         xlab = "Time (s)", ylab = "HRF",
         main = "Median HRF across voxels",
         ...)
    abline(h = 0, lty = 2, col = "gray")
    
    # Add parameter text
    param_text <- paste(x$parameter_names, "=", round(median_params, 2), collapse = ", ")
    mtext(param_text, side = 3, line = 0.5, cex = 0.8)
    
  } else {
    # Plot multiple voxels
    n_vox <- length(voxel_indices)
    if (n_vox > 20) {
      warning("Plotting only first 20 voxels")
      voxel_indices <- voxel_indices[1:20]
      n_vox <- 20
    }
    
    # Set up colors
    cols <- rainbow(n_vox)
    
    # Initialize plot
    plot(NULL, xlim = c(0, hrf_time_max), 
         ylim = c(-0.5, 1.5) * ifelse(add_amplitude, max(abs(x$amplitudes[voxel_indices])), 1),
         xlab = "Time (s)", ylab = "HRF",
         main = paste("HRF shapes for", n_vox, "voxels"),
         ...)
    abline(h = 0, lty = 2, col = "gray")
    
    # Plot each voxel
    for (i in seq_along(voxel_indices)) {
      v <- voxel_indices[i]
      params_v <- x$estimated_parameters[v, ]
      
      # Generate HRF
      hrf_shape <- hrf_function(t_plot, params_v)
      if (add_amplitude) {
        hrf_shape <- x$amplitudes[v] * hrf_shape
      }
      
      lines(t_plot, hrf_shape, col = cols[i], lwd = 1.5)
    }
    
    # Add legend if few voxels
    if (n_vox <= 10) {
      legend("topright", 
             legend = paste("Voxel", voxel_indices),
             col = cols, lty = 1, lwd = 1.5, cex = 0.8)
    }
  }
}

#' Plot parameter distributions
#' @keywords internal  
plot_parameter_distributions <- function(x, ...) {
  n_params <- length(x$parameter_names)
  
  # Set up multi-panel plot
  if (n_params <= 3) {
    par(mfrow = c(1, n_params))
  } else {
    par(mfrow = c(2, ceiling(n_params/2)))
  }
  
  # Plot each parameter
  for (i in seq_len(n_params)) {
    param_vals <- x$estimated_parameters[, i]
    param_name <- x$parameter_names[i]
    
    # Basic histogram
    hist(param_vals, 
         main = paste("Distribution of", param_name),
         xlab = param_name,
         col = "lightblue",
         border = "darkblue",
         ...)
    
    # Add median line
    abline(v = median(param_vals, na.rm = TRUE), 
           col = "red", lwd = 2, lty = 2)
    
    # Add R² information if available
    if (!is.null(x$r_squared)) {
      # Color by R²
      good_vox <- which(x$r_squared > 0.2)
      if (length(good_vox) > 0) {
        hist(param_vals[good_vox], 
             add = TRUE,
             col = rgb(0, 0, 1, 0.3),
             border = NA)
      }
    }
  }
  
  # Add amplitude distribution
  if (n_params < 4) {
    hist(x$amplitudes,
         main = "Distribution of Amplitudes",
         xlab = "Amplitude",
         col = "lightgreen",
         border = "darkgreen",
         ...)
    abline(v = median(x$amplitudes, na.rm = TRUE),
           col = "red", lwd = 2, lty = 2)
  }
  
  # Reset par
  par(mfrow = c(1, 1))
}