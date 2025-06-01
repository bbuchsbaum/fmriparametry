#' Enhanced S3 methods for parametric_hrf_fit (Sprint 2)
#'
#' These methods provide comprehensive functionality for working with
#' parametric HRF fit objects, including support for Sprint 2 enhancements.

#' Print a parametric_hrf_fit object
#'
#' Displays a concise summary of the parametric HRF fit including the HRF model 
#' used, number of voxels analyzed, basic parameter statistics, and R-squared
#' distribution if available.
#'
#' @param x An object of class \code{parametric_hrf_fit}
#' @param ... Additional arguments (currently ignored)
#' @return The input object \code{x} invisibly
#' @export
print.parametric_hrf_fit <- function(x, ...) {
  cat("Parametric HRF Fit\n")
  cat("Model:", x$hrf_model, "\n")
  cat("Voxels:", nrow(x$estimated_parameters), "\n")
  
  # Check if Sprint 2 enhancements present
  if (!is.null(x$r_squared)) {
    cat("Mean R²:", round(mean(x$r_squared, na.rm = TRUE), 3), "\n")
  }
  
  if (!is.null(x$convergence_info) && length(x$convergence_info$trajectory) > 1) {
    cat("Global iterations:", x$convergence_info$n_iterations, "\n")
  }
  
  cat("\nParameter Summary:\n")
  print(summary(x$estimated_parameters))
  invisible(x)
}

#' Summarize a parametric_hrf_fit object
#'
#' Produces comprehensive summary statistics of the estimated HRF parameters,
#' amplitudes, R-squared values, and convergence information.
#'
#' @param object An object of class \code{parametric_hrf_fit}
#' @param ... Additional arguments (currently ignored)
#' @return A list containing summary information
#' @export
summary.parametric_hrf_fit <- function(object, ...) {
  # Basic summaries
  param_sum <- apply(object$estimated_parameters, 2, summary)
  amp_sum <- summary(object$amplitudes)
  
  # Sprint 2 additions
  r2_sum <- NULL
  if (!is.null(object$r_squared)) {
    r2_sum <- summary(object$r_squared)
  }
  
  # Parameter standard errors summary
  se_sum <- NULL
  if (!is.null(object$parameter_ses)) {
    se_sum <- apply(object$parameter_ses, 2, function(x) {
      summary(x[!is.na(x)])
    })
  }
  
  # Convergence information
  conv_info <- NULL
  if (!is.null(object$convergence_info) && length(object$convergence_info$trajectory) > 1) {
    conv_info <- list(
      n_iterations = object$convergence_info$n_iterations,
      converged = object$convergence_info$converged,
      final_global_theta = object$convergence_info$final_global_theta
    )
    
    # Delta theta from start to finish
    if (length(object$convergence_info$trajectory) >= 2) {
      theta_start <- object$convergence_info$trajectory[[1]]
      theta_end <- object$convergence_info$trajectory[[length(object$convergence_info$trajectory)]]
      conv_info$delta_theta <- theta_end - theta_start
    }
  }
  
  structure(
    list(
      parameter_summary = param_sum,
      amplitude_summary = amp_sum,
      r_squared_summary = r2_sum,
      parameter_se_summary = se_sum,
      convergence_info = conv_info,
      hrf_model = object$hrf_model,
      n_voxels = n_voxels(object),
      n_timepoints = n_timepoints(object),
      call = object$metadata$call
    ),
    class = "summary.parametric_hrf_fit"
  )
}

#' Print method for summary.parametric_hrf_fit
#' @param x summary.parametric_hrf_fit object
#' @param ... Additional arguments
#' @export
print.summary.parametric_hrf_fit <- function(x, ...) {
  cat("\nCall:\n")
  print(x$call)
  
  cat("\nModel:", x$hrf_model, "\n")
  cat("Voxels:", x$n_voxels, "\n")
  cat("Timepoints:", x$n_timepoints, "\n")
  
  cat("\nParameter Estimates:\n")
  print(x$parameter_summary)
  
  cat("\nAmplitude Summary:\n")
  print(x$amplitude_summary)
  
  if (!is.null(x$r_squared_summary)) {
    cat("\nR-squared Summary:\n")
    print(x$r_squared_summary)
  }
  
  if (!is.null(x$parameter_se_summary)) {
    cat("\nParameter Standard Errors:\n")
    print(x$parameter_se_summary)
  }
  
  if (!is.null(x$convergence_info)) {
    cat("\nConvergence Information:\n")
    cat("  Iterations:", x$convergence_info$n_iterations, "\n")
    cat("  Converged:", x$convergence_info$converged, "\n")
    if (!is.null(x$convergence_info$delta_theta)) {
      cat("  Parameter change:", 
          paste(round(x$convergence_info$delta_theta, 3), collapse = ", "), "\n")
    }
  }
  
  invisible(x)
}

#' Extract coefficients from a parametric_hrf_fit object
#'
#' Returns the matrix of estimated HRF parameters or amplitudes. Can optionally
#' return standard errors.
#'
#' @param object An object of class \code{parametric_hrf_fit}
#' @param type Character string: "parameters" (default), "amplitude", or "se"
#' @param ... Additional arguments (currently ignored)
#' @return A numeric matrix or vector depending on type
#' @export
coef.parametric_hrf_fit <- function(object, type = c("parameters", "amplitude", "se"), ...) {
  type <- match.arg(type)
  
  switch(type,
    parameters = object$estimated_parameters,
    amplitude = object$amplitudes,
    se = {
      if (is.null(object$parameter_ses)) {
        warning("Standard errors not computed for this fit")
        NULL
      } else {
        object$parameter_ses
      }
    }
  )
}

#' Extract fitted values from a parametric_hrf_fit object
#'
#' Returns the matrix of fitted values. Requires either stored residuals
#' or the original data to reconstruct fitted values.
#'
#' @param object An object of class \code{parametric_hrf_fit}
#' @param Y_proj Original projected Y data (required if residuals not stored)
#' @param ... Additional arguments (currently ignored)
#' @return Numeric matrix of fitted values (timepoints x voxels)
#' @export
fitted.parametric_hrf_fit <- function(object, Y_proj = NULL, ...) {
  if (!is.null(object$residuals)) {
    if (is.null(Y_proj)) {
      stop("Y_proj required to compute fitted values from residuals")
    }
    # Fitted = Y - residuals
    return(Y_proj - object$residuals)
  } else {
    stop("Cannot compute fitted values without residuals. ",
         "Re-fit with compute_residuals = TRUE or provide design matrix.")
  }
}

#' Extract residuals from a parametric_hrf_fit object
#'
#' Returns the matrix of residuals if computed during fitting.
#'
#' @param object An object of class \code{parametric_hrf_fit}
#' @param ... Additional arguments (currently ignored)
#' @return Numeric matrix of residuals (timepoints x voxels) or NULL
#' @export
residuals.parametric_hrf_fit <- function(object, ...) {
  if (is.null(object$residuals)) {
    warning("Residuals not computed for this fit. ",
            "Re-fit with compute_residuals = TRUE.")
  }
  object$residuals
}

#' Extract variance-covariance matrix for a specific voxel
#'
#' Returns the variance-covariance matrix of parameter estimates for a
#' specified voxel. Only available if standard errors were computed.
#'
#' @param object An object of class \code{parametric_hrf_fit}
#' @param voxel_index Integer index of the voxel
#' @param ... Additional arguments (currently ignored)
#' @return Numeric matrix (parameters x parameters) or NULL
#' @export
vcov.parametric_hrf_fit <- function(object, voxel_index = 1, ...) {
  if (is.null(object$parameter_ses)) {
    warning("Standard errors not computed for this fit")
    return(NULL)
  }
  
  assertthat::assert_that(
    is.numeric(voxel_index),
    length(voxel_index) == 1,
    voxel_index >= 1,
    voxel_index <= nrow(object$estimated_parameters)
  )
  
  # For now, return diagonal matrix with SE^2
  # Full covariance would require storing more information
  se_vox <- object$parameter_ses[voxel_index, ]
  if (any(is.na(se_vox))) {
    warning("Standard errors not available for voxel ", voxel_index)
    return(NULL)
  }
  
  diag(se_vox^2)
}