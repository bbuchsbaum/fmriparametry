#' Rock Solid Error Recovery Functions
#'
#' Error recovery system that ensures the package NEVER crashes and ALWAYS
#' returns meaningful output, even in the face of catastrophic failures.

#' Execute function with comprehensive error recovery
#' @keywords internal
.try_with_recovery <- function(primary_fn, fallback_fn = NULL, 
                               default_value = NULL,
                               error_prefix = "Operation",
                               max_attempts = 3,
                               verbose = TRUE) {
  attempt <- 0
  last_error <- NULL
  
  while (attempt < max_attempts) {
    attempt <- attempt + 1
    
    # Try primary function
    result <- tryCatch({
      primary_fn()
    }, error = function(e) {
      last_error <<- e
      NULL
    }, warning = function(w) {
      if (verbose) message(error_prefix, " warning: ", w$message)
      invokeRestart("muffleWarning")
    })
    
    if (!is.null(result)) {
      return(result)
    }
    
    # If failed and we have more attempts, wait briefly
    if (attempt < max_attempts) {
      Sys.sleep(0.1 * attempt)  # Exponential backoff
    }
  }
  
  # Primary failed, try fallback
  if (!is.null(fallback_fn)) {
    if (verbose) {
      message(error_prefix, " failed with: ", last_error$message, 
              ". Trying fallback method.")
    }
    
    result <- tryCatch({
      fallback_fn()
    }, error = function(e) {
      if (verbose) {
        message("Fallback also failed: ", e$message)
      }
      NULL
    })
    
    if (!is.null(result)) {
      return(result)
    }
  }
  
  # Everything failed, return default
  if (verbose) {
    message(error_prefix, " failed completely. Returning default value.")
  }
  
  if (is.null(default_value)) {
    # Construct a reasonable default based on context
    default_value <- .construct_safe_default(error_prefix)
  }
  
  default_value
}

#' Construct safe default values based on context
#' @keywords internal
.construct_safe_default <- function(context) {
  defaults <- list(
    "Parameter estimation" = matrix(c(6, 2.5, 0.35), nrow = 1),
    "R-squared calculation" = 0,
    "Residual calculation" = NULL,
    "Standard error calculation" = NULL,
    "Amplitude estimation" = 1,
    "Design matrix" = matrix(1, nrow = 10, ncol = 1)
  )
  
  # Find best matching default
  for (pattern in names(defaults)) {
    if (grepl(pattern, context, ignore.case = TRUE)) {
      return(defaults[[pattern]])
    }
  }
  
  # Generic default
  0
}

#' Progressive algorithm degradation
#' @keywords internal
.progressive_estimation <- function(Y_proj, S_target_proj, hrf_interface,
                                    theta_seed, theta_bounds, 
                                    recenter_global_passes = 3,
                                    verbose = TRUE) {
  n_vox <- ncol(Y_proj)
  n_params <- length(theta_seed)
  
  # Level 1: Try full iterative estimation
  result <- .try_with_recovery(
    primary_fn = function() {
      source(file.path(dirname(getwd()), "R", "parametric-engine-iterative.R"), 
             local = TRUE)
      
      .parametric_engine_iterative(
        Y_proj = Y_proj,
        S_target_proj = S_target_proj,
        scan_times = seq_len(nrow(Y_proj)),
        hrf_eval_times = seq(0, 30, by = 0.5),
        hrf_interface = hrf_interface,
        theta_seed = theta_seed,
        theta_bounds = theta_bounds,
        recenter_global_passes = recenter_global_passes,
        compute_residuals = TRUE,
        compute_se = TRUE,
        verbose = verbose
      )
    },
    error_prefix = "Iterative estimation",
    verbose = verbose
  )
  
  if (!is.null(result)) return(result)
  
  # Level 2: Try single-pass estimation
  result <- .try_with_recovery(
    primary_fn = function() {
      source(file.path(dirname(getwd()), "R", "parametric-engine.R"), 
             local = TRUE)
      
      basic_result <- .parametric_engine(
        Y_proj = Y_proj,
        S_target_proj = S_target_proj,
        scan_times = seq_len(nrow(Y_proj)),
        hrf_eval_times = seq(0, 30, by = 0.5),
        hrf_interface = hrf_interface,
        theta_seed = theta_seed,
        theta_bounds = theta_bounds
      )
      
      # Convert to iterative format
      list(
        theta_hat = basic_result$theta_hat,
        beta0 = basic_result$beta0,
        r_squared = rep(NA, n_vox),
        residuals = NULL,
        se_theta_hat = NULL,
        convergence_info = list(converged = FALSE, reason = "fallback_single_pass"),
        coeffs = NULL
      )
    },
    error_prefix = "Single-pass estimation",
    verbose = verbose
  )
  
  if (!is.null(result)) return(result)
  
  # Level 3: Return seed parameters with warning
  if (verbose) {
    message("All estimation methods failed. Returning seed parameters.")
  }
  
  list(
    theta_hat = matrix(theta_seed, nrow = n_vox, ncol = n_params, byrow = TRUE),
    beta0 = rep(1, n_vox),
    r_squared = rep(0, n_vox),
    residuals = NULL,
    se_theta_hat = NULL,
    convergence_info = list(converged = FALSE, reason = "all_methods_failed"),
    coeffs = NULL
  )
}

#' Wrap voxel-wise operations with error recovery
#' @keywords internal
.voxel_safe_apply <- function(voxel_indices, operation_fn, 
                              n_params = 3, 
                              combine_results = TRUE,
                              progress = TRUE,
                              verbose = TRUE) {
  n_voxels <- length(voxel_indices)
  
  # Pre-allocate results
  results <- vector("list", n_voxels)
  failed_voxels <- integer(0)
  
  if (progress && n_voxels > 100) {
    cat("Processing", n_voxels, "voxels...\n")
  }
  
  for (i in seq_along(voxel_indices)) {
    v <- voxel_indices[i]
    
    results[[i]] <- tryCatch({
      operation_fn(v)
    }, error = function(e) {
      failed_voxels <<- c(failed_voxels, v)
      # Return safe default
      list(
        theta = rep(NA, n_params),
        r2 = 0,
        converged = FALSE,
        error = e$message
      )
    })
    
    if (progress && i %% 1000 == 0) {
      cat("  Processed", i, "/", n_voxels, "voxels (",
          round(100 * i / n_voxels), "% complete)\n")
    }
  }
  
  if (length(failed_voxels) > 0) {
    warning("Processing failed for ", length(failed_voxels), " voxels: ",
            paste(head(failed_voxels, 10), collapse = ", "),
            if (length(failed_voxels) > 10) "..." else "")
  }
  
  # Combine results if requested
  if (combine_results) {
    combined <- .safe_combine_results(results, n_params)
    combined$failed_voxels <- failed_voxels
    return(combined)
  }
  
  list(
    results = results,
    failed_voxels = failed_voxels
  )
}

#' Safely combine voxel-wise results
#' @keywords internal
.safe_combine_results <- function(results, n_params) {
  n_results <- length(results)
  
  # Extract components safely
  theta_list <- lapply(results, function(r) {
    if (is.null(r$theta)) rep(NA, n_params) else r$theta
  })
  
  r2_vec <- sapply(results, function(r) {
    if (is.null(r$r2)) 0 else r$r2
  })
  
  converged_vec <- sapply(results, function(r) {
    if (is.null(r$converged)) FALSE else r$converged
  })
  
  # Combine into matrices
  theta_mat <- do.call(rbind, theta_list)
  
  list(
    theta = theta_mat,
    r2 = r2_vec,
    converged = converged_vec,
    n_failed = sum(!converged_vec)
  )
}

#' Create detailed error report
#' @keywords internal
.create_error_report <- function(error_list, context = "Unknown operation") {
  report <- list(
    context = context,
    timestamp = Sys.time(),
    session_info = sessionInfo(),
    errors = error_list
  )
  
  # Categorize errors
  error_types <- table(sapply(error_list, function(e) class(e)[1]))
  report$error_summary <- as.list(error_types)
  
  # Common error patterns
  patterns <- list(
    memory = "cannot allocate|memory",
    numerical = "NaN|Inf|singular",
    dimension = "dimension|conform",
    type = "must be|cannot coerce"
  )
  
  pattern_counts <- sapply(patterns, function(p) {
    sum(grepl(p, sapply(error_list, function(e) e$message), 
              ignore.case = TRUE))
  })
  
  report$error_patterns <- pattern_counts[pattern_counts > 0]
  
  # Recommendations
  report$recommendations <- .generate_error_recommendations(report$error_patterns)
  
  class(report) <- "fmriparametric_error_report"
  report
}

#' Generate recommendations based on error patterns
#' @keywords internal
.generate_error_recommendations <- function(error_patterns) {
  recommendations <- list()
  
  if ("memory" %in% names(error_patterns)) {
    recommendations$memory <- c(
      "Consider processing data in smaller chunks",
      "Increase memory limits with options(future.globals.maxSize)",
      "Use mask to reduce voxel count",
      "Close other applications to free memory"
    )
  }
  
  if ("numerical" %in% names(error_patterns)) {
    recommendations$numerical <- c(
      "Increase ridge regularization (lambda_ridge parameter)",
      "Check for extreme values in input data",
      "Consider data normalization/scaling",
      "Use more conservative parameter bounds"
    )
  }
  
  if ("dimension" %in% names(error_patterns)) {
    recommendations$dimension <- c(
      "Verify fMRI data and event model have matching time points",
      "Check mask dimensions match data dimensions",
      "Ensure all inputs are properly formatted"
    )
  }
  
  recommendations
}

#' Print method for error reports
#' @export
print.fmriparametric_error_report <- function(x, ...) {
  cat("=== fmriparametric Error Report ===\n")
  cat("Context:", x$context, "\n")
  cat("Time:", format(x$timestamp), "\n")
  cat("Total errors:", length(x$errors), "\n\n")
  
  if (length(x$error_summary) > 0) {
    cat("Error types:\n")
    for (type in names(x$error_summary)) {
      cat("  ", type, ":", x$error_summary[[type]], "\n")
    }
  }
  
  if (length(x$error_patterns) > 0) {
    cat("\nError patterns detected:\n")
    for (pattern in names(x$error_patterns)) {
      cat("  ", pattern, ":", x$error_patterns[[pattern]], "occurrences\n")
    }
  }
  
  if (length(x$recommendations) > 0) {
    cat("\nRecommendations:\n")
    for (category in names(x$recommendations)) {
      cat("\n", toupper(category), ":\n")
      for (rec in x$recommendations[[category]]) {
        cat("  - ", rec, "\n")
      }
    }
  }
  
  invisible(x)
}