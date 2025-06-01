#' Engineering Standards for fmriparametric
#'
#' This file defines the engineering standards and utilities that ensure
#' code quality throughout the package.

# Validation Standards --------------------------------------------------------

#' Validate function inputs with standardized error messages
#'
#' @param x Object to validate
#' @param arg_name Character name of argument for error messages
#' @param type Expected type(s)
#' @param dims Expected dimensions (NULL to skip)
#' @param constraints List of constraints (range, properties, etc.)
#' @keywords internal
.validate_input <- function(x, arg_name, type = NULL, dims = NULL, 
                           constraints = NULL, null_ok = FALSE) {
  # NULL handling
  if (is.null(x)) {
    if (null_ok) return(invisible(TRUE))
    stop(sprintf("Argument '%s' cannot be NULL", arg_name), call. = FALSE)
  }
  
  # Type validation
  if (!is.null(type)) {
    if (!inherits(x, type)) {
      stop(sprintf(
        "Argument '%s' must be of type %s, got %s",
        arg_name, 
        paste(type, collapse = " or "),
        paste(class(x), collapse = ", ")
      ), call. = FALSE)
    }
  }
  
  # Dimension validation
  if (!is.null(dims)) {
    actual_dims <- if (is.matrix(x) || is.array(x)) dim(x) else length(x)
    
    if (length(dims) != length(actual_dims) || !all(dims == actual_dims | dims == -1)) {
      stop(sprintf(
        "Argument '%s' must have dimensions %s, got %s",
        arg_name,
        paste(dims, collapse = " x "),
        paste(actual_dims, collapse = " x ")
      ), call. = FALSE)
    }
  }
  
  # Constraint validation
  if (!is.null(constraints)) {
    .validate_constraints(x, arg_name, constraints)
  }
  
  invisible(TRUE)
}

#' Validate numerical constraints
#' @keywords internal
.validate_constraints <- function(x, arg_name, constraints) {
  # Range constraints
  if (!is.null(constraints$range)) {
    if (any(x < constraints$range[1] | x > constraints$range[2], na.rm = TRUE)) {
      stop(sprintf(
        "Argument '%s' contains values outside range [%g, %g]",
        arg_name, constraints$range[1], constraints$range[2]
      ), call. = FALSE)
    }
  }
  
  # Finite values
  if (isTRUE(constraints$finite)) {
    if (any(!is.finite(x))) {
      stop(sprintf(
        "Argument '%s' must contain only finite values",
        arg_name
      ), call. = FALSE)
    }
  }
  
  # Positive values
  if (isTRUE(constraints$positive)) {
    if (any(x <= 0, na.rm = TRUE)) {
      stop(sprintf(
        "Argument '%s' must contain only positive values",
        arg_name
      ), call. = FALSE)
    }
  }
}

# Performance Utilities -------------------------------------------------------

# Global timing storage
.timing_data <- new.env(parent = emptyenv())

#' Time operation execution
#'
#' @param expr Expression to time
#' @param name Operation name
#' @param verbose Whether to print timing
#' @return Result of expression
#' @keywords internal
.with_timing <- function(expr, name = "operation", verbose = getOption("fmriparametric.verbose", FALSE)) {
  start_time <- Sys.time()
  result <- expr
  end_time <- Sys.time()
  
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  if (verbose) {
    cat(sprintf("[%s] %s: %.3f seconds\n", format(Sys.time(), "%H:%M:%S"), name, elapsed))
  }
  
  # Record timing
  .record_timing(name, elapsed)
  
  result
}

#' Record timing for an operation
#'
#' @param operation Character name of operation
#' @param duration Numeric duration in seconds
#' @keywords internal
.record_timing <- function(operation, duration) {
  if (!exists("timings", envir = .timing_data)) {
    .timing_data$timings <- data.frame(
      operation = character(0),
      duration = numeric(0),
      timestamp = as.POSIXct(character(0))
    )
  }
  
  .timing_data$timings <- rbind(.timing_data$timings, data.frame(
    operation = operation,
    duration = duration,
    timestamp = Sys.time()
  ))
}

#' Get timing report
#'
#' @return Data frame with timing information
#' @export
get_timing_report <- function() {
  if (!exists("timings", envir = .timing_data)) {
    return(data.frame(
      operation = character(0),
      duration = numeric(0),
      timestamp = as.POSIXct(character(0))
    ))
  }
  
  .timing_data$timings
}


# Numerical Robustness --------------------------------------------------------

#' Safe division avoiding numerical issues
#'
#' @param numerator Numeric vector/matrix
#' @param denominator Numeric vector/matrix  
#' @param epsilon Small value to prevent division by zero
#' @return Result of division with safety checks
#' @keywords internal
.safe_divide <- function(numerator, denominator, epsilon = .Machine$double.eps) {
  # Check for zero denominator
  small_denom <- abs(denominator) < epsilon
  
  if (any(small_denom, na.rm = TRUE)) {
    warning("Near-zero denominator detected, applying epsilon correction")
    denominator[small_denom] <- sign(denominator[small_denom]) * epsilon
    denominator[denominator == 0] <- epsilon
  }
  
  result <- numerator / denominator
  
  # Check for non-finite results (including NaN from Inf/Inf)
  if (any(!is.finite(result), na.rm = TRUE)) {
    n_invalid <- sum(!is.finite(result), na.rm = TRUE)
    warning(sprintf(
      "Division produced %d non-finite values, setting to zero",
      n_invalid
    ))
    result[!is.finite(result)] <- 0
  }
  
  result
}

#' Numerically stable matrix inversion
#'
#' @param X Matrix to invert
#' @param method Method to use: "qr", "svd", "chol"
#' @param tol Tolerance for singular values
#' @param regularize Whether to add ridge regularization if needed
#' @return Inverse or pseudo-inverse of X
#' @keywords internal
.safe_solve <- function(X, method = "auto", tol = .Machine$double.eps^0.5,
                       regularize = TRUE) {
  n <- nrow(X)
  p <- ncol(X)
  
  # Auto-select method
  if (method == "auto") {
    if (n == p && all(X == t(X))) {
      method <- "chol"  # Symmetric
    } else if (n >= p) {
      method <- "qr"    # Overdetermined
    } else {
      method <- "svd"   # Underdetermined
    }
  }
  
  # Check condition number
  condition <- kappa(X)
  if (condition > 1e10) {
    warning(sprintf(
      "Matrix is poorly conditioned (kappa = %.2e)",
      condition
    ))
    
    if (regularize) {
      # Adaptive regularization
      lambda <- tol * max(diag(crossprod(X)))
      X <- X + lambda * diag(n)
      message(sprintf("Applied regularization lambda = %.2e", lambda))
    }
  }
  
  # Method-specific inversion
  result <- switch(method,
    qr = {
      qr_decomp <- qr(X)
      if (qr_decomp$rank < p) {
        warning("Matrix is rank-deficient, using pseudo-inverse")
        .svd_pinv(X, tol)
      } else {
        solve.qr(qr_decomp)
      }
    },
    
    chol = {
      tryCatch({
        chol_decomp <- chol(X)
        chol2inv(chol_decomp)
      }, error = function(e) {
        warning("Cholesky failed, falling back to SVD: ", e$message)
        .svd_pinv(X, tol)
      })
    },
    
    svd = .svd_pinv(X, tol),
    
    stop("Unknown method: ", method)
  )
  
  # Verify result
  if (any(!is.finite(result))) {
    stop("Matrix inversion produced non-finite values")
  }
  
  result
}

#' SVD-based pseudo-inverse
#' @keywords internal
.svd_pinv <- function(X, tol = .Machine$double.eps^0.5) {
  svd_decomp <- svd(X)
  d <- svd_decomp$d
  
  # Threshold small singular values
  d_inv <- ifelse(d > tol * max(d), 1/d, 0)
  
  # Reconstruct
  svd_decomp$v %*% (d_inv * t(svd_decomp$u))
}

# Performance Utilities -------------------------------------------------------

#' Time and profile a code block
#'
#' @param expr Expression to evaluate
#' @param label Character label for profiling
#' @return Result of expression
#' @keywords internal
.with_timing <- function(expr, label = NULL) {
  if (getOption("fmriparametric.verbose", FALSE)) {
    start_time <- Sys.time()
    
    result <- force(expr)
    
    elapsed <- as.numeric(Sys.time() - start_time, units = "secs")
    
    if (!is.null(label)) {
      message(sprintf("[%s] Elapsed time: %.3f seconds", label, elapsed))
    }
    
    # Record in global profile
    .record_timing(label, elapsed)
    
    result
  } else {
    expr
  }
}

#' Record timing information
#' @keywords internal
.record_timing <- local({
  timings <- new.env(parent = emptyenv())
  
  function(label, elapsed) {
    if (exists(label, envir = timings)) {
      current <- get(label, envir = timings)
      current$count <- current$count + 1
      current$total <- current$total + elapsed
      current$mean <- current$total / current$count
      assign(label, current, envir = timings)
    } else {
      assign(label, list(
        count = 1,
        total = elapsed,
        mean = elapsed
      ), envir = timings)
    }
  }
})

#' Get timing report
#' @export
get_timing_report <- function() {
  timings_env <- environment(.record_timing)$timings
  if (length(ls(timings_env)) == 0) {
    message("No timing data recorded")
    return(invisible(NULL))
  }
  
  # Convert to data frame
  timing_list <- as.list(timings_env)
  timing_df <- data.frame(
    operation = names(timing_list),
    count = sapply(timing_list, `[[`, "count"),
    total_time = sapply(timing_list, `[[`, "total"),
    mean_time = sapply(timing_list, `[[`, "mean"),
    row.names = NULL
  )
  
  # Sort by total time
  timing_df[order(timing_df$total_time, decreasing = TRUE), ]
}

# Error Handling --------------------------------------------------------------

#  Consolidated utility functions are defined in `internal-utils.R`.


# Quality Assertions ----------------------------------------------------------

#' Assert output quality
#'
#' @param result Result to check
#' @param checks List of checks to perform
#' @keywords internal
.assert_output_quality <- function(result, checks = list()) {
  # Default checks
  if ("finite" %in% names(checks) && isTRUE(checks$finite)) {
    if (any(!is.finite(unlist(result)))) {
      stop("Output contains non-finite values")
    }
  }
  
  if ("positive_r2" %in% names(checks) && isTRUE(checks$positive_r2)) {
    if (any(result$r_squared < 0 | result$r_squared > 1, na.rm = TRUE)) {
      stop("R-squared values outside [0, 1] range")
    }
  }
  
  if ("bounded_params" %in% names(checks) && !is.null(checks$bounded_params)) {
    bounds <- checks$bounded_params
    params <- result$parameters
    
    for (i in seq_len(ncol(params))) {
      if (any(params[, i] < bounds$lower[i] | params[, i] > bounds$upper[i])) {
        stop(sprintf(
          "Parameter %d outside bounds [%g, %g]",
          i, bounds$lower[i], bounds$upper[i]
        ))
      }
    }
  }
  
  invisible(TRUE)
}

# Engineering Standards Options -----------------------------------------------

#' Set engineering standards options
#' @export
set_engineering_options <- function(
  verbose = FALSE,
  validate = TRUE,
  profile = FALSE,
  precision = "double",
  debug = FALSE
) {
  options(
    fmriparametric.verbose = verbose,
    fmriparametric.validate = validate,
    fmriparametric.profile = profile,
    fmriparametric.precision = precision,
    fmriparametric.debug = debug
  )
  
  if (debug) {
    message("Engineering standards: DEBUG mode enabled")
    message("  - All validations active")
    message("  - Profiling enabled")
    message("  - Verbose output")
    options(
      fmriparametric.verbose = TRUE,
      fmriparametric.validate = TRUE,
      fmriparametric.profile = TRUE
    )
  }
  
  invisible(NULL)
}