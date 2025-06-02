# Internal utility helpers
# Consolidated versions of .try_with_context(), .get_available_memory(),
# and .check_memory_available() stored in a dedicated environment to
# avoid accidental masking.

.fmriparametric_internal <- new.env(parent = emptyenv())

.fmriparametric_internal$try_with_context <- function(expr, context = "", fallback = NULL) {
  tryCatch(
    expr,
    error = function(e) {
      enhanced_msg <- sprintf(
        "Error in %s:\n  %s\n  Call stack: %s",
        context,
        e$message,
        paste(deparse(sys.calls()), collapse = " > ")
      )

      if (!is.null(fallback)) {
        warning(enhanced_msg, "\n  Attempting fallback...")
        tryCatch(
          fallback,
          error = function(e2) {
            stop(sprintf(
              "%s\n  Fallback also failed: %s",
              enhanced_msg, e2$message
            ), call. = FALSE)
          }
        )
      } else {
        stop(enhanced_msg, call. = FALSE)
      }
    }
  )
}

.fmriparametric_internal$get_available_memory <- function() {
  if (.Platform$OS.type == "windows") {
    memory.limit() * 1024^2
  } else if (Sys.info()["sysname"] == "Darwin") {
    # macOS
    tryCatch({
      vm_stat <- system("vm_stat", intern = TRUE)
      free_line <- grep("Pages free:", vm_stat, value = TRUE)
      inactive_line <- grep("Pages inactive:", vm_stat, value = TRUE)
      
      free_pages <- as.numeric(gsub(".*:\\s*(\\d+).*", "\\1", free_line))
      inactive_pages <- as.numeric(gsub(".*:\\s*(\\d+).*", "\\1", inactive_line))
      
      # Page size is typically 4096 bytes on macOS
      (free_pages + inactive_pages) * 4096
    }, error = function(e) {
      4e9  # Default 4GB
    })
  } else {
    # Linux
    tryCatch({
      as.numeric(system("awk '/MemAvailable/ {print $2}' /proc/meminfo", intern = TRUE)) * 1024
    }, error = function(e) {
      4e9  # Default 4GB
    })
  }
}

.fmriparametric_internal$check_memory_available <- function(required_bytes, operation = "") {
  available <- .fmriparametric_internal$get_available_memory()
  
  # Handle case where available memory cannot be determined
  if (is.na(available) || is.null(available) || !is.numeric(available)) {
    # Can't determine memory, assume it's OK
    return(TRUE)
  }

  if (required_bytes > available * 0.8) {
    warning(sprintf(
      "Operation '%s' requires %.1f GB but only %.1f GB available",
      operation,
      required_bytes / 1e9,
      available / 1e9
    ))

    message("Consider:")
    message("  - Using smaller chunks")
    message("  - Increasing memory limit: memory.limit(size = ...)")
    message("  - Closing other applications")

    return(FALSE)
  }

  TRUE
}

# Backwards compatible wrappers
.try_with_context <- function(...) .fmriparametric_internal$try_with_context(...)
.get_available_memory <- function(...) .fmriparametric_internal$get_available_memory(...)
.check_memory_available <- function(...) .fmriparametric_internal$check_memory_available(...)

# ---------------------------------------------------------------------------
# Timing Utilities -----------------------------------------------------------

#' Time and profile a code block
#'
#' Executes `expr` while optionally recording and reporting the elapsed time.
#' When `fmriparametric.verbose` is TRUE and a `label` is provided, the timing
#' information is stored for later inspection via `get_timing_report()`.
#'
#' @param expr  Expression to evaluate
#' @param label Optional character label for the operation
#' @return Result of evaluating `expr`
#' @keywords internal
.with_timing <- function(expr, label = NULL) {
  if (getOption("fmriparametric.verbose", FALSE)) {
    start_time <- Sys.time()

    result <- force(expr)

    elapsed <- as.numeric(Sys.time() - start_time, units = "secs")

    if (!is.null(label)) {
      message(sprintf("[%s] Elapsed time: %.3f seconds", label, elapsed))
      .record_timing(label, elapsed)
    }

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
    return(data.frame(
      operation = character(0),
      count = numeric(0),
      total_time = numeric(0),
      mean_time = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  timing_list <- as.list(timings_env)
  timing_df <- data.frame(
    operation = names(timing_list),
    count = sapply(timing_list, `[[`, "count"),
    total_time = sapply(timing_list, `[[`, "total"),
    mean_time = sapply(timing_list, `[[`, "mean"),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  timing_df[order(timing_df$total_time, decreasing = TRUE), ]
}

# ---------------------------------------------------------------------------
# Numerical Helpers ---------------------------------------------------------

#' Safe division avoiding numerical issues
#'
#' @param numerator Numeric vector/matrix
#' @param denominator Numeric vector/matrix
#' @param epsilon Small value to prevent division by zero
#' @return Result of division with safety checks
#' @keywords internal
.safe_divide <- function(numerator, denominator, epsilon = .Machine$double.eps) {
  inf_num <- is.infinite(numerator)
  inf_denom <- is.infinite(denominator)
  small_denom <- abs(denominator) < epsilon

  inf_inf_case <- inf_num & inf_denom
  if (any(inf_inf_case, na.rm = TRUE)) {
    warning("Division produced non-finite values")
  }

  if (any(small_denom, na.rm = TRUE)) {
    warning("Near-zero denominator detected, applying epsilon correction")
    denominator[small_denom] <- sign(denominator[small_denom]) * epsilon
    denominator[denominator == 0] <- epsilon
  }

  result <- numerator / denominator

  result[inf_inf_case] <- 0

  remaining_invalid <- !is.finite(result) & !inf_inf_case
  if (any(remaining_invalid, na.rm = TRUE)) {
    n_invalid <- sum(remaining_invalid, na.rm = TRUE)
    warning(sprintf(
      "Division produced %d non-finite values, setting to zero",
      n_invalid
    ))
    result[remaining_invalid] <- 0
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

  if (method == "auto") {
    if (n == p && all(X == t(X))) {
      method <- "chol"
    } else if (n >= p) {
      method <- "qr"
    } else {
      method <- "svd"
    }
  }

  condition <- kappa(X)
  if (condition > 1e10) {
    warning(sprintf(
      "Matrix is poorly conditioned (kappa = %.2e)",
      condition
    ))

    if (regularize) {
      lambda <- condition * .Machine$double.eps * max(diag(crossprod(X)))
      X <- X + lambda * diag(n)
      message(sprintf("Applied regularization lambda = %.2e", lambda))
    }
  }

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

  d_inv <- ifelse(d > tol * max(d), 1 / d, 0)

  svd_decomp$v %*% (d_inv * t(svd_decomp$u))
}

#' Assert output quality
#'
#' Performs a series of simple checks on a result object to ensure it does not
#' contain obviously problematic values.
#'
#' @param result Result to check
#' @param checks List of checks to perform
#' @keywords internal
.assert_output_quality <- function(result, checks = list()) {
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

