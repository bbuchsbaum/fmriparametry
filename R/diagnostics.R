#' Diagnostics and messaging utilities
#'
#' These helper functions centralize error, warning, and information
#' messages as well as optional progress reporting. Verbosity can be
#' controlled via the option `fmriparametric.verbose`.
#'
#' @name diagnostics
#' @keywords internal
NULL

#' Emit an informational message respecting package verbosity
#' @keywords internal
.diag_inform <- function(msg, ...) {
  if (isTRUE(getOption("fmriparametric.verbose", TRUE))) {
    rlang::inform(sprintf(msg, ...))
  }
}

#' Emit a warning message using rlang
#' @keywords internal
.diag_warn <- function(msg, ...) {
  rlang::warn(sprintf(msg, ...))
}

#' Emit an error message using rlang
#' @keywords internal
.diag_abort <- function(msg, ..., caller = NULL) {
  # Only use sprintf if there are additional arguments
  dots <- list(...)
  if (length(dots) > 0 && !("caller" %in% names(dots))) {
    msg <- sprintf(msg, ...)
  }
  rlang::abort(msg, call = if (!is.null(caller)) rlang::caller_env() else NULL)
}

#' Create a progressr progressor if available
#' @keywords internal
.diag_progressor <- function(steps) {
  if (requireNamespace("progressr", quietly = TRUE)) {
    progressr::progressor(steps)
  } else {
    function(...) NULL
  }
}

# Memory Usage Tracking -------------------------------------------------------

#' Track memory usage for an operation
#'
#' Records memory usage before and after an operation, optionally
#' logging the information.
#'
#' @param expr Expression to evaluate
#' @param label Label for the operation
#' @param log Logical whether to log the usage
#' @return Result of evaluating expr
#' @keywords internal
.diag_memory_track <- function(expr, label = "Operation", log = TRUE) {
  # Get initial memory state
  gc(verbose = FALSE, reset = TRUE)
  mem_before <- gc()[, "used"]
  
  # Execute the operation
  result <- force(expr)
  
  # Get final memory state
  gc(verbose = FALSE)
  mem_after <- gc()[, "used"]
  
  # Calculate usage - gc() returns Ncells and Vcells
  # Vcells are 8 bytes each on 64-bit systems
  mem_delta <- mem_after - mem_before
  # Convert Vcells (second row) to MB: Vcells * 8 bytes / 1024^2
  mem_mb <- round(mem_delta[2] * 8 / (1024^2), 2)
  
  if (log && getOption("fmriparametric.verbose", TRUE)) {
    .diag_inform("[Memory] %s: %.1f MB", label, mem_mb)
  }
  
  # Store in diagnostics environment
  .record_memory_usage(label, mem_mb)
  
  result
}

#' Record memory usage in diagnostics environment
#' @keywords internal
.record_memory_usage <- local({
  memory_log <- new.env(parent = emptyenv())
  
  function(label, mb_used) {
    if (!exists(label, envir = memory_log)) {
      memory_log[[label]] <- list(
        count = 0,
        total_mb = 0,
        max_mb = 0,
        measurements = numeric()
      )
    }
    
    entry <- memory_log[[label]]
    entry$count <- entry$count + 1
    entry$total_mb <- entry$total_mb + mb_used
    entry$max_mb <- max(entry$max_mb, mb_used)
    entry$measurements <- c(entry$measurements, mb_used)
    
    memory_log[[label]] <- entry
  }
})

#' Get memory usage report
#' @export
get_memory_report <- function() {
  memory_env <- environment(.record_memory_usage)$memory_log
  if (length(ls(memory_env)) == 0) {
    return(data.frame(
      operation = character(0),
      count = numeric(0),
      total_mb = numeric(0),
      mean_mb = numeric(0),
      max_mb = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  
  memory_list <- as.list(memory_env)
  memory_df <- data.frame(
    operation = names(memory_list),
    count = sapply(memory_list, `[[`, "count"),
    total_mb = sapply(memory_list, `[[`, "total_mb"),
    mean_mb = sapply(memory_list, function(x) x$total_mb / x$count),
    max_mb = sapply(memory_list, `[[`, "max_mb"),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  memory_df[order(memory_df$total_mb, decreasing = TRUE), ]
}

# Convergence Monitoring ------------------------------------------------------

#' Monitor convergence of iterative algorithms
#'
#' Creates a convergence monitor that tracks objective values,
#' parameter changes, and convergence criteria.
#'
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @param verbose Whether to print progress
#' @return List with monitor functions
#' @keywords internal
.diag_convergence_monitor <- function(max_iter = 100, tol = 1e-6, verbose = TRUE) {
  # Storage for convergence history
  history <- list(
    objective = numeric(),
    param_change = numeric(),
    gradient_norm = numeric(),
    converged = FALSE,
    iter = 0
  )
  
  list(
    update = function(objective, param_change = NULL, gradient_norm = NULL) {
      history$iter <<- history$iter + 1
      history$objective <<- c(history$objective, objective)
      
      if (!is.null(param_change)) {
        history$param_change <<- c(history$param_change, param_change)
      }
      
      if (!is.null(gradient_norm)) {
        history$gradient_norm <<- c(history$gradient_norm, gradient_norm)
      }
      
      # Check convergence
      if (history$iter > 1) {
        obj_change <- abs(objective - history$objective[history$iter - 1])
        if (obj_change < tol) {
          history$converged <<- TRUE
          if (verbose) {
            .diag_inform("[Convergence] Achieved at iteration %d (objective change: %.2e)",
                        history$iter, obj_change)
          }
        }
      }
      
      # Check iteration limit
      if (history$iter >= max_iter) {
        history$converged <<- TRUE
        if (verbose) {
          .diag_warn("[Convergence] Maximum iterations (%d) reached", max_iter)
        }
      }
      
      history$converged
    },
    
    get_history = function() {
      history
    },
    
    plot = function() {
      if (length(history$objective) > 0) {
        par(mfrow = c(2, 2))
        
        # Objective
        plot(history$objective, type = "l", 
             xlab = "Iteration", ylab = "Objective",
             main = "Objective Function")
        
        # Log objective (for better visualization)
        if (all(history$objective > 0)) {
          plot(log(history$objective), type = "l",
               xlab = "Iteration", ylab = "Log Objective",
               main = "Log Objective Function")
        }
        
        # Parameter change
        if (length(history$param_change) > 0) {
          plot(history$param_change, type = "l",
               xlab = "Iteration", ylab = "Parameter Change",
               main = "Parameter Change Norm")
        }
        
        # Gradient norm
        if (length(history$gradient_norm) > 0) {
          plot(history$gradient_norm, type = "l",
               xlab = "Iteration", ylab = "Gradient Norm",
               main = "Gradient Norm")
        }
        
        par(mfrow = c(1, 1))
      }
    }
  )
}

# Numerical Stability Checks --------------------------------------------------

#' Check numerical stability of a matrix
#'
#' Performs comprehensive stability checks including condition number,
#' rank deficiency, and near-singularity detection.
#'
#' @param X Matrix to check
#' @param operation Description of the operation
#' @param tol Tolerance for rank detection
#' @return Invisible NULL, issues warnings as needed
#' @keywords internal
.diag_check_stability <- function(X, operation = "Matrix operation", 
                                 tol = .Machine$double.eps^0.5) {
  # Condition number
  cond <- kappa(X)
  if (cond > 1e10) {
    .diag_warn("[Stability] %s: High condition number (%.2e)", operation, cond)
  }
  
  # Rank check
  qr_decomp <- qr(X)
  rank_def <- qr_decomp$rank < min(nrow(X), ncol(X))
  if (rank_def) {
    .diag_warn("[Stability] %s: Rank deficient (rank %d < %d)", 
              operation, qr_decomp$rank, min(nrow(X), ncol(X)))
  }
  
  # Check for near-zero singular values
  if (min(dim(X)) <= 100) {  # Only for smaller matrices
    svd_decomp <- svd(X)
    small_sv <- sum(svd_decomp$d < tol * max(svd_decomp$d))
    if (small_sv > 0) {
      .diag_warn("[Stability] %s: %d near-zero singular values", 
                operation, small_sv)
    }
  }
  
  # Check for NaN/Inf
  if (any(!is.finite(X))) {
    .diag_abort("[Stability] %s: Matrix contains non-finite values", operation)
  }
  
  invisible(NULL)
}

# Stage-wise Performance Profiling --------------------------------------------

#' Profile performance across estimation stages
#'
#' Records detailed timing and resource usage for each stage
#' of the estimation pipeline.
#'
#' @return Environment with profiling functions
#' @keywords internal
.diag_stage_profiler <- function() {
  profile_data <- new.env(parent = emptyenv())
  profile_data$stages <- list()
  
  list(
    start_stage = function(stage_name) {
      profile_data$stages[[stage_name]] <- list(
        start_time = Sys.time(),
        start_memory = gc()[2, "used"],
        substages = list()
      )
    },
    
    end_stage = function(stage_name, metrics = NULL) {
      if (!stage_name %in% names(profile_data$stages)) {
        warning("Stage ", stage_name, " was not started")
        return(invisible(NULL))
      }
      
      stage <- profile_data$stages[[stage_name]]
      stage$end_time <- Sys.time()
      stage$end_memory <- gc()[2, "used"]
      stage$elapsed <- as.numeric(stage$end_time - stage$start_time, units = "secs")
      stage$memory_mb <- stage$end_memory - stage$start_memory
      
      if (!is.null(metrics)) {
        stage$metrics <- metrics
      }
      
      profile_data$stages[[stage_name]] <- stage
      
      if (getOption("fmriparametric.verbose", TRUE)) {
        .diag_inform("[Profile] Stage '%s': %.2f sec, %.1f MB",
                    stage_name, stage$elapsed, stage$memory_mb)
      }
    },
    
    add_substage = function(stage_name, substage_name, elapsed, metrics = NULL) {
      if (!stage_name %in% names(profile_data$stages)) {
        return(invisible(NULL))
      }
      
      profile_data$stages[[stage_name]]$substages[[substage_name]] <- list(
        elapsed = elapsed,
        metrics = metrics
      )
    },
    
    get_report = function() {
      stages <- profile_data$stages
      if (length(stages) == 0) {
        return(NULL)
      }
      
      # Summary data frame
      stage_df <- data.frame(
        stage = names(stages),
        elapsed_sec = sapply(stages, function(s) s$elapsed %||% NA),
        memory_mb = sapply(stages, function(s) s$memory_mb %||% NA),
        stringsAsFactors = FALSE
      )
      
      list(
        summary = stage_df[order(stage_df$elapsed_sec, decreasing = TRUE), ],
        details = stages
      )
    }
  )
}

# Parameter Identifiability Diagnostics ---------------------------------------

#' Check parameter identifiability
#'
#' Analyzes the Fisher information matrix and parameter correlations
#' to assess identifiability issues.
#'
#' @param hessian Hessian matrix
#' @param param_names Parameter names
#' @param threshold Correlation threshold for warnings
#' @return List with identifiability diagnostics
#' @keywords internal
.diag_check_identifiability <- function(hessian, param_names = NULL, 
                                       threshold = 0.95) {
  n_params <- nrow(hessian)
  if (is.null(param_names)) {
    param_names <- paste0("param", seq_len(n_params))
  }
  
  # Try to invert for Fisher information
  fisher_info <- tryCatch({
    solve(hessian)
  }, error = function(e) {
    .diag_warn("[Identifiability] Cannot compute Fisher information: %s", e$message)
    NULL
  })
  
  if (is.null(fisher_info)) {
    return(list(
      identifiable = FALSE,
      correlations = NULL,
      standard_errors = rep(NA, n_params)
    ))
  }
  
  # Standard errors
  se <- sqrt(diag(fisher_info))
  
  # Parameter correlations
  D <- diag(1 / se)
  corr_matrix <- D %*% fisher_info %*% D
  
  # Check for high correlations
  high_corr <- which(abs(corr_matrix) > threshold & 
                    upper.tri(corr_matrix), arr.ind = TRUE)
  
  if (nrow(high_corr) > 0) {
    for (i in seq_len(nrow(high_corr))) {
      .diag_warn("[Identifiability] High correlation (%.3f) between %s and %s",
                abs(corr_matrix[high_corr[i, 1], high_corr[i, 2]]),
                param_names[high_corr[i, 1]],
                param_names[high_corr[i, 2]])
    }
  }
  
  list(
    identifiable = nrow(high_corr) == 0,
    correlations = corr_matrix,
    standard_errors = se,
    condition_number = kappa(hessian)
  )
}

# Master Diagnostic Report ----------------------------------------------------

#' Generate comprehensive diagnostic report
#'
#' Compiles all diagnostic information into a single report.
#'
#' @return List with all diagnostic data
#' @export
get_diagnostic_report <- function() {
  list(
    timing = get_timing_report(),
    memory = get_memory_report(),
    session_info = sessionInfo(),
    timestamp = Sys.time()
  )
}
