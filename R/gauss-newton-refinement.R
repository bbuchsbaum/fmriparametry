#' Gauss-Newton refinement for hard voxels
#'
#' Implements full nonlinear Gauss-Newton optimization for the most challenging
#' voxels that don't respond well to Taylor approximation methods. This is the
#' most computationally intensive refinement but can recover good fits for
#' difficult cases.
#'
#' @param theta_hat_voxel Matrix of current parameter estimates (voxels x parameters)
#' @param r2_voxel Numeric vector of current R-squared values
#' @param Y_proj Numeric matrix of projected BOLD data (timepoints x voxels)
#' @param S_target_proj Numeric matrix of projected stimulus design
#' @param scan_times Numeric vector of scan acquisition times
#' @param hrf_eval_times Numeric vector of HRF evaluation time points
#' @param hrf_interface List with HRF interface functions
#' @param theta_bounds List with elements `lower` and `upper`
#' @param queue_labels Character vector of refinement queue assignments
#' @param max_iter_gn Maximum iterations for Gauss-Newton
#' @param tol_gn Convergence tolerance
#' @param lambda_ridge Ridge penalty
#' @param step_size Initial step size for line search
#' @param verbose Logical whether to print progress
#' @param convergence_config Convergence configuration list
#' @param use_conv_cache Logical; if TRUE, reuse precomputed FFT signal context
#'   across objective/Jacobian calls in this refinement pass.
#' @param max_line_search_steps Maximum backtracking attempts per GN iteration.
#' @param min_line_search_alpha Smallest admissible line-search step.
#' @param parallel Logical; if TRUE, run the hard-voxel loop in parallel when
#'   feasible on supported platforms.
#' @param n_cores Optional number of worker cores for parallel hard-voxel
#'   refinement (`NULL` = auto-detect).
#' @param parallel_min_voxels Minimum hard-voxel count required before enabling
#'   parallel execution.
#'
.gauss_newton_refinement <- function(
  theta_hat_voxel,
  r2_voxel,
  Y_proj,
  S_target_proj,
  scan_times,
  hrf_eval_times,
  hrf_interface,
  theta_bounds,
  queue_labels,
  max_iter_gn = 5,
  tol_gn = 1e-4,
  lambda_ridge = 0.01,
  step_size = 1.0,
  verbose = FALSE,
  convergence_config = NULL,
  use_conv_cache = TRUE,
  max_line_search_steps = 4L,
  min_line_search_alpha = 1e-6,
  parallel = FALSE,
  n_cores = NULL,
  parallel_min_voxels = 200L
) {
  n_params <- ncol(theta_hat_voxel)
  n_time <- nrow(Y_proj)
  
  # Use unified convergence config if provided, otherwise use defaults
  if (is.null(convergence_config)) {
    convergence_config <- .create_convergence_config(
      param_tol = tol_gn,
      objective_tol = tol_gn,
      max_iter = max_iter_gn
    )
  }
  
  # Identify hard voxels
  idx_hard <- which(queue_labels == "hard_GN")
  n_hard <- length(idx_hard)
  
  if (n_hard == 0) {
    if (verbose) cat("No hard voxels to refine\n")
    return(list(
      theta_hat = theta_hat_voxel,
      r2 = r2_voxel,
      n_refined = 0,
      n_converged = 0,
      n_improved = 0,
      convergence_status = character(0)
    ))
  }
  
  if (verbose) {
    cat("\nGauss-Newton refinement for", n_hard, "hard voxels\n")
    cat("  Max iterations:", max_iter_gn, "\n")
    cat("  Convergence tolerance:", tol_gn, "\n")
  }
  
  # Store original values
  r2_orig <- r2_voxel

  # Shared stimulus/convolution context for all hard voxels.
  stim_signal <- .extract_primary_stimulus(S_target_proj)
  conv_context <- NULL
  if (isTRUE(use_conv_cache)) {
    conv_context <- .prepare_fft_convolution_context(
      signal = stim_signal,
      output_length = n_time,
      kernel_length = length(hrf_eval_times)
    )
  }
  parallel_cfg <- .resolve_gn_parallel_config(
    parallel = parallel,
    n_cores = n_cores,
    n_hard = n_hard,
    parallel_min_voxels = parallel_min_voxels
  )

  if (verbose) {
    if (isTRUE(parallel_cfg$use_parallel)) {
      cat("  Hard-voxel loop: parallel (", parallel_cfg$n_cores, " cores)\n", sep = "")
    } else if (isTRUE(parallel)) {
      cat("  Hard-voxel loop: serial (", parallel_cfg$reason, ")\n", sep = "")
    }
  }

  refine_one <- function(v) {
    .gauss_newton_refine_single_voxel(
      v = v,
      theta_initial = theta_hat_voxel[v, ],
      r2_initial = r2_voxel[v],
      y_v = Y_proj[, v],
      stim_signal = stim_signal,
      hrf_eval_times = hrf_eval_times,
      hrf_interface = hrf_interface,
      theta_bounds = theta_bounds,
      n_params = n_params,
      n_time = n_time,
      convergence_config = convergence_config,
      lambda_ridge = lambda_ridge,
      step_size = step_size,
      max_line_search_steps = max_line_search_steps,
      min_line_search_alpha = min_line_search_alpha,
      conv_context = conv_context
    )
  }

  hard_results <- if (isTRUE(parallel_cfg$use_parallel)) {
    parallel::mclapply(
      idx_hard,
      refine_one,
      mc.cores = parallel_cfg$n_cores,
      mc.preschedule = TRUE
    )
  } else {
    lapply(idx_hard, refine_one)
  }

  convergence_status <- vapply(hard_results, `[[`, character(1), "status")
  iteration_counts <- vapply(hard_results, `[[`, numeric(1), "iterations")
  converged_flags <- vapply(hard_results, `[[`, logical(1), "converged")
  improved_flags <- vapply(hard_results, `[[`, logical(1), "improved")
  n_converged <- sum(converged_flags)
  n_improved <- sum(improved_flags)

  improved_idx <- which(improved_flags)
  if (length(improved_idx) > 0L) {
    improved_vox <- idx_hard[improved_idx]
    theta_updates <- do.call(rbind, lapply(hard_results[improved_idx], `[[`, "theta"))
    if (!is.matrix(theta_updates)) {
      theta_updates <- matrix(theta_updates, nrow = 1L)
    }
    theta_hat_voxel[improved_vox, ] <- theta_updates
    r2_voxel[improved_vox] <- vapply(hard_results[improved_idx], `[[`, numeric(1), "r2")
  }
  
  if (verbose) {
    cat("  Gauss-Newton refinement complete:\n")
    cat("    Refined:", n_hard, "voxels\n")
    cat("    Converged:", n_converged, "(", 
        round(100 * n_converged / n_hard, 1), "%)\n")
    cat("    Improved:", n_improved, "(", 
        round(100 * n_improved / n_hard, 1), "%)\n")
    cat("    Mean iterations:", round(mean(iteration_counts), 1), "\n")
    
    # Convergence summary
    conv_table <- table(convergence_status)
    cat("    Convergence status:\n")
    for (status in names(conv_table)) {
      cat("      ", status, ":", conv_table[status], "\n")
    }
  }
  
  # Update queue labels for successfully refined voxels
  successfully_refined <- idx_hard[which(r2_voxel[idx_hard] > r2_orig[idx_hard])]
  if (length(successfully_refined) > 0) {
    queue_labels[successfully_refined] <- "easy"
  }
  
  # Return results
  list(
    theta_hat = theta_hat_voxel,
    r2 = r2_voxel,
    queue_labels = queue_labels,
    n_refined = n_hard,
    n_converged = n_converged,
    n_improved = n_improved,
    convergence_status = convergence_status,
    iteration_counts = iteration_counts,
    improvement_summary = list(
      mean_r2_improvement = mean(r2_voxel[idx_hard] - r2_orig[idx_hard], na.rm = TRUE),
      max_r2_improvement = max(r2_voxel[idx_hard] - r2_orig[idx_hard], na.rm = TRUE)
    )
  )
}

.resolve_gn_parallel_config <- function(parallel, n_cores, n_hard, parallel_min_voxels) {
  min_voxels <- suppressWarnings(as.integer(parallel_min_voxels[1]))
  if (!is.finite(min_voxels) || min_voxels < 1L) {
    min_voxels <- 1L
  }
  if (!isTRUE(parallel)) {
    return(list(use_parallel = FALSE, n_cores = 1L, reason = "parallel=FALSE"))
  }
  if (n_hard < min_voxels) {
    return(list(use_parallel = FALSE, n_cores = 1L, reason = "too_few_voxels"))
  }
  if (.Platform$OS.type == "windows") {
    return(list(use_parallel = FALSE, n_cores = 1L, reason = "mclapply_unavailable_on_windows"))
  }

  detected <- suppressWarnings(tryCatch(
    parallel::detectCores(logical = FALSE),
    error = function(...) NA_integer_
  ))
  if (!is.finite(detected) || detected < 1L) {
    detected <- 1L
  }

  requested <- if (is.null(n_cores) || length(n_cores) == 0L) {
    detected
  } else {
    suppressWarnings(as.integer(n_cores[1]))
  }
  if (length(requested) != 1L || !is.finite(requested) || requested < 2L) {
    return(list(use_parallel = FALSE, n_cores = 1L, reason = "n_cores<2"))
  }

  cores <- min(requested, n_hard)
  if (cores < 2L) {
    return(list(use_parallel = FALSE, n_cores = 1L, reason = "single_core_after_capping"))
  }
  list(use_parallel = TRUE, n_cores = cores, reason = "ok")
}

.gauss_newton_refine_single_voxel <- function(
  v,
  theta_initial,
  r2_initial,
  y_v,
  stim_signal,
  hrf_eval_times,
  hrf_interface,
  theta_bounds,
  n_params,
  n_time,
  convergence_config,
  lambda_ridge,
  step_size,
  max_line_search_steps,
  min_line_search_alpha,
  conv_context
) {
  theta_initial <- as.numeric(theta_initial)
  if (length(theta_initial) != n_params || any(!is.finite(theta_initial))) {
    return(list(
      voxel = v,
      theta = theta_initial,
      r2 = r2_initial,
      improved = FALSE,
      converged = FALSE,
      status = "invalid_input",
      iterations = 0
    ))
  }

  if (is.finite(r2_initial) && r2_initial > 0.9) {
    return(list(
      voxel = v,
      theta = theta_initial,
      r2 = r2_initial,
      improved = FALSE,
      converged = FALSE,
      status = "already_good",
      iterations = 0
    ))
  }

  y_centered_ss <- sum((y_v - mean(y_v))^2)
  theta_current <- theta_initial
  theta_best <- theta_current
  r2_best <- r2_initial

  obj_current <- .calculate_objective_gn_fast(
    theta_current,
    y_v,
    stim_signal,
    hrf_eval_times,
    hrf_interface,
    n_time,
    conv_context = conv_context,
    theta_bounds = theta_bounds
  )
  if (!is.finite(obj_current)) {
    return(list(
      voxel = v,
      theta = theta_initial,
      r2 = r2_initial,
      improved = FALSE,
      converged = FALSE,
      status = "singular_system",
      iterations = 0
    ))
  }

  alpha_base <- step_size
  identity_n <- diag(n_params)
  converged <- FALSE
  status <- "max_iterations"
  iter_used <- 0L

  for (iter in seq_len(convergence_config$max_iter)) {
    iter_used <- iter
    jacob_info <- .get_jacobian_and_residuals_fast(
      theta_current,
      y_v,
      stim_signal,
      hrf_eval_times,
      hrf_interface,
      n_time,
      conv_context = conv_context,
      theta_bounds = theta_bounds
    )
    if (is.null(jacob_info)) {
      status <- "singular_system"
      break
    }

    J <- jacob_info$jacobian
    residuals <- jacob_info$residuals
    JtJ_ridge <- crossprod(J) + lambda_ridge * identity_n
    Jtr <- crossprod(J, residuals)
    delta <- .solve_gn_direction(JtJ_ridge, Jtr)
    if (is.null(delta) || length(delta) != n_params) {
      status <- "singular_system"
      break
    }

    alpha <- alpha_base
    accepted <- FALSE
    accepted_alpha <- alpha_base
    theta_new <- theta_current
    obj_new <- obj_current

    for (ls_iter in seq_len(max_line_search_steps)) {
      theta_proposal <- .apply_bounds(theta_current + alpha * delta, theta_bounds)
      obj_proposal <- .calculate_objective_gn_fast(
        theta_proposal,
        y_v,
        stim_signal,
        hrf_eval_times,
        hrf_interface,
        n_time,
        conv_context = conv_context,
        theta_bounds = theta_bounds
      )

      if (is.finite(obj_proposal) && obj_proposal < obj_current) {
        theta_new <- theta_proposal
        obj_new <- obj_proposal
        accepted <- TRUE
        accepted_alpha <- alpha
        break
      }

      alpha <- alpha * 0.5
      if (alpha < min_line_search_alpha) {
        break
      }
    }

    if (!accepted) {
      status <- "line_search_failed"
      break
    }

    if (accepted_alpha >= alpha_base * 0.99) {
      alpha_base <- min(1.0, alpha_base * 1.2)
    } else {
      alpha_base <- max(min_line_search_alpha, accepted_alpha)
    }

    conv_check <- .check_convergence(
      current = theta_new,
      previous = theta_current,
      obj_current = obj_new,
      obj_previous = obj_current,
      config = convergence_config
    )

    theta_current <- theta_new
    obj_current <- obj_new

    r2_current <- if (y_centered_ss > 1e-12) {
      1 - obj_current / y_centered_ss
    } else {
      0
    }
    if (is.finite(r2_current) && r2_current > r2_best) {
      theta_best <- theta_current
      r2_best <- r2_current
    }

    if (conv_check$converged) {
      converged <- TRUE
      status <- conv_check$reason
      break
    }
  }

  improved <- is.finite(r2_best) && is.finite(r2_initial) && (r2_best > r2_initial)
  list(
    voxel = v,
    theta = if (improved) theta_best else theta_initial,
    r2 = if (improved) r2_best else r2_initial,
    improved = improved,
    converged = converged,
    status = status,
    iterations = iter_used
  )
}

#' Calculate objective function for Gauss-Newton
#'
#' Returns `Inf` if the HRF predictor has near-zero magnitude, which
#' indicates that the system is singular and the amplitude cannot be
#' estimated reliably.
#'
#' @param theta Numeric vector of HRF parameters
#' @param y Numeric vector of observed BOLD data
#' @param S Numeric matrix of stimulus design
#' @param t_hrf Numeric vector of HRF evaluation time points
#' @param hrf_interface List with HRF interface functions
#' @param n_time Integer number of time points
#' @return Numeric scalar sum of squared residuals
#' @noRd
.extract_primary_stimulus <- function(S) {
  if (is.numeric(S) && is.null(dim(S))) {
    return(S)
  }
  if (is.null(dim(S))) {
    return(as.numeric(S))
  }
  as.numeric(S[, 1])
}

.solve_gn_direction <- function(JtJ_ridge, Jtr) {
  if (!is.matrix(JtJ_ridge)) {
    return(NULL)
  }
  if (!all(is.finite(JtJ_ridge)) || !all(is.finite(Jtr))) {
    return(NULL)
  }

  n_params <- nrow(JtJ_ridge)
  if (n_params == 0L || ncol(JtJ_ridge) != n_params || length(Jtr) != n_params) {
    return(NULL)
  }

  # Symmetrize to damp minor numerical asymmetry before eigensolve.
  A <- 0.5 * (JtJ_ridge + t(JtJ_ridge))
  eig <- eigen(A, symmetric = TRUE)
  values <- eig$values
  vectors <- eig$vectors

  if (!all(is.finite(values)) || !all(is.finite(vectors))) {
    return(NULL)
  }

  scale <- max(abs(values))
  if (!is.finite(scale) || scale <= 0) {
    return(NULL)
  }

  tol <- scale * 1e-12
  inv_values <- ifelse(values > tol, 1 / values, 0)
  rhs <- as.numeric(crossprod(vectors, Jtr))
  delta <- -as.numeric(vectors %*% (inv_values * rhs))

  if (length(delta) != n_params || !all(is.finite(delta))) {
    return(NULL)
  }

  delta
}

.resolve_lwu_cpp_bounds <- function(hrf_interface, theta_bounds = NULL) {
  bounds <- theta_bounds
  if (is.null(bounds) && !is.null(hrf_interface$active_bounds)) {
    bounds <- hrf_interface$active_bounds
  }
  if (is.null(bounds) && is.function(hrf_interface$default_bounds)) {
    bounds <- hrf_interface$default_bounds()
  }
  if (is.null(bounds) || !is.list(bounds) || is.null(bounds$lower) || is.null(bounds$upper)) {
    return(NULL)
  }

  lower <- as.numeric(bounds$lower)
  upper <- as.numeric(bounds$upper)
  if (length(lower) != 3L || length(upper) != 3L) {
    return(NULL)
  }
  if (any(!is.finite(lower)) || any(!is.finite(upper))) {
    return(NULL)
  }

  list(lower = lower, upper = upper)
}

.can_use_lwu_cpp_gn <- function(theta, y, stim_signal, t_hrf, hrf_interface, lwu_bounds) {
  has_cpp <- exists("lwu_gn_objective_cpp", mode = "function") &&
    exists("lwu_gn_jacobian_cpp", mode = "function")

  is_lwu <- !is.null(hrf_interface$parameter_names) &&
    identical(unname(as.character(hrf_interface$parameter_names)), c("tau", "sigma", "rho"))

  has_cpp &&
    is_lwu &&
    !is.null(lwu_bounds) &&
    is.numeric(theta) &&
    length(theta) == 3L &&
    is.numeric(y) &&
    is.numeric(stim_signal) &&
    is.numeric(t_hrf)
}

.calculate_objective_gn <- function(theta, y, S, t_hrf, hrf_interface, n_time) {
  # Validate theta
  if (!is.numeric(theta) || length(theta) != length(hrf_interface$parameter_names)) {
    stop(sprintf("Invalid theta: expected numeric vector of length %d, got %s of length %d",
                 length(hrf_interface$parameter_names),
                 class(theta)[1],
                 length(theta)))
  }

  stim_signal <- .extract_primary_stimulus(S)
  .calculate_objective_gn_fast(
    theta, y, stim_signal, t_hrf, hrf_interface, n_time,
    conv_context = NULL,
    theta_bounds = NULL
  )
}

.calculate_objective_gn_fast <- function(theta, y, stim_signal, t_hrf, hrf_interface, n_time,
                                         conv_context = NULL,
                                         theta_bounds = NULL) {
  if (!is.numeric(theta) || length(theta) != length(hrf_interface$parameter_names)) {
    return(Inf)
  }
  if (!is.numeric(y) || length(y) != n_time || any(!is.finite(y))) {
    return(Inf)
  }
  if (!is.numeric(stim_signal) || any(!is.finite(stim_signal))) {
    return(Inf)
  }

  lwu_bounds <- .resolve_lwu_cpp_bounds(hrf_interface, theta_bounds = theta_bounds)
  if (.can_use_lwu_cpp_gn(theta, y, stim_signal, t_hrf, hrf_interface, lwu_bounds)) {
    cpp_obj <- tryCatch(
      lwu_gn_objective_cpp(
        theta = as.numeric(theta),
        y = as.numeric(y),
        signal = as.numeric(stim_signal),
        t_hrf_eval = as.numeric(t_hrf),
        lower = lwu_bounds$lower,
        upper = lwu_bounds$upper
      ),
      error = function(...) NA_real_
    )
    if (is.numeric(cpp_obj) && length(cpp_obj) == 1L) {
      if (is.finite(cpp_obj)) {
        return(as.numeric(cpp_obj))
      }
      if (is.infinite(cpp_obj)) {
        return(Inf)
      }
    }
  }

  # Fallback: R implementation
  hrf_vals <- hrf_interface$hrf_function(t_hrf, theta)
  if (!is.numeric(hrf_vals) || any(!is.finite(hrf_vals))) {
    return(Inf)
  }

  x_pred_raw <- .convolve_signal_with_kernels(
    stim_signal, matrix(hrf_vals, ncol = 1), n_time, conv_context = conv_context
  )[, 1]
  if (!is.numeric(x_pred_raw) || any(!is.finite(x_pred_raw))) {
    return(Inf)
  }

  denom <- sum(x_pred_raw^2)
  if (!is.finite(denom) || denom < 1e-8) {
    return(Inf)
  }

  beta <- as.numeric(crossprod(x_pred_raw, y)) / denom
  if (!is.finite(beta)) {
    return(Inf)
  }
  x_pred <- beta * x_pred_raw
  if (any(!is.finite(x_pred))) {
    return(Inf)
  }

  obj <- sum((y - x_pred)^2)
  if (!is.finite(obj)) {
    return(Inf)
  }
  obj
}

#' Get Jacobian matrix and residuals for Gauss-Newton
#'
#' Returns `NULL` if the HRF predictor has near-zero magnitude, which
#' indicates a singular system.
#'
#' @param theta Numeric vector of HRF parameters
#' @param y Numeric vector of observed BOLD data
#' @param S Numeric matrix of stimulus design
#' @param t_hrf Numeric vector of HRF evaluation time points
#' @param hrf_interface List with HRF interface functions
#' @param n_time Integer number of time points
#' @return List with jacobian matrix, residuals, and amplitude, or NULL if singular
#' @noRd
.get_jacobian_and_residuals <- function(theta, y, S, t_hrf, hrf_interface, n_time) {
  stim_signal <- .extract_primary_stimulus(S)
  .get_jacobian_and_residuals_fast(
    theta, y, stim_signal, t_hrf, hrf_interface, n_time,
    conv_context = NULL,
    theta_bounds = NULL
  )
}

.get_jacobian_and_residuals_fast <- function(theta, y, stim_signal, t_hrf, hrf_interface, n_time,
                                             conv_context = NULL,
                                             theta_bounds = NULL) {
  n_params <- length(theta)
  if (n_params == 0L || !is.numeric(theta) || any(!is.finite(theta))) {
    return(NULL)
  }
  if (!is.numeric(y) || length(y) != n_time || any(!is.finite(y))) {
    return(NULL)
  }
  if (!is.numeric(stim_signal) || any(!is.finite(stim_signal))) {
    return(NULL)
  }
  
  lwu_bounds <- .resolve_lwu_cpp_bounds(hrf_interface, theta_bounds = theta_bounds)
  if (.can_use_lwu_cpp_gn(theta, y, stim_signal, t_hrf, hrf_interface, lwu_bounds)) {
    cpp_info <- tryCatch(
      lwu_gn_jacobian_cpp(
        theta = as.numeric(theta),
        y = as.numeric(y),
        signal = as.numeric(stim_signal),
        t_hrf_eval = as.numeric(t_hrf),
        lower = lwu_bounds$lower,
        upper = lwu_bounds$upper,
        rel_step = 1e-4,
        min_step = 1e-5,
        bound_eps = 1e-8
      ),
      error = function(...) NULL
    )

    if (is.list(cpp_info) && isTRUE(cpp_info$ok)) {
      jacobian <- cpp_info$jacobian
      residuals <- cpp_info$residuals
      amplitude <- as.numeric(cpp_info$amplitude)
      if (!is.matrix(jacobian) || nrow(jacobian) != n_time || ncol(jacobian) != 3L) {
        return(NULL)
      }
      if (!is.numeric(residuals) || length(residuals) != n_time) {
        return(NULL)
      }
      if (length(amplitude) != 1L || !is.finite(amplitude)) {
        return(NULL)
      }
      if (any(!is.finite(jacobian)) || any(!is.finite(residuals))) {
        return(NULL)
      }
      return(list(
        jacobian = jacobian,
        residuals = residuals,
        amplitude = amplitude
      ))
    }

    # C++ path ran and marked singular/invalid
    if (is.list(cpp_info) && isFALSE(cpp_info$ok)) {
      return(NULL)
    }
  }

  # Fallback: R implementation
  # Get Taylor basis (HRF and derivatives)
  taylor_basis <- hrf_interface$taylor_basis(theta, t_hrf)
  if (!is.matrix(taylor_basis)) {
    taylor_basis <- matrix(taylor_basis, ncol = n_params + 1)
  }
  if (ncol(taylor_basis) != n_params + 1 || nrow(taylor_basis) == 0L) {
    return(NULL)
  }
  if (any(!is.finite(taylor_basis))) {
    return(NULL)
  }
  
  # Convolve HRF and derivative basis in a single batched call.
  X_conv <- .fast_batch_convolution(
    stim_signal, taylor_basis, n_time, conv_context = conv_context
  )
  if (!is.matrix(X_conv) || any(!is.finite(X_conv))) {
    return(NULL)
  }
  
  # Fit amplitude for current HRF
  x_hrf <- X_conv[, 1]
  denom <- sum(x_hrf^2)
  if (!is.finite(denom) || denom < 1e-8) {
    return(NULL)
  }
  beta <- as.numeric(crossprod(x_hrf, y)) / denom
  if (!is.finite(beta)) {
    return(NULL)
  }
  
  # Residuals
  residuals <- y - beta * x_hrf
  if (any(!is.finite(residuals))) {
    return(NULL)
  }
  
  deriv <- X_conv[, -1, drop = FALSE]
  dot_dy <- as.numeric(crossprod(deriv, y))
  dot_xd <- as.numeric(crossprod(x_hrf, deriv))
  dbeta <- (dot_dy - 2 * beta * dot_xd) / denom
  jacobian <- -beta * deriv - tcrossprod(x_hrf, dbeta)
  if (any(!is.finite(jacobian))) {
    return(NULL)
  }
  
  list(
    jacobian = jacobian,
    residuals = residuals,
    amplitude = beta
  )
}
