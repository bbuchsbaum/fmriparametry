# Shared utilities for HRF estimation
# These functions are used by both legacy and refactored implementations


#' Generate a grid of seed points for multi-seed estimation
#'
#' Creates a set of seed points spanning the LWU parameter space.
#' Each seed is clamped to the active bounds.
#'
#' @param theta_bounds List with `lower` and `upper` bound vectors
#' @param user_grid Optional matrix (n_seeds x n_params) of user-specified seeds
#' @return Matrix of seed points (n_seeds x n_params)
#' @keywords internal
.generate_seed_grid <- function(theta_bounds, user_grid = NULL) {
  n_params <- length(theta_bounds$lower)

  if (!is.null(user_grid)) {
    if (!is.matrix(user_grid)) {
      user_grid <- matrix(user_grid, ncol = n_params, byrow = TRUE)
    } else if (ncol(user_grid) != n_params) {
      stop("user_grid must have one column per parameter", call. = FALSE)
    }
    # Clamp each seed to bounds
    for (i in seq_len(nrow(user_grid))) {
      user_grid[i, ] <- pmax(theta_bounds$lower, pmin(theta_bounds$upper, user_grid[i, ]))
    }
    return(user_grid)
  }

  # Fallback for non-LWU models: center point + lower/upper corners.
  if (n_params != 3L) {
    center <- (theta_bounds$lower + theta_bounds$upper) / 2
    grid <- rbind(center, theta_bounds$lower, theta_bounds$upper)
    return(unique(grid))
  }

  # Default 5-point grid covering the LWU space
  # [tau, sigma, rho]
  grid <- rbind(
    c(6.0, 2.5, 0.35),   # default center
    c(3.0, 1.5, 0.25),   # early, narrow HRF
    c(9.0, 3.0, 0.50),   # late, wide HRF
    c(5.0, 1.0, 0.15),   # narrow with minimal undershoot
    c(7.0, 4.0, 0.60)    # wide with strong undershoot
  )

  # Clamp to active bounds
  for (i in seq_len(nrow(grid))) {
    grid[i, ] <- pmax(theta_bounds$lower, pmin(theta_bounds$upper, grid[i, ]))
  }

  # Remove duplicate rows after clamping
  grid <- unique(grid)
  grid
}

#' Compute data-driven initial HRF parameters
#'
#' Estimates latency and width by cross-correlating the average high-variance
#' voxels with the stimulus.
#'
#' @param Y BOLD signal matrix (time x voxels)
#' @param S Stimulus design matrix
#' @param hrf_interface HRF model interface
#' @param theta_bounds List with parameter lower and upper bounds
#'
#' @return Numeric vector of initial parameter estimates
#' @keywords internal
.compute_data_driven_seed <- function(Y, S, hrf_interface, theta_bounds) {
  default_seed <- hrf_interface$default_seed()
  bounds <- theta_bounds

  # Average high variance voxels
  vars <- apply(Y, 2, var)
  idx <- which(vars >= quantile(vars, 0.75))
  y_mean <- rowMeans(Y[, idx, drop = FALSE])

  # Cross correlation with stimulus
  lags <- seq(-8, 12, by = 1)
  
  # Extract stimulus column as vector and ensure proper length
  stim_vec <- if (is.matrix(S)) S[, 1] else as.vector(S)
  n_time <- length(y_mean)
  
  # Ensure stimulus vector has same length as y_mean
  if (length(stim_vec) != n_time) {
    stim_vec <- stim_vec[1:n_time]
  }
  
  cors <- vapply(lags, function(l) {
    if (l >= 0) {
      s_shift <- c(rep(0, l), stim_vec)[1:n_time]
    } else {
      s_shift <- c(stim_vec, rep(0, -l))[(-l + 1):(n_time - l)]
    }
    # Ensure both are vectors of same length
    if (length(s_shift) != n_time) {
      s_shift <- s_shift[1:n_time]
    }
    stats::cor(as.vector(y_mean), as.vector(s_shift), use = "complete.obs")
  }, numeric(1))

  best_lag <- lags[which.max(abs(cors))]
  tau_est <- default_seed[1] + best_lag
  tau_est <- max(bounds$lower[1], min(bounds$upper[1], tau_est))

  cors_abs <- abs(cors)
  half <- max(cors_abs) * 0.5
  width_idx <- which(cors_abs >= half)
  sigma_est <- if (length(width_idx) >= 2) diff(range(lags[width_idx]))/2 else default_seed[2]
  sigma_est <- max(bounds$lower[2], min(bounds$upper[2], sigma_est))

  c(tau_est, sigma_est, default_seed[3])
}

#' Initialize HRF parameters with K-means clustering
#'
#' Voxels are clustered by their estimated latency from cross-correlation and
#' cluster-specific seeds are derived.
#'
#' @param Y BOLD signal matrix
#' @param S Stimulus design matrix
#' @param k Number of clusters
#' @param hrf_interface HRF model interface
#' @param theta_bounds List with parameter lower and upper bounds
#'
#' @return List with cluster assignments and center parameter seeds
#' @keywords internal
.perform_kmeans_initialization <- function(Y, S, k, hrf_interface, theta_bounds) {
  n_vox <- ncol(Y)
  if (k <= 1 || n_vox <= k) {
    return(list(cluster = rep(1, n_vox),
                centers = matrix(hrf_interface$default_seed(), nrow = 1)))
  }

  bounds <- theta_bounds
  delays <- seq(-8, 12, by = 1)
  n_time <- nrow(Y)

  # Matrix-based cross-correlation for all voxels simultaneously
  stim <- if (is.matrix(S)) S[, 1] else S
  stim_shifts <- vapply(delays, function(l) {
    if (l >= 0) {
      c(rep(0, l), stim)[1:n_time]
    } else {
      c(stim, rep(0, -l))[(-l + 1):n_time]
    }
  }, numeric(n_time))

  Y_cent <- scale(Y, center = TRUE, scale = FALSE)
  stim_cent <- scale(stim_shifts, center = TRUE, scale = FALSE)
  Y_sd <- sqrt(colSums(Y_cent^2) / (n_time - 1))
  stim_sd <- sqrt(colSums(stim_cent^2) / (n_time - 1))

  cors_mat <- (t(Y_cent) %*% stim_cent) / (n_time - 1)
  cors_mat <- sweep(cors_mat, 1, Y_sd, "/")
  cors_mat <- sweep(cors_mat, 2, stim_sd, "/")

  tau_est <- delays[max.col(abs(cors_mat), ties.method = "first")]

  features <- matrix(tau_est, ncol = 1)
  km <- kmeans(features, centers = k, nstart = 20, iter.max = 50)

  centers <- matrix(rep(hrf_interface$default_seed(), k), nrow = k, byrow = TRUE)
  for (cl in seq_len(k)) {
    idx <- which(km$cluster == cl)
    if (length(idx) > 0) {
      centers[cl, 1] <- median(tau_est[idx])
    }
  }
  centers[, 1] <- pmax(bounds$lower[1], pmin(bounds$upper[1], centers[, 1]))

  list(cluster = km$cluster, centers = centers, iterations = km$iter)
}

#' Classify voxels for refinement
#' 
#' Classifies voxels into easy, moderate, and hard categories based on
#' their R-squared values and, optionally, standard errors.
#' 
#' @param r_squared Vector of R-squared values
#' @param se_theta Matrix of standard errors (voxels x parameters), or NULL for R^2-only mode
#' @param thresholds List with classification thresholds (r2_easy, r2_hard, se_low, se_high)
#' @return Character vector of voxel classifications
#' @keywords internal
.classify_voxels_for_refinement <- function(r_squared, se_theta = NULL, thresholds) {
  classes <- rep("moderate", length(r_squared))
  use_se <- !is.null(se_theta)

  # 1. Define base masks using R-squared
  easy_mask_r2 <- r_squared > thresholds$r2_easy
  hard_mask_r2 <- r_squared < thresholds$r2_hard

  # 2. Initialize final masks
  easy_mask <- easy_mask_r2
  hard_mask <- hard_mask_r2

  # 3. Refine masks if standard errors are provided
  if (use_se) {
    # Validate se_theta dimensions
    if (nrow(se_theta) != length(r_squared)) {
      stop("nrow(se_theta) must match length(r_squared).")
    }
    mean_se <- rowMeans(se_theta, na.rm = TRUE)
    
    # Easy: high R^2 AND low SE
    easy_mask <- easy_mask_r2 & (mean_se < thresholds$se_low)
    
    # Hard: low R^2 OR high SE
    hard_mask <- hard_mask_r2 | (mean_se > thresholds$se_high)
  }
  
  # 4. Apply classifications (order matters: hard can override easy)
  classes[easy_mask] <- "easy"
  classes[hard_mask] <- "hard"
  
  return(classes)
}

#' Compute standard errors via Delta method
#' 
#' Computes standard errors for HRF parameters and amplitudes using
#' the Delta method approximation
#' 
#' @param theta_hat Matrix of estimated parameters (voxels x parameters)
#' @param beta0 Vector of amplitudes
#' @param Y_proj Projected BOLD data
#' @param S_target_proj Projected stimulus matrix
#' @param hrf_interface HRF model interface
#' @param hrf_eval_times HRF evaluation times
#' @return List with se_theta_hat and se_beta0
#' @keywords internal
.compute_standard_errors_delta <- function(theta_hat, beta0, Y_proj, S_target_proj,
                                          hrf_interface, hrf_eval_times) {
  n_vox <- ncol(Y_proj)
  n_params <- length(hrf_interface$parameter_names)
  n_time <- nrow(Y_proj)

  # Ensure theta_hat is a matrix
  if (!is.matrix(theta_hat) || length(dim(theta_hat)) != 2) {
    theta_hat <- matrix(theta_hat, nrow = n_vox, ncol = n_params,
                        byrow = FALSE)
    colnames(theta_hat) <- hrf_interface$parameter_names
  }
  
  basis_list <- vector("list", n_vox)
  hrf_list <- vector("list", n_vox)
  for (v in seq_len(n_vox)) {
    theta_v <- theta_hat[v, ]
    taylor_basis <- hrf_interface$taylor_basis(theta_v, hrf_eval_times)
    if (!is.matrix(taylor_basis)) {
      taylor_basis <- matrix(taylor_basis, ncol = n_params + 1)
    }
    # Ensure numeric storage mode for C++
    derivatives <- taylor_basis[, -1, drop = FALSE]
    storage.mode(derivatives) <- "double"
    basis_list[[v]] <- derivatives
    
    hrf_vals <- hrf_interface$hrf_function(hrf_eval_times, theta_v)
    storage.mode(hrf_vals) <- "double"
    hrf_list[[v]] <- hrf_vals
  }

  # Ensure Y_proj and S_target_proj have the right storage mode
  if (storage.mode(Y_proj) != "double") {
    storage.mode(Y_proj) <- "double"
  }
  if (storage.mode(S_target_proj) != "double") {
    storage.mode(S_target_proj) <- "double"
  }

  res <- compute_standard_errors_bulk_cpp(
    basis_list, hrf_list, Y_proj, S_target_proj, as.numeric(beta0)
  )

  return(res)
}

.invert_psd_safe <- function(A, rel_tol = 1e-10) {
  if (!is.matrix(A) || nrow(A) != ncol(A) || nrow(A) == 0) {
    return(NULL)
  }
  if (any(!is.finite(A))) {
    return(NULL)
  }

  A <- 0.5 * (A + t(A))
  eig <- tryCatch(eigen(A, symmetric = TRUE), error = function(...) NULL)
  if (is.null(eig) || any(!is.finite(eig$values)) || any(!is.finite(eig$vectors))) {
    return(NULL)
  }

  scale <- max(abs(eig$values))
  if (!is.finite(scale) || scale <= 0) {
    return(NULL)
  }
  keep <- eig$values > (scale * rel_tol)
  if (!any(keep)) {
    return(NULL)
  }

  V <- eig$vectors[, keep, drop = FALSE]
  inv_vals <- 1 / eig$values[keep]
  inv_A <- V %*% (diag(inv_vals, nrow = length(inv_vals)) %*% t(V))
  if (any(!is.finite(inv_A))) {
    return(NULL)
  }
  inv_A
}

#' Compute heteroskedasticity-robust standard errors
#'
#' Uses a sandwich covariance estimator (HC0) for profiled Gauss-Newton Jacobians.
#' This is slower than the delta approximation but more robust under
#' heteroskedastic residual variance.
#'
#' @keywords internal
.compute_standard_errors_sandwich <- function(
  theta_hat, beta0, Y_proj, S_target_proj,
  hrf_interface, hrf_eval_times, baseline_model = NULL
) {
  n_time <- nrow(Y_proj)
  n_vox <- ncol(Y_proj)
  n_params <- length(hrf_interface$parameter_names)
  has_intercept <- .is_intercept_baseline(baseline_model)

  if (!is.matrix(theta_hat) || length(dim(theta_hat)) != 2) {
    theta_hat <- matrix(theta_hat, nrow = n_vox, ncol = n_params, byrow = FALSE)
  }

  se_theta_hat <- matrix(NA_real_, nrow = n_vox, ncol = n_params)
  colnames(se_theta_hat) <- hrf_interface$parameter_names
  se_beta0 <- rep(NA_real_, n_vox)

  stim_signal <- .extract_primary_stimulus(S_target_proj)

  for (v in seq_len(n_vox)) {
    theta_v <- as.numeric(theta_hat[v, ])
    y_v <- as.numeric(Y_proj[, v])

    jac_info <- .get_jacobian_and_residuals_fast(
      theta = theta_v,
      y = y_v,
      stim_signal = stim_signal,
      t_hrf = hrf_eval_times,
      hrf_interface = hrf_interface,
      n_time = n_time,
      conv_context = NULL,
      theta_bounds = NULL,
      has_intercept = has_intercept
    )
    if (is.null(jac_info)) {
      next
    }

    J <- jac_info$jacobian
    e <- as.numeric(jac_info$residuals)
    if (!is.matrix(J) || nrow(J) != n_time || ncol(J) != n_params) {
      next
    }
    if (length(e) != n_time || any(!is.finite(e)) || any(!is.finite(J))) {
      next
    }

    bread_inv <- .invert_psd_safe(crossprod(J))
    if (is.null(bread_inv)) {
      next
    }
    Je <- J * e
    meat <- crossprod(Je)
    cov_theta <- bread_inv %*% meat %*% bread_inv
    d_theta <- diag(cov_theta)
    d_theta[!is.finite(d_theta) | d_theta < 0] <- 0
    se_theta_hat[v, ] <- sqrt(d_theta)

    hrf_vals <- hrf_interface$hrf_function(hrf_eval_times, theta_v)
    if (!is.numeric(hrf_vals) || any(!is.finite(hrf_vals))) {
      next
    }
    x_hrf <- .convolve_signal_with_kernels(
      signal = stim_signal,
      kernels = matrix(hrf_vals, ncol = 1),
      output_length = n_time
    )[, 1]
    if (any(!is.finite(x_hrf))) {
      next
    }

    X <- if (has_intercept) cbind(1, x_hrf) else matrix(x_hrf, ncol = 1)
    XtX_inv <- .invert_psd_safe(crossprod(X))
    if (is.null(XtX_inv)) {
      next
    }
    beta_hat <- XtX_inv %*% crossprod(X, y_v)
    resid_lm <- as.numeric(y_v - X %*% beta_hat)
    if (any(!is.finite(resid_lm))) {
      next
    }
    Xe <- X * resid_lm
    cov_beta <- XtX_inv %*% crossprod(Xe) %*% XtX_inv
    beta_idx <- if (has_intercept) 2L else 1L
    beta_var <- cov_beta[beta_idx, beta_idx]
    if (is.finite(beta_var) && beta_var >= 0) {
      se_beta0[v] <- sqrt(beta_var)
    }
  }

  list(
    se_theta_hat = se_theta_hat,
    se_beta0 = se_beta0
  )
}
