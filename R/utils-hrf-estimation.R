# Shared utilities for HRF estimation
# These functions are used by both legacy and refactored implementations


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