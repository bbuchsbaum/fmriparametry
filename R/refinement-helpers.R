# Refinement helper functions shared by Stage 4 refinement paths.

#' Refine moderate voxels using grouped local recentering
#'
#' Groups voxels with identical parameters and refines them together
#' using local Taylor expansion. This is more efficient than individual
#' voxel refinement.
#'
#' @noRd
.refine_moderate_voxels_grouped <- function(
  voxel_idx, Y_proj, S_target_proj, theta_current, amplitudes_current, r_squared,
  hrf_interface, hrf_eval_times, theta_bounds = NULL,
  lambda_ridge = 0.01, baseline_model = NULL
) {

  if (length(voxel_idx) == 0) {
    return(list(
      theta_refined = theta_current[voxel_idx, , drop = FALSE],
      amplitudes = amplitudes_current[voxel_idx],
      r_squared = r_squared[voxel_idx],
      n_improved = 0
    ))
  }

  bounds <- theta_bounds %||% hrf_interface$default_bounds()
  n_time <- nrow(Y_proj)
  n_params <- length(hrf_interface$parameter_names)

  theta_block <- theta_current[voxel_idx, , drop = FALSE]
  id <- apply(theta_block, 1, paste, collapse = ",")
  groups <- split(seq_along(id), id)

  theta_out <- theta_block
  amps_out <- amplitudes_current[voxel_idx]
  r2_out <- r_squared[voxel_idx]
  n_improved <- 0

  for (g in groups) {
    theta_v <- theta_block[g[1], ]

    basis <- hrf_interface$taylor_basis(theta_v, hrf_eval_times)
    if (!is.matrix(basis)) {
      basis <- matrix(basis, ncol = n_params + 1)
    }

    X <- .convolve_signal_with_kernels(S_target_proj[, 1], basis, n_time)
    has_intercept <- .is_intercept_baseline(baseline_model)
    if (has_intercept) {
      X <- cbind(1, X)
    }

    Y_block <- Y_proj[, voxel_idx[g], drop = FALSE]

    # Use centralized C++ ridge solver for consistency/performance.
    coeffs <- .ridge_linear_solve(X, Y_block, lambda_ridge)

    if (has_intercept) {
      beta0_new <- as.numeric(coeffs[2, ])
      coeff_start_idx <- 3
    } else {
      beta0_new <- as.numeric(coeffs[1, ])
      coeff_start_idx <- 2
    }

    delta <- matrix(0, length(beta0_new), n_params)
    valid <- abs(beta0_new) >= 1e-6
    if (any(valid)) {
      delta[valid, ] <- t(
        coeffs[coeff_start_idx:(coeff_start_idx + n_params - 1), valid, drop = FALSE]
      ) / beta0_new[valid]
    }

    theta_new <- sweep(delta, 2, theta_v, FUN = "+")
    theta_new <- .clamp_theta_to_bounds(theta_new, bounds)
    if (!is.matrix(theta_new)) {
      theta_new <- matrix(theta_new, nrow = length(beta0_new), ncol = n_params, byrow = TRUE)
    }

    fitted <- X %*% coeffs
    r2_new <- .compute_r_squared(Y_block, fitted, has_intercept = has_intercept)

    keep <- !is.na(r2_new) & r2_new > r_squared[voxel_idx[g]] & valid
    if (any(keep)) {
      theta_out[g[keep], ] <- theta_new[keep, , drop = FALSE]
      amps_out[g[keep]] <- beta0_new[keep]
      r2_out[g[keep]] <- r2_new[keep]
      n_improved <- n_improved + sum(keep)
    }
  }

  list(
    theta_refined = theta_out,
    amplitudes = amps_out,
    r_squared = r2_out,
    n_improved = n_improved
  )
}

#' Compute amplitudes for specific voxels
#'
#' Computes optimal amplitudes given HRF parameters.
#'
#' @noRd
.compute_amplitudes_for_voxels <- function(
  voxel_idx, theta_hat, Y_proj, S_target_proj,
  hrf_interface, hrf_eval_times
) {
  if (length(voxel_idx) == 0) {
    return(numeric(0))
  }

  n_time <- nrow(Y_proj)
  n_eval <- length(hrf_eval_times)
  signal <- S_target_proj[, 1]

  kernels <- vapply(
    voxel_idx,
    function(v) hrf_interface$hrf_function(hrf_eval_times, theta_hat[v, ]),
    FUN.VALUE = numeric(n_eval)
  )
  if (!is.matrix(kernels)) {
    kernels <- matrix(kernels, ncol = 1)
  }

  x_pred <- .convolve_signal_with_kernels(signal, kernels, n_time)
  y_block <- Y_proj[, voxel_idx, drop = FALSE]

  numer <- colSums(x_pred * y_block)
  denom <- colSums(x_pred^2)
  amps <- rep(0, length(voxel_idx))
  valid <- is.finite(denom) & denom >= 1e-10
  amps[valid] <- numer[valid] / denom[valid]
  amps
}
