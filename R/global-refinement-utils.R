# Utilities for Stage 3 global refinement and shared bounds handling.

# Resolve bounds from interface/config without forcing defaults.
# Returns NULL when no explicit bounds are available.
# @noRd
.resolve_optional_theta_bounds <- function(hrf_interface, config = NULL) {
  if (!is.null(hrf_interface$active_bounds)) {
    return(hrf_interface$active_bounds)
  }
  if (!is.null(config) && !is.null(config$theta_bounds)) {
    return(config$theta_bounds)
  }
  NULL
}

# Apply bounds with optional interior epsilon padding.
# @noRd
.clamp_theta_to_bounds <- function(theta, theta_bounds, epsilon = NULL) {
  if (is.null(theta_bounds)) {
    return(theta)
  }
  .apply_bounds(theta, theta_bounds, epsilon = epsilon)
}

# Identify rows that are not pinned at parameter bounds.
# @noRd
.rows_strictly_inside_bounds <- function(theta, theta_bounds, eps = 1e-6) {
  if (is.null(theta_bounds)) {
    return(rep(TRUE, nrow(theta)))
  }

  lower <- matrix(
    theta_bounds$lower + eps,
    nrow = nrow(theta),
    ncol = ncol(theta),
    byrow = TRUE
  )
  upper <- matrix(
    theta_bounds$upper - eps,
    nrow = nrow(theta),
    ncol = ncol(theta),
    byrow = TRUE
  )

  rowSums(theta <= lower | theta >= upper) == 0
}

# Compute a stable global recentering seed.
# @noRd
.compute_global_refinement_center <- function(
  theta_current,
  default_seed,
  theta_bounds = NULL,
  min_default_pull = 0.25,
  boundary_eps = 1e-6
) {
  theta_median <- apply(theta_current, 2, stats::median)
  if (is.null(theta_bounds)) {
    return(theta_median)
  }

  lower <- matrix(
    theta_bounds$lower + boundary_eps,
    nrow = nrow(theta_current),
    ncol = ncol(theta_current),
    byrow = TRUE
  )
  upper <- matrix(
    theta_bounds$upper - boundary_eps,
    nrow = nrow(theta_current),
    ncol = ncol(theta_current),
    byrow = TRUE
  )

  on_boundary <- theta_current <= lower | theta_current >= upper
  boundary_fraction <- mean(on_boundary)
  weight_default <- min(1, max(min_default_pull, 2 * boundary_fraction))
  theta_center <- (1 - weight_default) * theta_median + weight_default * default_seed

  near_lower <- colMeans(theta_current <= lower)
  near_upper <- colMeans(theta_current >= upper)
  boundary_locked <- (near_lower > 0.5) | (near_upper > 0.5)
  if (any(boundary_locked)) {
    theta_center[boundary_locked] <- default_seed[boundary_locked]
  }

  .clamp_theta_to_bounds(theta_center, theta_bounds, epsilon = boundary_eps)
}

# Select voxel-wise updates for global refinement.
# @noRd
.select_global_refinement_updates <- function(
  theta_prev,
  theta_candidate,
  r_squared_prev,
  r_squared_candidate,
  default_seed,
  theta_bounds = NULL,
  r2_tol = 1e-5
) {
  seed_mat <- matrix(
    default_seed,
    nrow = nrow(theta_prev),
    ncol = ncol(theta_prev),
    byrow = TRUE
  )
  closer_to_default <- rowSums((theta_candidate - seed_mat)^2) < rowSums((theta_prev - seed_mat)^2)
  interior_candidate <- .rows_strictly_inside_bounds(theta_candidate, theta_bounds)

  improved <- is.finite(r_squared_candidate) & (
    (r_squared_candidate > (r_squared_prev + r2_tol)) |
      ((r_squared_candidate >= (r_squared_prev - r2_tol)) & closer_to_default & interior_candidate)
  )

  which(improved)
}
