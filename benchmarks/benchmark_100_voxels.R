#!/usr/bin/env Rscript

# Optional library path override for benchmarking installed builds:
#   FMRIPARAMETRY_LIB=/tmp/Rlib_bench Rscript benchmarks/benchmark_100_voxels.R
lib_override <- Sys.getenv("FMRIPARAMETRY_LIB", "")
if (nzchar(lib_override)) {
  .libPaths(c(lib_override, .libPaths()))
}

if (!requireNamespace("fmriparametric", quietly = TRUE)) {
  stop(
    "Package 'fmriparametric' is not installed in current library paths. ",
    "Install it first (e.g., `R CMD INSTALL -l /tmp/Rlib_bench .`) or set ",
    "FMRIPARAMETRY_LIB to a library path containing the package."
  )
}

suppressPackageStartupMessages(library(fmriparametric))

set.seed(42)

n_time <- 240L
n_vox <- 100L
t_hrf <- seq(0, 30, length.out = 61)

hrf_interface <- fmriparametric:::.create_hrf_interface("lwu")
bounds <- hrf_interface$active_bounds

theta <- matrix(
  rep(fmriparametric:::.lwu_hrf_default_seed(), n_vox),
  nrow = n_vox,
  byrow = TRUE
)

S <- matrix(0, nrow = n_time, ncol = 1)
S[seq(10, n_time - 20, by = 20), 1] <- 1
signal <- S[, 1]

Y <- matrix(rnorm(n_time * n_vox, sd = 0.2), nrow = n_time)
for (v in seq_len(n_vox)) {
  h <- hrf_interface$hrf_function(t_hrf, theta[v, ])
  x <- fmriparametric:::.fast_batch_convolution(signal, matrix(h, ncol = 1), n_time)[, 1]
  Y[, v] <- Y[, v] + 0.8 * x
}

y1 <- Y[, 1]
th1 <- theta[1, ]

legacy_taylor_basis <- function(theta0, t_eval, bounds) {
  theta0 <- fmriparametric:::.apply_bounds(theta0, bounds, epsilon = c(0.01, 0.01, 0.01))
  names(theta0) <- c("tau", "sigma", "rho")
  basis <- fmrihrf::hrf_basis_lwu(theta0 = theta0, t = t_eval, normalize_primary = "none")
  if (!is.matrix(basis)) {
    basis <- matrix(basis, ncol = 4)
  }
  storage.mode(basis) <- "double"
  basis
}

legacy_objective <- function(theta, y, signal, t_hrf, hrf_interface, n_time) {
  hrf_vals <- hrf_interface$hrf_function(t_hrf, theta)
  x_pred_raw <- fmriparametric:::.convolve_signal_with_kernels(
    signal, matrix(hrf_vals, ncol = 1), n_time
  )[, 1]
  denom <- sum(x_pred_raw^2)
  if (denom < 1e-8) {
    return(Inf)
  }
  beta <- as.numeric(crossprod(x_pred_raw, y)) / denom
  sum((y - beta * x_pred_raw)^2)
}

legacy_jacobian <- function(theta, y, signal, t_hrf, hrf_interface, n_time, bounds) {
  n_params <- length(theta)
  taylor_basis <- legacy_taylor_basis(theta, t_hrf, bounds)
  X_conv <- fmriparametric:::.convolve_signal_with_kernels(signal, taylor_basis, n_time)

  x_hrf <- X_conv[, 1]
  denom <- sum(x_hrf^2)
  if (denom < 1e-8) {
    return(NULL)
  }

  beta <- as.numeric(crossprod(x_hrf, y)) / denom
  residuals <- y - beta * x_hrf

  jacobian <- matrix(0, nrow = n_time, ncol = n_params)
  for (k in seq_len(n_params)) {
    dx <- X_conv[, k + 1]
    dbeta <- (as.numeric(crossprod(dx, y)) - 2 * beta * as.numeric(crossprod(x_hrf, dx))) / denom
    jacobian[, k] <- -beta * dx - dbeta * x_hrf
  }

  list(jacobian = jacobian, residuals = residuals, amplitude = beta)
}

legacy_amplitudes <- function(voxel_idx, theta_hat, Y_proj, signal, hrf_interface, t_hrf, n_time) {
  out <- numeric(length(voxel_idx))
  for (i in seq_along(voxel_idx)) {
    v <- voxel_idx[i]
    h <- hrf_interface$hrf_function(t_hrf, theta_hat[v, ])
    x <- fmriparametric:::.convolve_signal_with_kernels(signal, matrix(h, ncol = 1), n_time)[, 1]
    denom <- sum(x^2)
    if (is.finite(denom) && denom >= 1e-10) {
      out[i] <- sum(x * Y_proj[, v]) / denom
    }
  }
  out
}

bench <- function(fn, n) {
  gc()
  as.numeric(system.time(for (i in seq_len(n)) fn())[3]) / n
}

obj_old <- function() legacy_objective(th1, y1, signal, t_hrf, hrf_interface, n_time)
obj_new <- function() fmriparametric:::.calculate_objective_gn(th1, y1, signal, t_hrf, hrf_interface, n_time)
jac_old <- function() legacy_jacobian(th1, y1, signal, t_hrf, hrf_interface, n_time, bounds)
jac_new <- function() fmriparametric:::.get_jacobian_and_residuals(th1, y1, signal, t_hrf, hrf_interface, n_time)
amp_old <- function() legacy_amplitudes(seq_len(n_vox), theta, Y, signal, hrf_interface, t_hrf, n_time)
amp_new <- function() fmriparametric:::.compute_amplitudes_for_voxels(seq_len(n_vox), theta, Y, S, hrf_interface, t_hrf)

obj_old_t <- bench(obj_old, 120L)
obj_new_t <- bench(obj_new, 120L)
jac_old_t <- bench(jac_old, 20L)
jac_new_t <- bench(jac_new, 120L)
amp_old_t <- bench(amp_old, 40L)
amp_new_t <- bench(amp_new, 40L)

obj_old_v <- obj_old()
obj_new_v <- obj_new()
jac_old_v <- jac_old()
jac_new_v <- jac_new()
amp_old_v <- amp_old()
amp_new_v <- amp_new()

cat("100-voxel benchmark\n")
cat("===================\n")
cat(sprintf(
  "Objective: old=%.6f s new=%.6f s speedup=%.2fx diff=%.3e\n",
  obj_old_t, obj_new_t, obj_old_t / obj_new_t, abs(obj_old_v - obj_new_v)
))
cat(sprintf(
  "Jacobian:  old=%.6f s new=%.6f s speedup=%.2fx jac_diff=%.3e resid_diff=%.3e amp_diff=%.3e\n",
  jac_old_t, jac_new_t, jac_old_t / jac_new_t,
  max(abs(jac_old_v$jacobian - jac_new_v$jacobian)),
  max(abs(jac_old_v$residuals - jac_new_v$residuals)),
  abs(jac_old_v$amplitude - jac_new_v$amplitude)
))
cat(sprintf(
  "Amplitude: old=%.6f s new=%.6f s speedup=%.2fx diff=%.3e\n",
  amp_old_t, amp_new_t, amp_old_t / amp_new_t, max(abs(amp_old_v - amp_new_v))
))

# -------------------------------------------------------------------------
# GN cache benchmark (new path only): cached vs uncached convolution context

theta_start <- theta + matrix(
  rep(c(3.0, 1.2, -0.15), n_vox),
  nrow = n_vox,
  byrow = TRUE
)
r2_start <- rep(0.05, n_vox)
queue_labels <- rep("hard_GN", n_vox)

run_gn <- function(use_cache, use_parallel = FALSE, n_cores = NULL) {
  fmriparametric:::.gauss_newton_refinement(
    theta_hat_voxel = theta_start,
    r2_voxel = r2_start,
    Y_proj = Y,
    S_target_proj = S,
    scan_times = seq_len(n_time),
    hrf_eval_times = t_hrf,
    hrf_interface = hrf_interface,
    theta_bounds = bounds,
    queue_labels = queue_labels,
    max_iter_gn = 6,
    tol_gn = 1e-4,
    lambda_ridge = 0.01,
    step_size = 1.0,
    verbose = FALSE,
    use_conv_cache = use_cache,
    parallel = use_parallel,
    n_cores = n_cores,
    parallel_min_voxels = 1L
  )
}

gn_uncached_t <- bench(function() run_gn(FALSE), 6L)
gn_cached_t <- bench(function() run_gn(TRUE), 6L)

gn_uncached <- run_gn(FALSE)
gn_cached <- run_gn(TRUE)

cat(sprintf(
  "GN refine: uncached=%.6f s cached=%.6f s speedup=%.2fx theta_diff=%.3e r2_diff=%.3e\n",
  gn_uncached_t,
  gn_cached_t,
  gn_uncached_t / gn_cached_t,
  max(abs(gn_uncached$theta_hat - gn_cached$theta_hat)),
  max(abs(gn_uncached$r2 - gn_cached$r2))
))

gn_parallel_cores <- suppressWarnings(tryCatch(parallel::detectCores(logical = FALSE), error = function(...) 1L))
if (.Platform$OS.type != "windows") {
  if (!is.finite(gn_parallel_cores) || gn_parallel_cores < 2L) {
    gn_parallel_cores <- 2L
  } else {
    gn_parallel_cores <- min(gn_parallel_cores, 4L)
  }
  gn_parallel_ok <- TRUE
  gn_cached_parallel_t <- tryCatch(
    bench(function() run_gn(TRUE, use_parallel = TRUE, n_cores = gn_parallel_cores), 6L),
    error = function(...) {
      gn_parallel_ok <<- FALSE
      NA_real_
    }
  )
  if (gn_parallel_ok) {
    gn_cached_parallel <- run_gn(TRUE, use_parallel = TRUE, n_cores = gn_parallel_cores)
    cat(sprintf(
      "GN parallel: cached_serial=%.6f s cached_parallel=%.6f s speedup=%.2fx theta_diff=%.3e r2_diff=%.3e (cores=%d)\n",
      gn_cached_t,
      gn_cached_parallel_t,
      gn_cached_t / gn_cached_parallel_t,
      max(abs(gn_cached$theta_hat - gn_cached_parallel$theta_hat)),
      max(abs(gn_cached$r2 - gn_cached_parallel$r2)),
      gn_parallel_cores
    ))
  } else {
    cat("GN parallel: attempted but unavailable in this runtime\n")
  }
} else {
  cat("GN parallel: skipped (unsupported platform)\n")
}

# -------------------------------------------------------------------------
# Public interface benchmark: direct fit vs precomputed design fit

direct_fit <- function() {
  estimate_parametric_hrf(
    fmri_data = Y,
    event_model = S,
    hrf_eval_times = t_hrf,
    global_refinement = FALSE,
    compute_se = FALSE,
    progress = FALSE,
    verbose = FALSE
  )
}

design_obj <- create_parametric_design(
  fmri_data = Y,
  event_model = S,
  hrf_eval_times = t_hrf
)

from_design_fit <- function() {
  estimate_parametric_hrf_from_design(
    design = design_obj,
    global_refinement = FALSE,
    compute_se = FALSE,
    progress = FALSE,
    verbose = FALSE
  )
}

# Warm-up and repeated timings for stable medians
invisible(direct_fit())
invisible(from_design_fit())

reps <- 10L
direct_times <- numeric(reps)
from_design_times <- numeric(reps)
for (i in seq_len(reps)) {
  direct_times[i] <- as.numeric(system.time(invisible(direct_fit()))[3])
  from_design_times[i] <- as.numeric(system.time(invisible(from_design_fit()))[3])
}

fit_direct <- direct_fit()
fit_from_design <- from_design_fit()

cat(sprintf(
  "Interface: direct_med=%.6f s from_design_med=%.6f s speedup=%.2fx param_diff=%.3e amp_diff=%.3e r2_diff=%.3e\n",
  median(direct_times),
  median(from_design_times),
  median(direct_times) / median(from_design_times),
  max(abs(get_parameters(fit_direct) - get_parameters(fit_from_design))),
  max(abs(fit_direct$amplitudes - fit_from_design$amplitudes)),
  max(abs(get_gof_per_voxel(fit_direct) - get_gof_per_voxel(fit_from_design)))
))
