#!/usr/bin/env Rscript

# Compact kernel benchmark for:
# 1) LWU Taylor basis
# 2) Gauss-Newton objective + Jacobian
# 3) Fit-metrics kernel
#
# Usage:
#   Rscript benchmarks/benchmark_kernel_triplet.R
# Optional env:
#   FMRIPARAMETRY_LIB=/tmp/Rlib_bench
#   BENCH_REPS=40
#   BENCH_OUT=benchmarks/results

lib_override <- Sys.getenv("FMRIPARAMETRY_LIB", "")
if (nzchar(lib_override)) {
  .libPaths(c(lib_override, .libPaths()))
}

if (!requireNamespace("fmriparametric", quietly = TRUE)) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop(
      "Package 'fmriparametric' is not installed and 'pkgload' is unavailable. ",
      "Install fmriparametric or install pkgload for source fallback."
    )
  }
  message("fmriparametric not installed in current library; using pkgload::load_all() fallback")
  pkgload::load_all(path = ".", quiet = TRUE, export_all = FALSE)
}

suppressPackageStartupMessages(library(fmriparametric))

bench_reps <- suppressWarnings(as.integer(Sys.getenv("BENCH_REPS", unset = "40")))
if (!is.finite(bench_reps) || bench_reps < 5L) {
  bench_reps <- 40L
}

set.seed(10101)

bench_elapsed <- function(fn, reps = 20L, warmup = 3L) {
  for (i in seq_len(warmup)) fn()
  n_inner <- calibrate_inner_loops(fn)
  gc()
  times <- numeric(reps)
  for (i in seq_len(reps)) {
    times[i] <- as.numeric(system.time({
      for (j in seq_len(n_inner)) {
        invisible(fn())
      }
    })[["elapsed"]]) / n_inner
  }
  c(
    inner_loops = n_inner,
    median_sec = stats::median(times),
    mean_sec = mean(times),
    sd_sec = stats::sd(times)
  )
}

calibrate_inner_loops <- function(fn, target_elapsed = 0.03, max_loops = 16384L) {
  loops <- 1L
  repeat {
    elapsed <- as.numeric(system.time({
      for (j in seq_len(loops)) {
        invisible(fn())
      }
    })[["elapsed"]])
    if (!is.finite(elapsed)) {
      return(1L)
    }
    if (elapsed >= target_elapsed || loops >= max_loops) {
      return(loops)
    }
    loops <- loops * 2L
  }
}

safe_speedup <- function(ref, cur) {
  if (!is.finite(ref) || !is.finite(cur) || cur <= 0) {
    return(NA_real_)
  }
  ref / cur
}

or_else <- function(x, y) {
  if (is.null(x)) y else x
}

hrf_interface <- fmriparametric:::.create_hrf_interface("lwu")
bounds <- hrf_interface$active_bounds
t_hrf <- seq(0, 30, length.out = 61)
theta <- c(6.2, 2.4, 0.35)

# --- LWU Taylor basis ----------------------------------------------------------

fd_taylor_reference <- function(theta0, t_eval, b) {
  theta0 <- fmriparametric:::.apply_bounds(theta0, b, epsilon = c(1e-6, 1e-6, 1e-6))
  theta0[2] <- max(theta0[2], 0.051)
  basis <- matrix(0, nrow = length(t_eval), ncol = 4)

  hrf_fn <- function(p) {
    exp(-((t_eval - p[1])^2) / (2 * p[2]^2)) -
      p[3] * exp(-((t_eval - (p[1] + 2 * p[2]))^2) / (2 * (1.6 * p[2])^2))
  }

  basis[, 1] <- hrf_fn(theta0)
  step_scale <- pmax(abs(theta0), c(1, 0.2, 0.1))
  step_size <- pmax(1e-4 * step_scale, 1e-5)

  for (k in seq_len(3)) {
    plus <- theta0
    minus <- theta0
    plus[k] <- min(theta0[k] + step_size[k], b$upper[k] - 1e-8)
    minus[k] <- max(theta0[k] - step_size[k], b$lower[k] + 1e-8)
    if (k == 2) {
      plus[2] <- max(plus[2], 0.051)
      minus[2] <- max(minus[2], 0.051)
    }
    denom <- plus[k] - minus[k]
    if (denom > 1e-10) {
      basis[, k + 1] <- (hrf_fn(plus) - hrf_fn(minus)) / denom
    }
  }

  basis
}

basis_fast_exists <- exists(
  ".lwu_hrf_taylor_basis_fast",
  envir = asNamespace("fmriparametric"),
  mode = "function",
  inherits = FALSE
)
basis_new_fn <- if (basis_fast_exists) {
  function() fmriparametric:::.lwu_hrf_taylor_basis_fast(theta, t_hrf, bounds)
} else {
  function() fmriparametric:::.lwu_hrf_taylor_basis_function(theta, t_hrf, bounds)
}
basis_ref_fn <- function() {
  fd_taylor_reference(theta, t_hrf, bounds)
}
basis_new <- basis_new_fn()
basis_ref <- basis_ref_fn()

basis_new_t <- bench_elapsed(basis_new_fn, reps = bench_reps)
basis_ref_t <- bench_elapsed(basis_ref_fn, reps = bench_reps)
basis_err <- max(abs(basis_new - basis_ref))

# --- GN objective + Jacobian ---------------------------------------------------

n_time <- 180L
signal <- numeric(n_time)
signal[seq(10, n_time - 10, by = 18)] <- 1
h <- hrf_interface$hrf_function(t_hrf, theta)
x <- fmriparametric:::.convolve_signal_with_kernels(signal, matrix(h, ncol = 1), n_time)[, 1]
y <- 0.8 * x + rnorm(n_time, sd = 0.05)

legacy_objective <- function(theta0, y0, s0, t0, iface, n_t) {
  hrf_vals <- iface$hrf_function(t0, theta0)
  x_pred <- fmriparametric:::.convolve_signal_with_kernels(s0, matrix(hrf_vals, ncol = 1), n_t)[, 1]
  denom <- sum(x_pred^2)
  if (!is.finite(denom) || denom < 1e-8) {
    return(Inf)
  }
  beta <- as.numeric(crossprod(x_pred, y0)) / denom
  sum((y0 - beta * x_pred)^2)
}

legacy_jacobian <- function(theta0, y0, s0, t0, iface, n_t, b) {
  basis <- fd_taylor_reference(theta0, t0, b)
  X_conv <- fmriparametric:::.convolve_signal_with_kernels(s0, basis, n_t)
  x_hrf <- X_conv[, 1]
  denom <- sum(x_hrf^2)
  if (!is.finite(denom) || denom < 1e-8) {
    return(NULL)
  }
  beta <- as.numeric(crossprod(x_hrf, y0)) / denom
  residuals <- y0 - beta * x_hrf

  jac <- matrix(0, nrow = n_t, ncol = 3)
  for (k in seq_len(3)) {
    dx <- X_conv[, k + 1]
    dbeta <- (as.numeric(crossprod(dx, y0)) - 2 * beta * as.numeric(crossprod(x_hrf, dx))) / denom
    jac[, k] <- -beta * dx - dbeta * x_hrf
  }
  list(jacobian = jac, residuals = residuals, amplitude = beta)
}

gn_obj_new_fn <- function() {
  fmriparametric:::.calculate_objective_gn(theta, y, signal, t_hrf, hrf_interface, n_time)
}
gn_obj_ref_fn <- function() {
  legacy_objective(theta, y, signal, t_hrf, hrf_interface, n_time)
}

gn_jac_new_fn <- function() {
  fmriparametric:::.get_jacobian_and_residuals(theta, y, signal, t_hrf, hrf_interface, n_time)
}
gn_jac_ref_fn <- function() {
  legacy_jacobian(theta, y, signal, t_hrf, hrf_interface, n_time, bounds)
}

gn_obj_new <- gn_obj_new_fn()
gn_obj_ref <- gn_obj_ref_fn()
gn_jac_new <- gn_jac_new_fn()
gn_jac_ref <- gn_jac_ref_fn()

gn_obj_new_t <- bench_elapsed(gn_obj_new_fn, reps = bench_reps)
gn_obj_ref_t <- bench_elapsed(gn_obj_ref_fn, reps = bench_reps)
gn_jac_new_t <- bench_elapsed(gn_jac_new_fn, reps = bench_reps)
gn_jac_ref_t <- bench_elapsed(gn_jac_ref_fn, reps = bench_reps)

gn_obj_err <- abs(gn_obj_new - gn_obj_ref)
gn_jac_err <- if (is.null(gn_jac_new) || is.null(gn_jac_ref)) {
  NA_real_
} else {
  max(abs(gn_jac_new$jacobian - gn_jac_ref$jacobian))
}
gn_resid_err <- if (is.null(gn_jac_new) || is.null(gn_jac_ref)) {
  NA_real_
} else {
  max(abs(gn_jac_new$residuals - gn_jac_ref$residuals))
}

# --- Fit metrics kernel --------------------------------------------------------

n_obs <- 240L
n_vox <- 300L
y_true <- matrix(rnorm(n_obs * n_vox), nrow = n_obs, ncol = n_vox)
y_pred <- y_true + matrix(rnorm(n_obs * n_vox, sd = 0.2), nrow = n_obs, ncol = n_vox)

fit_metrics_reference <- function(y_t, y_p, has_intercept = TRUE, tol = 1e-10) {
  n_t <- nrow(y_t)
  n_v <- ncol(y_t)
  r2 <- numeric(n_v)
  r2_raw <- numeric(n_v)
  rss <- numeric(n_v)
  tss <- numeric(n_v)
  mse <- numeric(n_v)
  rmse <- numeric(n_v)
  mae <- numeric(n_v)

  for (v in seq_len(n_v)) {
    yv <- y_t[, v]
    yp <- y_p[, v]
    res <- yv - yp
    rss_v <- sum(res^2)
    tss_v <- if (has_intercept) {
      sum((yv - mean(yv))^2)
    } else {
      sum(yv^2)
    }
    if (abs(tss_v) < tol) {
      r2_raw_v <- if (abs(rss_v) < tol) 1 else 0
    } else {
      r2_raw_v <- 1 - rss_v / tss_v
    }
    r2[v] <- min(1, max(0, r2_raw_v))
    r2_raw[v] <- r2_raw_v
    rss[v] <- rss_v
    tss[v] <- tss_v
    mse[v] <- rss_v / n_t
    rmse[v] <- sqrt(mse[v])
    mae[v] <- mean(abs(res))
  }

  list(r_squared = r2, r_squared_raw = r2_raw, rss = rss, tss = tss, mse = mse, rmse = rmse, mae = mae)
}

fit_new_fn <- function() {
  fmriparametric:::.calculate_fit_metrics(
    y_true = y_true,
    y_pred = y_pred,
    n_predictors = 2L,
    has_intercept = TRUE
  )
}
fit_ref_fn <- function() {
  fit_metrics_reference(y_true, y_pred, has_intercept = TRUE)
}

fit_new <- fit_new_fn()
fit_ref <- fit_ref_fn()
fit_new_t <- bench_elapsed(fit_new_fn, reps = max(10L, bench_reps %/% 2L))
fit_ref_t <- bench_elapsed(fit_ref_fn, reps = max(5L, bench_reps %/% 4L))

fit_r2_err <- max(abs(fit_new$r_squared - fit_ref$r_squared))
fit_new_r2_raw <- or_else(fit_new$r_squared_raw, fit_new$r_squared)
fit_ref_r2_raw <- or_else(fit_ref$r_squared_raw, fit_ref$r_squared)
fit_r2_raw_err <- max(abs(fit_new_r2_raw - fit_ref_r2_raw))

# --- Report -------------------------------------------------------------------

benchmark_rows <- rbind(
  data.frame(
    kernel = "lwu_taylor_basis",
    impl = c("current", "reference_fd"),
    inner_loops = c(basis_new_t[["inner_loops"]], basis_ref_t[["inner_loops"]]),
    median_sec = c(basis_new_t[["median_sec"]], basis_ref_t[["median_sec"]]),
    mean_sec = c(basis_new_t[["mean_sec"]], basis_ref_t[["mean_sec"]]),
    sd_sec = c(basis_new_t[["sd_sec"]], basis_ref_t[["sd_sec"]]),
    speedup_vs_reference = c(safe_speedup(basis_ref_t[["median_sec"]], basis_new_t[["median_sec"]]), 1),
    max_abs_error = c(basis_err, basis_err),
    stringsAsFactors = FALSE
  ),
  data.frame(
    kernel = "gn_objective",
    impl = c("current", "reference_r"),
    inner_loops = c(gn_obj_new_t[["inner_loops"]], gn_obj_ref_t[["inner_loops"]]),
    median_sec = c(gn_obj_new_t[["median_sec"]], gn_obj_ref_t[["median_sec"]]),
    mean_sec = c(gn_obj_new_t[["mean_sec"]], gn_obj_ref_t[["mean_sec"]]),
    sd_sec = c(gn_obj_new_t[["sd_sec"]], gn_obj_ref_t[["sd_sec"]]),
    speedup_vs_reference = c(safe_speedup(gn_obj_ref_t[["median_sec"]], gn_obj_new_t[["median_sec"]]), 1),
    max_abs_error = c(gn_obj_err, gn_obj_err),
    stringsAsFactors = FALSE
  ),
  data.frame(
    kernel = "gn_jacobian",
    impl = c("current", "reference_r"),
    inner_loops = c(gn_jac_new_t[["inner_loops"]], gn_jac_ref_t[["inner_loops"]]),
    median_sec = c(gn_jac_new_t[["median_sec"]], gn_jac_ref_t[["median_sec"]]),
    mean_sec = c(gn_jac_new_t[["mean_sec"]], gn_jac_ref_t[["mean_sec"]]),
    sd_sec = c(gn_jac_new_t[["sd_sec"]], gn_jac_ref_t[["sd_sec"]]),
    speedup_vs_reference = c(safe_speedup(gn_jac_ref_t[["median_sec"]], gn_jac_new_t[["median_sec"]]), 1),
    max_abs_error = c(gn_jac_err, gn_jac_err),
    stringsAsFactors = FALSE
  ),
  data.frame(
    kernel = "fit_metrics",
    impl = c("current", "reference_r"),
    inner_loops = c(fit_new_t[["inner_loops"]], fit_ref_t[["inner_loops"]]),
    median_sec = c(fit_new_t[["median_sec"]], fit_ref_t[["median_sec"]]),
    mean_sec = c(fit_new_t[["mean_sec"]], fit_ref_t[["mean_sec"]]),
    sd_sec = c(fit_new_t[["sd_sec"]], fit_ref_t[["sd_sec"]]),
    speedup_vs_reference = c(safe_speedup(fit_ref_t[["median_sec"]], fit_new_t[["median_sec"]]), 1),
    max_abs_error = c(max(fit_r2_err, fit_r2_raw_err), max(fit_r2_err, fit_r2_raw_err)),
    stringsAsFactors = FALSE
  )
)

git_sha <- tryCatch(
  system("git rev-parse --short HEAD", intern = TRUE),
  error = function(...) NA_character_
)
if (length(git_sha) == 0L) {
  git_sha <- NA_character_
}
timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
out_dir <- Sys.getenv("BENCH_OUT", unset = "benchmarks/results")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

meta <- data.frame(
  timestamp = timestamp,
  git_sha = git_sha[1],
  reps = bench_reps,
  stringsAsFactors = FALSE
)

out_file <- file.path(out_dir, paste0("kernel_triplet_", timestamp, ".csv"))
latest_file <- file.path(out_dir, "kernel_triplet_latest.csv")
utils::write.csv(cbind(meta[rep(1, nrow(benchmark_rows)), ], benchmark_rows), out_file, row.names = FALSE)
utils::write.csv(cbind(meta[rep(1, nrow(benchmark_rows)), ], benchmark_rows), latest_file, row.names = FALSE)

cat("Kernel Triplet Benchmark\n")
cat("========================\n")
cat(sprintf("timestamp: %s\n", timestamp))
cat(sprintf("git_sha:   %s\n", git_sha[1]))
cat(sprintf("reps:      %d\n\n", bench_reps))

print(benchmark_rows, row.names = FALSE, digits = 6)
cat("\nAdditional accuracy checks:\n")
cat(sprintf("  gn_residual_max_abs_diff: %0.6e\n", gn_resid_err))
cat(sprintf("  fit_r2_max_abs_diff:      %0.6e\n", fit_r2_err))
cat(sprintf("  fit_r2_raw_max_abs_diff:  %0.6e\n", fit_r2_raw_err))
cat(sprintf("\nSaved: %s\n", out_file))
cat(sprintf("Saved: %s\n", latest_file))
