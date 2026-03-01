library(fmriparametric)

# Test that the Gauss-Newton refinement step can recover parameters
# for a simple single voxel toy example.

test_that(".gauss_newton_refinement improves fit on toy data", {
  set.seed(123)

  n_time <- 60
  onsets <- rep(0, n_time)
  onsets[c(10, 30, 50)] <- 1

  true_theta <- c(tau = 6, sigma = 2.5, rho = 0.4)
  hrf_times <- seq(0, 30, by = 0.5)
  bounds <- fmriparametric:::.lwu_hrf_default_bounds()
  true_hrf <- fmriparametric:::.lwu_hrf_function(hrf_times, true_theta, bounds)

  conv_full <- stats::convolve(onsets, rev(true_hrf), type = "open")
  signal <- conv_full[seq_len(n_time)]
  y <- signal + rnorm(n_time, sd = 0.01)

  Y_proj <- matrix(y, ncol = 1)
  S_target_proj <- matrix(onsets, ncol = 1)

  # Create proper HRF interface with bounds
  # Use the factory function which properly handles bounds
  interface <- fmriparametric:::.create_hrf_interface("lwu", user_bounds = bounds)
  
  # Use the actual default seed function
  theta_start <- fmriparametric:::.lwu_hrf_default_seed() + c(1, 0.5, 0.1)

  res <- fmriparametric:::.gauss_newton_refinement(
    theta_hat_voxel = matrix(theta_start, nrow = 1),
    r2_voxel = 0,
    Y_proj = Y_proj,
    S_target_proj = S_target_proj,
    scan_times = seq_len(n_time),
    hrf_eval_times = hrf_times,
    hrf_interface = interface,
    theta_bounds = bounds,
    queue_labels = "hard_GN",
    max_iter_gn = 10,
    verbose = FALSE
  )

  # Check that refinement improved the fit
  expect_true(res$r2 > 0.9)
  
  # Queue labels are only changed to "easy" if there was improvement
  # Since we started with r2=0 and should get high r2, it should change
  if (res$r2 > 0) {
    expect_equal(res$queue_labels, "easy")
  } else {
    expect_equal(res$queue_labels, "hard_GN")
  }
  
  # Check parameter recovery (allowing some tolerance)
  expect_true(all(abs(res$theta_hat - true_theta) < c(0.5, 0.5, 0.2)))
})

test_that(".gauss_newton_refinement parallel path matches serial results", {
  if (.Platform$OS.type == "windows") {
    testthat::skip("Parallel mclapply path is not supported on Windows")
  }
  cores_raw <- suppressWarnings(tryCatch(parallel::detectCores(logical = FALSE), error = function(...) 1L))
  cores <- suppressWarnings(as.integer(cores_raw[1]))
  if (!is.finite(cores) || cores < 2L) {
    testthat::skip("Need at least 2 physical cores for parallel comparison")
  }

  set.seed(42)
  n_time <- 80
  n_vox <- 6
  hrf_times <- seq(0, 30, by = 0.5)
  bounds <- fmriparametric:::.lwu_hrf_default_bounds()
  interface <- fmriparametric:::.create_hrf_interface("lwu", user_bounds = bounds)

  S_target_proj <- matrix(0, nrow = n_time, ncol = 1)
  S_target_proj[c(8, 22, 40, 58, 72), 1] <- 1
  signal <- S_target_proj[, 1]

  true_theta <- matrix(
    rep(c(6, 2.5, 0.35), n_vox),
    nrow = n_vox,
    byrow = TRUE
  )
  Y_proj <- matrix(0, nrow = n_time, ncol = n_vox)
  for (v in seq_len(n_vox)) {
    h <- interface$hrf_function(hrf_times, true_theta[v, ])
    x <- fmriparametric:::.fast_batch_convolution(signal, matrix(h, ncol = 1), n_time)[, 1]
    Y_proj[, v] <- x + rnorm(n_time, sd = 0.05)
  }

  theta_start <- matrix(
    rep(fmriparametric:::.lwu_hrf_default_seed() + c(2, 0.8, -0.1), n_vox),
    nrow = n_vox,
    byrow = TRUE
  )
  queue_labels <- rep("hard_GN", n_vox)
  r2_start <- rep(0.05, n_vox)

  serial <- fmriparametric:::.gauss_newton_refinement(
    theta_hat_voxel = theta_start,
    r2_voxel = r2_start,
    Y_proj = Y_proj,
    S_target_proj = S_target_proj,
    scan_times = seq_len(n_time),
    hrf_eval_times = hrf_times,
    hrf_interface = interface,
    theta_bounds = bounds,
    queue_labels = queue_labels,
    max_iter_gn = 6,
    use_conv_cache = TRUE,
    parallel = FALSE,
    parallel_min_voxels = 1L,
    verbose = FALSE
  )

  parallel_res <- fmriparametric:::.gauss_newton_refinement(
    theta_hat_voxel = theta_start,
    r2_voxel = r2_start,
    Y_proj = Y_proj,
    S_target_proj = S_target_proj,
    scan_times = seq_len(n_time),
    hrf_eval_times = hrf_times,
    hrf_interface = interface,
    theta_bounds = bounds,
    queue_labels = queue_labels,
    max_iter_gn = 6,
    use_conv_cache = TRUE,
    parallel = TRUE,
    n_cores = min(cores, 2L),
    parallel_min_voxels = 1L,
    verbose = FALSE
  )

  expect_equal(parallel_res$theta_hat, serial$theta_hat, tolerance = 1e-10)
  expect_equal(parallel_res$r2, serial$r2, tolerance = 1e-12)
  expect_equal(parallel_res$n_converged, serial$n_converged)
  expect_equal(parallel_res$n_improved, serial$n_improved)
  expect_equal(parallel_res$convergence_status, serial$convergence_status)
})
