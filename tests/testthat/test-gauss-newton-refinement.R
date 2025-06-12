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
  base_interface <- fmriparametric:::.get_hrf_interface("lwu")
  theta_bounds <- base_interface$default_bounds()
  interface <- list(
    hrf_function = function(t, params_vector, ...) {
      base_interface$hrf_function(t, params_vector, bounds = theta_bounds, ...)
    },
    taylor_basis = function(params_vector0, t_hrf_eval, ...) {
      base_interface$taylor_basis(params_vector0, t_hrf_eval, bounds = theta_bounds, ...)
    },
    parameter_names = base_interface$parameter_names,
    default_seed = base_interface$default_seed,
    default_bounds = base_interface$default_bounds
  )
  
  theta_start <- interface$default_seed() + c(1, 0.5, 0.1)

  res <- fmriparametric:::.gauss_newton_refinement(
    theta_hat_voxel = matrix(theta_start, nrow = 1),
    r2_voxel = 0,
    Y_proj = Y_proj,
    S_target_proj = S_target_proj,
    scan_times = seq_len(n_time),
    hrf_eval_times = hrf_times,
    hrf_interface = interface,
    theta_bounds = interface$default_bounds(),
    queue_labels = "hard_GN",
    max_iter_gn = 10,
    verbose = FALSE
  )

  expect_true(res$r2 > 0.9)
  expect_equal(res$queue_labels, "easy")
  expect_true(all(abs(res$theta_hat - true_theta) < c(0.5, 0.5, 0.2)))
})
