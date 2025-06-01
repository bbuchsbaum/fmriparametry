context("SIMD Taylor rho derivative")

test_that("rho parameter can be recovered with SIMD Taylor basis", {
  set.seed(123)
  n_time <- 60
  t_hrf <- seq(0, 30, length.out = 61)
  true_theta <- c(6, 2.5, 0.35)

  S <- matrix(0, nrow = n_time, ncol = 1)
  S[seq(10, n_time, by = 20), 1] <- 1

  true_hrf <- exp(-(t_hrf - true_theta[1])^2 / (2 * true_theta[2]^2)) -
    true_theta[3] * exp(-(t_hrf - true_theta[1] - 2 * true_theta[2])^2 /
                         (2 * (1.6 * true_theta[2])^2))
  conv_full <- convolve(S[, 1], rev(true_hrf), type = "open")
  Y <- matrix(conv_full[1:n_time] + rnorm(n_time, sd = 0.05), ncol = 1)

  taylor_basis_simd <- function(theta0, t) {
    z1 <- (t - theta0[1]) / theta0[2]
    h0 <- exp(-0.5 * z1 * z1)
    dh_dt1 <- h0 * z1 / theta0[2]
    dh_dt2 <- h0 * z1 * z1 / theta0[2]
    z_u <- (t - theta0[1] - 2 * theta0[2]) / (1.6 * theta0[2])
    dh_dt3 <- -exp(-0.5 * z_u * z_u)
    cbind(h0, dh_dt1, dh_dt2, dh_dt3)
  }

  hrf_iface <- list(
    taylor_basis = taylor_basis_simd,
    parameter_names = c("tau", "sigma", "rho")
  )

  fit <- .parametric_engine(
    Y_proj = Y,
    S_target_proj = S,
    scan_times = seq_len(n_time),
    hrf_eval_times = t_hrf,
    hrf_interface = hrf_iface,
    theta_seed = c(6, 2.5, 0.1),
    theta_bounds = list(lower = c(0, 0.05, 0), upper = c(20, 10, 1.5))
  )

  expect_lt(abs(fit$theta_hat[1, 3] - true_theta[3]), 0.1)
})
