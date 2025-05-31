library(testthat)
library(fmriparametric)

context("parametric engine")

# Dummy interface that returns a simple Taylor basis
hrf_iface <- list(
  taylor_basis = function(theta0, t_hrf) {
    matrix(rep(c(1, 0.1, 0.2, 0.05), each = length(t_hrf)), nrow = length(t_hrf), byrow = FALSE)
  }
)

Y <- matrix(rnorm(20), nrow = 10, ncol = 2)
S <- matrix(rbinom(10, 1, 0.2), ncol = 1)
scan_t <- seq_len(10)
t_hrf <- c(0, 1)

res <- .parametric_engine(
  Y_proj = Y,
  S_target_proj = S,
  scan_times = scan_t,
  hrf_eval_times = t_hrf,
  hrf_interface = hrf_iface,
  theta_seed = c(1, 1, 1),
  theta_bounds = list(lower = c(0, 0, 0), upper = c(2, 2, 2))
)

test_that("engine returns correct structure", {
  expect_type(res, "list")
  expect_true(all(c("theta_hat", "beta0") %in% names(res)))
  expect_equal(nrow(res$theta_hat), ncol(Y))
  expect_equal(ncol(res$theta_hat), 3)
  expect_length(res$beta0, ncol(Y))
})

# near-zero amplitude handling
Y_zero <- matrix(0, nrow = 10, ncol = 1)
res2 <- .parametric_engine(
  Y_proj = Y_zero,
  S_target_proj = S,
  scan_times = scan_t,
  hrf_eval_times = t_hrf,
  hrf_interface = hrf_iface,
  theta_seed = c(1, 1, 1),
  theta_bounds = list(lower = c(0, 0, 0), upper = c(2, 2, 2))
)

test_that("no NaNs when amplitude is zero", {
  expect_false(any(is.na(res2$theta_hat)))
  expect_false(any(is.na(res2$beta0)))
})
