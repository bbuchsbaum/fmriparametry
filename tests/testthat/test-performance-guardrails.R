library(testthat)
library(fmriparametric)

test_that("analytic LWU Taylor basis remains faster than finite-difference reference", {
  skip_if(
    !identical(Sys.getenv("FMRIPARAMETRIC_RUN_PERF", unset = "0"), "1"),
    "Set FMRIPARAMETRIC_RUN_PERF=1 to run performance guardrails"
  )

  t <- seq(0, 30, by = 0.1)
  bounds <- fmriparametric:::.lwu_hrf_default_bounds()
  theta <- c(6.2, 2.3, 0.4)

  fd_reference <- function(theta0, t_eval, b) {
    theta0 <- fmriparametric:::.apply_bounds(theta0, b, epsilon = c(1e-6, 1e-6, 1e-6))
    theta0[2] <- max(theta0[2], 0.051)

    hrf <- function(p) {
      exp(-((t_eval - p[1])^2) / (2 * p[2]^2)) -
        p[3] * exp(-((t_eval - (p[1] + 2 * p[2]))^2) / (2 * (1.6 * p[2])^2))
    }

    basis <- matrix(0, nrow = length(t_eval), ncol = 4)
    basis[, 1] <- hrf(theta0)
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
        basis[, k + 1] <- (hrf(plus) - hrf(minus)) / denom
      } else {
        basis[, k + 1] <- 0
      }
    }
    basis
  }

  n_iter <- 250L
  fast_elapsed <- system.time({
    for (i in seq_len(n_iter)) {
      fmriparametric:::.lwu_hrf_taylor_basis_fast(theta, t, bounds)
    }
  })[["elapsed"]]
  ref_elapsed <- system.time({
    for (i in seq_len(n_iter)) {
      fd_reference(theta, t, bounds)
    }
  })[["elapsed"]]

  expect_lt(fast_elapsed, ref_elapsed * 0.9)
})
