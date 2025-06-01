library(testthat)

context("performance regression")

test_that("C++ ridge solver not slower than R implementation", {
  set.seed(1)
  X <- matrix(rnorm(200 * 4), 200, 4)
  Y <- matrix(rnorm(200 * 500), 200, 500)

  t_r <- system.time({
    qr_decomp <- qr(X)
    Q <- qr.Q(qr_decomp)
    R <- qr.R(qr_decomp)
    coeffs_r <- solve(R + 0.01 * diag(ncol(R)), t(Q) %*% Y)
  })["elapsed"]

  t_cpp <- system.time({
    coeffs_cpp <- .ridge_linear_solve(X, Y, 0.01)
  })["elapsed"]

  expect_true(all.equal(coeffs_r, coeffs_cpp, tolerance = 1e-8))
  expect_lte(t_cpp, t_r * 1.1)
})
