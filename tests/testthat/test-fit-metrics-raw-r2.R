library(testthat)
library(fmriparametric)

test_that("fit metrics expose unclipped raw R-squared", {
  set.seed(41)
  y <- matrix(rnorm(120), nrow = 40, ncol = 3)
  y_bad <- -5 * y

  metrics <- fmriparametric:::.calculate_fit_metrics(
    y_true = y,
    y_pred = y_bad,
    n_predictors = 1,
    has_intercept = TRUE
  )

  expect_true(all(metrics$r_squared == 0))
  expect_true(all(metrics$r_squared_raw < 0))
})

test_that("raw and clipped R-squared agree on well-fitted models", {
  set.seed(42)
  y <- matrix(rnorm(180), nrow = 60, ncol = 3)
  y_hat <- y + matrix(rnorm(180, sd = 0.05), nrow = 60, ncol = 3)

  metrics <- fmriparametric:::.calculate_fit_metrics(
    y_true = y,
    y_pred = y_hat,
    n_predictors = 2,
    has_intercept = TRUE
  )

  expect_true(all(metrics$r_squared_raw <= 1 + 1e-10))
  expect_equal(metrics$r_squared, pmin(1, pmax(0, metrics$r_squared_raw)), tolerance = 1e-12)
})
