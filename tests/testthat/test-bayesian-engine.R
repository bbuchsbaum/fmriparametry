context(".bayesian_engine")

test_that("posterior mean and variance match analytic results", {
  Y <- c(1, 2)
  X <- c(1, 1)
  priors <- list(mean = 0, var = 1, sigma2 = 1)

  result <- .bayesian_engine(Y, X, priors)

  expected_var <- 1 / (1/priors$var + sum(X * X) / priors$sigma2)
  expected_mean <- expected_var * (priors$mean / priors$var + sum(X * Y) / priors$sigma2)

  expect_equal(result$posterior_mean, expected_mean, tolerance = 1e-8)
  expect_equal(result$posterior_sd^2, expected_var, tolerance = 1e-8)
  expect_length(result$samples, 1000)
})

test_that("validation error for mismatched input lengths", {
  Y <- c(1, 2, 3)
  X <- c(1, 1)
  priors <- list(mean = 0, var = 1, sigma2 = 1)

  expect_error(.bayesian_engine(Y, X, priors))
})

test_that("validation error for invalid prior parameters", {
  Y <- c(1, 2)
  X <- c(1, 1)
  priors_bad_var <- list(mean = 0, var = -1, sigma2 = 1)
  priors_bad_sigma2 <- list(mean = 0, var = 1, sigma2 = 0)

  expect_error(.bayesian_engine(Y, X, priors_bad_var))
  expect_error(.bayesian_engine(Y, X, priors_bad_sigma2))
})
