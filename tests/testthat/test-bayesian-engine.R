
test_that("posterior mean and variance match analytic results", {
  Y <- c(1, 2)
  X <- c(1, 1)
  priors <- list(mean = 0, var = 1, sigma2 = 1)

  result <- fmriparametric:::.bayesian_engine(Y, X, priors)

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

  expect_error(fmriparametric:::.bayesian_engine(Y, X, priors))
})

test_that("validation error for invalid prior parameters", {
  Y <- c(1, 2)
  X <- c(1, 1)
  priors_bad_var <- list(mean = 0, var = -1, sigma2 = 1)
  priors_bad_sigma2 <- list(mean = 0, var = 1, sigma2 = 0)

  expect_error(fmriparametric:::.bayesian_engine(Y, X, priors_bad_var))
  expect_error(fmriparametric:::.bayesian_engine(Y, X, priors_bad_sigma2))
})

test_that(".bayesian_engine validates inputs", {
  y <- 1:10
  x <- 1:10
  priors <- list(mean = 0, var = 1, sigma2 = 1)

  expect_error(fmriparametric:::.bayesian_engine("a", x, priors), "Y must be a numeric vector")
  expect_error(fmriparametric:::.bayesian_engine(y, c(1,2), priors), "equal length")
  expect_error(fmriparametric:::.bayesian_engine(c(1,2,Inf), c(1,2,3), priors), "finite values")
  expect_error(fmriparametric:::.bayesian_engine(y, x, list(mean=0, var=-1, sigma2=1)), "priors\\$var must be positive")
  expect_error(fmriparametric:::.bayesian_engine(y, x, priors, mcmc_samples = 0), "positive integer")

})
