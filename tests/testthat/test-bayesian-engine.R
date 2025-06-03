test_that(".bayesian_engine validates inputs", {
  y <- 1:10
  x <- 1:10
  priors <- list(mean = 0, var = 1, sigma2 = 1)

  expect_error(.bayesian_engine("a", x, priors), "Y must be a numeric vector")
  expect_error(.bayesian_engine(y, c(1,2), priors), "equal length")
  expect_error(.bayesian_engine(c(1,2,Inf), c(1,2,3), priors), "finite values")
  expect_error(.bayesian_engine(y, x, list(mean=0, var=-1, sigma2=1)), "priors\$var must be positive")
  expect_error(.bayesian_engine(y, x, priors, mcmc_samples = 0), "positive integer")
})
