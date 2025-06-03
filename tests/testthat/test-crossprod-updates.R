test_that(".bayesian_engine computes posterior quantities correctly", {
  X <- 1:5
  Y <- 2:6
  priors <- list(mean = 0, var = 1, sigma2 = 1)

  result <- fmriparametric:::.bayesian_engine(Y, X, priors, mcmc_samples = 0)

  XtX <- as.numeric(crossprod(X))
  Xty <- as.numeric(crossprod(X, Y))
  post_var <- 1 / (1/priors$var + XtX / priors$sigma2)
  post_mean <- post_var * (priors$mean / priors$var + Xty / priors$sigma2)

  expect_equal(result$posterior_mean, post_mean)
  expect_equal(result$posterior_sd, sqrt(post_var))
  expect_length(result$samples, 0)
})

test_that(".calculate_objective_gn returns correct sum of squares", {
  n_time <- 5
  y <- seq_len(n_time)
  S <- matrix(c(1, rep(0, n_time - 1)), ncol = 1)
  t_hrf <- 1
  hrf_interface <- list(hrf_function = function(t, theta) rep(1, length(t)))

  obj <- fmriparametric:::.calculate_objective_gn(numeric(), y, S, t_hrf,
                                                  hrf_interface, n_time)

  expect_equal(obj, sum(y[-1]^2))
})
