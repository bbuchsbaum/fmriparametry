#' Simple Bayesian HRF estimation engine
#'
#' This is a lightweight Bayesian alternative used mainly for testing the
#' plugin architecture. It performs conjugate Bayesian inference on the
#' amplitude parameter while treating the HRF parameters as fixed.
#'
#' @param Y Numeric vector of observed BOLD data
#' @param X Design matrix for the predicted response
#' @param priors List with elements `mean`, `var` and `sigma2`
#' @param mcmc_samples Number of posterior samples to draw
#' @return List with posterior mean, sd and samples
#' @keywords internal
.bayesian_engine <- function(Y, X, priors, mcmc_samples = 1000) {
  XtX <- as.numeric(crossprod(X))
  Xty <- as.numeric(crossprod(X, Y))
  post_var <- 1 / (1/priors$var + XtX / priors$sigma2)
  post_mean <- post_var * (priors$mean / priors$var + Xty / priors$sigma2)
  samples <- rnorm(mcmc_samples, post_mean, sqrt(post_var))
  list(posterior_mean = post_mean,
       posterior_sd = sqrt(post_var),
       samples = samples)
}
