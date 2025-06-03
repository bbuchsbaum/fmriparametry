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
  # Validate Y and X
  if (!is.numeric(Y) || !is.vector(Y)) {
    stop(".bayesian_engine: Y must be a numeric vector", call. = FALSE)
  }
  if (!is.numeric(X) || !is.vector(X)) {
    stop(".bayesian_engine: X must be a numeric vector", call. = FALSE)
  }
  if (length(Y) != length(X)) {
    stop(".bayesian_engine: Y and X must have equal length", call. = FALSE)
  }
  if (any(!is.finite(Y)) || any(!is.finite(X))) {
    stop(".bayesian_engine: Y and X must contain only finite values", call. = FALSE)
  }

  # Validate priors
  if (!is.list(priors) || !all(c("mean", "var", "sigma2") %in% names(priors))) {
    stop(".bayesian_engine: priors must be a list with numeric scalars 'mean', 'var', and 'sigma2'", call. = FALSE)
  }
  if (!is.numeric(priors$mean) || length(priors$mean) != 1 || !is.finite(priors$mean) ||
      !is.numeric(priors$var) || length(priors$var) != 1 || !is.finite(priors$var) ||
      !is.numeric(priors$sigma2) || length(priors$sigma2) != 1 || !is.finite(priors$sigma2)) {
    stop(".bayesian_engine: priors must contain numeric scalars 'mean', 'var', and 'sigma2'", call. = FALSE)
  }
  if (priors$var <= 0) {
    stop(".bayesian_engine: priors$var must be positive", call. = FALSE)
  }
  if (priors$sigma2 <= 0) {
    stop(".bayesian_engine: priors$sigma2 must be positive", call. = FALSE)
  }

  # Validate mcmc_samples
  if (length(mcmc_samples) != 1 || !is.numeric(mcmc_samples) ||
      !is.finite(mcmc_samples) || mcmc_samples <= 0 ||
      mcmc_samples != as.integer(mcmc_samples)) {
    stop(".bayesian_engine: mcmc_samples must be a positive integer", call. = FALSE)
  }

  XtX <- as.numeric(crossprod(X))
  Xty <- as.numeric(crossprod(X, Y))
  post_var <- 1 / (1/priors$var + XtX / priors$sigma2)
  post_mean <- post_var * (priors$mean / priors$var + Xty / priors$sigma2)
  samples <- rnorm(mcmc_samples, post_mean, sqrt(post_var))
  list(posterior_mean = post_mean,
       posterior_sd = sqrt(post_var),
       samples = samples)
}
