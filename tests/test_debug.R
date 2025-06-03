# Debug test
library(fmrireg)
source("R/parametric-engine.R")
source("R/hrf-interface-lwu.R")

# Minimal test
set.seed(123)
n_time <- 50
n_vox <- 3

# Events
S <- matrix(0, n_time, 1)
S[c(10, 20, 30, 40), 1] <- 1

# Simple data
Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)

# HRF times
t_hrf <- seq(0, 30, length.out = 61)

# Interface
hrf_interface <- list(
  hrf_function = .lwu_hrf_function,
  taylor_basis = .lwu_hrf_taylor_basis_function,
  parameter_names = c("tau", "sigma", "rho"),
  default_seed = c(6, 2.5, 0.35),
  default_bounds = list(lower = c(2, 1, 0), upper = c(12, 5, 1))
)

cat("Testing parametric engine...\n")
cat("Y dimensions:", dim(Y), "\n")
cat("S dimensions:", dim(S), "\n")

# Try the engine
tryCatch({
  fit <- .parametric_engine(
    Y_proj = Y,
    S_target_proj = S,
    hrf_eval_times = t_hrf,
    hrf_interface = hrf_interface,
    theta_seed = c(6, 2.5, 0.35),
    theta_bounds = list(lower = c(2, 1, 0), upper = c(12, 5, 1)),
    lambda_ridge = 0.01
  )
  
  cat("\nFit object names:", names(fit), "\n")
  cat("theta_hat class:", class(fit$theta_hat), "\n")
  cat("theta_hat dimensions:", dim(fit$theta_hat), "\n")
  
  if (!is.null(dim(fit$theta_hat))) {
    cat("theta_hat first row:", fit$theta_hat[1,], "\n")
    cat("Mean parameters:", colMeans(fit$theta_hat), "\n")
  }
  
  cat("beta0 length:", length(fit$beta0), "\n")
  cat("r_squared length:", length(fit$r_squared), "\n")
  
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  cat("Traceback:\n")
  traceback()
})