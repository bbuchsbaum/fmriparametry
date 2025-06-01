# Mock functions to test without fmrireg dependency
# ===================================================

# Mock LWU HRF functions
.lwu_hrf_function <- function(t, params_vector, ...) {
  tau <- params_vector[1]
  sigma <- params_vector[2]
  rho <- params_vector[3]
  
  # Simple LWU implementation without fmrireg
  main <- exp(-(t - tau)^2 / (2 * sigma^2))
  undershoot <- rho * exp(-(t - tau - 2*sigma)^2 / (2 * (1.6*sigma)^2))
  hrf <- main - undershoot
  hrf[t < 0] <- 0
  return(hrf)
}

.lwu_hrf_taylor_basis_function <- function(params_vector0, t_hrf_eval, ...) {
  n_t <- length(t_hrf_eval)
  basis <- matrix(0, n_t, 4)
  
  # Base function
  basis[, 1] <- .lwu_hrf_function(t_hrf_eval, params_vector0)
  
  # Numerical derivatives
  delta <- 1e-4
  for (i in 1:3) {
    params_plus <- params_minus <- params_vector0
    params_plus[i] <- params_vector0[i] + delta
    params_minus[i] <- params_vector0[i] - delta
    
    h_plus <- .lwu_hrf_function(t_hrf_eval, params_plus)
    h_minus <- .lwu_hrf_function(t_hrf_eval, params_minus)
    
    basis[, i + 1] <- (h_plus - h_minus) / (2 * delta)
  }
  
  return(basis)
}

.lwu_hrf_parameter_names <- function() {
  c("tau", "sigma", "rho")
}

.lwu_hrf_default_seed <- function() {
  c(6, 2.5, 0.35)
}

.lwu_hrf_default_bounds <- function() {
  list(
    lower = c(0, 0.05, 0),
    upper = c(20, 10, 1.5)
  )
}