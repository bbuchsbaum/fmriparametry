# Debug the scale mismatch issue

library(fmriparametric)
library(fmrihrf)

set.seed(42)
n_time <- 50

# Simple impulse at time 10
S <- matrix(0, n_time, 1)
S[10, 1] <- 1

# Simple data - constant with noise
Y <- matrix(100 + rnorm(n_time, 0, 5), n_time, 1)

# Add a scaled HRF response at event time
hrf_vals <- fmrihrf::hrf_spmg1(seq(0, 20, by = 1))
Y[10:(10+length(hrf_vals)-1), 1] <- Y[10:(10+length(hrf_vals)-1), 1] + hrf_vals * 20

cat("Data setup:\n")
cat("Y range:", range(Y), "\n")
cat("Y mean:", mean(Y), "\n")
cat("Y variance:", var(Y), "\n\n")

# Get HRF interface and create Taylor basis
hrf_interface <- fmriparametric:::.get_hrf_interface("lwu")
theta <- c(6, 1, 0.5)
hrf_times <- seq(0, 30, by = 2)
basis <- hrf_interface$taylor_basis(theta, hrf_times)

cat("HRF basis:\n")
cat("Basis dim:", dim(basis), "\n")
cat("HRF (col 1) range:", range(basis[,1]), "\n")
cat("HRF sum:", sum(basis[,1]), "\n\n")

# Manual convolution
X <- matrix(0, n_time, ncol(basis))
for (j in 1:ncol(basis)) {
  conv_result <- stats::convolve(S[,1], rev(basis[,j]), type = "open")
  X[,j] <- conv_result[1:n_time]
}

cat("Design matrix X:\n")
cat("X dim:", dim(X), "\n")
cat("X[,1] (HRF) non-zeros at:", which(X[,1] != 0), "\n")
cat("X[,1] range:", range(X[,1]), "\n\n")

# Manual ridge regression
XtX <- crossprod(X)
XtY <- crossprod(X, Y)
lambda <- 0.01
coeffs <- solve(XtX + lambda * diag(ncol(X)), XtY)

cat("Ridge regression:\n")
cat("Coefficients:", as.vector(coeffs), "\n")
cat("Beta0 (amplitude):", coeffs[1], "\n\n")

# Fitted values
Y_fitted <- X %*% coeffs

cat("Fit quality:\n")
cat("Y_fitted range:", range(Y_fitted), "\n")
cat("Y_fitted at event (rows 10-15):\n")
print(Y_fitted[10:15,1])
cat("\nY actual at event (rows 10-15):\n") 
print(Y[10:15,1])

# Compute R-squared properly
Y_centered <- Y - mean(Y)
Y_fitted_centered <- Y_fitted - mean(Y)
SS_tot <- sum(Y_centered^2)
SS_res <- sum((Y - Y_fitted)^2)
R2 <- 1 - SS_res/SS_tot

cat("\nR-squared calculation:\n")
cat("SS_tot:", SS_tot, "\n")
cat("SS_res:", SS_res, "\n")
cat("R²:", R2, "\n")

# Check centering issue
cat("\nMeans:\n")
cat("mean(Y):", mean(Y), "\n")
cat("mean(Y_fitted):", mean(Y_fitted), "\n")
cat("Difference in means:", mean(Y) - mean(Y_fitted), "\n")