# Debug script for convolution issue
# Following Gemini's hypothesis 1

library(fmriparametric)
library(fmrihrf)

cat("=== Testing fast_batch_convolution ===\n\n")

# 1. Create a simple HRF (using fmrihrf functions)
time_points <- seq(0, 32, by = 2)
hrf_values <- fmrihrf::hrf_spmg1(time_points)  # Use a simple SPM HRF
cat("HRF values:\n")
print(hrf_values)
cat("Sum of HRF:", sum(hrf_values), "\n\n")

# 2. Create a delta function stimulus (event at time t=10)
n_time <- 100
stim <- numeric(n_time)
stim[10] <- 1  # Single spike

cat("Stimulus vector (showing indices 5-15):\n")
print(stim[5:15])
cat("Sum of stimulus:", sum(stim), "\n\n")

# 3. Test convolution
# Convert to matrix format as expected by fast_batch_convolution
stim_matrix <- matrix(stim, ncol = 1)
hrf_matrix <- matrix(hrf_values, ncol = 1)

cat("Testing .fast_batch_convolution...\n")
result <- .fast_batch_convolution(
  signal = stim,
  kernels = hrf_matrix,
  output_length = n_time
)

cat("\nConvolution result:\n")
cat("Dimensions:", dim(result), "\n")
cat("Result sum:", sum(result), "\n")
cat("Max value:", max(result), "\n")
cat("Non-zero elements:", sum(result != 0), "\n")

# Show the peak region (should see HRF shape starting at index 10)
cat("\nResult around stimulus time (indices 8-25):\n")
print(round(result[8:25, 1], 4))

# Plot if we can
if (interactive()) {
  plot(1:n_time, result[,1], type = 'l', 
       main = "Convolution Result", 
       xlab = "Time", ylab = "Signal")
  abline(v = 10, col = "red", lty = 2)
}

# Also test with the internal convolution used in parametric engine
cat("\n\n=== Testing stats::convolve (as comparison) ===\n")
conv_stats <- stats::convolve(stim, rev(hrf_values), type = "open")
conv_stats <- conv_stats[1:n_time]  # Trim to original length
cat("stats::convolve sum:", sum(conv_stats), "\n")
cat("stats::convolve max:", max(conv_stats), "\n")

# Compare
cat("\nDifference between methods:", sum(abs(result[,1] - conv_stats)), "\n")