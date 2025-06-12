library(devtools)
load_all(quiet = TRUE)

# Test what S_target looks like
set.seed(123)
Y <- matrix(rnorm(50*2), 50, 2)
events <- matrix(0, 50, 1)
events[c(10, 30), 1] <- 1

# Prepare inputs like the pipeline does
inputs <- .prepare_parametric_inputs(
  fmri_data = Y,
  event_model = events,
  hrf_eval_times = seq(0, 30, by = 0.1),
  hrf_span = 30
)

cat("Input dimensions:\n")
cat("Y_proj:", dim(inputs$Y_proj), "\n")
cat("S_target:", dim(inputs$S_target), "\n")
cat("S_target_proj:", dim(inputs$S_target_proj), "\n")

cat("\nS_target content (first 20 rows):\n")
print(head(inputs$S_target, 20))

cat("\nS_target_proj content (first 20 rows):\n")
print(head(inputs$S_target_proj, 20))

# Check if the issue is with how GN processes the stimulus
cat("\nChecking convolution with S_target[,1]:\n")
hrf <- c(0, 0.5, 1, 0.5, 0.2, 0.1, 0.05)
conv_result <- stats::convolve(inputs$S_target[,1], rev(hrf), type = "open")
cat("Convolution result length:", length(conv_result), "\n")
cat("Non-zero elements:", sum(conv_result != 0), "\n")

# Now test if this works in calculate_objective_gn
hrf_interface <- .get_hrf_interface("lwu")
theta_test <- c(5, 2, 0.5)
n_time <- nrow(inputs$Y_proj)

cat("\nTesting objective calculation:\n")
obj <- tryCatch({
  .calculate_objective_gn(
    theta = theta_test,
    y = inputs$Y_proj[,1],
    S = inputs$S_target,
    t_hrf = inputs$hrf_eval_times,
    hrf_interface = hrf_interface,
    n_time = n_time
  )
}, error = function(e) {
  cat("Error in objective:", e$message, "\n")
  NULL
})

if (!is.null(obj)) {
  cat("Objective value:", obj, "\n")
}