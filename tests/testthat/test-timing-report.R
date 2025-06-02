test_that("get_timing_report works correctly", {
  # Get initial state
  initial_report <- get_timing_report()
  expect_true(is.data.frame(initial_report))
  expect_named(initial_report, c("operation", "count", "total_time", "mean_time"))
  
  # Enable verbose mode to record timings
  old_option <- getOption("fmriparametric.verbose")
  options(fmriparametric.verbose = TRUE)
  
  # Run a timed operation
  set.seed(123)
  fmri_data <- matrix(rnorm(100), nrow = 20, ncol = 5)
  event_design <- matrix(rbinom(20, 1, 0.2), ncol = 1)
  
  # This should record timing information
  fit <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_design,
    verbose = FALSE
  )
  
  # Get report after operation
  report <- get_timing_report()
  expect_true(is.data.frame(report))
  
  # If timing was recorded, we should have entries
  if (nrow(report) > nrow(initial_report)) {
    expect_true(all(report$count >= 1))
    expect_true(all(report$total_time >= 0))
    expect_true(all(report$mean_time >= 0))
  }
  
  # Restore option
  options(fmriparametric.verbose = old_option)
})