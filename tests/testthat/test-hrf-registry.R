test_that("register_hrf_model works correctly", {
  # Define a dummy HRF model
  dummy_model <- list(
    hrf_function = function(t, p) rep(0, length(t)),
    taylor_basis = function(p, t) matrix(0, length(t), 4),
    parameter_names = c("a", "b", "c"),
    default_seed = function() c(1, 1, 1),
    default_bounds = function() list(lower = c(0, 0, 0), upper = c(2, 2, 2))
  )
  
  # Register the model
  expect_silent(register_hrf_model("dummy", dummy_model))
  
  # Test that the model can be used with estimate_parametric_hrf
  set.seed(123)
  fmri_data <- matrix(rnorm(100), nrow = 20, ncol = 5)
  event_design <- matrix(rbinom(20, 1, 0.2), ncol = 1)
  
  # This should work with the registered model
  expect_error(
    estimate_parametric_hrf(
      fmri_data = fmri_data,
      event_model = event_design,
      parametric_hrf = "dummy",
      verbose = FALSE
    ),
    NA  # Expect no error
  )
})