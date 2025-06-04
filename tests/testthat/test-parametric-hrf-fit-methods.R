library(fmriparametric)
library(testthat)

pars <- matrix(rep(1:3, each = 2), nrow = 2, byrow = TRUE)
amps <- c(1, 2)
obj <- fmriparametric:::new_parametric_hrf_fit(
  estimated_parameters = pars,
  amplitudes = amps,
  parameter_names = c("tau", "sigma", "rho"),
  metadata = list(n_timepoints = 10)
)

test_that("print method produces output", {
  expect_output(print(obj), "Parametric HRF Fit")
})

test_that("coef method returns matrix", {
  coef_result <- coef(obj)
  expect_equal(unname(coef_result), unname(pars))
  expect_equal(dim(coef_result), dim(pars))
})

test_that("summary method returns list with summaries", {
  s <- summary(obj)
  expect_true(is.list(s))
  expect_true("parameter_summary" %in% names(s))
  expect_true("amplitude_summary" %in% names(s))
})
