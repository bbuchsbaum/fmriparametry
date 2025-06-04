library(fmriparametric)
library(testthat)

set.seed(123)

data_obj <- matrix(rnorm(20), nrow = 10, ncol = 2)
event_obj <- matrix(rbinom(10, 1, 0.2), ncol = 1)

test_that("global refinement respects package option", {
  options(fmriparametric.refine_global = FALSE)
  fit_off <- estimate_parametric_hrf(
    fmri_data = data_obj,
    event_model = event_obj,
    global_refinement = TRUE,
    global_passes = 1,
    verbose = FALSE
  )
  expect_equal(length(fit_off$convergence), 0)

  options(fmriparametric.refine_global = TRUE)
  fit_on <- estimate_parametric_hrf(
    fmri_data = data_obj,
    event_model = event_obj,
    global_refinement = TRUE,
    global_passes = 1,
    verbose = FALSE
  )
  expect_true(length(fit_on$convergence) > 0)

  options(fmriparametric.refine_global = NULL)
})
