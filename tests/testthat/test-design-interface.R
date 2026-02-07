# Tests for design-first interface harmonization

test_that("create_parametric_design builds projected inputs", {
  set.seed(123)
  n_time <- 50
  n_vox <- 4

  fmri_data <- matrix(rnorm(n_time * n_vox), nrow = n_time, ncol = n_vox)
  event_model <- matrix(0, nrow = n_time, ncol = 1)
  event_model[c(8, 20, 35, 45), 1] <- 1

  design <- create_parametric_design(
    fmri_data = fmri_data,
    event_model = event_model,
    hrf_eval_times = seq(0, 24, by = 0.5)
  )

  expect_s3_class(design, "parametric_design")
  expect_true(is.matrix(design$Y_proj))
  expect_true(is.matrix(design$S_target_proj))
  expect_equal(dim(design$Y_proj), c(n_time, n_vox))
  expect_equal(nrow(design$S_target_proj), n_time)
  expect_equal(length(design$hrf_eval_times), length(seq(0, 24, by = 0.5)))
  expect_equal(design$metadata$n_vox, n_vox)
  expect_equal(design$metadata$n_cond, 1)
})

test_that("estimate_parametric_hrf_from_design matches direct fit", {
  set.seed(321)
  n_time <- 60
  n_vox <- 5

  fmri_data <- matrix(rnorm(n_time * n_vox), nrow = n_time, ncol = n_vox)
  event_model <- matrix(0, nrow = n_time, ncol = 1)
  event_model[c(10, 25, 40, 55), 1] <- 1

  hrf_grid <- seq(0, 24, by = 0.5)

  fit_direct <- estimate_parametric_hrf(
    fmri_data = fmri_data,
    event_model = event_model,
    hrf_eval_times = hrf_grid,
    global_refinement = FALSE,
    compute_se = FALSE,
    progress = FALSE,
    verbose = FALSE
  )

  design <- create_parametric_design(
    fmri_data = fmri_data,
    event_model = event_model,
    hrf_eval_times = hrf_grid
  )

  fit_design <- estimate_parametric_hrf_from_design(
    design = design,
    global_refinement = FALSE,
    compute_se = FALSE,
    progress = FALSE,
    verbose = FALSE
  )

  expect_equal(get_parameters(fit_design), get_parameters(fit_direct), tolerance = 1e-8)
  expect_equal(fit_design$amplitudes, fit_direct$amplitudes, tolerance = 1e-8)
  expect_equal(get_gof_per_voxel(fit_design), get_gof_per_voxel(fit_direct), tolerance = 1e-8)

  info <- get_design_info(fit_design)
  expect_equal(info$n_time, n_time)
  expect_equal(info$n_vox, n_vox)
  expect_equal(info$n_cond, ncol(event_model))
  expect_equal(get_method_used(fit_design), "parametric_taylor")
})

test_that("estimate_parametric_hrf_from_design validates design", {
  y <- matrix(rnorm(100), nrow = 20)

  bad_design <- list(S_target_proj = matrix(0, nrow = 20, ncol = 1))
  expect_error(
    estimate_parametric_hrf_from_design(y = y, design = bad_design, verbose = FALSE),
    "must contain"
  )
})
