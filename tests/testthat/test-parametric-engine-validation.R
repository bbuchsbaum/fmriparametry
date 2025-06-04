test_that(".parametric_engine validates Y_proj input", {
  hrf_interface <- list(
    hrf_function = fmriparametric:::.lwu_hrf_function,
    taylor_basis = fmriparametric:::.lwu_hrf_taylor_basis_function,
    parameter_names = fmriparametric:::.lwu_hrf_parameter_names(),
    default_seed = fmriparametric:::.lwu_hrf_default_seed,
    default_bounds = fmriparametric:::.lwu_hrf_default_bounds
  )
  bad_Y <- 1:10
  good_S <- matrix(0, nrow = 10, ncol = 1)
  expect_error(
    fmriparametric:::.parametric_engine(
      Y_proj = bad_Y,
      S_target_proj = good_S,
      hrf_eval_times = seq(0, 30, length.out = 10),
      hrf_interface = hrf_interface,
      theta_seed = c(6, 2.5, 0.35),
      theta_bounds = list(lower = c(2,1,0), upper = c(12,5,1)),
      lambda_ridge = 0.01,
      verbose = FALSE
    ),
    "Y_proj"
  )
})

test_that(".parametric_engine validates S_target_proj input", {
  hrf_interface <- list(
    hrf_function = fmriparametric:::.lwu_hrf_function,
    taylor_basis = fmriparametric:::.lwu_hrf_taylor_basis_function,
    parameter_names = fmriparametric:::.lwu_hrf_parameter_names(),
    default_seed = fmriparametric:::.lwu_hrf_default_seed,
    default_bounds = fmriparametric:::.lwu_hrf_default_bounds
  )
  Y <- matrix(rnorm(20), nrow = 10, ncol = 2)
  bad_S <- 1:10
  expect_error(
    fmriparametric:::.parametric_engine(
      Y_proj = Y,
      S_target_proj = bad_S,
      hrf_eval_times = seq(0, 30, length.out = 10),
      hrf_interface = hrf_interface,
      theta_seed = c(6, 2.5, 0.35),
      theta_bounds = list(lower = c(2,1,0), upper = c(12,5,1)),
      lambda_ridge = 0.01,
      verbose = FALSE
    ),
    "S_target_proj"
  )
})

test_that(".parametric_engine checks row mismatch", {
  hrf_interface <- list(
    hrf_function = fmriparametric:::.lwu_hrf_function,
    taylor_basis = fmriparametric:::.lwu_hrf_taylor_basis_function,
    parameter_names = fmriparametric:::.lwu_hrf_parameter_names(),
    default_seed = fmriparametric:::.lwu_hrf_default_seed,
    default_bounds = fmriparametric:::.lwu_hrf_default_bounds
  )
  Y <- matrix(rnorm(20), nrow = 10, ncol = 2)
  S <- matrix(0, nrow = 8, ncol = 1)
  expect_error(
    fmriparametric:::.parametric_engine(
      Y_proj = Y,
      S_target_proj = S,
      hrf_eval_times = seq(0, 30, length.out = 10),
      hrf_interface = hrf_interface,
      theta_seed = c(6, 2.5, 0.35),
      theta_bounds = list(lower = c(2,1,0), upper = c(12,5,1)),
      lambda_ridge = 0.01,
      verbose = FALSE
    ),
    "S_target_proj must have"
  )
})
