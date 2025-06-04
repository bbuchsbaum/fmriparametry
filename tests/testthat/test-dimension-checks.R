test_that(".parametric_engine errors on mismatched input sizes", {
  Y <- matrix(rnorm(20), nrow = 10, ncol = 2)
  S <- matrix(0, nrow = 9, ncol = 1)

  hrf_interface <- list(
    hrf_function = .lwu_hrf_function,
    taylor_basis = .lwu_hrf_taylor_basis_function,
    parameter_names = .lwu_hrf_parameter_names(),
    default_seed = .lwu_hrf_default_seed(),
    default_bounds = .lwu_hrf_default_bounds()
  )

  expect_error(
    .parametric_engine(
      Y_proj = Y,
      S_target_proj = S,
      hrf_eval_times = seq(0, 30, length.out = 61),
      hrf_interface = hrf_interface,
      theta_seed = hrf_interface$default_seed(),
      theta_bounds = hrf_interface$default_bounds(),
      lambda_ridge = 0.01,
      verbose = FALSE
    ),
    "S_target_proj must have"
  )
})
