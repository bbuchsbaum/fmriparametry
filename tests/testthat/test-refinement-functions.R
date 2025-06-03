test_that("refinement helper functions run without errors", {
  set.seed(123)
  n_time <- 20
  n_vox <- 3
  Y_proj <- matrix(rnorm(n_time * n_vox), nrow = n_time)
  S_target_proj <- matrix(rbinom(n_time, 1, 0.2), ncol = 1)
  hrf_eval_times <- seq(0, 30, by = 0.5)
  hrf_interface <- fmriparametric:::.create_hrf_interface("lwu")
  theta_bounds <- hrf_interface$default_bounds()
  theta_current <- matrix(rep(hrf_interface$default_seed(), n_vox),
                          nrow = n_vox, byrow = TRUE)
  r_squared <- rep(0.5, n_vox)

  expect_silent(
    res_mod <- fmriparametric:::.refine_moderate_voxels(
      voxel_idx = 1:n_vox,
      Y_proj = Y_proj,
      S_target_proj = S_target_proj,
      theta_current = theta_current,
      r_squared = r_squared,
      hrf_interface = hrf_interface,
      hrf_eval_times = hrf_eval_times,
      theta_bounds = theta_bounds,
      parallel = FALSE,
      n_cores = 1
    )
  )
  expect_type(res_mod, "list")

  expect_silent(
    res_hard <- fmriparametric:::.refine_hard_voxels(
      voxel_idx = 1:n_vox,
      Y_proj = Y_proj,
      S_target_proj = S_target_proj,
      theta_current = theta_current,
      r_squared = r_squared,
      hrf_interface = hrf_interface,
      hrf_eval_times = hrf_eval_times,
      theta_bounds = theta_bounds,
      max_iter = 2,
      parallel = FALSE,
      n_cores = 1
    )
  )
  expect_type(res_hard, "list")
})
