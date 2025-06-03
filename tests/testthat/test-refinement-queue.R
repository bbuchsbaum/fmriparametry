library(fmriparametric)

 test_that(".classify_refinement_queue errors on mismatched SE rows", {
   r2_vals <- c(0.5, 0.8)
   se_mat <- matrix(1:6, nrow = 3, ncol = 2)
   expect_error(
     fmriparametric:::.classify_refinement_queue(r2_voxel = r2_vals,
                                                 se_theta_hat_voxel = se_mat),
     "se_theta_hat_voxel must have one row per voxel"
   )
 })

