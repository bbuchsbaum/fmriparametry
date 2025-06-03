library(fmriparametric)

test_that(".classify_refinement_queue returns full structure when apply_refinement = FALSE", {
  r2 <- runif(5)
  se <- matrix(runif(15), nrow = 5)
  res <- fmriparametric:::.classify_refinement_queue(
    r2_voxel = r2,
    se_theta_hat_voxel = se,
    refinement_opts = list(apply_refinement = FALSE)
  )
  expect_named(res, c(
    "queue_labels",
    "queue_summary",
    "queue_proportions",
    "queue_details",
    "refinement_needed",
    "classification_criteria"
  ))
  expect_false(res$refinement_needed)
  expect_equal(as.character(res$queue_labels), rep("easy", length(r2)))
})
