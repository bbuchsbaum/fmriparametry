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



# Representative inputs
r2_vals <- c(0.2, 0.5, 0.8)
se_mat <- matrix(c(0.6, 0.6,
                   0.35, 0.35,
                   0.2, 0.2), ncol = 2, byrow = TRUE)

# helper function to capture print output
capture_summary <- function(res) {
  paste(capture.output(fmriparametric:::.print_refinement_summary(res)), collapse = "\n")
}

test_that("queue classification works with SE data", {
  res <- fmriparametric:::.classify_refinement_queue(r2_vals, se_mat)
  expect_equal(res$queue_labels,
               c("hard_GN", "moderate_local_recenter", "easy"))
  expect_true(res$refinement_needed)
  expect_equal(as.numeric(res$queue_summary["hard_GN"]), 1)
  expect_equal(as.numeric(res$queue_summary["moderate_local_recenter"]), 1)
  expect_equal(as.numeric(res$queue_summary["easy"]), 1)

  out <- capture_summary(res)
  expect_match(out, "Refinement Queue Classification")
  expect_match(out, "hard_GN")
  expect_match(out, "SE thresholds: hard")
})

test_that("queue classification works without SE data", {
  res <- fmriparametric:::.classify_refinement_queue(r2_vals)
  expect_equal(res$queue_labels,
               c("hard_GN", "moderate_local_recenter", "easy"))
  expect_true(res$refinement_needed)
  expect_true(!res$classification_criteria$se_available)
})

  out <- capture_summary(res)
  expect_match(out, "SE thresholds: not used")
})

test_that("queue classification can be turned off", {
  res <- fmriparametric:::.classify_refinement_queue(
    r2_vals, se_mat, refinement_opts = list(apply_refinement = FALSE))
  expect_equal(res$queue_labels, rep("easy", length(r2_vals)))
  expect_false(res$refinement_needed)
  expect_equal(as.numeric(res$queue_summary["easy"]), length(r2_vals))

})

