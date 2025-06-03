library(fmriparametric)

context("refinement queue classification and summary")

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
