
# Test that .ridge_linear_solve errors on mismatched row counts

test_that(".ridge_linear_solve errors when X and Y have different rows", {
  X <- matrix(rnorm(6), nrow = 2, ncol = 3)
  Y <- matrix(rnorm(9), nrow = 3, ncol = 3)
  expect_error(fmriparametric:::.ridge_linear_solve(X, Y),
               "same number of rows")
})
