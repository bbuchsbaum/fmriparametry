library(testthat)
library(fmriparametric)

context("HRF registry")

dummy <- list(
  hrf_function = function(t,p) rep(0,length(t)),
  taylor_basis = function(p,t) matrix(0,length(t),4),
  parameter_names = c("a","b","c"),
  default_seed = function() c(1,1,1),
  default_bounds = function() list(lower=c(0,0,0), upper=c(2,2,2))
)
register_hrf_model("dummy", dummy)

int <- .create_hrf_interface("dummy")

test_that("custom model is retrievable", {
  expect_equal(int$parameter_names, c("a","b","c"))
})
