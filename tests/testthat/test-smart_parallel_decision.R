test_that(".smart_parallel_decision handles detectCores edge cases", {
  # simulate detectCores returning 1
  testthat::with_mock(
    `parallel::detectCores` = function(...) 1L,
    {
      res <- .smart_parallel_decision(n_voxels = 100)
      expect_false(res$use_parallel)
    }
  )

  # simulate detectCores returning NA
  testthat::with_mock(
    `parallel::detectCores` = function(...) NA_integer_,
    {
      res <- .smart_parallel_decision(n_voxels = 100)
      expect_false(res$use_parallel)
    }
  )
})
