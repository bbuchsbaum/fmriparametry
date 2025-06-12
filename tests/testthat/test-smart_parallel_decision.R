# Skip all smart parallel decision tests - function is in attic for a reason
skip("Smart parallel decision function is in attic - skipping all tests")

test_that(".smart_parallel_decision respects system configuration", {
  # Test with small n_voxels
  res <- fmriparametric:::.smart_parallel_decision(n_voxels = 10)
  expect_false(res$use_parallel)
  expect_match(res$reason, "Too few voxels")
  
  # Test with large n_voxels
  res <- fmriparametric:::.smart_parallel_decision(n_voxels = 10000)
  # This will depend on the system, but we can check the structure
  expect_true(is.list(res))
  expect_true("use_parallel" %in% names(res))
  expect_true("recommended_cores" %in% names(res))
  expect_true("reason" %in% names(res))
  expect_true("estimated_speedup" %in% names(res))
  
  # Test with specified n_cores
  res <- fmriparametric:::.smart_parallel_decision(n_voxels = 10000, n_cores = 1)
  expect_false(res$use_parallel)
  expect_equal(res$reason, "Single core system")
  
  # Test with high overhead
  res <- fmriparametric:::.smart_parallel_decision(n_voxels = 100, overhead_per_voxel = 1)
  expect_false(res$use_parallel)
})