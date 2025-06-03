test_that("memory_limit NA triggers fallback to default", {
  decision <- .smart_memory_decision(100, 100, memory_limit = NA)
  expect_equal(decision$available_memory_gb, 8)
  expect_false(decision$use_chunking)
})
