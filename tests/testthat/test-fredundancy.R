test_that("fredundancy computes redundancy measures", {

  # Test data
  functions <- c(0.8, 0.1, 0.05, 0.05, 0)
  abundance <- c(0.2, 0.1, 0.05, 0.05, 0.6)
  n_reference <- 7

  # Expected results (manually obtained or precomputed)
  expected <- list(
    sample_based = -0.901091330,
    reference_based = -1.23756357,
    abundance_based = -1.10903549,
    interdependency = 0.192744760
  )

  # Call the function
  actual <- fredundancy(functions = functions, abundance = abundance, n_reference = 7)

  # Round both actual and expected values to 8 decimal places for consistency
  actual_rounded <- lapply(actual, function(x) round(x, 8))  # Round to 8 decimal places
  expected_rounded <- lapply(expected, function(x) round(x, 8))  # Round to 8 decimal places

  # Use expect_equal with a larger tolerance for floating-point comparison
  expect_equal(actual_rounded, expected_rounded, tolerance = 1e-5)  # Allow even finer differences
})
