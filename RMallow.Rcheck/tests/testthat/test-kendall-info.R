# Tests for KendallInfo function

test_that("KendallInfo computes correct pairwise comparisons", {
  # Simple increasing sequence: all comparisons should be 0 (increase)
  r <- matrix(1:5, nrow = 1)
  info <- KendallInfo(r)
  
  expect_equal(nrow(info), 1)
  expect_equal(ncol(info), 10)  # 5 choose 2 = 10 pairs
  expect_true(all(info == 0))
})

test_that("KendallInfo handles decreasing sequence", {

# Decreasing sequence: all comparisons should be 1 (decrease)
  r <- matrix(5:1, nrow = 1)
  info <- KendallInfo(r)
  
  expect_true(all(info == 1))
})

test_that("KendallInfo handles ties correctly", {
  # Sequence with ties
  r <- matrix(c(1, 1, 2, 3, 3), nrow = 1)
  info <- KendallInfo(r)
  
  # Tied positions should produce NA
  expect_true(any(is.na(info)))
  
  # Count NAs - ties at positions (1,2) and (4,5)
  # Pairs involving ties: (1,2), (4,5) = 2 NAs
  expect_equal(sum(is.na(info)), 2)
})

test_that("KendallInfo works with multiple rows", {
  r <- rbind(1:4, 4:1)
  info <- KendallInfo(r)
  
  expect_equal(nrow(info), 2)
  expect_equal(ncol(info), 6)  # 4 choose 2 = 6
  
  # First row: all increases (0)
  expect_true(all(info[1, ] == 0))
  
  # Second row: all decreases (1)
  expect_true(all(info[2, ] == 1))
})

test_that("KendallInfo validates input", {
  # Empty matrix should error
  expect_error(KendallInfo(matrix(nrow = 0, ncol = 0)))
})

test_that("KendallInfo handles single column", {
  r <- matrix(c(1, 2, 3), ncol = 1)
  info <- KendallInfo(r)
  
  # Single column means no pairs to compare
  expect_equal(ncol(info), 0)
})

test_that("KendallInfo accepts custom indices", {
  r <- matrix(1:5, nrow = 1)
  # Only compare columns 1-2 and 1-3
  custom_inds <- matrix(c(1, 1, 2, 3), nrow = 2)
  info <- KendallInfo(r, inds = custom_inds)
  
  expect_equal(ncol(info), 2)
})
