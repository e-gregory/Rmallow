# Tests for SimplifySequences function

test_that("SimplifySequences handles already simplified sequences", {
  rankings <- matrix(1:5, nrow = 1)
  result <- SimplifySequences(rankings)
  
  expect_equal(as.numeric(result), 1:5)
})

test_that("SimplifySequences removes gaps", {
  # (1, 1, 2, 4, 4, 5) -> (1, 1, 2, 3, 3, 4)
  rankings <- matrix(c(1, 1, 2, 4, 4, 5), nrow = 1)
  result <- SimplifySequences(rankings)
  
  expect_equal(as.numeric(result), c(1, 1, 2, 3, 3, 4))
})

test_that("SimplifySequences handles multiple rows", {
  rankings <- rbind(
    c(1, 3, 3, 5),
    c(2, 2, 4, 6)
  )
  result <- SimplifySequences(rankings)
  
  expect_equal(as.numeric(result[1, ]), c(1, 2, 2, 3))
  expect_equal(as.numeric(result[2, ]), c(1, 1, 2, 3))
})

test_that("SimplifySequences preserves ties", {
  rankings <- matrix(c(1, 1, 1, 2, 2), nrow = 1)
  result <- SimplifySequences(rankings)
  
  expect_equal(as.numeric(result), c(1, 1, 1, 2, 2))
})

test_that("SimplifySequences handles large gaps", {
  rankings <- matrix(c(10, 20, 30), nrow = 1)
  result <- SimplifySequences(rankings)
  
  expect_equal(as.numeric(result), c(1, 2, 3))
})

test_that("SimplifySequences handles empty input", {
  rankings <- matrix(nrow = 0, ncol = 5)
  result <- SimplifySequences(rankings)
  
  expect_equal(nrow(result), 0)
})

test_that("SimplifySequences handles vector input", {
  rankings <- c(1, 3, 3, 5)
  result <- SimplifySequences(rankings)
  
  expect_equal(as.numeric(result), c(1, 2, 2, 3))
})
