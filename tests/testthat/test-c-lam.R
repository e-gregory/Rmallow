# Tests for C_lam function (normalizing constant)

test_that("C_lam returns 1/n! when lambda = 0", {
  # When lambda = 0, all permutations are equally likely
  # So C = 1 / (sum of all 1s) = 1 / n!
  
  dist4 <- DistanceDistribution(4)
  c_val <- C_lam(0, dists_table = dist4)
  
  expect_equal(c_val, 1 / factorial(4))
})

test_that("C_lam increases as lambda increases", {
  dist5 <- DistanceDistribution(5)
  
  c0 <- C_lam(0, dists_table = dist5)
  c1 <- C_lam(1, dists_table = dist5)
  c2 <- C_lam(2, dists_table = dist5)
  
  expect_true(c1 > c0)
  expect_true(c2 > c1)
})

test_that("C_lam approaches 1 as lambda approaches infinity", {
  dist4 <- DistanceDistribution(4)
  
  c_large <- C_lam(100, dists_table = dist4)
  
  # Should be very close to 1
  expect_true(c_large > 0.99)
})

test_that("C_lam validates lambda", {
  dist4 <- DistanceDistribution(4)
  
  expect_error(C_lam(-1, dists_table = dist4), "non-negative")
})

test_that("C_lam requires dists_table or dists", {
  expect_error(C_lam(1), "must be provided")
})

test_that("C_lam computes from raw distances", {
  # Can provide raw distances instead of table
  dists <- c(0, 1, 1, 2, 2, 3)  # Some example distances
  c_val <- C_lam(1, dists = dists)
  
  expect_true(is.numeric(c_val))
  expect_true(c_val > 0)
})
