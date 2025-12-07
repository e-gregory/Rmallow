# Tests for Lambda and UpdateLambda functions

test_that("Lambda objective function is zero at correct value", {
  dist_table <- DistanceDistribution(4)
  
  # When rhs equals the expected distance, Lambda should be ~0
  # At lambda=0, expected distance is mean of all distances
  distances <- as.numeric(names(dist_table))
  counts <- as.numeric(dist_table)
  expected_dist_at_zero <- sum(distances * counts) / sum(counts)
  
  result <- Lambda(0, rhs = expected_dist_at_zero, dists_table = dist_table)
  expect_equal(result, 0, tolerance = 1e-10)
})

test_that("Lambda decreases as lambda increases", {
  dist_table <- DistanceDistribution(4)
  rhs <- 2  # Some target distance
  
  l0 <- Lambda(0, rhs = rhs, dists_table = dist_table)
  l1 <- Lambda(1, rhs = rhs, dists_table = dist_table)
  l2 <- Lambda(2, rhs = rhs, dists_table = dist_table)
  
  # Lambda output should decrease as lambda increases
  # (expected distance decreases, so difference from rhs changes)
  expect_true(l1 < l0)
  expect_true(l2 < l1)
})

test_that("Lambda validates inputs", {
  expect_error(Lambda(1, rhs = 1), "must be provided")
  
  dist_table <- DistanceDistribution(4)
  expect_error(Lambda(-1, rhs = 1, dists_table = dist_table), "non-negative")
})

test_that("UpdateLambda finds correct lambda for simple case", {
  set.seed(42)
  
  # Create data all at distance 0 from modal sequence
  r <- rbind(1:4, 1:4, 1:4, 1:4)
  R <- list(1:4)
  z <- matrix(1, nrow = 4, ncol = 1)
  dist_table <- DistanceDistribution(4)
  dists_to_Rg <- AllKendall(r, do.call(rbind, R))
  
  lambda <- UpdateLambda(r, R, z, G = 1, dists_to_Rg, dist_table)
  
  # When all observations are at modal sequence (distance 0),
  # lambda should be very large (concentrated distribution)
  expect_true(lambda[1] > 100)
})

test_that("UpdateLambda handles multiple clusters", {
  set.seed(42)
  
  r <- rbind(
    1:4, 1:4, c(1, 2, 4, 3),  # Close to 1:4
    4:1, 4:1, c(4, 3, 1, 2)   # Close to 4:1
  )
  R <- list(1:4, 4:1)
  z <- rbind(
    c(0.9, 0.1), c(0.9, 0.1), c(0.8, 0.2),
    c(0.1, 0.9), c(0.1, 0.9), c(0.2, 0.8)
  )
  dist_table <- DistanceDistribution(4)
  dists_to_Rg <- AllKendall(r, do.call(rbind, R))
  
  lambda <- UpdateLambda(r, R, z, G = 2, dists_to_Rg, dist_table)
  
  expect_length(lambda, 2)
  expect_true(all(lambda > 0))
})

test_that("UpdateLambda respects top_bound", {
  set.seed(42)
  
  # Data all at modal sequence should push lambda to bound
  r <- rbind(1:4, 1:4)
  R <- list(1:4)
  z <- matrix(1, nrow = 2, ncol = 1)
  dist_table <- DistanceDistribution(4)
  dists_to_Rg <- AllKendall(r, do.call(rbind, R))
  
  lambda <- UpdateLambda(r, R, z, G = 1, dists_to_Rg, dist_table, top_bound = 50)
  
  expect_true(lambda[1] <= 50)
})
