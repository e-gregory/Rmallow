# Tests for EM algorithm components

test_that("UpdateP computes correct proportions", {
  z <- matrix(c(
    0.9, 0.1,
    0.8, 0.2,
    0.2, 0.8,
    0.1, 0.9
  ), nrow = 4, byrow = TRUE)
  
  p <- UpdateP(z)
  
  expect_length(p, 2)
  expect_equal(sum(p), 1)
  expect_equal(p[1], mean(z[, 1]))
  expect_equal(p[2], mean(z[, 2]))
})

test_that("UpdateP validates input", {
  expect_error(UpdateP(c(0.5, 0.5)), "must be a matrix")
})

test_that("Rgen generates correct number of clusters", {
  R <- Rgen(3, abils = 5)
  
  expect_length(R, 3)
  for (i in 1:3) {
    expect_length(R[[i]], 5)
    expect_setequal(R[[i]], 1:5)
  }
})

test_that("Rgen uses hypothesis sequence", {
  hyp <- c(5, 4, 3, 2, 1)
  R <- Rgen(3, hyp = hyp, abils = 5)
  
  expect_equal(R[[1]], hyp)
})

test_that("Rgen validates inputs", {
  expect_error(Rgen(0, abils = 5), "positive integer")
  expect_error(Rgen(3, abils = 0), "positive integer")
  expect_error(Rgen(3, hyp = 1:4, abils = 5), "length equal to")
  expect_error(Rgen(3, hyp = c(1, 1, 3, 4, 5), abils = 5), "permutation")
})

test_that("EStep returns valid probabilities", {
  set.seed(42)
  r <- rbind(1:5, 5:1, c(1, 3, 2, 4, 5))
  R <- list(1:5, 5:1)
  p <- c(0.5, 0.5)
  lambda <- c(1, 1)
  dist_table <- DistanceDistribution(5)
  C <- sapply(lambda, function(l) C_lam(l, dists_table = dist_table))
  
  z <- EStep(R, r, p, lambda, G = 2, N = 3, C)
  
  # All probabilities should be between 0 and 1
  expect_true(all(z >= 0))
  expect_true(all(z <= 1))
  
  # Rows should sum to 1
  expect_equal(rowSums(z), rep(1, 3), tolerance = 1e-10)
  
  # First observation (1:5) should be closer to first cluster
  expect_true(z[1, 1] > z[1, 2])
  
  # Second observation (5:1) should be closer to second cluster
  expect_true(z[2, 2] > z[2, 1])
})

test_that("UpdateR finds correct modal sequences", {
  # Create data strongly favoring two modes
  r <- rbind(
    1:5, 1:5, 1:5, 1:5, 1:5,  # 5 copies of 1:5
    5:1, 5:1, 5:1, 5:1, 5:1   # 5 copies of 5:1
  )
  
  # Membership probabilities strongly assign to correct clusters
  z <- rbind(
    matrix(c(1, 0), nrow = 5, ncol = 2, byrow = TRUE),
    matrix(c(0, 1), nrow = 5, ncol = 2, byrow = TRUE)
  )
  
  R <- UpdateR(r, z)
  
  expect_length(R, 2)
  expect_equal(R[[1]], 1:5)
  expect_equal(R[[2]], 5:1)
})

test_that("Likelihood returns finite value", {
  set.seed(42)
  r <- rbind(1:5, 5:1)
  R <- list(1:5, 5:1)
  p <- c(0.5, 0.5)
  lambda <- c(1, 1)
  dist_table <- DistanceDistribution(5)
  C_lam_vals <- sapply(lambda, function(l) C_lam(l, dists_table = dist_table))
  
  z <- matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2, byrow = TRUE)
  all_dists <- AllKendall(r, do.call(rbind, R))
  
  like <- Likelihood(z, p, C_lam_vals, lambda, all_dists)
  
  expect_true(is.finite(like))
  expect_true(like < 0)  # Log-likelihood should be negative
})
