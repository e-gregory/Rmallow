# Tests for FormatOut function

test_that("FormatOut creates correct output structure", {
  R <- list(1:4, 4:1)
  p <- c(0.6, 0.4)
  lambda <- c(1.5, 2.0)
  z <- matrix(c(0.9, 0.1, 0.2, 0.8, 0.7, 0.3), nrow = 3, byrow = TRUE)
  datas <- rbind(1:4, 4:1, c(1, 2, 4, 3))
  likelihood <- c(-100, -80, -70)
  
  result <- FormatOut(R, p, lambda, z, datas, likelihood)
  
  # Check structure
  expect_type(result, "list")
  expect_named(result, c("R", "p", "lambda", "datas", "min.like"))
  
  # Check values
  expect_equal(result$R, R)
  expect_equal(result$p, p)
  expect_equal(result$lambda, lambda)
  expect_equal(result$min.like, likelihood)
})

test_that("FormatOut adds cluster assignments", {
  R <- list(1:4, 4:1)
  p <- c(0.5, 0.5)
  lambda <- c(1, 1)
  z <- matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2, byrow = TRUE)
  datas <- rbind(1:4, 4:1)
  likelihood <- c(-50)
  
  result <- FormatOut(R, p, lambda, z, datas, likelihood)
  
  # Should have cluster column
  expect_true("clust" %in% names(result$datas))
  expect_equal(result$datas$clust, c(1, 2))
})

test_that("FormatOut adds probability columns", {
  R <- list(1:3, 3:1)
  p <- c(0.5, 0.5)
  lambda <- c(1, 1)
  z <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
  datas <- rbind(1:3, 3:1)
  likelihood <- c(-30)
  
  result <- FormatOut(R, p, lambda, z, datas, likelihood)
  
  expect_true("pvals.1" %in% names(result$datas))
  expect_true("pvals.2" %in% names(result$datas))
  expect_equal(result$datas$pvals.1, c(0.8, 0.3))
  expect_equal(result$datas$pvals.2, c(0.2, 0.7))
})

test_that("FormatOut adds distance columns", {
  R <- list(1:3, 3:1)
  p <- c(0.5, 0.5)
  lambda <- c(1, 1)
  z <- matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2, byrow = TRUE)
  datas <- rbind(1:3, 3:1)
  likelihood <- c(-30)
  
  result <- FormatOut(R, p, lambda, z, datas, likelihood)
  
  expect_true("dists.1" %in% names(result$datas))
  expect_true("dists.2" %in% names(result$datas))
  
  # Distance from 1:3 to 1:3 should be 0
  expect_equal(result$datas$dists.1[1], 0)
  
  # Distance from 3:1 to 1:3 should be 3 (max for n=3)
  expect_equal(result$datas$dists.1[2], 3)
})

test_that("FormatOut adds sequence factor", {
  R <- list(1:3, 3:1)
  p <- c(0.5, 0.5)
  lambda <- c(1, 1)
  z <- matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2, byrow = TRUE)
  datas <- rbind(1:3, 3:1)
  likelihood <- c(-30)
  
  result <- FormatOut(R, p, lambda, z, datas, likelihood)
  
  expect_true("seq" %in% names(result$datas))
  expect_s3_class(result$datas$seq, "factor")
})
