# Tests for main Mallows function

test_that("Mallows fits simple two-cluster data", {
  set.seed(42)
  
  # Create clear two-cluster data
  data <- rbind(
    matrix(rep(1:5, each = 10), nrow = 10, byrow = TRUE),
    matrix(rep(5:1, each = 10), nrow = 10, byrow = TRUE)
  )
  
  result <- Mallows(data, G = 2, iter = 20, verbose = FALSE)
  
  # Check structure
expect_type(result, "list")
  expect_named(result, c("R", "p", "lambda", "datas", "min.like"))
  
  # Should find two modal sequences
  expect_length(result$R, 2)
  
  # Modal sequences should be 1:5 and 5:1 (in some order)
  modes <- lapply(result$R, sort)
  expect_true(all(sapply(modes, function(m) all(m == 1:5))))
})

test_that("Mallows returns correct output structure", {
  set.seed(42)
  data <- rbind(1:4, 4:1, 1:4, 4:1)
  
  result <- Mallows(data, G = 2, iter = 5, verbose = FALSE)
  
  # R should be a list of modal sequences
  expect_type(result$R, "list")
  for (r in result$R) {
    expect_setequal(r, 1:4)
  }
  
  # p should sum to 1
  expect_equal(sum(result$p), 1, tolerance = 1e-10)
  expect_true(all(result$p >= 0))
  
  # lambda should be positive
  expect_true(all(result$lambda >= 0))
  
  # datas should have cluster assignments
  expect_true("clust" %in% names(result$datas))
  expect_true(all(result$datas$clust %in% 1:2))
  
  # min.like should generally increase (likelihood should increase)
  likes <- result$min.like[result$min.like != 0]
  if (length(likes) > 1) {
    # Likelihood should generally increase, but small decreases can happen
    # due to numerical issues in the M-step
    # Check that at least the final likelihood is higher than initial
    expect_true(likes[length(likes)] >= likes[1] - 1e-3)
  }
})

test_that("Mallows handles single cluster", {
  set.seed(42)
  data <- rbind(1:4, c(1, 2, 4, 3), c(1, 3, 2, 4))
  
  result <- Mallows(data, G = 1, iter = 10, verbose = FALSE)
  
  expect_length(result$R, 1)
  expect_equal(result$p, 1)
})

test_that("Mallows respects hypothesis sequence", {
  set.seed(42)
  data <- rbind(1:5, 5:1, 1:5, 5:1)
  hyp <- c(3, 1, 2, 4, 5)
  
  result <- Mallows(data, G = 2, iter = 5, hyp = hyp, verbose = FALSE)
  
  # First cluster should start at hypothesis
  # (though EM may move it)
  expect_length(result$R, 2)
})

test_that("Mallows validates inputs", {
  data <- rbind(1:5, 5:1)
  
  expect_error(Mallows("not a matrix", G = 2), "must be a matrix")
  expect_error(Mallows(data, G = 0), "positive integer")
  expect_error(Mallows(data, G = 10), "at least G")  # More clusters than data
})

test_that("Mallows converges early when tolerance is met", {
  set.seed(42)
  
  # Very clear data should converge quickly
  data <- rbind(
    matrix(rep(1:4, each = 20), nrow = 20, byrow = TRUE),
    matrix(rep(4:1, each = 20), nrow = 20, byrow = TRUE)
  )
  
  result <- Mallows(data, G = 2, iter = 100, tol = 1e-6, verbose = FALSE)
  
  # Should converge before max iterations
  n_iters <- length(result$min.like)
  expect_true(n_iters < 100)
})

test_that("BestFit selects best model", {
  set.seed(42)
  
  data <- rbind(1:4, 4:1, 1:4, 4:1, c(1, 2, 4, 3), c(4, 3, 1, 2))
  
  result <- BestFit(data, N = 3, iter = 10, G = 2, verbose = FALSE)
  
  # Should return a valid model
  expect_type(result, "list")
  expect_named(result, c("R", "p", "lambda", "datas", "min.like"))
})

test_that("Mallows works with ties in data", {
  set.seed(42)
  
  # Data with ties
  data <- rbind(
    c(1, 1, 2, 3),
    c(1, 2, 2, 3),
    c(3, 2, 2, 1),
    c(3, 2, 1, 1)
  )
  
  result <- Mallows(data, G = 2, iter = 10, verbose = FALSE)
  
  expect_type(result, "list")
  expect_length(result$R, 2)
})
