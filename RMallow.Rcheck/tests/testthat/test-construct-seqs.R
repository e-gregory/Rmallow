# Tests for ConstructSeqs function

test_that("ConstructSeqs reconstructs identity from all zeros", {
  # All zeros means all increases -> identity permutation
  prefs <- matrix(rep(0, 6), nrow = 1)  # 4 items -> 6 pairs
  result <- ConstructSeqs(prefs, 4)
  
  expect_equal(result[[1]], 1:4)
})

test_that("ConstructSeqs reconstructs reverse from all ones", {
  # All ones means all decreases -> reverse permutation
  prefs <- matrix(rep(1, 6), nrow = 1)
  result <- ConstructSeqs(prefs, 4)
  
  expect_equal(result[[1]], 4:1)
})

test_that("ConstructSeqs handles multiple sequences", {
  prefs <- rbind(
    rep(0, 6),  # identity
    rep(1, 6)   # reverse
  )
  result <- ConstructSeqs(prefs, 4)
  
  expect_length(result, 2)
  expect_equal(result[[1]], 1:4)
  expect_equal(result[[2]], 4:1)
})

test_that("ConstructSeqs validates column count", {
  # 4 items needs 6 columns, not 5
  prefs <- matrix(rep(0, 5), nrow = 1)
  expect_error(ConstructSeqs(prefs, 4), "should have 6 columns")
})

test_that("ConstructSeqs validates n_abils", {
  expect_error(ConstructSeqs(matrix(0), 1), "must be an integer >= 2")
  expect_error(ConstructSeqs(matrix(0), 0), "must be an integer >= 2")
})

test_that("ConstructSeqs handles 3 items correctly", {
  # 3 items -> 3 pairs
  # All zeros: 1 2 3
  prefs <- matrix(c(0, 0, 0), nrow = 1)
  result <- ConstructSeqs(prefs, 3)
  expect_equal(result[[1]], 1:3)
  
  # All ones: 3 2 1
  prefs <- matrix(c(1, 1, 1), nrow = 1)
  result <- ConstructSeqs(prefs, 3)
  expect_equal(result[[1]], 3:1)
})

test_that("ConstructSeqs is inverse of KendallInfo for full rankings", {
  # Generate a random permutation
  set.seed(42)
  perm <- sample(5)
  
  # Get Kendall info
  info <- KendallInfo(matrix(perm, nrow = 1))
  
  # Convert to binary preferences (info uses 0/1 already for full rankings)
  prefs <- matrix(as.numeric(info), nrow = 1)
  
  # Reconstruct
  result <- ConstructSeqs(prefs, 5)
  
  expect_equal(result[[1]], perm)
})
