# Tests for AllSeqDists function

test_that("AllSeqDists returns 0 for canonical sequence", {
  seqs <- matrix(1:5, nrow = 1)
  dists <- AllSeqDists(seqs)
  
  expect_equal(dists[1], 0)
})

test_that("AllSeqDists returns maximum for reversed sequence", {
  seqs <- matrix(5:1, nrow = 1)
  dists <- AllSeqDists(seqs)
  
  # Maximum distance is n*(n-1)/2 = 10 for n=5
  expect_equal(dists[1], 10)
})

test_that("AllSeqDists handles multiple sequences", {
  seqs <- rbind(1:5, 5:1, c(1, 3, 2, 4, 5))
  dists <- AllSeqDists(seqs)
  
  expect_length(dists, 3)
  expect_equal(dists[1], 0)
  expect_equal(dists[2], 10)
  expect_equal(dists[3], 1)
})

test_that("AllSeqDists returns integers", {
  seqs <- rbind(1:4, 4:1)
  dists <- AllSeqDists(seqs)
  
  expect_type(dists, "integer")
})

test_that("AllSeqDists handles empty input", {
  seqs <- matrix(nrow = 0, ncol = 5)
  dists <- AllSeqDists(seqs)
  
  expect_length(dists, 0)
})
