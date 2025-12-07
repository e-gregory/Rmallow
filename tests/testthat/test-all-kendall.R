# Tests for AllKendall function

test_that("AllKendall computes correct distances", {
  # Distance from 1:5 to itself should be 0
  r <- matrix(1:5, nrow = 1)
  seqs <- matrix(1:5, nrow = 1)
  
  dists <- AllKendall(r, seqs)
  expect_equal(dists[1, 1], 0)
})

test_that("AllKendall computes maximum distance correctly", {
  # Distance from 1:5 to 5:1 should be n*(n-1)/2 = 10
  r <- matrix(1:5, nrow = 1)
  seqs <- matrix(5:1, nrow = 1)
  
  dists <- AllKendall(r, seqs)
  expect_equal(dists[1, 1], 10)
})

test_that("AllKendall handles multiple observations", {
  r <- rbind(1:5, 5:1, c(1, 3, 2, 4, 5))
  seqs <- matrix(1:5, nrow = 1)
  
  dists <- AllKendall(r, seqs)
  
  expect_equal(nrow(dists), 3)
  expect_equal(ncol(dists), 1)
  expect_equal(dists[1, 1], 0)   # 1:5 to 1:5
  expect_equal(dists[2, 1], 10)  # 5:1 to 1:5
  expect_equal(dists[3, 1], 1)   # one swap from 1:5
})

test_that("AllKendall handles multiple reference sequences", {
  r <- matrix(1:5, nrow = 1)
  seqs <- rbind(1:5, 5:1)
  
  dists <- AllKendall(r, seqs)
  
  expect_equal(nrow(dists), 1)
  expect_equal(ncol(dists), 2)
  expect_equal(dists[1, 1], 0)
  expect_equal(dists[1, 2], 10)
})

test_that("AllKendall is symmetric in a sense", {
  # d(a, b) should equal d(b, a)
  a <- matrix(c(1, 3, 2, 4, 5), nrow = 1)
  b <- matrix(c(2, 1, 3, 5, 4), nrow = 1)
  
  d_ab <- AllKendall(a, b)[1, 1]
  d_ba <- AllKendall(b, a)[1, 1]
  
  expect_equal(d_ab, d_ba)
})

test_that("AllKendall validates dimension mismatch", {
  r <- matrix(1:5, nrow = 1)
  seqs <- matrix(1:4, nrow = 1)
  
  expect_error(AllKendall(r, seqs), "same number of columns")
})

test_that("AllKendall accepts precomputed data_info", {
  r <- rbind(1:5, 5:1)
  seqs <- matrix(1:5, nrow = 1)
  
  # Precompute info
  data_info <- KendallInfo(r)
  
  # Should give same results
  dists1 <- AllKendall(r, seqs)
  dists2 <- AllKendall(r, seqs, data_info = data_info)
  
  expect_equal(dists1, dists2)
})

test_that("AllKendall handles single reference sequence correctly", {
  # This tests the bug fix for n_seqs == 1
  r <- rbind(1:4, 4:1, c(1, 2, 4, 3))
  seqs <- matrix(1:4, nrow = 1)
  
  dists <- AllKendall(r, seqs)
  
  expect_equal(dists[1, 1], 0)
  expect_equal(dists[2, 1], 6)  # max distance for n=4
  expect_equal(dists[3, 1], 1)
})
