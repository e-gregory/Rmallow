# Tests for DistanceDistribution and related functions

test_that("DistanceDistribution sums to n! for small n", {
  # n=3: 3! = 6
  dist3 <- DistanceDistribution(3)
  expect_equal(sum(dist3), 6)
  
  # n=4: 4! = 24
  dist4 <- DistanceDistribution(4)
  expect_equal(sum(dist4), 24)
  
  # n=5: 5! = 120
  dist5 <- DistanceDistribution(5)
  expect_equal(sum(dist5), 120)
})

test_that("DistanceDistribution has correct length", {
  # Maximum distance for n items is n*(n-1)/2
  # So distribution has n*(n-1)/2 + 1 entries (including 0)
  
  for (n in 3:6) {
    dist <- DistanceDistribution(n)
    expected_len <- n * (n - 1) / 2 + 1
    expect_equal(length(dist), expected_len)
  }
})

test_that("DistanceDistribution is symmetric", {
  # The number of permutations at distance k equals
  # the number at distance max - k
  dist5 <- DistanceDistribution(5)
  
  expect_equal(as.numeric(dist5), rev(as.numeric(dist5)))
})

test_that("DistanceDistribution starts with 1", {
  # There's exactly one permutation at distance 0 (the identity)
  for (n in 2:6) {
    dist <- DistanceDistribution(n)
    expect_equal(as.numeric(dist[1]), 1)
  }
})

test_that("DistanceDistribution validates input", {
  expect_error(DistanceDistribution(0))
  expect_error(DistanceDistribution(-1))
  expect_error(DistanceDistribution(1.5))
  expect_error(DistanceDistribution("a"))
})

test_that("DistanceDistribution handles n=1", {
  dist1 <- DistanceDistribution(1)
  expect_equal(sum(dist1), 1)
  expect_equal(length(dist1), 1)
})

test_that("DistanceDistribution handles n=2", {
  dist2 <- DistanceDistribution(2)
  expect_equal(sum(dist2), 2)
  expect_equal(length(dist2), 2)
  expect_equal(as.numeric(dist2), c(1, 1))
})

test_that("NextTable produces correct results", {
  # Test that NextTable gives same results as direct computation
  dist4 <- DistanceDistribution(4)
  dist5_direct <- DistanceDistribution(5)
  dist5_next <- NextTable(dist4, 4)
  
  expect_equal(as.numeric(dist5_next), as.numeric(dist5_direct))
})

test_that("SeqDistribution matches DistanceDistribution for small n", {
  for (n in 3:5) {
    seq_dist <- SeqDistribution(n)
    dist_dist <- DistanceDistribution(n)
    
    # Tables should match
    expect_equal(as.numeric(table(seq_dist)), as.numeric(dist_dist))
  }
})
