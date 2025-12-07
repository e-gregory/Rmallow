#!/usr/bin/env Rscript
#' Quick Benchmark: R vs Rcpp Performance in RMallow Package

library(RMallow)

cat("\n")
cat("======================================================================\n")
cat("  RMallow Package: R vs Rcpp Performance Comparison\n")
cat("======================================================================\n\n")

benchmark <- function(name, fn, args_cpp, args_r, times = 5) {
  cat(sprintf("\n%s:\n", name))
  
  # Warm up
  invisible(do.call(fn, args_cpp))
  invisible(do.call(fn, args_r))
  
  # Benchmark C++
  cpp_time <- system.time({
    for(i in 1:times) do.call(fn, args_cpp)
  })[["elapsed"]] / times
  
  # Benchmark R
  r_time <- system.time({
    for(i in 1:times) do.call(fn, args_r)
  })[["elapsed"]] / times
  
  speedup <- r_time / cpp_time
  
  cat(sprintf("  C++ (Rcpp): %.4f sec\n", cpp_time))
  cat(sprintf("  Pure R:     %.4f sec\n", r_time))
  cat(sprintf("  Speedup:    %.1fx faster\n", speedup))
  
  return(speedup)
}

speedups <- c()

# Test 1: KendallInfo
cat("\n--- KendallInfo (1000 obs x 15 items) ---")
set.seed(42)
r <- t(replicate(1000, sample(15)))
s <- benchmark("KendallInfo", KendallInfo,
               list(r = r, use_cpp = TRUE),
               list(r = r, use_cpp = FALSE))
speedups <- c(speedups, s)

# Test 2: AllKendall
cat("\n--- AllKendall (500 obs x 15 items, 10 seqs) ---")
set.seed(42)
r <- t(replicate(500, sample(15)))
seqs <- t(replicate(10, sample(15)))
s <- benchmark("AllKendall", AllKendall,
               list(r = r, seqs = seqs, use_cpp = TRUE),
               list(r = r, seqs = seqs, use_cpp = FALSE))
speedups <- c(speedups, s)

# Test 3: AllSeqDists
cat("\n--- AllSeqDists (1000 obs x 15 items) ---")
set.seed(42)
seqs <- t(replicate(1000, sample(15)))
s <- benchmark("AllSeqDists", AllSeqDists,
               list(seqs = seqs, use_cpp = TRUE),
               list(seqs = seqs, use_cpp = FALSE))
speedups <- c(speedups, s)

# Test 4: UpdateR
cat("\n--- UpdateR (500 obs x 15 items, 3 clusters) ---")
set.seed(42)
r <- t(replicate(500, sample(15)))
z <- matrix(runif(500 * 3), nrow = 500, ncol = 3)
z <- z / rowSums(z)
s <- benchmark("UpdateR", UpdateR,
               list(r = r, z = z, use_cpp = TRUE),
               list(r = r, z = z, use_cpp = FALSE))
speedups <- c(speedups, s)

# Test 5: EStep
cat("\n--- EStep (500 obs x 10 items, 3 clusters) ---")
set.seed(42)
n_items <- 10
r <- t(replicate(500, sample(n_items)))
R <- lapply(1:3, function(i) sample(n_items))
p <- rep(1/3, 3)
lambda <- runif(3, 0.5, 2)
dists_table <- DistanceDistribution(n_items)
C <- sapply(lambda, function(l) C_lam(l, dists_table = dists_table))
s <- benchmark("EStep", EStep,
               list(R = R, r = r, p = p, lambda = lambda, G = 3, N = 500, C = C, use_cpp = TRUE),
               list(R = R, r = r, p = p, lambda = lambda, G = 3, N = 500, C = C, use_cpp = FALSE))
speedups <- c(speedups, s)

# Test 6: SimplifySequences
cat("\n--- SimplifySequences (1000 obs x 15 items) ---")
set.seed(42)
rankings <- t(replicate(1000, sample(30, 15)))
s <- benchmark("SimplifySequences", SimplifySequences,
               list(rankings = rankings, use_cpp = TRUE),
               list(rankings = rankings, use_cpp = FALSE))
speedups <- c(speedups, s)

# Test 7: Full Mallows model
cat("\n--- Mallows Model (200 obs x 8 items, 2 clusters, 20 iter) ---")
set.seed(42)
data_matrix <- t(replicate(200, sample(8)))
s <- benchmark("Mallows (full model)", Mallows,
               list(datas = data_matrix, G = 2, iter = 20, verbose = FALSE, use_cpp = TRUE),
               list(datas = data_matrix, G = 2, iter = 20, verbose = FALSE, use_cpp = FALSE),
               times = 3)
speedups <- c(speedups, s)

# Summary
cat("\n\n")
cat("======================================================================\n")
cat("  SUMMARY\n")
cat("======================================================================\n")
cat(sprintf("\nAverage speedup across all functions: %.1fx faster\n", mean(speedups)))
cat(sprintf("Range of speedups: %.1fx - %.1fx\n", min(speedups), max(speedups)))
cat("\n")
cat("The C++ implementations via Rcpp provide significant performance\n")
cat("improvements, especially for larger datasets and more complex models.\n")
cat("\n")
