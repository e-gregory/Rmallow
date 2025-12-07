#
#' Compute distance distribution for (N+1)! from N! space
#'
#' Given the Kendall distance distribution in N! space, computes the
#' distribution in (N+1)! space. This uses the recurrence relation for
#' Mahonian numbers (number of permutations with k inversions).
#'
#' @param last_table Table or numeric vector of distance counts in N! space.
#' @param n_last The value of N (current space is N!).
#' @return Numeric vector of distance counts in (N+1)! space.
#' @author Erik Gregory
#' @keywords bubblesort Kendall
#' @export
#' @examples
#' # Get distribution for 5! from 4!
#' tab4 <- DistanceDistribution(4)
#' tab5 <- NextTable(tab4, 4)
#' # Verify: should have 5! = 120 total permutations
#' sum(tab5)
NextTable <- function(last_table, n_last) {
  # Input validation
  if (n_last < 1L) {
    stop("'n_last' must be at least 1.")
  }

  last_table <- as.numeric(last_table)

  # Maximum distance in (N+1)! space is N*(N+1)/2
  # (adding element N+1 can add at most N inversions to any permutation)
  new_len <- ((n_last + 1L) * n_last) %/% 2L + 1L

  # Length of the middle section (not easily computed from boundaries)
  middle_len <- new_len - 2L * (n_last + 1L)

  # The recurrence relation for Mahonian numbers:
  # T(n+1, k) = sum of T(n, k-j) for j = 0 to min(k, n)
  #
  # This can be computed efficiently using cumulative sums:
  # - First (n+1) entries: cumsum of first (n+1) entries of last_table
  # - Middle section: sliding window sums of width (n+1)

  # - Last (n+1) entries: reverse cumsum (symmetric)

  cumsum_last <- cumsum(last_table)

  # First n+1 entries
  first_part <- cumsum_last[seq_len(n_last + 1L)]

  # Middle section (sliding window of width n+1)
  if (middle_len > 0L) {
    end_idx <- (n_last + 2L):(n_last + 1L + middle_len)
    start_idx <- seq_len(middle_len)
    middle_part <- cumsum_last[end_idx] - cumsum_last[start_idx]
  } else {
    middle_part <- numeric(0L)
  }

  # Last n+1 entries (symmetric to first part)
  last_part <- rev(first_part)

  c(first_part, middle_part, last_part)
}

# Backward compatibility alias
#' @rdname NextTable
#' @param last.table Deprecated. Use \code{last_table} instead.
#' @param N.last Deprecated. Use \code{n_last} instead.
NextTable.legacy <- function(last.table, N.last) {
  NextTable(last_table = last.table, n_last = N.last)
}
