#
#' Compute distance distribution for (N+1)! from N! space
#'
#' Given the Kendall distance distribution in N! space, computes the
#' distribution in (N+1)! space. This uses the recurrence relation for
#' Mahonian numbers (number of permutations with k inversions).
#'
#' @param last_table Table or numeric vector of distance counts in N! space.
#' @param n_last The value of N (current space is N!).
#' @param use_cpp Logical; if TRUE (default), use the fast C++ implementation.
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
NextTable <- function(last_table, n_last, use_cpp = TRUE) {
  # Input validation
  if (n_last < 1L) {
    stop("'n_last' must be at least 1.")
  }

  last_table <- as.numeric(last_table)

  # Use C++ implementation
  if (use_cpp) {
    return(next_table_cpp(last_table, n_last))
  }

  # Fallback to R implementation
  NextTable_R(last_table, n_last)
}

#' Pure R implementation of NextTable (for benchmarking)
#' @keywords internal
NextTable_R <- function(last_table, n_last) {
  if (n_last < 1L) {
    stop("'n_last' must be at least 1.")
  }

  last_table <- as.numeric(last_table)

  new_len <- ((n_last + 1L) * n_last) %/% 2L + 1L
  middle_len <- new_len - 2L * (n_last + 1L)

  cumsum_last <- cumsum(last_table)

  first_part <- cumsum_last[seq_len(n_last + 1L)]

  if (middle_len > 0L) {
    end_idx <- (n_last + 2L):(n_last + 1L + middle_len)
    start_idx <- seq_len(middle_len)
    middle_part <- cumsum_last[end_idx] - cumsum_last[start_idx]
  } else {
    middle_part <- numeric(0L)
  }

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
