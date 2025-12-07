#' Construct sequences from Kendall preference vectors
#'
#' Reconstructs ranking sequences from their Kendall information representation.
#' Each fully-ordered sequence has a unique Kendall information vector.
#'
#' @param prefs Matrix of pairwise preferences, where each row represents
#'   one sequence's preferences. A value of 1 indicates a decrease (first
#'   element > second), 0 indicates an increase.
#' @param n_abils Number of items in the original ranking.
#' @return List of integer vectors, each representing a reconstructed ranking.
#' @author Erik Gregory
#' @keywords Sequences
#' @export
#' @examples
#' # Reconstruct a sequence where item 4 is ranked first
#' prefs <- matrix(c(1, 1, 1, 0, 0, 0), nrow = 1)
#' ConstructSeqs(prefs, 4)
#' # Returns list(c(4, 1, 2, 3))
#'
#' # Multiple sequences at once
#' prefs2 <- rbind(c(0, 0, 0, 0, 0, 0), c(1, 1, 1, 1, 1, 1))
#' ConstructSeqs(prefs2, 4)
#' # Returns list(c(1,2,3,4), c(4,3,2,1))
ConstructSeqs <- function(prefs, n_abils) {
  # Input validation
  if (!is.numeric(n_abils) || n_abils < 2L) {
    stop("'n_abils' must be an integer >= 2.")
  }
  n_abils <- as.integer(n_abils)

  # Handle single sequence (vector) input
  if (!is.matrix(prefs)) {
    prefs <- matrix(prefs, nrow = 1L)
  }

  n_seqs <- nrow(prefs)

  # Expected number of pairwise comparisons: n*(n-1)/2
  expected_cols <- (n_abils * (n_abils - 1L)) %/% 2L
  if (ncol(prefs) != expected_cols) {
    stop(sprintf("'prefs' should have %d columns for n_abils = %d", expected_cols, n_abils))
  }

  # Initialize result list
  R <- vector("list", n_seqs)

  # For each sequence, maintain available items
  nums <- lapply(seq_len(n_seqs), function(j) seq_len(n_abils))

  # Positions where each item's comparisons start in the preference vector
  # Item i is compared with items (i+1), (i+2), ..., n_abils
  # cumulative positions: 0, n-1, n-1 + n-2, ...
  tops <- c(0L, cumsum(rev(seq(2L, n_abils)) - 1L))

  # Build sequences item by item
  for (i in seq_len(n_abils - 1L)) {
    start_col <- tops[i] + 1L
    end_col <- tops[i + 1L]

    # Sum of preferences determines position in remaining items
    if (start_col == end_col) {
      sumz <- prefs[, start_col]
    } else {
      sumz <- rowSums(prefs[, start_col:end_col, drop = FALSE])
    }

    # For each sequence, pick the appropriate item and remove it
    for (j in seq_len(n_seqs)) {
      # Position in remaining items (1-indexed)
      pos <- sumz[j] + 1L
      item <- nums[[j]][pos]

      # Add to result
      if (i == 1L) {
        R[[j]] <- item
      } else {
        R[[j]] <- c(R[[j]], item)
      }

      # Remove used item
      nums[[j]] <- nums[[j]][-pos]

      # Add final item if this is the second-to-last iteration
      if (length(R[[j]]) == n_abils - 1L) {
        R[[j]] <- c(R[[j]], nums[[j]])
      }
    }
  }

  R
}

# Backward compatibility alias
#' @rdname ConstructSeqs
#' @param n.abils Deprecated. Use \code{n_abils} instead.
ConstructSeqs.legacy <- function(prefs, n.abils) {
  ConstructSeqs(prefs, n_abils = n.abils)
}
