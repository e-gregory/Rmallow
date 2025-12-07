#' Update modal sequences (M-step)
#'
#' Updates the cluster centers (modal sequences) to maximize the likelihood
#' given current membership probabilities.
#'
#' @param r Matrix of rankings (N x n_items).
#' @param z Matrix of membership probabilities (N x G).
#' @param infos Optional precomputed Kendall information matrix for \code{r}.
#' @return List of G updated modal sequences.
#' @author Erik Gregory
#' @references Murphy, T.B. and Martin, D. (2003). Mixtures of distance-based
#'   models for ranking data. Computational Statistics & Data Analysis,
#'   41, 645-655.
#' @keywords cluster center
#' @export
UpdateR <- function(r, z, infos = NULL) {
  # Ensure inputs are matrices
  if (!is.matrix(z)) {
    stop("'z' must be a matrix.")
  }
  if (!is.matrix(r)) {
    r <- as.matrix(r)
  }

  G <- ncol(z)  # Number of clusters
  N <- nrow(r)  # Number of observations
  n_abils <- ncol(r)  # Number of items

  # Compute Kendall information if not provided
  if (is.null(infos)) {
    infos <- KendallInfo(r)
  }

  n_pairs <- ncol(infos)  # Number of pairwise comparisons

  # Initialize matrices for weighted sums
  zero_sums <- matrix(0, nrow = G, ncol = n_pairs)
  one_sums <- matrix(0, nrow = G, ncol = n_pairs)

  # For each pair comparison, compute weighted counts of 0s and 1s
  for (i in seq_len(n_pairs)) {
    zero_idx <- which(infos[, i] == 0)
    one_idx <- which(infos[, i] == 1)

    if (length(zero_idx) > 0L) {
      if (G > 1L) {
        zero_sums[, i] <- colSums(z[zero_idx, , drop = FALSE])
      } else {
        zero_sums[, i] <- sum(z[zero_idx, ])
      }
    }

    if (length(one_idx) > 0L) {
      if (G > 1L) {
        one_sums[, i] <- colSums(z[one_idx, , drop = FALSE])
      } else {
        one_sums[, i] <- sum(z[one_idx, ])
      }
    }
  }

  # Determine preference direction for each cluster and pair
  # Positive diff = more 0s (increases), prefer 0
  # Negative diff = more 1s (decreases), prefer 1
  prefs <- zero_sums - one_sums

  # Convert to binary preferences (0 or 1)
  # When tied (prefs == 0), default to 0 (arbitrary but consistent)
  prefs <- ifelse(prefs < 0, 1, 0)

  # Reconstruct sequences from preferences
  ConstructSeqs(prefs, n_abils)
}
