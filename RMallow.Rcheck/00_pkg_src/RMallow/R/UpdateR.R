#' Update modal sequences (M-step)
#'
#' Updates the cluster centers (modal sequences) to maximize the likelihood
#' given current membership probabilities.
#'
#' @param r Matrix of rankings (N x n_items).
#' @param z Matrix of membership probabilities (N x G).
#' @param infos Optional precomputed Kendall information matrix for \code{r}.
#' @param use_cpp Logical; if TRUE (default), use the fast C++ implementation.
#' @return List of G updated modal sequences.
#' @author Erik Gregory
#' @references Murphy, T.B. and Martin, D. (2003). Mixtures of distance-based
#'   models for ranking data. Computational Statistics & Data Analysis,
#'   41, 645-655.
#' @keywords cluster center
#' @export
UpdateR <- function(r, z, infos = NULL, use_cpp = TRUE) {
  # Ensure inputs are matrices
  if (!is.matrix(z)) {
    stop("'z' must be a matrix.")
  }
  if (!is.matrix(r)) {
    r <- as.matrix(r)
  }

  G <- ncol(z)  # Number of clusters
  n_abils <- ncol(r)  # Number of items

  # Compute Kendall information if not provided
  if (is.null(infos)) {
    infos <- KendallInfo(r, use_cpp = use_cpp)
  }

  # Use C++ implementation for the weighted sums
  if (use_cpp) {
    sums <- update_r_sums_cpp(infos, z)
    zero_sums <- sums$zero_sums
    one_sums <- sums$one_sums

    # Determine preference direction
    prefs <- zero_sums - one_sums
    prefs <- ifelse(prefs < 0, 1, 0)

    # Reconstruct sequences from preferences
    return(construct_seqs_cpp(prefs, n_abils))
  }

  # Fallback to R implementation
  UpdateR_R(r, z, infos)
}

#' Pure R implementation of UpdateR (for benchmarking)
#' @keywords internal
UpdateR_R <- function(r, z, infos = NULL) {
  if (!is.matrix(z)) {
    stop("'z' must be a matrix.")
  }
  if (!is.matrix(r)) {
    r <- as.matrix(r)
  }

  G <- ncol(z)
  N <- nrow(r)
  n_abils <- ncol(r)

  if (is.null(infos)) {
    infos <- KendallInfo_R(r)
  }

  n_pairs <- ncol(infos)

  zero_sums <- matrix(0, nrow = G, ncol = n_pairs)
  one_sums <- matrix(0, nrow = G, ncol = n_pairs)

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

  prefs <- zero_sums - one_sums
  prefs <- ifelse(prefs < 0, 1, 0)

  ConstructSeqs(prefs, n_abils, use_cpp = FALSE)
}
