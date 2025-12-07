#' Compute Kendall's distances between two sets of rankings
#'
#' Calculates all Kendall's distances between each ranking in \code{r} and
#' each sequence in \code{seqs}.
#'
#' @param r Matrix of rankings (each row is a ranking).
#' @param seqs Matrix of reference sequences (each row is a sequence).
#' @param data_info Optional precomputed Kendall information matrix for \code{r}.
#'   If NULL, it will be computed internally. Providing this can improve
#'   efficiency when calling the function repeatedly with the same data.
#' @return Matrix where entry \code{[i, j]} represents the Kendall distance
#'   from ranking \code{i} in \code{r} to sequence \code{j} in \code{seqs}.
#' @author Erik Gregory
#' @keywords Kendall distance
#' @export
#' @examples
#' data1 <- rbind(1:5, 5:1, c(3, 2, 1, 4, 5))
#' data2 <- rbind(1:5, 5:1)
#' AllKendall(data1, data2)
AllKendall <- function(r, seqs, data_info = NULL) {
  # Ensure inputs are matrices
  if (!is.matrix(r)) {
    r <- as.matrix(r)
  }
  if (!is.matrix(seqs)) {
    seqs <- as.matrix(seqs)
  }

  # Validate dimensions
  if (ncol(r) != ncol(seqs)) {
    stop("'r' and 'seqs' must have the same number of columns.")
  }

  n_obs <- nrow(r)
  n_seqs <- nrow(seqs)

  # Initialize distance matrix
  dists <- matrix(0, nrow = n_obs, ncol = n_seqs)

  # Generate column pair indices once
  inds <- combn(ncol(r), 2L)

  # Compute Kendall information for data if not provided
  if (is.null(data_info)) {
    data_info <- KendallInfo(r, inds)
  }

  # Compute Kendall information for reference sequences
  seqs_info <- KendallInfo(seqs, inds)

  # Compute distances using vectorized operations
  # Distance is sum of absolute differences, ignoring NAs (ties)
  if (n_seqs == 1L) {
    # Single reference sequence - seqs_info is a vector
    dists[, 1L] <- rowSums(abs(sweep(data_info, 2L, seqs_info, "-")), na.rm = TRUE)
  } else {
    # Multiple reference sequences
    for (j in seq_len(n_seqs)) {
      dists[, j] <- rowSums(abs(sweep(data_info, 2L, seqs_info[j, ], "-")), na.rm = TRUE)
    }
  }

  dists
}

# Backward compatibility alias (deprecated parameter name)
#' @rdname AllKendall
#' @param data.info Deprecated. Use \code{data_info} instead.
AllKendall.legacy <- function(r, seqs, data.info = NULL) {
  AllKendall(r, seqs, data_info = data.info)
}
