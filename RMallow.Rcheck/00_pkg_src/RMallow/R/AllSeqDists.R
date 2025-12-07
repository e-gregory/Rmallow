#' Compute distances from sequences to the canonical ordering
#'
#' Calculates the Kendall distance from each sequence to the canonical
#' ordering (1, 2, ..., N).
#'
#' @param seqs Matrix or data frame of sequences, where each row is a ranking.
#' @param use_cpp Logical; if TRUE (default), use the fast C++ implementation.
#' @return Integer vector of distances from each sequence to \code{1:N}.
#' @author Erik Gregory
#' @keywords Kendall Distance
#' @export
#' @examples
#' seqs <- rbind(1:5, 5:1, c(1, 3, 2, 4, 5))
#' AllSeqDists(seqs)
#' # Returns: 0, 10, 1
AllSeqDists <- function(seqs, use_cpp = TRUE) {
  # Ensure seqs is a matrix
  if (!is.matrix(seqs)) {
    seqs <- as.matrix(seqs)
  }

  if (nrow(seqs) == 0L) {
    return(integer(0L))
  }

  # Use C++ implementation
  if (use_cpp) {
    infos <- KendallInfo(seqs, use_cpp = TRUE)
    if (is.matrix(infos) && ncol(infos) > 0L) {
      return(all_seq_dists_cpp(infos))
    } else {
      return(rep(0L, nrow(seqs)))
    }
  }

  # Fallback to R implementation
  AllSeqDists_R(seqs)
}

#' Pure R implementation of AllSeqDists (for benchmarking)
#' @keywords internal
AllSeqDists_R <- function(seqs) {
  if (!is.matrix(seqs)) {
    seqs <- as.matrix(seqs)
  }

  if (nrow(seqs) == 0L) {
    return(integer(0L))
  }

  # Compute Kendall information using R version
  infos <- KendallInfo_R(seqs)

  # Distance = number of inversions = count of 1s (decreases)
  if (is.matrix(infos) && ncol(infos) > 0L) {
    dists <- rowSums(infos == 1, na.rm = TRUE)
  } else if (length(infos) > 0L) {
    dists <- sum(infos == 1, na.rm = TRUE)
  } else {
    dists <- rep(0L, nrow(seqs))
  }

  as.integer(dists)
}
