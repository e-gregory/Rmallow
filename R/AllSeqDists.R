#' Compute distances from sequences to the canonical ordering
#'
#' Calculates the Kendall distance from each sequence to the canonical
#' ordering (1, 2, ..., N).
#'
#' @param seqs Matrix or data frame of sequences, where each row is a ranking.
#' @return Integer vector of distances from each sequence to \code{1:N}.
#' @author Erik Gregory
#' @keywords Kendall Distance
#' @export
#' @examples
#' seqs <- rbind(1:5, 5:1, c(1, 3, 2, 4, 5))
#' AllSeqDists(seqs)
#' # Returns: 0, 10, 1
AllSeqDists <- function(seqs) {
  # Ensure seqs is a matrix
  if (!is.matrix(seqs)) {
    seqs <- as.matrix(seqs)
  }

  if (nrow(seqs) == 0L) {
    return(integer(0L))
  }

  # Compute Kendall information
  infos <- KendallInfo(seqs)

  # Distance = number of inversions = count of 1s (decreases)
  # Using sum instead of which() for efficiency
  if (is.matrix(infos) && ncol(infos) > 0L) {
    dists <- rowSums(infos == 1, na.rm = TRUE)
  } else if (length(infos) > 0L) {
    dists <- sum(infos == 1, na.rm = TRUE)
  } else {
    dists <- rep(0L, nrow(seqs))
  }

  as.integer(dists)
}
