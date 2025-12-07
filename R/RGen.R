#' Initialize cluster modal sequences
#'
#' Generates initial modal sequences for the Mallows mixture model clustering.
#' If a hypothesis sequence is provided, it becomes the first cluster center.
#' Remaining cluster centers are initialized randomly.
#'
#' @param G Number of cluster centers (modes).
#' @param hyp Optional hypothesis sequence to use as the first cluster center.
#'   Must be a permutation of 1:abils.
#' @param abils Number of items being ranked.
#' @return List of G integer vectors, each a permutation of 1:abils.
#' @author Erik Gregory
#' @keywords initialization cluster
#' @export
#' @examples
#' # Three random cluster centers for 5-item rankings
#' set.seed(42)
#' Rgen(3, abils = 5)
#'
#' # With a specific hypothesis sequence
#' Rgen(3, hyp = c(5, 4, 3, 2, 1), abils = 5)
Rgen <- function(G, hyp = NULL, abils) {
  # Input validation
  if (!is.numeric(G) || G < 1L) {
    stop("'G' must be a positive integer.")
  }
  G <- as.integer(G)

  if (!is.numeric(abils) || abils < 1L) {
    stop("'abils' must be a positive integer.")
  }
  abils <- as.integer(abils)

  if (!is.null(hyp)) {
    if (length(hyp) != abils) {
      stop("'hyp' must have length equal to 'abils'.")
    }
    if (!setequal(hyp, seq_len(abils))) {
      stop("'hyp' must be a permutation of 1:abils.")
    }
  }

  # Initialize result list
  R <- vector("list", G)

  # Set first cluster center
  if (!is.null(hyp)) {
    R[[1L]] <- as.integer(hyp)
  } else {
    R[[1L]] <- sample.int(abils)
  }

  # Generate random cluster centers for remaining clusters
  if (G > 1L) {
    for (i in 2L:G) {
      R[[i]] <- sample.int(abils)
    }
  }

  R
}
