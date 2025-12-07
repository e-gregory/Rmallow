#' Compute pairwise comparison information for Kendall's distance
#'
#' Performs each column-wise comparison on a matrix of sequences. A 0 value
#' denotes an increase between two columns, 1 denotes a decrease, and NA
#' indicates tied values.
#'
#' @param r Matrix of sequences, where each row is a ranking.
#' @param inds Optional matrix of column index pairs for comparisons. If NULL,
#'   all pairs are computed using \code{combn(ncol(r), 2)}.
#' @return Matrix of 0s, 1s, and NAs representing pairwise comparisons. Each
#'   column corresponds to a pair of positions in the original ranking.
#' @author Erik Gregory
#' @references \url{https://en.wikipedia.org/wiki/Kendall_tau_distance}
#' @keywords Kendall Distance
#' @export
#' @examples
#' # Full ranking
#' r <- matrix(c(1, 2, 3, 4, 5, 5, 4, 3, 2, 1), nrow = 2, byrow = TRUE)
#' KendallInfo(r)
#'
#' # Ranking with ties
#' r_ties <- matrix(c(1, 1, 2, 3, 3), nrow = 1)
#' KendallInfo(r_ties)
KendallInfo <- function(r, inds = NULL) {
  # Ensure r is a matrix

if (!is.matrix(r)) {
    r <- as.matrix(r)
    attr(r, "dimnames") <- NULL
  }

  # Input validation
  if (nrow(r) == 0L || ncol(r) == 0L) {
    stop("Input matrix 'r' must have at least one row and one column.")
  }

  # Generate column pairs if not provided
  if (is.null(inds)) {
    if (ncol(r) < 2L) {
      return(matrix(NA_real_, nrow = nrow(r), ncol = 0L))
    }
    inds <- combn(ncol(r), 2L)
  }

  # Compute pairwise differences
  diffs <- r[, inds[1L, ], drop = FALSE] - r[, inds[2L, ], drop = FALSE]

  # Convert to Kendall information: 1 = decrease, 0 = increase, NA = tie
  # Using vectorized operations for efficiency
  infos <- matrix(NA_real_, nrow = nrow(diffs), ncol = ncol(diffs))
  infos[diffs > 0] <- 1
  infos[diffs < 0] <- 0

  infos
}
