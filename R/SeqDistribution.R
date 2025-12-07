#' Compute all Kendall distances in N! space by enumeration
#'
#' Computes Kendall's distance from 1:N to each permutation in N! space by
#' enumerating all permutations. This is computationally expensive for N >= 8.
#' For larger N, use \code{\link{DistanceDistribution}} instead.
#'
#' @param N Length of the ranking. Should be less than 9 for practical use.
#' @return Integer vector of Kendall distances from 1:N to each permutation.
#' @author Erik Gregory
#' @keywords distance bubblesort
#' @seealso \code{\link{DistanceDistribution}} for an efficient alternative
#' @export
#' @examples
#' # All distances in 4! space
#' dists <- SeqDistribution(4)
#' table(dists)
#'
#' # Warning: this is slow for N >= 8
#' \dontrun{
#' SeqDistribution(10)  # Will take a very long time!
#' }
SeqDistribution <- function(N) {
  # Input validation
  if (!is.numeric(N) || length(N) != 1L || N < 1L || N != as.integer(N)) {
    stop("'N' must be a positive integer.")
  }
  N <- as.integer(N)

  if (N >= 8L) {
    message("Note: For N >= 8, consider using DistanceDistribution() for efficiency.")
  }

  # Generate all permutations
  seqs <- do.call(rbind, combinat::permn(N))

  # Compute distances to canonical ordering
  AllSeqDists(seqs)
}
