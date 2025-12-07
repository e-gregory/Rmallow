#' Compute the Kendall distance distribution in N! space
#'
#' Efficiently computes the count of permutations at each Kendall distance
#' from the canonical ordering (1, 2, ..., N). Uses a recurrence relation
#' rather than enumeration, making it efficient even for large N.
#'
#' @param N Integer value, the length of rankings (must be >= 1).
#' @return Named numeric vector where names are distances (0 to N*(N-1)/2) and
#'   values are counts of permutations at that distance.
#' @author Erik Gregory
#' @keywords bubblesort Kendall
#' @seealso \code{\link{SeqDistribution}} for the brute-force approach
#' @export
#' @examples
#' # Distribution for 5-item rankings
#' dist_tab <- DistanceDistribution(5)
#' dist_tab
#'
#' # Verify total equals 5!
#' sum(dist_tab)
#'
#' # Works efficiently even for large N
#' dist_tab_20 <- DistanceDistribution(20)
#' length(dist_tab_20)  # Maximum distance is 20*19/2 = 190
DistanceDistribution <- function(N = 3L) {
  # Input validation
  if (!is.numeric(N) || length(N) != 1L || N < 1L || N != as.integer(N)) {
    stop("'N' must be a positive integer.")
  }
  N <- as.integer(N)

  # Base cases
  if (N == 1L) {
    out <- 1L
    names(out) <- "0"
    return(out)
  }

  if (N == 2L) {
    out <- c(1L, 1L)
    names(out) <- c("0", "1")
    return(out)
  }

  # For small N, direct enumeration is fine
  if (N <= 4L) {
    out <- table(SeqDistribution(N))
    return(out)
  }

  # For larger N, use recurrence relation
  # Start with N=4 and build up
  out <- table(SeqDistribution(4L))
  out <- as.numeric(out)

  for (i in 4L:(N - 1L)) {
    out <- NextTable(out, i)
  }

  # Add names (distances from 0 to max)
  names(out) <- as.character(seq(0L, length(out) - 1L))

  out
}
