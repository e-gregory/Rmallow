#' Compute the normalizing constant for Mallows' model
#'
#' Calculates the normalizing constant C(lambda) for Mallows' model in a
#' sequence space of size N!.
#'
#' @param lambda Spread parameter for Mallows' model (non-negative).
#' @param dists Optional vector of all distances from each sequence to the
#'   modal sequence 1:N. If NULL, \code{dists_table} must be provided.
#' @param dists_table Optional table of distance counts. If NULL, computed
#'   from \code{dists}. This is more efficient for repeated calls.
#' @return The normalizing constant C(lambda).
#' @author Erik Gregory
#' @keywords normalize
#' @export
#' @examples
#' # For a 4-item ranking space
#' dist_tab <- DistanceDistribution(4)
#' C_lam(0.5, dists_table = dist_tab)
#' C_lam(1.0, dists_table = dist_tab)
C_lam <- function(lambda, dists = NULL, dists_table = NULL) {
  # Input validation
  if (lambda < 0) {
    stop("'lambda' must be non-negative.")
  }

  # Compute distance table if not provided
  if (is.null(dists_table)) {
    if (is.null(dists)) {
      stop("Either 'dists' or 'dists_table' must be provided.")
    }
    dists_table <- table(dists)
  }

  # Extract distances as numeric (convert from table names)
  distances <- as.numeric(names(dists_table))
  counts <- as.numeric(dists_table)

  # Compute unnormalized probabilities: count * exp(-lambda * distance)
  unnorm_probs <- counts * exp(-lambda * distances)

  # Return normalizing constant (inverse of sum)
  1 / sum(unnorm_probs)
}

# Backward compatibility alias
#' @rdname C_lam
#' @param dists.table Deprecated. Use \code{dists_table} instead.
C_lam.legacy <- function(lambda, dists = NULL, dists.table = NULL) {
  C_lam(lambda, dists = dists, dists_table = dists.table)
}
