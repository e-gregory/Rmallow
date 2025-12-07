#' Objective function for lambda estimation
#'
#' Computes the objective function whose root gives the MLE of lambda
#' in Mallows' model. Used internally by \code{\link{UpdateLambda}}.
#'
#' @param lambda Lambda value to evaluate the function at.
#' @param rhs Right-hand side of the estimating equation (weighted mean distance).
#' @param dists Not used (kept for backward compatibility).
#' @param dists_table Table of distance distribution in N! space.
#' @return Value of the objective function. The goal is to find lambda where
#'   this equals zero.
#' @author Erik Gregory
#' @references Murphy, T.B. and Martin, D. (2003). Mixtures of distance-based
#'   models for ranking data. Computational Statistics & Data Analysis,
#'   41, 645-655.
#' @keywords Mallow lambda
#' @export
Lambda <- function(lambda, rhs, dists = NULL, dists_table = NULL) {
  # Input validation
  if (is.null(dists_table)) {
    stop("'dists_table' must be provided.")
  }

  if (lambda < 0) {
    stop("'lambda' must be non-negative.")
  }

  # Extract distances and counts
  distances <- as.numeric(names(dists_table))
  counts <- as.numeric(dists_table)

  # Compute expected distance under Mallows model with this lambda
  # E[d] = sum(d * p(d)) where p(d) proportional to count(d) * exp(-lambda * d)
  weights <- counts * exp(-lambda * distances)
  expected_dist <- sum(distances * weights) / sum(weights)

  # Return difference from observed weighted mean distance
  expected_dist - rhs
}

# Backward compatibility alias
#' @rdname Lambda
#' @param dists.table Deprecated. Use \code{dists_table} instead.
Lambda.legacy <- function(lambda, rhs, dists = NULL, dists.table = NULL) {
  Lambda(lambda, rhs, dists = dists, dists_table = dists.table)
}
