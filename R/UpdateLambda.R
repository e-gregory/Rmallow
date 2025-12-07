#' Update spread parameters (M-step)
#'
#' Updates the lambda (spread) parameters for each cluster to maximize
#' the likelihood given current membership probabilities.
#'
#' @param r Matrix of rankings (N x n_items).
#' @param R List of current modal sequences.
#' @param z Matrix of membership probabilities (N x G).
#' @param G Number of clusters.
#' @param dists_to_Rg Matrix of distances from data to modal sequences (N x G).
#' @param dists_table Table of distance distribution in N! space.
#' @param top_bound Maximum allowed value for lambda (default 1000).
#' @return Vector of updated lambda parameters (length G).
#' @author Erik Gregory
#' @references Murphy, T.B. and Martin, D. (2003). Mixtures of distance-based
#'   models for ranking data. Computational Statistics & Data Analysis,
#'   41, 645-655.
#' @keywords lambda maximization
#' @export
UpdateLambda <- function(r, R, z, G, dists_to_Rg,
                         dists_table, top_bound = 1000) {
  # Initialize output
  lambda <- numeric(G)

  # Precompute values for root-finding bounds
  distances <- as.numeric(names(dists_table))
  counts <- as.numeric(dists_table)

  # Function value at lambda = 0
  # E[d | lambda=0] = mean distance over uniform distribution
  f_at_zero <- sum(distances * counts) / sum(counts)

  # Function value at lambda = top_bound
  # E[d | lambda=top_bound] approaches min distance (0) as lambda -> inf
  weights_upper <- counts * exp(-top_bound * distances)
  f_at_upper <- sum(distances * weights_upper) / sum(weights_upper)

  # Update lambda for each cluster
  for (i in seq_len(G)) {
    # Compute weighted mean distance (RHS of estimating equation)
    total_weight <- sum(z[, i])
    if (total_weight < .Machine$double.eps) {
      # Empty cluster - use default value
      lambda[i] <- 0
      next
    }

    rhs <- sum(z[, i] * dists_to_Rg[, i]) / total_weight

    # Compute function values at bounds for this RHS
    f_lower <- f_at_zero - rhs
    f_upper <- f_at_upper - rhs

    # Check if root exists in interval
    if (sign(f_lower) != sign(f_upper)) {
      # Find root using uniroot
      result <- uniroot(
        Lambda,
        interval = c(0, top_bound),
        rhs = rhs,
        dists_table = dists_table,
        f.lower = f_lower,
        f.upper = f_upper,
        tol = 1e-5
      )
      lambda[i] <- result$root
    } else if (rhs < .Machine$double.eps) {
      # All observations exactly at modal sequence
      lambda[i] <- top_bound
    } else {
      # Solution exceeds bound
      message(sprintf("Lambda for cluster %d exceeds bound, using top_bound = %g", i, top_bound))
      lambda[i] <- top_bound
    }
  }

  lambda
}

# Backward compatibility alias
#' @rdname UpdateLambda
#' @param dists.to.Rg Deprecated. Use \code{dists_to_Rg} instead.
#' @param dists.table Deprecated. Use \code{dists_table} instead.
#' @param top.bound Deprecated. Use \code{top_bound} instead.
UpdateLambda.legacy <- function(r, R, z, G, dists.to.Rg,
                                dists.table, top.bound = 1000) {
  UpdateLambda(r, R, z, G,
               dists_to_Rg = dists.to.Rg,
               dists_table = dists.table,
               top_bound = top.bound)
}
