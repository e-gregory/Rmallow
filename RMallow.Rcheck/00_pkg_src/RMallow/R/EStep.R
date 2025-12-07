#' E-Step of the EM algorithm
#'
#' Computes the probability that each ranking belongs to each cluster,
#' given the current model parameters.
#'
#' @param R List of current cluster modal sequences.
#' @param r Matrix of rankings (N x n_items).
#' @param p Vector of cluster proportions (length G).
#' @param lambda Vector of spread parameters (length G).
#' @param G Number of clusters.
#' @param N Number of observations.
#' @param C Vector of normalizing constants (length G).
#' @param all_dists Optional precomputed distance matrix (N x G). If NULL,
#'   computed from \code{r} and \code{R}.
#' @param use_cpp Logical; if TRUE (default), use the fast C++ implementation.
#' @return Matrix (N x G) where entry \code{[i, j]} is the probability that
#'   observation i belongs to cluster j.
#' @author Erik Gregory
#' @references Murphy, T.B. and Martin, D. (2003). Mixtures of distance-based
#'   models for ranking data. Computational Statistics & Data Analysis,
#'   41, 645-655.
#' @keywords expectation maximization
#' @export
EStep <- function(R, r, p, lambda, G, N, C, all_dists = NULL, use_cpp = TRUE) {
  # Compute distances if not provided
  if (is.null(all_dists)) {
    all_dists <- AllKendall(r, do.call(rbind, R), use_cpp = use_cpp)
  }

  # Use C++ implementation
  if (use_cpp) {
    z <- estep_cpp(all_dists, lambda, C, p)
    # Handle any NaN from 0/0
    if (any(is.nan(z))) {
      z[is.nan(z)] <- 1 / G
    }
    return(z)
  }

  # Fallback to R implementation
  EStep_R(R, r, p, lambda, G, N, C, all_dists)
}

#' Pure R implementation of EStep (for benchmarking)
#' @keywords internal
EStep_R <- function(R, r, p, lambda, G, N, C, all_dists = NULL) {
  if (is.null(all_dists)) {
    all_dists <- AllKendall_R(r, do.call(rbind, R))
  }

  # Compute exp(-lambda * distance) for each obs and cluster
  log_probs <- sweep(-all_dists, 2L, lambda, `*`)

  # Add log(C * p) for numerical stability
  log_probs <- sweep(log_probs, 2L, log(C * p), `+`)

  # Convert to probabilities using softmax
  log_probs <- log_probs - apply(log_probs, 1L, max)
  z <- exp(log_probs)

  # Normalize rows to sum to 1
  z <- z / rowSums(z)

  # Handle any NaN from 0/0
  if (any(is.nan(z))) {
    z[is.nan(z)] <- 1 / G
  }

  z
}

# Backward compatibility alias
#' @rdname EStep
#' @param all.dists Deprecated. Use \code{all_dists} instead.
EStep.legacy <- function(R, r, p, lambda, G, N, C, all.dists = NULL) {
  EStep(R, r, p, lambda, G, N, C, all_dists = all.dists)
}
