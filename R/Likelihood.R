#' Compute log-likelihood of Mallows mixture model
#'
#' Calculates the complete-data log-likelihood of the rankings given the
#' current model parameters.
#'
#' @param z Matrix of cluster membership probabilities (N x G).
#' @param p Vector of cluster proportions (length G).
#' @param C_lam Vector of normalizing constants (length G).
#' @param lambda Vector of spread parameters (length G).
#' @param all_dists_data Matrix of distances from data to modal sequences (N x G).
#' @return Numeric log-likelihood value.
#' @author Erik Gregory
#' @references Murphy, T.B. and Martin, D. (2003). Mixtures of distance-based
#'   models for ranking data. Computational Statistics & Data Analysis,
#'   41, 645-655.
#' @keywords likelihood Mallow
#' @export
Likelihood <- function(z, p, C_lam, lambda, all_dists_data) {
  # Input validation
  G <- length(p)
  if (length(C_lam) != G || length(lambda) != G) {
    stop("'p', 'C_lam', and 'lambda' must have the same length.")
  }
  if (ncol(z) != G || ncol(all_dists_data) != G) {
    stop("'z' and 'all_dists_data' must have G columns.")
  }

  # Compute log(C * p) for each cluster (additive constants)
  log_const <- log(C_lam * p)

  # Compute -lambda * distance for each observation and cluster
  # t(-lambda * t(all_dists_data)) is equivalent to sweep but less readable
  neg_lambda_dists <- sweep(all_dists_data, 2L, lambda, `*`)
  neg_lambda_dists <- -neg_lambda_dists

  # Add log constants
  log_probs <- sweep(neg_lambda_dists, 2L, log_const, `+`)

  # Weight by membership probabilities and sum
  sum(log_probs * z)
}

# Backward compatibility alias
#' @rdname Likelihood
#' @param C.lam Deprecated. Use \code{C_lam} instead.
#' @param all.dists.data Deprecated. Use \code{all_dists_data} instead.
Likelihood.legacy <- function(z, p, C.lam, lambda, all.dists.data) {
  Likelihood(z, p, C_lam = C.lam, lambda, all_dists_data = all.dists.data)
}
