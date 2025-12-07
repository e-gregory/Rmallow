#' Fit Mallows model multiple times and select best
#'
#' The EM algorithm for Mallows' model can converge to local maxima.
#' This function runs the algorithm multiple times with different random
#' initializations and returns the model with the highest likelihood.
#'
#' @param datas Matrix of rankings to fit (N x n_items).
#' @param N Number of times to run the model with different initializations.
#' @param iter Maximum number of EM iterations for each run.
#' @param G Number of cluster centers (modal sequences).
#' @param verbose Logical; if TRUE, print progress messages.
#' @return The model fit with the highest log-likelihood. See \code{\link{Mallows}}
#'   for details of the return value.
#' @author Erik Gregory
#' @seealso \code{\link{Mallows}}
#' @keywords Mallows mixture model
#' @export
#' @examples
#' \dontrun{
#' data(datas)
#' # Run 5 times, keeping best result
#' best_model <- BestFit(datas, N = 5, iter = 100, G = 3)
#' }
BestFit <- function(datas, N, iter, G, verbose = TRUE) {
  # Input validation
  if (!is.numeric(N) || N < 1L) {
    stop("'N' must be a positive integer.")
  }
  N <- as.integer(N)

  if (verbose) {
    message(sprintf("Running %d model fits with G=%d clusters...", N, G))
  }

  # Store all models
  models <- vector("list", N)

  for (i in seq_len(N)) {
    if (verbose) {
      message(sprintf("  Fit %d of %d...", i, N))
    }
    models[[i]] <- Mallows(datas, iter = iter, G = G, plot_like = FALSE)
  }

  # Extract final log-likelihoods
  # Use the last non-zero likelihood value from each model
  final_likes <- vapply(models, function(mod) {
    likes <- mod$min.like
    non_zero <- likes[likes != 0]
    if (length(non_zero) > 0L) {
      non_zero[length(non_zero)]
    } else {
      -Inf
    }
  }, numeric(1L))

  # Select best model
  best_idx <- which.max(final_likes)

  if (verbose) {
    message(sprintf("Best model: fit %d (log-likelihood: %.2f)", best_idx, final_likes[best_idx]))
  }

  models[[best_idx]]
}
