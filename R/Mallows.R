#' Fit a Multi-Modal Mallows' Model to Ranking Data
#'
#' Fits a mixture of Mallows' models to partial or full ranking data using
#' an EM algorithm with Kendall's distance metric.
#'
#' @param datas Matrix of partial or fully-ranked data (N x n_items).
#' @param G Number of modes (clusters), must be >= 1.
#' @param iter Maximum number of EM iterations (default 100).
#' @param hyp Optional hypothesis sequence to initialize one cluster center.
#' @param plot_like Logical; if TRUE, plot likelihood at each iteration.
#' @param verbose Logical; if TRUE, print progress messages.
#' @param tol Convergence tolerance. Algorithm stops when likelihood change
#'   is less than this value.
#' @return List with components:
#'   \describe{
#'     \item{R}{List of modal sequences for each cluster}
#'     \item{p}{Vector of cluster proportions}
#'     \item{lambda}{Vector of spread parameters}
#'     \item{datas}{Data frame with cluster assignments and diagnostics}
#'     \item{min.like}{Vector of log-likelihood at each iteration}
#'   }
#' @author Erik Gregory
#' @references
#'   Murphy, T.B. and Martin, D. (2003). Mixtures of distance-based models
#'   for ranking data. Computational Statistics & Data Analysis, 41, 645-655.
#'
#'   Adkins, L. and Fligner, M. (1998). A Non-iterative procedure for maximum
#'   likelihood estimation of the parameters of Mallows' Model Based on
#'   Partial Rankings. Communications in Statistics - Theory and Methods,
#'   27:9, 2199-2220.
#' @keywords cluster Mallow
#' @export
#' @examples
#' \dontrun{
#' # Load example data
#' data(datas)
#'
#' # Fit a 2-mode model
#' result <- Mallows(datas, G = 2, iter = 50)
#'
#' # View modal sequences
#' result$R
#'
#' # View cluster proportions
#' result$p
#'
#' # View spread parameters
#' result$lambda
#' }
Mallows <- function(datas, G, iter = 100L, hyp = NULL, plot_like = FALSE,
                    verbose = TRUE, tol = 1e-6) {
  # Input validation
  if (!is.matrix(datas) && !is.data.frame(datas)) {
    stop("'datas' must be a matrix or data frame.")
  }
  datas <- as.matrix(datas)

  if (!is.numeric(G) || G < 1L) {
    stop("'G' must be a positive integer.")
  }
  G <- as.integer(G)

  if (!is.numeric(iter) || iter < 1L) {
    stop("'iter' must be a positive integer.")
  }
  iter <- as.integer(iter)

  # Dimensions
  N <- nrow(datas)
  n_items <- ncol(datas)

  if (N < G) {
    stop("Number of observations must be at least G.")
  }

  # Precompute distance distribution table
  dists_table <- DistanceDistribution(n_items)

  # Initialize parameters
  p <- rep(1 / G, G)  # Equal cluster proportions
  R <- Rgen(G, hyp, n_items)  # Random modal sequences
  lambda <- runif(G, min = 0.1, max = 2)  # Random spread parameters

  if (verbose) {
    message("Fitting Mallows mixture model...")
    message(sprintf("  N = %d observations, %d items, G = %d clusters", N, n_items, G))
  }

  # Precompute Kendall information for data (used repeatedly)
  infos <- KendallInfo(datas)

  # Initialize likelihood tracking
  likelihood <- numeric(iter)
  all_dists_data <- NULL

  # Set up plotting if requested
  if (plot_like) {
    if (interactive()) {
      dev.new()
      on.exit(dev.off(), add = TRUE)
    } else {
      message("Note: plot_like=TRUE ignored in non-interactive mode")
      plot_like <- FALSE
    }
  }

  # EM iterations
  for (i in seq_len(iter)) {
    # Compute normalizing constants
    C_lam <- vapply(lambda, function(lam) C_lam(lam, dists_table = dists_table),
                    numeric(1L))

    # E-Step: update membership probabilities
    z <- EStep(R, datas, p, lambda, G, N, C_lam, all_dists_data)

    # M-Step: update parameters
    R <- UpdateR(datas, z, infos)
    p <- UpdateP(z)

    # Compute distances for lambda update and likelihood
    all_dists_data <- AllKendall(datas, do.call(rbind, R), infos)

    lambda <- UpdateLambda(datas, R, z, G, all_dists_data, dists_table = dists_table)

    # Compute likelihood
    likelihood[i] <- Likelihood(z, p, C_lam, lambda, all_dists_data)

    # Plot progress
    if (plot_like) {
      plot(
        seq_len(i), likelihood[seq_len(i)],
        type = "l", col = "red",
        main = "EM Algorithm Convergence",
        xlab = "Iteration",
        ylab = "Log-Likelihood",
        xlim = c(1, iter)
      )
    }

    # Check convergence
    if (i > 1L) {
      delta <- abs(likelihood[i] - likelihood[i - 1L])
      if (delta < tol) {
        if (verbose) {
          message(sprintf("  Converged at iteration %d (delta = %.2e)", i, delta))
        }
        likelihood <- likelihood[seq_len(i)]
        break
      }
    }

    # Progress update every 10 iterations
    if (verbose && i %% 10L == 0L) {
      message(sprintf("  Iteration %d: log-likelihood = %.2f", i, likelihood[i]))
    }
  }

  if (verbose) {
    message("Formatting output...")
  }

  FormatOut(R, p, lambda, z, datas, likelihood)
}

# Backward compatibility alias
#' @rdname Mallows
#' @param plot.like Deprecated. Use \code{plot_like} instead.
Mallows.legacy <- function(datas, G, iter = 100L, hyp = NULL, plot.like = FALSE) {
  Mallows(datas, G, iter = iter, hyp = hyp, plot_like = plot.like)
}
