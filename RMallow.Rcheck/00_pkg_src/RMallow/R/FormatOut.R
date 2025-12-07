#' Format Mallows model output
#'
#' Formats the results of the Mallows model fitting procedure into a
#' structured list with cluster assignments and diagnostics.
#'
#' @param R List of modal sequences (length G).
#' @param p Vector of cluster proportions (length G).
#' @param lambda Vector of spread parameters (length G).
#' @param z Matrix of membership probabilities (N x G).
#' @param datas Original data matrix of rankings (N x n_items).
#' @param likelihood Vector of log-likelihood values at each iteration.
#' @return List with components:
#'   \describe{
#'     \item{R}{Modal sequences for each cluster}
#'     \item{p}{Proportion of data in each cluster}
#'     \item{lambda}{Spread parameters for each cluster}
#'     \item{datas}{Data frame with original rankings plus:
#'       \describe{
#'         \item{clust}{Assigned cluster (hard assignment)}
#'         \item{pvals.*}{Membership probabilities for each cluster}
#'         \item{seq}{Character representation of assigned modal sequence}
#'         \item{dists.*}{Distance to each cluster center}
#'       }
#'     }
#'     \item{min.like}{Log-likelihood at each iteration}
#'   }
#' @author Erik Gregory
#' @keywords BubbleSort Kendall
#' @export
FormatOut <- function(R, p, lambda, z, datas, likelihood) {
  # Convert datas to matrix if needed
  if (!is.matrix(datas)) {
    datas <- as.matrix(datas)
  }

  G <- length(R)
  N <- ncol(datas)

  # Hard cluster assignment (maximum probability)
  clust <- apply(z, 1L, which.max)

  # Create string representation of modal sequences
  seqs <- vapply(R, function(seq) paste(seq, collapse = " "), character(1L))
  
  # Handle case where multiple clusters converge to same sequence
  unique_seqs <- unique(seqs)

  # Compute distances to all cluster centers
  R_matrix <- do.call(rbind, R)
  dists <- AllKendall(datas, R_matrix)

  # Build output data frame
  out_df <- as.data.frame(datas)
  out_df$clust <- clust

  # Add membership probabilities
  for (g in seq_len(G)) {
    out_df[[paste0("pvals.", g)]] <- z[, g]
  }

  # Add assigned sequence (use unique levels to avoid duplicate factor level error)
  out_df$seq <- factor(seqs[clust], levels = unique_seqs)

  # Add distances
  for (g in seq_len(G)) {
    out_df[[paste0("dists.", g)]] <- dists[, g]
  }

  # Return structured list
  list(
    R = R,
    p = p,
    lambda = lambda,
    datas = out_df,
    min.like = likelihood
  )
}
