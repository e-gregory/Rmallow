#' Update cluster proportions (M-step)
#'
#' Updates the proportion of data assigned to each cluster based on
#' current membership probabilities.
#'
#' @param z Matrix of membership probabilities (N x G).
#' @return Vector of updated cluster proportions (length G).
#' @author Erik Gregory
#' @references Murphy, T.B. and Martin, D. (2003). Mixtures of distance-based
#'   models for ranking data. Computational Statistics & Data Analysis,
#'   41, 645-655.
#' @keywords proportion
#' @export
#' @examples
#' # Example membership matrix for 10 observations, 3 clusters
#' z <- matrix(c(
#'   0.9, 0.05, 0.05,
#'   0.8, 0.1, 0.1,
#'   0.1, 0.8, 0.1,
#'   0.1, 0.1, 0.8,
#'   0.33, 0.33, 0.34,
#'   0.9, 0.05, 0.05,
#'   0.05, 0.9, 0.05,
#'   0.05, 0.05, 0.9,
#'   0.7, 0.2, 0.1,
#'   0.2, 0.3, 0.5
#' ), nrow = 10, byrow = TRUE)
#' UpdateP(z)
UpdateP <- function(z) {
  # Validate input
  if (!is.matrix(z)) {
    stop("'z' must be a matrix.")
  }

  # MLE of mixing proportions = column means
  colMeans(z)
}
