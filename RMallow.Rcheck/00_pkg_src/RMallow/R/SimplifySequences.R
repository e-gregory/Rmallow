#' Simplify tied rankings to consecutive integers
#'
#' Transforms rankings so that tie groups use consecutive integers.
#' For example, (1, 1, 2, 4, 4, 5) becomes (1, 1, 2, 3, 3, 4).
#'
#' @param rankings Matrix of rankings to simplify, where each row is a ranking.
#' @param use_cpp Logical; if TRUE (default), use the fast C++ implementation.
#' @return Matrix of simplified rankings with the same dimensions.
#' @author Erik Gregory
#' @keywords simplify sequence
#' @export
#' @examples
#' # Single ranking with gaps
#' rankings <- matrix(c(1, 1, 2, 4, 4, 5), nrow = 1)
#' SimplifySequences(rankings)
#' # Returns: 1 1 2 3 3 4
#'
#' # Multiple rankings
#' rankings <- rbind(
#'   c(1, 3, 3, 5),
#'   c(2, 2, 4, 6)
#' )
#' SimplifySequences(rankings)
SimplifySequences <- function(rankings, use_cpp = TRUE) {
  # Handle vector input
  if (!is.matrix(rankings)) {
    rankings <- matrix(rankings, nrow = 1L)
    was_vector <- TRUE
  } else {
    was_vector <- FALSE
  }

  if (nrow(rankings) == 0L) {
    return(rankings)
  }

  # Use C++ implementation
  if (use_cpp) {
    result <- simplify_sequences_cpp(rankings)

    # Ensure result is a matrix with correct dimensions
    if (was_vector || nrow(rankings) == 1L) {
      result <- matrix(result, nrow = 1L)
    }

    dimnames(result) <- dimnames(rankings)
    return(result)
  }

  # Fallback to R implementation
  SimplifySequences_R(rankings)
}

#' Pure R implementation of SimplifySequences (for benchmarking)
#' @keywords internal
SimplifySequences_R <- function(rankings) {
  if (!is.matrix(rankings)) {
    rankings <- matrix(rankings, nrow = 1L)
    was_vector <- TRUE
  } else {
    was_vector <- FALSE
  }

  if (nrow(rankings) == 0L) {
    return(rankings)
  }

  result <- t(apply(rankings, 1L, function(row) {
    unique_vals <- sort(unique(row))
    mapping <- match(row, unique_vals)
    mapping
  }))

  if (was_vector || nrow(rankings) == 1L) {
    result <- matrix(result, nrow = 1L)
  }

  dimnames(result) <- dimnames(rankings)

  result
}

# Backward compatibility alias
#' @rdname SimplifySequences
#' @param loss.time Deprecated. Use \code{rankings} instead.
SimplifySequences.legacy <- function(loss.time) {
  SimplifySequences(rankings = loss.time)
}
