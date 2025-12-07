\name{KendallInfo}
\alias{KendallInfo}
\title{Compute pairwise comparison information for Kendall's distance}
\usage{
KendallInfo(r, inds = NULL, use_cpp = TRUE)
}
\arguments{
\item{r}{Matrix of sequences, where each row is a ranking.}

\item{inds}{Optional matrix of column index pairs for comparisons. If NULL,
all pairs are computed using \code{combn(ncol(r), 2)}.}

\item{use_cpp}{Logical; if TRUE (default), use the fast C++ implementation.}
}
\value{
Matrix of 0s, 1s, and NAs representing pairwise comparisons. Each
column corresponds to a pair of positions in the original ranking.
}
\description{
Performs each column-wise comparison on a matrix of sequences. A 0 value
denotes an increase between two columns, 1 denotes a decrease, and NA
indicates tied values.
}
\examples{
# Full ranking
r <- matrix(c(1, 2, 3, 4, 5, 5, 4, 3, 2, 1), nrow = 2, byrow = TRUE)
KendallInfo(r)

# Ranking with ties
r_ties <- matrix(c(1, 1, 2, 3, 3), nrow = 1)
KendallInfo(r_ties)
}
\author{
Erik Gregory
}
\references{
\url{https://en.wikipedia.org/wiki/Kendall_tau_distance}
}
\keyword{Distance}
\keyword{Kendall}
