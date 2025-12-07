\name{C_lam}
\alias{C_lam}
\title{Compute the normalizing constant for Mallows' model}
\usage{
C_lam(lambda, dists = NULL, dists_table = NULL)
}
\arguments{
\item{lambda}{Spread parameter for Mallows' model (non-negative).}

\item{dists}{Optional vector of all distances from each sequence to the
modal sequence 1:N. If NULL, \code{dists_table} must be provided.}

\item{dists_table}{Optional table of distance counts. If NULL, computed
from \code{dists}. This is more efficient for repeated calls.}
}
\value{
The normalizing constant C(lambda).
}
\description{
Calculates the normalizing constant C(lambda) for Mallows' model in a
sequence space of size N!.
}
\examples{
# For a 4-item ranking space
dist_tab <- DistanceDistribution(4)
C_lam(0.5, dists_table = dist_tab)
C_lam(1.0, dists_table = dist_tab)
}
\author{
Erik Gregory
}
\keyword{normalize}
