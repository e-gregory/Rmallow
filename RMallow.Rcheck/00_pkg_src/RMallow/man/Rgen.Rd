\name{Rgen}
\alias{Rgen}
\title{Initialize cluster modal sequences}
\usage{
Rgen(G, hyp = NULL, abils)
}
\arguments{
\item{G}{Number of cluster centers (modes).}

\item{hyp}{Optional hypothesis sequence to use as the first cluster center.
Must be a permutation of 1:abils.}

\item{abils}{Number of items being ranked.}
}
\value{
List of G integer vectors, each a permutation of 1:abils.
}
\description{
Generates initial modal sequences for the Mallows mixture model clustering.
If a hypothesis sequence is provided, it becomes the first cluster center.
Remaining cluster centers are initialized randomly.
}
\examples{
# Three random cluster centers for 5-item rankings
set.seed(42)
Rgen(3, abils = 5)

# With a specific hypothesis sequence
Rgen(3, hyp = c(5, 4, 3, 2, 1), abils = 5)
}
\author{
Erik Gregory
}
\keyword{cluster}
\keyword{initialization}
