\name{SeqDistribution}
\alias{SeqDistribution}
\title{Compute all Kendall distances in N! space by enumeration}
\usage{
SeqDistribution(N)
}
\arguments{
\item{N}{Length of the ranking. Should be less than 9 for practical use.}
}
\value{
Integer vector of Kendall distances from 1:N to each permutation.
}
\description{
Computes Kendall's distance from 1:N to each permutation in N! space by
enumerating all permutations. This is computationally expensive for N >= 8.
For larger N, use \code{\link{DistanceDistribution}} instead.
}
\examples{
# All distances in 4! space
dists <- SeqDistribution(4)
table(dists)
}
\seealso{
\code{\link{DistanceDistribution}} for an efficient alternative
}
\author{
Erik Gregory
}
\keyword{bubblesort}
\keyword{distance}
