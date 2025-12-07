\name{DistanceDistribution}
\alias{DistanceDistribution}
\title{Compute the Kendall distance distribution in N! space}
\usage{
DistanceDistribution(N = 3, use_cpp = TRUE)
}
\arguments{
\item{N}{Integer value, the length of rankings (must be >= 1).}

\item{use_cpp}{Logical; if TRUE (default), use the fast C++ implementation.}
}
\value{
Named numeric vector where names are distances (0 to N*(N-1)/2) and
values are counts of permutations at that distance.
}
\description{
Efficiently computes the count of permutations at each Kendall distance
from the canonical ordering (1, 2, ..., N). Uses a recurrence relation
rather than enumeration, making it efficient even for large N.
}
\examples{
# Distribution for 5-item rankings
dist_tab <- DistanceDistribution(5)
dist_tab

# Verify total equals 5!
sum(dist_tab)

# Works efficiently even for large N
dist_tab_20 <- DistanceDistribution(20)
length(dist_tab_20)  # Maximum distance is 20*19/2 = 190
}
\seealso{
\code{\link{SeqDistribution}} for the brute-force approach
}
\author{
Erik Gregory
}
\keyword{Kendall}
\keyword{bubblesort}
