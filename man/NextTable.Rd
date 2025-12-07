\name{NextTable}
\alias{NextTable}
\title{Compute distance distribution for (N+1)! from N! space}
\usage{
NextTable(last_table, n_last)
}
\arguments{
\item{last_table}{Table or numeric vector of distance counts in N! space.}

\item{n_last}{The value of N (current space is N!).}
}
\value{
Numeric vector of distance counts in (N+1)! space.
}
\description{
Given the Kendall distance distribution in N! space, computes the
distribution in (N+1)! space. This uses the recurrence relation for
Mahonian numbers (number of permutations with k inversions).
}
\examples{
# Get distribution for 5! from 4!
tab4 <- DistanceDistribution(4)
tab5 <- NextTable(tab4, 4)
# Verify: should have 5! = 120 total permutations
sum(tab5)
}
\author{
Erik Gregory
}
\keyword{Kendall}
\keyword{bubblesort}
