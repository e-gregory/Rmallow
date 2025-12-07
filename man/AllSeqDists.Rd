\name{AllSeqDists}
\alias{AllSeqDists}
\title{Compute distances from sequences to the canonical ordering}
\usage{
AllSeqDists(seqs)
}
\arguments{
\item{seqs}{Matrix or data frame of sequences, where each row is a ranking.}
}
\value{
Integer vector of distances from each sequence to \code{1:N}.
}
\description{
Calculates the Kendall distance from each sequence to the canonical
ordering (1, 2, ..., N).
}
\examples{
seqs <- rbind(1:5, 5:1, c(1, 3, 2, 4, 5))
AllSeqDists(seqs)
# Returns: 0, 10, 1
}
\author{
Erik Gregory
}
\keyword{Distance}
\keyword{Kendall}
