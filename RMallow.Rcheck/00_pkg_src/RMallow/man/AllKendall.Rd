\name{AllKendall}
\alias{AllKendall}
\title{Compute Kendall's distances between two sets of rankings}
\usage{
AllKendall(r, seqs, data_info = NULL, use_cpp = TRUE)
}
\arguments{
\item{r}{Matrix of rankings (each row is a ranking).}

\item{seqs}{Matrix of reference sequences (each row is a sequence).}

\item{data_info}{Optional precomputed Kendall information matrix for \code{r}.
If NULL, it will be computed internally. Providing this can improve
efficiency when calling the function repeatedly with the same data.}

\item{use_cpp}{Logical; if TRUE (default), use the fast C++ implementation.}
}
\value{
Matrix where entry \code{[i, j]} represents the Kendall distance
from ranking \code{i} in \code{r} to sequence \code{j} in \code{seqs}.
}
\description{
Calculates all Kendall's distances between each ranking in \code{r} and
each sequence in \code{seqs}.
}
\examples{
data1 <- rbind(1:5, 5:1, c(3, 2, 1, 4, 5))
data2 <- rbind(1:5, 5:1)
AllKendall(data1, data2)
}
\author{
Erik Gregory
}
\keyword{Kendall}
\keyword{distance}
