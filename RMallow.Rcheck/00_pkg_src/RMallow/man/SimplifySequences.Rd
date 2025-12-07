\name{SimplifySequences}
\alias{SimplifySequences}
\title{Simplify tied rankings to consecutive integers}
\usage{
SimplifySequences(rankings, use_cpp = TRUE)
}
\arguments{
\item{rankings}{Matrix of rankings to simplify, where each row is a ranking.}

\item{use_cpp}{Logical; if TRUE (default), use the fast C++ implementation.}
}
\value{
Matrix of simplified rankings with the same dimensions.
}
\description{
Transforms rankings so that tie groups use consecutive integers.
For example, (1, 1, 2, 4, 4, 5) becomes (1, 1, 2, 3, 3, 4).
}
\examples{
# Single ranking with gaps
rankings <- matrix(c(1, 1, 2, 4, 4, 5), nrow = 1)
SimplifySequences(rankings)
# Returns: 1 1 2 3 3 4

# Multiple rankings
rankings <- rbind(
  c(1, 3, 3, 5),
  c(2, 2, 4, 6)
)
SimplifySequences(rankings)
}
\author{
Erik Gregory
}
\keyword{sequence}
\keyword{simplify}
