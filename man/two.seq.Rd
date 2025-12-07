\docType{data}
\name{two.seq}
\alias{two.seq}
\title{Two-mode Mallows model fit to APA election data}
\format{
A list with components:
\describe{
  \item{R}{List of 2 modal sequences}
  \item{p}{Vector of cluster proportions}
  \item{lambda}{Vector of spread parameters}
  \item{datas}{Data frame with cluster assignments}
  \item{min.like}{Log-likelihood trajectory}
}
}
\source{
American National Election Studies
}
\description{
The \code{elect} data fitted with 2 modal sequences. The two modes
appear to correspond roughly to Democratic and Republican voting patterns.
}
\examples{
data(two.seq)
two.seq$R  # Modal sequences (voting patterns)
two.seq$p  # Proportion in each group
}
\seealso{
\code{\link{elect}}
}
\keyword{datasets}
