\docType{data}
\name{two.mode}
\alias{two.mode}
\title{Two-mode Mallows model fit to synthetic data}
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
\description{
The \code{datas} data set fitted with only 2 modal sequences.
Since the true data has 3 modes, this shows what happens with
model misspecification.
}
\examples{
data(two.mode)
two.mode$R  # Only 2 modal sequences
}
\seealso{
\code{\link{three.mode}}, \code{\link{datas}}
}
\keyword{datasets}
