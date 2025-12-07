\docType{data}
\name{three.mode}
\alias{three.mode}
\title{Three-mode Mallows model fit to synthetic data}
\format{
A list with components:
\describe{
  \item{R}{List of 3 modal sequences}
  \item{p}{Vector of cluster proportions}
  \item{lambda}{Vector of spread parameters}
  \item{datas}{Data frame with cluster assignments}
  \item{min.like}{Log-likelihood trajectory}
}
}
\description{
The \code{datas} data set fitted with 3 modal sequences.
Compare to \code{\link{two.mode}} to see the effect of model selection.
}
\examples{
data(three.mode)
three.mode$R  # Modal sequences
three.mode$p  # Cluster proportions
}
\seealso{
\code{\link{two.mode}}, \code{\link{datas}}
}
\keyword{datasets}
