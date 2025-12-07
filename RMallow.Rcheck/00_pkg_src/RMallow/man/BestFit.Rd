\name{BestFit}
\alias{BestFit}
\title{Fit Mallows model multiple times and select best}
\usage{
BestFit(datas, N, iter, G, verbose = TRUE)
}
\arguments{
\item{datas}{Matrix of rankings to fit (N x n_items).}

\item{N}{Number of times to run the model with different initializations.}

\item{iter}{Maximum number of EM iterations for each run.}

\item{G}{Number of cluster centers (modal sequences).}

\item{verbose}{Logical; if TRUE, print progress messages.}
}
\value{
The model fit with the highest log-likelihood. See \code{\link{Mallows}}
for details of the return value.
}
\description{
The EM algorithm for Mallows' model can converge to local maxima.
This function runs the algorithm multiple times with different random
initializations and returns the model with the highest likelihood.
}
\seealso{
\code{\link{Mallows}}
}
\author{
Erik Gregory
}
\keyword{Mallows}
\keyword{mixture}
\keyword{model}
