\name{FormatOut}
\alias{FormatOut}
\title{Format Mallows model output}
\usage{
FormatOut(R, p, lambda, z, datas, likelihood)
}
\arguments{
\item{R}{List of modal sequences (length G).}

\item{p}{Vector of cluster proportions (length G).}

\item{lambda}{Vector of spread parameters (length G).}

\item{z}{Matrix of membership probabilities (N x G).}

\item{datas}{Original data matrix of rankings (N x n_items).}

\item{likelihood}{Vector of log-likelihood values at each iteration.}
}
\value{
List with components:
\describe{
  \item{R}{Modal sequences for each cluster}
  \item{p}{Proportion of data in each cluster}
  \item{lambda}{Spread parameters for each cluster}
  \item{datas}{Data frame with original rankings plus cluster assignments,
    membership probabilities, assigned modal sequence, and distances to
    each cluster center}
  \item{min.like}{Log-likelihood at each iteration}
}
}
\description{
Formats the results of the Mallows model fitting procedure into a
structured list with cluster assignments and diagnostics.
}
\author{
Erik Gregory
}
\keyword{BubbleSort}
\keyword{Kendall}
