\name{Likelihood}
\alias{Likelihood}
\title{Compute log-likelihood of Mallows mixture model}
\usage{
Likelihood(z, p, C_lam, lambda, all_dists_data)
}
\arguments{
\item{z}{Matrix of cluster membership probabilities (N x G).}

\item{p}{Vector of cluster proportions (length G).}

\item{C_lam}{Vector of normalizing constants (length G).}

\item{lambda}{Vector of spread parameters (length G).}

\item{all_dists_data}{Matrix of distances from data to modal sequences (N x G).}
}
\value{
Numeric log-likelihood value.
}
\description{
Calculates the complete-data log-likelihood of the rankings given the
current model parameters.
}
\author{
Erik Gregory
}
\references{
Murphy, T.B. and Martin, D. (2003). Mixtures of distance-based models
for ranking data. Computational Statistics & Data Analysis, 41, 645-655.
}
\keyword{Mallow}
\keyword{likelihood}
