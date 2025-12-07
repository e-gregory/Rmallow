\name{EStep}
\alias{EStep}
\title{E-Step of the EM algorithm}
\usage{
EStep(R, r, p, lambda, G, N, C, all_dists = NULL, use_cpp = TRUE)
}
\arguments{
\item{R}{List of current cluster modal sequences.}

\item{r}{Matrix of rankings (N x n_items).}

\item{p}{Vector of cluster proportions (length G).}

\item{lambda}{Vector of spread parameters (length G).}

\item{G}{Number of clusters.}

\item{N}{Number of observations.}

\item{C}{Vector of normalizing constants (length G).}

\item{all_dists}{Optional precomputed distance matrix (N x G). If NULL,
computed from \code{r} and \code{R}.}

\item{use_cpp}{Logical; if TRUE (default), use the fast C++ implementation.}
}
\value{
Matrix (N x G) where entry \code{[i, j]} is the probability that
observation i belongs to cluster j.
}
\description{
Computes the probability that each ranking belongs to each cluster,
given the current model parameters.
}
\author{
Erik Gregory
}
\references{
Murphy, T.B. and Martin, D. (2003). Mixtures of distance-based models
for ranking data. Computational Statistics & Data Analysis, 41, 645-655.
}
\keyword{expectation}
\keyword{maximization}
