\name{UpdateLambda}
\alias{UpdateLambda}
\title{Update spread parameters (M-step)}
\usage{
UpdateLambda(r, R, z, G, dists_to_Rg, dists_table, top_bound = 1000)
}
\arguments{
\item{r}{Matrix of rankings (N x n_items).}

\item{R}{List of current modal sequences.}

\item{z}{Matrix of membership probabilities (N x G).}

\item{G}{Number of clusters.}

\item{dists_to_Rg}{Matrix of distances from data to modal sequences (N x G).}

\item{dists_table}{Table of distance distribution in N! space.}

\item{top_bound}{Maximum allowed value for lambda (default 1000).}
}
\value{
Vector of updated lambda parameters (length G).
}
\description{
Updates the lambda (spread) parameters for each cluster to maximize
the likelihood given current membership probabilities.
}
\author{
Erik Gregory
}
\references{
Murphy, T.B. and Martin, D. (2003). Mixtures of distance-based models
for ranking data. Computational Statistics & Data Analysis, 41, 645-655.
}
\keyword{lambda}
\keyword{maximization}
