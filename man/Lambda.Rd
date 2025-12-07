\name{Lambda}
\alias{Lambda}
\title{Objective function for lambda estimation}
\usage{
Lambda(lambda, rhs, dists = NULL, dists_table = NULL)
}
\arguments{
\item{lambda}{Lambda value to evaluate the function at.}

\item{rhs}{Right-hand side of the estimating equation (weighted mean distance).}

\item{dists}{Not used (kept for backward compatibility).}

\item{dists_table}{Table of distance distribution in N! space.}
}
\value{
Value of the objective function. The goal is to find lambda where
this equals zero.
}
\description{
Computes the objective function whose root gives the MLE of lambda
in Mallows' model. Used internally by \code{\link{UpdateLambda}}.
}
\author{
Erik Gregory
}
\references{
Murphy, T.B. and Martin, D. (2003). Mixtures of distance-based models
for ranking data. Computational Statistics & Data Analysis, 41, 645-655.
}
\keyword{Mallow}
\keyword{lambda}
