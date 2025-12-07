\docType{package}
\name{RMallow-package}
\alias{RMallow-package}
\alias{RMallow}
\title{RMallow: Fit Multi-Modal Mallows' Models to Ranking Data}
\description{
Fits the Mallows' model to ranking data using an EM algorithm.
Data can be partially or fully-ranked, with or without ties.
}
\details{
The main functions are:
\itemize{
  \item \code{\link{Mallows}}: Fit a Mallows mixture model
  \item \code{\link{BestFit}}: Fit multiple models and select the best
}

Distance computation functions:
\itemize{
  \item \code{\link{AllKendall}}: Compute Kendall distances between rankings
  \item \code{\link{KendallInfo}}: Compute pairwise comparison information
  \item \code{\link{DistanceDistribution}}: Distance distribution in N! space
}

Utility functions:
\itemize{
  \item \code{\link{SimplifySequences}}: Normalize tied rankings
  \item \code{\link{ConstructSeqs}}: Reconstruct sequences from preferences
}
}
\author{
Erik Gregory \email{egregory2011@gmail.com}
}
\references{
Murphy, T.B. and Martin, D. (2003). Mixtures of distance-based models
for ranking data. Computational Statistics & Data Analysis, 41, 645-655.
\doi{10.1016/S0167-9473(02)00165-2}

Beckett, L.A. and Evans, D.A. (1991). Estimating a Population Distribution
of Sequences of k Items from Cross-Sectional Data. Journal of the Royal
Statistical Society, Series C, 40(1), 31-42.

Adkins, L. and Fligner, M. (1998). A Non-iterative procedure for maximum
likelihood estimation of the parameters of Mallows' Model Based on
Partial Rankings. Communications in Statistics - Theory and Methods,
27(9), 2199-2220. \doi{10.1080/03610929808832223}
}
\keyword{package}
\keyword{ranking}
\keyword{Mallows}
\keyword{clustering}
