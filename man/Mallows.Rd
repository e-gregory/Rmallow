\name{Mallows}
\alias{Mallows}
\title{Fit a Multi-Modal Mallows' Model to Ranking Data}
\usage{
Mallows(datas, G, iter = 100L, hyp = NULL, plot_like = FALSE,
        verbose = TRUE, tol = 1e-06)
}
\arguments{
\item{datas}{Matrix of partial or fully-ranked data (N x n_items).}

\item{G}{Number of modes (clusters), must be >= 1.}

\item{iter}{Maximum number of EM iterations (default 100).}

\item{hyp}{Optional hypothesis sequence to initialize one cluster center.}

\item{plot_like}{Logical; if TRUE, plot likelihood at each iteration.}

\item{verbose}{Logical; if TRUE, print progress messages.}

\item{tol}{Convergence tolerance. Algorithm stops when likelihood change
is less than this value.}
}
\value{
List with components:
\describe{
  \item{R}{List of modal sequences for each cluster}
  \item{p}{Vector of cluster proportions}
  \item{lambda}{Vector of spread parameters}
  \item{datas}{Data frame with cluster assignments and diagnostics}
  \item{min.like}{Vector of log-likelihood at each iteration}
}
}
\description{
Fits a mixture of Mallows' models to partial or full ranking data using
an EM algorithm with Kendall's distance metric.
}
\author{
Erik Gregory
}
\references{
Murphy, T.B. and Martin, D. (2003). Mixtures of distance-based models
for ranking data. Computational Statistics & Data Analysis, 41, 645-655.

Adkins, L. and Fligner, M. (1998). A Non-iterative procedure for maximum
likelihood estimation of the parameters of Mallows' Model Based on
Partial Rankings. Communications in Statistics - Theory and Methods,
27:9, 2199-2220.
}
\keyword{Mallow}
\keyword{cluster}
