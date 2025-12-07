\name{UpdateR}
\alias{UpdateR}
\title{Update modal sequences (M-step)}
\usage{
UpdateR(r, z, infos = NULL, use_cpp = TRUE)
}
\arguments{
\item{r}{Matrix of rankings (N x n_items).}

\item{z}{Matrix of membership probabilities (N x G).}

\item{infos}{Optional precomputed Kendall information matrix for \code{r}.}

\item{use_cpp}{Logical; if TRUE (default), use the fast C++ implementation.}
}
\value{
List of G updated modal sequences.
}
\description{
Updates the cluster centers (modal sequences) to maximize the likelihood
given current membership probabilities.
}
\author{
Erik Gregory
}
\references{
Murphy, T.B. and Martin, D. (2003). Mixtures of distance-based models
for ranking data. Computational Statistics & Data Analysis, 41, 645-655.
}
\keyword{center}
\keyword{cluster}
