\name{ConstructSeqs}
\alias{ConstructSeqs}
\title{Construct sequences from Kendall preference vectors}
\usage{
ConstructSeqs(prefs, n_abils, use_cpp = TRUE)
}
\arguments{
\item{prefs}{Matrix of pairwise preferences, where each row represents
one sequence's preferences. A value of 1 indicates a decrease (first
element > second), 0 indicates an increase.}

\item{n_abils}{Number of items in the original ranking.}

\item{use_cpp}{Logical; if TRUE (default), use the fast C++ implementation.}
}
\value{
List of integer vectors, each representing a reconstructed ranking.
}
\description{
Reconstructs ranking sequences from their Kendall information representation.
Each fully-ordered sequence has a unique Kendall information vector.
}
\examples{
# Reconstruct a sequence where item 4 is ranked first
# For 4 items, there are 4*3/2 = 6 pairwise comparisons
prefs <- matrix(c(1, 1, 1, 0, 0, 0), nrow = 1)
ConstructSeqs(prefs, 4)
# Returns list(c(4, 1, 2, 3))

# Multiple sequences at once
prefs2 <- rbind(c(0, 0, 0, 0, 0, 0), c(1, 1, 1, 1, 1, 1))
ConstructSeqs(prefs2, 4)
# Returns list(c(1,2,3,4), c(4,3,2,1))
}
\author{
Erik Gregory
}
\keyword{Sequences}
