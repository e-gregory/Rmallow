\docType{data}
\name{elect}
\alias{elect}
\title{1980 APA Presidential Candidate Ranking Data}
\format{
An integer matrix with 1378 rows (voters) and 3 columns
(Carter, Reagan, Anderson).
}
\source{
American National Election Studies, \url{https://electionstudies.org/}
}
\description{
Pre-processed version of the 1980 American Psychological Association
Presidential candidate ranking data. Uninformative rankings have been
removed and values have been pre-simplified into partial rankings.
}
\examples{
data(elect)
head(elect)
table(elect[, "Carter"])  # Distribution of Carter rankings
}
\keyword{datasets}
