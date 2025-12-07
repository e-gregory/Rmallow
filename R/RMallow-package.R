#' RMallow: Fit Multi-Modal Mallows' Models to Ranking Data
#'
#' Fits the Mallows' model to ranking data using an EM algorithm.
#' Data can be partially or fully-ranked, with or without ties.
#'
#' @section Main Functions:
#' \describe{
#'   \item{\code{\link{Mallows}}}{Fit a Mallows mixture model}
#'   \item{\code{\link{BestFit}}}{Fit multiple models and select the best}
#' }
#'
#' @section Distance Functions:
#' \describe{
#'   \item{\code{\link{AllKendall}}}{Compute Kendall distances between rankings}
#'   \item{\code{\link{KendallInfo}}}{Compute pairwise comparison information}
#'   \item{\code{\link{DistanceDistribution}}}{Distance distribution in N! space}
#' }
#'
#' @section Utility Functions:
#' \describe{
#'   \item{\code{\link{SimplifySequences}}}{Normalize tied rankings}
#'   \item{\code{\link{ConstructSeqs}}}{Reconstruct sequences from preferences}
#' }
#'
#' @name RMallow-package
#' @aliases RMallow-package RMallow
#' @docType package
#' @author Erik Gregory \email{egregory2011@@gmail.com}
#'
#' @references
#'   Murphy, T.B. and Martin, D. (2003). Mixtures of distance-based models
#'   for ranking data. Computational Statistics & Data Analysis, 41, 645-655.
#'   \doi{10.1016/S0167-9473(02)00165-2}
#'
#'   Beckett, L.A. and Evans, D.A. (1991). Estimating a Population Distribution
#'   of Sequences of k Items from Cross-Sectional Data. Journal of the Royal
#'   Statistical Society, Series C, 40(1), 31-42.
#'
#'   Adkins, L. and Fligner, M. (1998). A Non-iterative procedure for maximum
#'   likelihood estimation of the parameters of Mallows' Model Based on
#'   Partial Rankings. Communications in Statistics - Theory and Methods,
#'   27(9), 2199-2220. \doi{10.1080/03610929808832223}
#'
#' @keywords package ranking Mallows mixture model clustering
#' @importFrom combinat permn combn
#' @importFrom grDevices dev.new dev.off
#' @importFrom graphics plot
#' @importFrom stats runif uniroot
"_PACKAGE"

#' Sample synthetic data set
#'
#' Synthetic data set containing 3 modal sequences in 15! space, with
#' some noise added. Useful for testing and demonstrating the package.
#'
#' @name datas
#' @docType data
#' @format A numeric matrix with 1700 rows (observations) and 15 columns (items).
#' @keywords datasets
#' @examples
#' data(datas)
#' head(datas)
#' dim(datas)
NULL

#' 1980 APA Presidential Candidate Ranking Data
#'
#' Pre-processed version of the 1980 American Psychological Association
#' Presidential candidate ranking data. Uninformative rankings have been
#' removed and values have been pre-simplified into partial rankings.
#'
#' @name elect
#' @docType data
#' @format An integer matrix with 1378 rows (voters) and 3 columns
#'   (Carter, Reagan, Anderson).
#' @source American National Election Studies,
#'   \url{https://electionstudies.org/}
#' @keywords datasets
#' @examples
#' data(elect)
#' head(elect)
#' table(elect[, "Carter"])  # Distribution of Carter rankings
NULL

#' Three-mode Mallows model fit to synthetic data
#'
#' The \code{datas} data set fitted with 3 modal sequences.
#' Compare to \code{\link{two.mode}} to see the effect of model selection.
#'
#' @name three.mode
#' @docType data
#' @format A list with components:
#' \describe{
#'   \item{R}{List of 3 modal sequences}
#'   \item{p}{Vector of cluster proportions}
#'   \item{lambda}{Vector of spread parameters}
#'   \item{datas}{Data frame with cluster assignments}
#'   \item{min.like}{Log-likelihood trajectory}
#' }
#' @keywords datasets
#' @seealso \code{\link{two.mode}}, \code{\link{datas}}
#' @examples
#' data(three.mode)
#' three.mode$R  # Modal sequences
#' three.mode$p  # Cluster proportions
NULL

#' Two-mode Mallows model fit to synthetic data
#'
#' The \code{datas} data set fitted with only 2 modal sequences.
#' Since the true data has 3 modes, this shows what happens with
#' model misspecification.
#'
#' @name two.mode
#' @docType data
#' @format A list with components:
#' \describe{
#'   \item{R}{List of 2 modal sequences}
#'   \item{p}{Vector of cluster proportions}
#'   \item{lambda}{Vector of spread parameters}
#'   \item{datas}{Data frame with cluster assignments}
#'   \item{min.like}{Log-likelihood trajectory}
#' }
#' @keywords datasets
#' @seealso \code{\link{three.mode}}, \code{\link{datas}}
#' @examples
#' data(two.mode)
#' two.mode$R  # Only 2 modal sequences
NULL

#' Two-mode Mallows model fit to APA election data
#'
#' The \code{elect} data fitted with 2 modal sequences. The two modes
#' appear to correspond roughly to Democratic and Republican voting patterns.
#'
#' @name two.seq
#' @docType data
#' @format A list with components:
#' \describe{
#'   \item{R}{List of 2 modal sequences}
#'   \item{p}{Vector of cluster proportions}
#'   \item{lambda}{Vector of spread parameters}
#'   \item{datas}{Data frame with cluster assignments}
#'   \item{min.like}{Log-likelihood trajectory}
#' }
#' @source American National Election Studies
#' @keywords datasets
#' @seealso \code{\link{elect}}
#' @examples
#' data(two.seq)
#' two.seq$R  # Modal sequences (voting patterns)
#' two.seq$p  # Proportion in each group
NULL
