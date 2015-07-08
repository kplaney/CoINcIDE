\name{computeAdjMatrices}
\alias{computeAdjMatrices}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
TEST
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
computeAdjMatrices(dataMatrixList, clustSampleIndexList, clustFeatureIndexList, edgeMethod = c("spearman", "pearson", "kendall", "Euclidean", "cosine", "Manhattan", "Minkowski", "Mahalanobis"), numParallelCores = 1, minTrueSimilThresh = -Inf, maxTrueSimilThresh = Inf, numSims = 500, outputFile = "./CoINcIDE_messages.txt", fractFeatIntersectThresh = 0, numFeatIntersectThresh = 0, clustSizeThresh = 0, clustSizeFractThresh = 0, checkNA = FALSE, centroidMethod = c("mean", "median"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataMatrixList}{
%%     ~~Describe \code{dataMatrixList} here~~
}
  \item{clustSampleIndexList}{
%%     ~~Describe \code{clustSampleIndexList} here~~
}
  \item{clustFeatureIndexList}{
%%     ~~Describe \code{clustFeatureIndexList} here~~
}
  \item{edgeMethod}{
%%     ~~Describe \code{edgeMethod} here~~
}
  \item{numParallelCores}{
%%     ~~Describe \code{numParallelCores} here~~
}
  \item{minTrueSimilThresh}{
%%     ~~Describe \code{minTrueSimilThresh} here~~
}
  \item{maxTrueSimilThresh}{
%%     ~~Describe \code{maxTrueSimilThresh} here~~
}
  \item{numSims}{
%%     ~~Describe \code{numSims} here~~
}
  \item{outputFile}{
%%     ~~Describe \code{outputFile} here~~
}
  \item{fractFeatIntersectThresh}{
%%     ~~Describe \code{fractFeatIntersectThresh} here~~
}
  \item{numFeatIntersectThresh}{
%%     ~~Describe \code{numFeatIntersectThresh} here~~
}
  \item{clustSizeThresh}{
%%     ~~Describe \code{clustSizeThresh} here~~
}
  \item{clustSizeFractThresh}{
%%     ~~Describe \code{clustSizeFractThresh} here~~
}
  \item{checkNA}{
%%     ~~Describe \code{checkNA} here~~
}
  \item{centroidMethod}{
%%     ~~Describe \code{centroidMethod} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
