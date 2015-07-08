\name{global_FDR}
\alias{global_FDR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
TEST
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
global_FDR(CoINcIDE_outputList, edgeMethod = c("distCor", "spearman", "pearson", "kendall", "Euclidean", "cosine", "Manhattan", "Minkowski", "Mahalanobis"), minTrueSimilThresh = -Inf, maxTrueSimilThresh = Inf, outputFile = "./CoINcIDE_messages.txt", fractFeatIntersectThresh = 0, numFeatIntersectThresh = 0, clustSizeThresh = 0, clustSizeFractThresh = 0, meanEdgePairPvalueThresh = 0.01, indEdgePvalueThresh = 0.05, saveDir = "./", experimentName = "nullTest", commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity = 3,  minFractNN = 0.8, findCommWithWeights = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CoINcIDE_outputList}{
%%     ~~Describe \code{CoINcIDE_outputList} here~~
}
  \item{edgeMethod}{
%%     ~~Describe \code{edgeMethod} here~~
}

  \item{minTrueSimilThresh}{
%%     ~~Describe \code{minTrueSimilThresh} here~~
}
  \item{maxTrueSimilThresh}{
%%     ~~Describe \code{maxTrueSimilThresh} here~~
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
  \item{meanEdgePairPvalueThresh}{
%%     ~~Describe \code{meanEdgePairPvalueThresh} here~~
}
  \item{indEdgePvalueThresh}{
%%     ~~Describe \code{indEdgePvalueThresh} here~~
}
  \item{saveDir}{
%%     ~~Describe \code{saveDir} here~~
}
  \item{experimentName}{
%%     ~~Describe \code{experimentName} here~~
}

  \item{commMethod}{
%%     ~~Describe \code{commMethod} here~~
}
  \item{minNumUniqueStudiesPerCommunity}{
%%     ~~Describe \code{minNumUniqueStudiesPerCommunity} here~~
}

  \item{minFractNN}{
%%     ~~Describe \code{minFractNN} here~~
}
  \item{findCommWithWeights}{
%%     ~~Describe \code{findCommWithWeights} here~~
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
