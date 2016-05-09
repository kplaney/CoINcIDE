\name{findCommunities}
\alias{findCommunities}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
TEST
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
findCommunities(edgeMatrix, edgeWeightMatrix, clustIndexMatrix, fileTag = "CoINcIDE_communityNodeAttributes_", saveDir = "./", minNumUniqueStudiesPerCommunity = 3, experimentName = "exp", commMethod = c("edgeBetween", "fastGreedy", "walktrap", "eigenvector", "optimal", "spinglass", "multilevel"), makePlots = TRUE, plotToScreen = FALSE, saveGraphData = TRUE, nodeFontSize = 0.7, nodePlotSize = 10, findCommWithWeights = FALSE, plotSimilEdgeWeight = TRUE,fractEdgesInVsOutComm=0,
                            fractEdgesInVsOutEdge=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{edgeMatrix}{
%%     ~~Describe \code{edgeMatrix} here~~
}
  \item{edgeWeightMatrix}{
%%     ~~Describe \code{edgeWeightMatrix} here~~
}
  \item{clustIndexMatrix}{
%%     ~~Describe \code{clustIndexMatrix} here~~
}
  \item{fileTag}{
%%     ~~Describe \code{fileTag} here~~
}
  \item{saveDir}{
%%     ~~Describe \code{saveDir} here~~
}
  \item{minNumUniqueStudiesPerCommunity}{
%%     ~~Describe \code{minNumUniqueStudiesPerCommunity} here~~
}
  \item{experimentName}{
%%     ~~Describe \code{experimentName} here~~
}
  \item{commMethod}{
%%     ~~Describe \code{commMethod} here~~
}
  \item{makePlots}{
%%     ~~Describe \code{makePlots} here~~
}
  \item{plotToScreen}{
%%     ~~Describe \code{plotToScreen} here~~
}
  \item{saveGraphData}{
%%     ~~Describe \code{saveGraphData} here~~
}
  \item{nodeFontSize}{
%%     ~~Describe \code{nodeFontSize} here~~
}
  \item{nodePlotSize}{
%%     ~~Describe \code{nodePlotSize} here~~
}
  \item{findCommWithWeights}{
%%     ~~Describe \code{findCommWithWeights} here~~
}
  \item{plotSimilEdgeWeight}{
%%     ~~Describe \code{plotSimilEdgeWeight} here~~
}
  \item{fractEdgesInVsOutComm}{
%%     ~~Describe \code{plotSimilEdgeWeight} here~~
}

  \item{fractEdgesInVsOutEdge}{
%%     ~~Describe \code{plotSimilEdgeWeight} here~~
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
