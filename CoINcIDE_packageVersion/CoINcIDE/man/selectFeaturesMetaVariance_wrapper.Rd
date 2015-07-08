\name{selectFeaturesMetaVariance_wrapper}
\alias{selectFeaturesMetaVariance_wrapper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
TEST
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
selectFeaturesMetaVariance_wrapper(dataMatrixList, rankMethod = c("sd,mad"), numNAstudiesAllowedPerFeat = ceiling(length(dataMatrixList)/10), numFeatSelectByGlobalRank = 1000, numTopFeatFromEachDataset = 20, fractNATopFeatAllowedPerDataset = 0.2, selectMethod = c("mean", "median"), outputFile = "./selectFeaturesMetaVarianceOut.txt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataMatrixList}{
%%     ~~Describe \code{dataMatrixList} here~~
}
  \item{rankMethod}{
%%     ~~Describe \code{rankMethod} here~~
}
  \item{numNAstudiesAllowedPerFeat}{
%%     ~~Describe \code{numNAstudiesAllowedPerFeat} here~~
}
  \item{numFeatSelectByGlobalRank}{
%%     ~~Describe \code{numFeatSelectByGlobalRank} here~~
}
  \item{numTopFeatFromEachDataset}{
%%     ~~Describe \code{numTopFeatFromEachDataset} here~~
}
  \item{fractNATopFeatAllowedPerDataset}{
%%     ~~Describe \code{fractNATopFeatAllowedPerDataset} here~~
}
  \item{selectMethod}{
%%     ~~Describe \code{selectMethod} here~~
}
  \item{outputFile}{
%%     ~~Describe \code{outputFile} here~~
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
