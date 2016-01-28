\name{clustMatrixListWrapper}
\alias{clustMatrixListWrapper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
TEST
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
clustMatrixListWrapper(dataMatrixList, clustFeaturesList, clustMethod = c("km", "hc"), pickKMethod = c("gap", "consensus"), numSims = 1000, maxNumClusters = 30, outputFile = "./cluster_output.txt", iter.max = 10, nstart = 1, distMethod = c("euclidean", "pearson", "spearman", "binary", "maximum", "canberra", "minkowski"), hclustAlgorithm = c("average", "complete", "ward.D", "ward.D2", "single", "mcquitty", "median", "centroid"), consensusHclustAlgorithm = c("average", "complete", "ward.D", "ward.D2", "single", "mcquitty", "median", "centroid"), minClustConsensus = 0.7, minMeanClustConsensus = 0.7, maxPAC = 0.1, corUse = "everything", pItem = 0.9)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataMatrixList}{
%%     ~~Describe \code{dataMatrixList} here~~
}
  \item{clustFeaturesList}{
%%     ~~Describe \code{clustFeaturesList} here~~
}
  \item{clustMethod}{
%%     ~~Describe \code{clustMethod} here~~
}
  \item{pickKMethod}{
%%     ~~Describe \code{pickKMethod} here~~
}
  \item{numSims}{
%%     ~~Describe \code{numSims} here~~
}
  \item{maxNumClusters}{
%%     ~~Describe \code{maxNumClusters} here~~
}
  \item{outputFile}{
%%     ~~Describe \code{outputFile} here~~
}
  \item{iter.max}{
%%     ~~Describe \code{iter.max} here~~
}
  \item{nstart}{
%%     ~~Describe \code{nstart} here~~
}
  \item{distMethod}{
%%     ~~Describe \code{distMethod} here~~
}
  \item{hclustAlgorithm}{
%%     ~~Describe \code{hclustAlgorithm} here~~
}
  \item{consensusHclustAlgorithm}{
%%     ~~Describe \code{consensusHclustAlgorithm} here~~
}
  \item{minClustConsensus}{
%%     ~~Describe \code{minClustConsensus} here~~
}
  \item{minMeanClustConsensus}{
%%     ~~Describe \code{minMeanClustConsensus} here~~
}
  \item{maxPAC}{
%%     ~~Describe \code{maxPAC} here~~
}
  \item{corUse}{
%%     ~~Describe \code{corUse} here~~
}
  \item{pItem}{
%%     ~~Describe \code{pItem} here~~
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
