\name{createTissueSimDatasets}
\alias{createTissueSimDatasets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
TEST
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
createTissueSimDatasets(numSimDatasets=10,
                                   eigenValueMin = -.001,simType=c("highQualityClust","mixedClustQualityClust","unevenSizeClust"),
                                   randNumClust=FALSE,
                                   numPerClust = c(50,50,50,50),
                                   stddevNoise=0,numRows=200,minRandSize=1,maxRandSize=100,minRandNumClust=2)
}
\arguments{
  \item{numSimDatasets}{
%%     ~~Describe \code{dataMatrixList} here~~
}
  \item{eigenValueMin}{
%%     ~~Describe \code{dataMatrixList} here~~
}
  \item{simType}{
%%     ~~Describe \code{dataMatrixList} here~~
}
  \item{randNumClust}{
%%     ~~Describe \code{dataMatrixList} here~~
}
  \item{numPerClust}{
%%     ~~Describe \code{dataMatrixList} here~~
}
  \item{stddevNoise}{
%%     ~~Describe \code{dataMatrixList} here~~
}
  \item{numRows}{
%%     ~~Describe \code{dataMatrixList} here~~
}
%- maybe also 'usage' for other objects documented here.

  \item{minRandSize}{
%%     ~~Describe \code{dataMatrixList} here~~
}

  \item{maxRandSize}{
%%     ~~Describe \code{dataMatrixList} here~~
}
  \item{minRandNumClust}{
%%     ~~Describe \code{dataMatrixList} here~~
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
