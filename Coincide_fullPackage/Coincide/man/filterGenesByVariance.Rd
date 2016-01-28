\name{filterGenesByVariance}
\alias{filterGenesByVariance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filter genes by variance
}
\description{
A function that filters genes by variance; it can simply threshold out genes that are above or below a certain magnitude of variance, filter out genes that fall outside of a minimum and maximum percentile, or simply select the top N varying genes.
}
\usage{
filterGenesByVariance(study, plotSaveDir = "~/", minVarPercentile, 
maxVarPercentile=1, maxVar, minVar, exprIndex = "expr", 
keysIndex = "keys", outputFile = "varCal.txt", plotVarianceHist = FALSE,
varMetric = c("everything", "all.obs", "complete.obs", 
"na.or.complete", "pairwise.complete.obs"), 
sampleCol = TRUE, numTopVarGenes)
}

\arguments{
  \item{study}{
A list, of minimally the gene expression or some molecular data matrix with keys (molecular features, such as genes) in the rows and patient samples in the columns and a keys list. In line with using this function for any type of molecular data, the exprIndex name, and also the keysIndex name, in the list can be altered.
}

  \item{plotSaveDir}{
If plotVarianceHist is TRUE, then the plotSaveDir is a character string specifying where this histogram plot should be saved.
}
  \item{minVarPercentile}{
Minimum variance percentile. Must be provided in conjunction with maxVarPercentile to use percentiles to threshold genes.
}
  \item{maxVarPercentile}{
Maximum variance percentile. Defaul is 1, i.e. 1\%. Must be provided in conjunction with minVarPercentile to use percentiles to threshold genes.
}
  \item{maxVar}{
If maxVar is provided, as opposed to minVarPercentile and maxVarPercentile, genes are removed that are above a certain variance magnitude. This may be useful if a user suspects very highly varying genes are actually technical noise/outliers. May be used in conjunction with minVar or in isolation.
}
  \item{minVar}{
If maxVar is provided, as opposed to minVarPercentile and maxVarPercentile, genes are removed that are below a certain variance magnitude. This is helpful before running certain algorithms, such as the popular Combat batch normalization technique, that can throw errors if genes with extremely low variances are in the data matrix. May be used in conjunction with maxVar or in isolation.
}
  \item{exprIndex}{
Character string. List slot name for the data matrix, presumably an expression matrix.
}
  \item{keysIndex}{
Character string. List slot name for the feature names, presumably probes or gene names.
}
  \item{outputFile}{
Output file for messages that print status of the filtering.  Include full directory if file should not be printed to current working directory.
}
  \item{plotVarianceHist}{
Plot the histogram of variances overall? Good for exploratory analyses to understand the distribution of variance across all data points. Default is FALSE to avoid saving a ggplot image for every function run.
}
  \item{varMetric}{
Standard options taken from the base var() function. May be important if you have NA values in your data matrix; otherwise, "everything" is usually fine.
}
  \item{sampleCol}{
Are samples in the columns of the expression matrix? If not, this function will first transpose the matrix, as the function assumes samples are in the columns features are in the rows.
}
  \item{numTopVarGenes}{
A numeric value indicating the number of genes (features) to select; the function will only take this number of genes that have the highest variance across all genes.
}
}

\value{
A list:  output <- list(study=study,filteredStudy=filteredStudy,p=p);
\item{study}{Original study list object}
\item{filteredStudy}{filteredStudy object, i.e. the gene expression and keys only for the desired filtered keys/features.}
}

\author{
Katie Planey <katie.planey@gmail.com>
}
\note{
Filtering by variance is equivalent to filtering on the coefficient of variation if data is logged.  Further work includes automatically allowing the user to use the coefficient of variation as opposed to baseline variation for a threshold.

It is highly suggested you use filterAndImputeSamples() beforehand to remove any NA values, to avoid -Inf or NA variance calculations.
}


\examples{
library("curatedBreastData")
#load up our datasets
data(curatedBreastDataExprSetList);

#just perform on one dataset as an example, GSE1379. 
#This dataset does not have NA values, which makes for a
#good example without extra pre-processing.
#highestVariance calculation make take a minute to run.
#create study list object. 
study <- list(expr=exprs(curatedBreastDataExprSetList[[1]]),
keys=curatedBreastDataExprSetList[[1]]@featureData$gene_symbol)
#take top 100 varying genes

filterGeneStudy <- filterGenesByVariance(study, exprIndex = "expr", 
keysIndex = "keys", outputFile = "./varCal.txt", 
plotVarianceHist = FALSE,
varMetric = c("everything"), sampleCol = TRUE, numTopVarGenes=100)

#names of output
names(filterGeneStudy)
}