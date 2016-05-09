\name{procExprSetList}
\alias{procExprSetList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Process a list of S4 expressionSet objects.
}
\description{
A wrapper function for the post-processing function processExpressionSet() on a list of S4 expressionSet objects. This function is run after initial dataset normalization, such as quantile normalization on microarray datasets.
}

\usage{
procExprSetList(exprSetList, outputFile = "./procExprSetMessages.txt",
                 numTopVarGenes, minVarPercentile, maxVarPercentile =
                 1, minVar, featureDataFieldName = "gene_symbol",
                 uniquePDataID = "unique_patient_ID")
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exprSetList}{
List of S4 expression sets.
}
  \item{outputFile}{
Full output file path for messages that print status of post-processing the ExpressionSets.
}
  \item{minVarPercentile}{
Minimum variance percentile. Must be provided in conjunction with maxVarPercentile to use percentiles to threshold genes.
}
  \item{maxVarPercentile}{
Maximum variance percentile. Defaul is 1, i.e. 1\%. Must be provided in conjunction with minVarPercentile to use percentiles to threshold genes.
}

  \item{minVar}{
If maxVar is provided, as opposed to minVarPercentile and maxVarPercentile, genes are removed that are below a certain variance magnitude. This is helpful before running certain algorithms, such as the popular Combat batch normalization technique, that can throw errors if genes with extremely low variances are in the data matrix. May be used in conjunction with maxVar or in isolation.
}
  \item{numTopVarGenes}{
A numeric value indicating the number of genes (features) to select; the function will only take this number of genes that have the highest variance across all genes.
}  \item{featureDataFieldName}{
%%     ~~Describe \code{featureDataFieldName} here~~
}
  \item{uniquePDataID}{
%%     ~~Describe \code{uniquePDataID} here~~
}
}
\name{procExpressionSetList}
\alias{procExpressionSetList}

\value{
A list of processed S4 ExpressionSet objects.
}

\author{
Katie Planey <katie.planey@gmail.com>
}


\seealso{
 \code{\link{procExprSet}} 
}
\examples{
#load up the small example dataset list 
outputPath <- "~/test.txt"
data(breastExData)
#first one has some NA values - imputing NA
#values is part of this post-processing function.
any(is.na(exprs(breastExData[[1]])))
post_procExprSet <- procExprSetList(exprSet=
breastExData, outputFile = outputPath,
minVarPercentile=.75, maxVarPercentile = 1, featureDataFieldName = "gene", 
uniquePDataID= "dbUniquePatientID")
#now NAs should be gone
any(is.na(exprs(post_procExprSet[[1]])))
}




