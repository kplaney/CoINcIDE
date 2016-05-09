\name{procExprSet}
\alias{procExprSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Post-process a normalized assayData in an ExpressionSet object
}

\usage{
procExprSet(exprSet, outputFile = "./procExprSetMessages.txt",
                 numTopVarGenes, minVarPercentile, maxVarPercentile =
                 1, minVar, featureDataFieldName = "gene_symbol",
                 uniquePDataID = "unique_patient_ID")
}

\arguments{
  \item{exprSet}{
expressionSet S4 object with expression (assay) data, featureData and phenoData.
}
  \item{outputFile}{
Full output file path for messages that print status of post-processing the ExpressionSet.  
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
}
  \item{featureDataFieldName}{
Name for feature (gene) field in the expression set slot.
}
  \item{uniquePDataID}{
Name for unique patient ID identifier in the pData slot.
}
}
\details{
This function performs several post-processing tasks: filtering out genes and samples with high NA rates, imputing missing values, collapsing duplicated features/genes to make a unique feature list, removing any samples for which there is already a sample with the sample patient ID (duplicated samples), and filtering genes by variance.  This function is a wrapper for the functions: filterAndImputeSamples(), collapseDupProbes(), removeDuplicatedPatients(), and filterGenesByVariance(). It is is run after initial dataset normalization, such as quantile normalization on microarray datasets.
}
\value{
A post-processed S4 expressionSet.  Tests are run to confirm the final S4 object is a valid ExpressionObject before it is returned.
}

\author{
Katie Planey <katie.planey@gmail.com>
}
\seealso{
 \code{\link{procExprSetList}} 
}
\examples{
#load up the small example dataset list and use the first one
#as it has some NA values
outputPath <- "~/test.txt"
data(breastExData)
any(is.na(exprs(breastExData[[1]])))
post_procExprSet <- procExprSet(exprSet=
breastExData[[1]], outputFile = outputPath,
minVarPercentile=.75, maxVarPercentile = 1, featureDataFieldName = "gene",
uniquePDataID= "dbUniquePatientID")
#now NAs should be gone
any(is.na(exprs(post_procExprSet)))
}
