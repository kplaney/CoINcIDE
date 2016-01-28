\name{removeDupPatients}
\alias{removeDupPatients}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Remove duplicated patient samples (samples from the same patient/column ID)
}
\description{
Function to keep only 1 sample per patient (column ID) in the data matrix. Keeps the sample that has the overall highest variance.
}
\usage{
removeDupPatients(exprMatrix, outputFile = "duplicatedPatientsOutput.txt", 
varMetric = c("everything", "all.obs", "complete.obs", "na.or.complete",
"pairwise.complete.obs"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exprMatrix}{
Expression matrix, with features in the rows and samples in the columns.
}
  \item{outputFile}{
Output file for messages that print status of removing duplicated samples.  Include full directory if file should not be printed to current working directory.
}
  \item{varMetric}{
Standard options taken from the base var() function. May be important if you have NA values in your data matrix; otherwise, "everything" is usually fine.
}
}

\value{
exprMatrix:  the final data matrix with only 1 sample per patient ID.
}

\author{
Katie Planey <katie.planey@gmail.com>
}
\note{
Suggestions are welcome for further ways to pick the best sample from samples from the same patient. No curatedBreastData matrices currently have samples that share the same patient ID, but this function is especially useful for say TCGA data, where this is often the case. 

It is suggseted one imputes missing values using the filterAndImpute function before running this function to avoid -Inf and NA values in the variance calculations.
}


\examples{
library("curatedBreastData")
#No curatedBreastData has duplicated samples, 
#but we can still run this function on one of the datasets:
#load up our datasets
data(curatedBreastDataExprSetList);

#This dataset does not have NA values, which makes for a good example without extra pre-processing.
outputMatrix <- removeDupPatients(exprMatrix=
exprs(curatedBreastDataExprSetList[[1]]), 
outputFile = "./duplicatedPatientsOutput.txt", varMetric = c("everything"))
#final dimensions - unchanged in this case with 
#no samples sharing the same patient ID.
dim(outputMatrix)
}