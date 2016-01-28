\name{filterAndImputeSamples}
\alias{filterAndImputeSamples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filter and Impute Samples
}
\description{
A method that removes samples or genes with high NA rates and then KNN imputes remaining missing values.
}
\usage{
filterAndImputeSamples(study, studyName = "study", 
outputFile = "createTestTrainSetsOutput.txt", 
impute = TRUE, knnFractionSize = 0.01, fractionSampleNAcutoff = 0.005,
fractionGeneNAcutoff = 0.01, exprIndex = "expr", classIndex, 
sampleCol = TRUE, returnErrorRate = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{study}{
A list, of minimally the gene expression or some molecular data matrix with keys (molecular features, such as genes) in the rows and patient samples in the columns and a keys list. It is assumed that the keys entity is named "keys", but in line with using this function for any type of molecular data, the exprIndex name in the list can be altered.
}
  \item{studyName}{
Character string to name the study. Useful in cases where looping over multiple datasets; output messages printed to the output file can then be identified by each individual study name. 
}
  \item{outputFile}{
Output File for printing progress and stats on gene/sample filtering and data imputation.  Include full directory if file should not be printed to current working directory.
}
  \item{impute}{
Impute data? A boolean TRUE or FALSE value. If FALSE, only genes and samples with high NA rates are removed, and the rest of the data is not imputed.
}
  \item{knnFractionSize}{
What is the fraction of neighbors out of the total dataset to be used for knn impute nearest neighbor? This is translated into the "k" numeric magnitude in impute.knn() from the impute package. Default is .01, or 1% of the data. 
}
  \item{fractionSampleNAcutoff}{
  
Max fraction of NAs allowed for a certain sample across all genes. Default is .005 (.005, or .5\%, still captures a large number of gnees for a sample if there are tens of thousands of genes in the data matrix.)

}
  \item{fractionGeneNAcutoff}{
  
Max fraction of NAs allowed for a certain gene across all samples. Default is .01. Thus, a certain gene cannot be missing in greater than 1\% of patients. It is recommended that this threshold be increased for smaller datasets unless a user wants a gene to be removed that is missing in only 1 sample.

}
  \item{exprIndex}{
Character string. List slot name for the data matrix, presumably an expression matrix.
}
  \item{classIndex}{
Optional character string giving the list slot name for a phenotype vector or matrix if available. If phenotype/class data such as survival is already in the list, filtering out samples with high NA rates will result in the need to remove these samples from the phenotype data matrix; filterAndImputeSamples will appropriately filter out these samples from the phenoteyp data.

}
  \item{sampleCol}{
Are samples in the columns of the expression matrix? If not, this function will first transpose the matrix to make sure impute.knn is running properly.
}
  \item{returnErrorRate}{
Boolean TRUE or FALSE. If TRUE, a small amount of real expression data points are held out, and knn.impute is performed. The accuracy rate of the imputed values vs. the real values is returned. THis is helpful in early data analysis stages to determine whether KNN imputation is appropriate for your type of data.  Default is FALSE to reduce computation time.
}
}

\value{
A list containing the following objects:
\item{expr}{original expression matrix}
\item{exprFilterImputed}{final filtered and imputed expression matrix}
\item{keys}{original keys}
\item{keys}{final filtered and imputed keys}
\item{classes}{original classes/phenotype data}
\item{classes}{final classes/phenotype data, removing any sample rows that were removed from the expression matrix after filtering.}

}

\author{

Katie Planey <katie.planey@gmail.com>

}



\examples{
library("curatedBreastData")
#load up our datasets
data(curatedBreastDataExprSetList);

#just perform on one dataset as an example, GSE9893. This dataset does have NA values.
#highestVariance calculation make take a minute to run.
#create study list object. 
study <- list(expr=exprs(curatedBreastDataExprSetList[[5]]),
keys=curatedBreastDataExprSetList[[2]]@featureData$gene_symbol,
phenoData=pData(curatedBreastDataExprSetList[[5]]))

filteredStudy <- filterAndImputeSamples(study, studyName = "study", 
outputFile = "createTestTrainSetsOutput.txt", impute = TRUE, 
knnFractionSize = 0.01, fractionSampleNAcutoff = 0.005, 
fractionGeneNAcutoff = 0.01, exprIndex = "expr", classIndex="phenoData",
sampleCol = TRUE, returnErrorRate = TRUE)

#see output list names 
names(filteredStudy)
#what is the imputation error fraction (rate)?
filteredStudy$errorRate

}