\name{createS4exprSet}
\alias{createS4exprSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
createS4exprSet
}
\description{
Create an S4 exprSet object from a standard R list containing gene expression, 
feature and phenotype tables.
}
\usage{
createS4exprSet(expr, phenoData, featureData, featureDataFieldName = "gene")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expr}{
A matrix of (presumably gene) expression data, with samples (patients) in the 
columns and features (genes) in the rows.
}
  \item{phenoData}{
A matrix of phenotypic data, i.e. covariates that are measured for each sample. 
Each row corresponds to one sample (patient), and each column is a covariate. 
For example, if this is clinical trial data, a covariate may be whether the 
patient survived or died.
}
  \item{featureData}{
A matrix of feature data, with feature type names in the columns and feature 
indices in the rows.  This is usually probe or gene symbol names in the rows, 
and just 
"probe" or "gene" symbol as the type/column name. There can be greater than one 
feature type/column, but one must contain no duplicates, as this is required for 
the featureDataFieldName to assign unique row names to the expr matrix.

}
  \item{featureDataFieldName}{
The unique feature key column name to properly index the featureData matrix to 
grab this key and make it the row names of the expr matrix. No duplicates are 
allowed, otherwise the eSet check function won't pass.
}
}
\value{
An exprSet; this is a widely used S4 object that is a subset of the broader eSet 
object class.  It allows you to use your expression data along with covariates 
for many Bioconductor packages that require your data to be in this structure. 
See references for "helper" functions you can then use to extract your data from 
this S4 object.
}
\references{
\url{https://www.bioconductor.org/packages/3.3/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf}
}
\author{
Katie Planey <katie.planey@gmail.com>
}

\examples{
#just put one of the small example breast cancer datasets into a regular list 
#form, and use this function to get it back into an exprSet.
data(breastExData)
exprData <- exprs(breastExData[[1]])
phenData <- pData(breastExData[[1]])
featData <- rownames(exprData)
featData <- data.frame(featData)
#pick the featureDataFieldName (just some numbers here for this toy dataset)
colnames(featData) <- "gene"
myExprSet = createS4exprSet(expr = exprData, phenoData = phenData, 
                            featureData = featData, featureDataFieldName = 
                              "gene")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~exprSet }
\keyword{ ~eSet }% __ONLY ONE__ keyword per line
