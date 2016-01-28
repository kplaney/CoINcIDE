\name{collapseDupProbes}
\alias{collapseDupProbes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Collapse/handle duplicated probes (genes) in a dataset
}
\description{
Used internally by processExpressionSet.  Code to either take the average across a set of duplicated "keys" (can be probes or genes, which correspond to the rows in the expression matrix "expr"), or take the keys that has the highest variance across the set of duplicated keys.
}
\usage{
collapseDupProbes(expr, sampleColNames=colnames(expr), 
keys, method = c("average", "highestVariance"), debug = TRUE, removeNA_keys = TRUE, 
varMetric = c("everything", "all.obs", "complete.obs", "na.or.complete",
"pairwise.complete.obs"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expr}{
An expression matrix with genes in the rows and samples in the columns.
}
  \item{sampleColNames}{
Sample column names. Needed for internal debugging; usually the default colnames(expr) is appropriate.
}
  \item{keys}{
Generally the list of gene symbols, or some molecular key, that needs to be "collapsed" because it contains duplicated names.
}
  \item{method}{
Method used to collapse probes: take the mean across all duplicated keys, or just pick the key with the highest variance?
}
  \item{debug}{
Use internal unit tests that will stop the code if it detects a bug?
}
  \item{removeNA_keys}{
Remove any NA keys?
}
  \item{varMetric}{
Standard options taken from the base var() function. May be important if you have NA values in your data matrix; otherwise, "everything" is usually fine.
}
}

\value{
Returns a processed list with the items "expr" and "keys", the expression matrix and final keys list.
}

\author{
Katie Planey <katie.planey@gmail.com>
}

\examples{
library("curatedBreastData")
#load up our datasets
data(curatedBreastDataExprSetList);

#just perform on second dataset, GSE2034, as an example.
#This dataset has no NAs already but does have duplicated genes
#highestVariance calculation make take a minute to run.
collapsedData <- collapseDupProbes(expr=exprs(curatedBreastDataExprSetList[[2]]),  
keys=curatedBreastDataExprSetList[[2]]@featureData$gene_symbol, 
method = c("highestVariance"), debug = TRUE, removeNA_keys = TRUE, 
varMetric = c("everything"))
#look at names of outputs
names(collapsedData)

}