curatedOvarianData meta-clustering with 500 (meta-rank) features using kmeans consensus and pearson mean matrix correlation
========================================================
**Data loading**


This is a report for the curatedOvarianData meta-clustering. Initial clustering was done using k-means and selecting k via the the rounded PACR metric with nstart=1.

First, I loaded up my data just past the CoINcIDE" "getAdjMatrices() function (that's the function that takes  a really long time to run and is hogging up server space, as this is the part that computes the edge p-values):


```r
#grab data matrix list, clust features list
CoINcIDE_rankFeatures <- "/home/kplaney/ovarian_analysis/metaFeatures_500.RData.gzip"
load(CoINcIDE_rankFeatures)
metaFeatures <- metaFeatures
CoINcIDE_clusterOutput <-  "/home/kplaney/ovarian_analysis/curatedOvarianData_kmeansConsensus_nstart1_500Features_2015-04-29.RData.gzip"
load(CoINcIDE_clusterOutput)
clusterOutput <- kmeansConsensus
load("/home/kplaney/ovarian_analysis/esets_proc_TCGAcombat.RData.gzip")


source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
```

```
## Loading required package: BiocGenerics
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## 
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## 
## The following object is masked from 'package:stats':
## 
##     xtabs
## 
## The following objects are masked from 'package:base':
## 
##     Filter, Find, Map, Position, Reduce, anyDuplicated, append,
##     as.data.frame, as.vector, cbind, colnames, do.call,
##     duplicated, eval, evalq, get, intersect, is.unsorted, lapply,
##     mapply, match, mget, order, paste, pmax, pmax.int, pmin,
##     pmin.int, rank, rbind, rep.int, rownames, sapply, setdiff,
##     sort, table, tapply, union, unique, unlist
## 
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```r
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_communityDetection.R")
```

```
## 
## Attaching package: 'limma'
## 
## The following object is masked from 'package:BiocGenerics':
## 
##     plotMA
```

```r
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_visualization.R")
```

```
## 
## Attaching package: 'gplots'
## 
## The following object is masked from 'package:stats':
## 
##     lowess
```

```r
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_metaClusterAnalysis.R")
```

```
## Loading required package: annotate
## Loading required package: AnnotationDbi
## Loading required package: GenomeInfoDb
## Loading required package: graph
## 
## Attaching package: 'graph'
## 
## The following objects are masked from 'package:igraph':
## 
##     degree, edges
## 
## 
## Attaching package: 'plyr'
## 
## The following object is masked from 'package:graph':
## 
##     join
## 
## Loading required package: bitops
## Loading required package: GOstats
## Loading required package: Category
## Loading required package: Matrix
## 
## Attaching package: 'Matrix'
## 
## The following objects are masked from 'package:base':
## 
##     crossprod, tcrossprod
## 
## Loading required package: GO.db
## Loading required package: DBI
## 
## 
## Attaching package: 'GOstats'
## 
## The following object is masked from 'package:AnnotationDbi':
## 
##     makeGOGraph
## 
## 
## Attaching package: 'RDAVIDWebService'
## 
## The following object is masked from 'package:GSEABase':
## 
##     ids
## 
## The following object is masked from 'package:AnnotationDbi':
## 
##     species
## 
## The following objects are masked from 'package:igraph':
## 
##     is.connected, membership
## 
## The following object is masked from 'package:BiocGenerics':
## 
##     counts
```

```r
#now format just as a list of data matrices.
dataMatrixList <- exprSetListToMatrixList(esets,featureDataFieldName="gene")

names(dataMatrixList) <- names(esets)
##do for each 200,500,1000,2000 (load different metaFeatures RData object each time.)
#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

origToNewIndexMap <- cbind(1:length(dataMatrixList),na.omit(match(names(dataMatrixList),names(esets))))

clustSampleIndexList <-  kmeansConsensus$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensus$clustFeatureIndexList_PACR


CoINcIDE_computeEdgesObject <- "/home/kplaney/ovarian_analysis/adjMatrices_500F_spearman_meanMatrix_2015-04-29RData.gzip"
load(CoINcIDE_computeEdgesObject)
output <- ov_500F_spearman_meanMatrix
inputVariablesDF <- output$inputVariablesDF
computeTrueSimilOutput <- output$computeTrueSimilOutput
pvalueMatrix <- output$pvalueMatrix
clustIndexMatrix <- output$clustIndexMatrix

###inputs for edge detection
meanEdgePairPvalueThresh <- .01
indEdgePvalueThresh <- .05
minTrueSimilThresh <- .5
maxTrueSimilThresh <- Inf
clustSizeFractThresh <- inputVariablesDF$clustSizeFractThresh
#KEEP: small cluster sizes: it turns out one TCGA cluster is mostly linked to all smaller clusters.
clustSizeThresh <- 5
fractFeatIntersectThresh <- inputVariablesDF$fractFeatIntersectThresh
  numFeatIntersectThresh <- inputVariablesDF$numFeatIntersectThresh
saveDir <- "/home/kplaney/ovarian_analysis/"
```

**Input variables used to derive adjacency matrix**



```r
message("Input variables used to derive the adjacency matrix:\n")
```

```
## Input variables used to derive the adjacency matrix:
```

```r
inputVariablesDF
```

```
##                  date edgeMethod numParallelCores minTrueSimilThresh
## 1 2015-04-29 07:06:19   spearman                3                0.3
##   maxTrueSimilThresh  sigMethod maxNullFractSize numSims
## 1                Inf meanMatrix              0.2     500
##   includeRefClustInNull fractFeatIntersectThresh numFeatIntersectThresh
## 1                  TRUE                      0.8                    425
##   clustSizeThresh clustSizeFractThresh
## 1               5                 0.05
```

```r
message("There were ",nrow(clustIndexMatrix), " total input clusters from ",length(unique(clustIndexMatrix[,2])), " studies")
```

```
## There were 74 total input clusters from 23 studies
```

```r
message("The total number of input features was ",length(metaFeatures$finalFeatures))
```

```
## The total number of input features was 527
```

```r
message("Across the entire square (nonsymmetric) p-value matrix, there are ",length(which(pvalueMatrix<=.1))," pvalues less than or equal to .1")
```

```
## Across the entire square (nonsymmetric) p-value matrix, there are 1113 pvalues less than or equal to .1
```

```r
message("Across the entire square (nonsymmetric) p-value matrix, there are ",length(which(pvalueMatrix<=.05))," pvalues less than or equal to .05")
```

```
## Across the entire square (nonsymmetric) p-value matrix, there are 1020 pvalues less than or equal to .05
```

```r
message("Across the entire square (symmetric) similarity matrix, there are ",length(which(computeTrueSimilOutput$similValueMatrix>=minTrueSimilThresh))," similarities greater than or equal to ",minTrueSimilThresh)
```

```
## Across the entire square (symmetric) similarity matrix, there are 1346 similarities greater than or equal to 0.5
```

**Network Analysis**

As we can see in these plots, 4 meta-clusters remained after edge filtering. 

```r
finalEdgeInfo <- assignFinalEdges(computeTrueSimilOutput=computeTrueSimilOutput,pvalueMatrix=pvalueMatrix,indEdgePvalueThresh=indEdgePvalueThresh,
                              meanEdgePairPvalueThresh=meanEdgePairPvalueThresh,
                              minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                              fractFeatIntersectThresh=fractFeatIntersectThresh,numFeatIntersectThresh=numFeatIntersectThresh ,
                              clustSizeThresh=clustSizeThresh, clustSizeFractThresh= clustSizeFractThresh,saveDir=saveDir,
                              fileTag="CoINcIDE_edges_",restrictEdges=FALSE
)
```

```
## 17 clusters dropped because they were below the clust size thresh clustSizeThresh
##           threshold of 5
```

```
## Warning in dir.create(saveDir): '/home/kplaney/ovarian_analysis' already
## exists
```

```
## A total of 30 clusters removed because they have no significant edges.
## A total of 44 clusters have significant edges.
```

```r
commInfo <- findCommunities(edgeMatrix=finalEdgeInfo$filterEdgeOutput$edgeMatrix,edgeWeightMatrix=finalEdgeInfo$filterEdgeOutput$edgeWeightMatrix,
                                        clustIndexMatrix=output$clustIndexMatrix,fileTag="autoReport",
                            saveDir=saveDir,minNumUniqueStudiesPerCommunity=3,clustMethodName="",
                            commMethod=c("edgeBetween"),
                            makePlots=TRUE,saveGraphData=FALSE,plotToScreen=TRUE)
```

```
## Warning in dir.create(saveDir, showWarnings = TRUE):
## '/home/kplaney/ovarian_analysis' already exists
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-2.png) 

```r
aggregateData <- returnSampleMemberMatrix(clustSampleIndexList=clustSampleIndexList,
                                          dataMatrixList=dataMatrixList,communityInfo=commInfo)
clustSizes <- table(aggregateData$sampleClustCommKey$globalClustNum)
#add clusterSizes to node attributes
commInfo$attrDF$size <- clustSizes[na.omit(match(commInfo$attrDF$clust,names(clustSizes)))]

networkStats <- advancedNetworkPlots(communityMembership=commInfo,
                                  brewPal = c("Set3"),
                                  saveDir=saveDir,saveName="network",
                      plotToScreen=TRUE)$network_stats
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-3.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-4.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-5.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-6.png) 

```r
message("Overall network stats (I believe the origNumClusters variable may be off-need to debug:")
```

```
## Overall network stats (I believe the origNumClusters variable may be off-need to debug:
```

```r
networkStats
```

```
##  numCommunities     numClusters origNumClusters        numEdges 
##               3              38              38             147 
##      numStudies 
##              18
```

```r
#save variable
sampleClustCommKey<-aggregateData$sampleClustCommKey
```

**Gene meta-rank Analysis**


I ranked genes within each meta-cluster for all samples, and then ran a Kruskal test to see which genes significantly stratified/differentiated patients across the 4 meta-clusters. I still need to implement GSEA; it turns out there's a base GSEA package in Biocondcutor so I've decided to just adapt my code and use their baseline functions.


In the heatmap: red means that gene was ranked high in terms of expression level for patients in that meta-cluster (I took the median rank across all samples in a meta-cluster to create the heatmap. Only significant genes are shown in the heatmap but at over 80 significant genes, of course the gene names are illegible...I print out the top 20 genes below.) The actual meta-cluster numbers are not 1:4 because meta-clusters with only 2 studies were thresholded out earlier.


```r
binInfo <- binarizeMetaclustStudyStatus(aggregateData$sampleClustCommKey)
rankInfo <- computeRankMatrix(metaClustSampleNames=binInfo$metaClustSampleNames,
                              featureNames=metaFeatures$finalFeatures,dataMatrixList,
                              sampleClustCommKey=aggregateData$sampleClustCommKey,onlyIntersectingFeat=TRUE)
pvalueInfo <- computeFeaturePvalues(rankMatrix=rankInfo$rankMatrix,featureNames=rankInfo$filteredFeatures,groupings=rankInfo$groupings)
message("There are ",length(which(pvalueInfo$fdr.qvalue<=.05)), " genes with an FDR corrected p-value below .05 whose ranks significantly differed among all of the meta-clusters (using Kruskal's test)")
```

```
## There are 153 genes with an FDR corrected p-value below .05 whose ranks significantly differed among all of the meta-clusters (using Kruskal's test)
```

```r
#message("Top 20 significant genes:\n")
#rownames(pvalueInfo[which(pvalueInfo$fdr.qvalue<=.05)[1:20],])
#cat("\n")

metaMatrix <- commMedianRank(rankInfo$rankMatrix[which(pvalueInfo$fdr.qvalue<=.05),],rankInfo$groupings)
message("Red in heatmap means genes were ranked higher across all samples in that meta-cluster.")
```

```
## Red in heatmap means genes were ranked higher across all samples in that meta-cluster.
```

```r
plotMetaFeatureRank(metaMatrix,saveFile=FALSE,plotToScreen=TRUE,
                             saveDir=saveDir,fileTag="test",
                             plotTitle="Median rank\nacross all samples/studies",
                             key.xlab="")
```

```
## Warning in heatmap.2(metaMatrix, Rowv = FALSE, Colv = FALSE, cexRow = 1, :
## Discrepancy: Rowv is FALSE, while dendrogram is `none'. Omitting row
## dendogram.
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png) 

```r
message("Running effect size analysis:")
```

```
## Running effect size analysis:
```

```r
ES_out <- computeMetaclustEffectSizes(metaClustSampleNames=binInfo$metaClustSampleNames,dataMatrixList=dataMatrixList,featureNames=metaFeatures$finalFeatures,minOtherClass=5)
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(t(matrix1[featureNames[g], ]),
## t(matrix2[featureNames[g], : cannot compute exact p-value with ties
```

```r
sigGenes <- selectMetaclustSigGenes(computeMetaclustEffectSizesOutput=ES_out,qvalueThresh=.1,
                                    ESthresh=.2)
summSigGenes <- summarizePosMetaclustGenes(selectMetaclustSigGenesOut=sigGenes,computeMetaclustEffectSizesOutput)
```

```
## Error in colnames(computeMetaclustEffectSizesOutput$summWilcoxon_qvalue): error in evaluating the argument 'x' in selecting a method for function 'colnames': Error: object 'computeMetaclustEffectSizesOutput' not found
```

**GSEA Analysis**

It turns out "GSEABase" isn't all that great..so just using some Broad gene sets I already downloaded and the standard hypergeometric tests:


```r
#this doesn't work too well - but that is perhaps because I'm calculating this wrong...
#this is filtering on wilcoxon p-value; by ES alone gives too many results for certain meta-clusters.
#(could filter by higher ES in future.)
DAVID_out <- list()
for(c in 1:length(sigGenes$sigMetaclustGenes_pos)){
  
  testGeneSet <- sigGenes$sigMetaclustGenes_pos[[c]]
  cat("\n",length(testGeneSet), " genes for commmunity number ",c, "(name) ",names(sigGenes$sigMetaclustGenes_pos)[c],
      " have a positive effect size above .2 with a wilcoxon rank q-value below .1")
  #write to feed into David for now?
  write.table(testGeneSet,file=paste0("/home/kplaney/ovarian_analysis/DAVID/upGenes_comm",
              names(sigGenes$sigMetaclustGenes_pos)[c],"_spearman500F.txt"),quote=FALSE,
              row.names=FALSE,col.names=FALSE)
  
  DAVID_out[[c]] <- GSEA_DAVID(testGeneVector=testGeneSet,idType="GENE_SYMBOL",yourEmail="katie.planey@stanford.edu",
                       EASE_thresh=0.1,count_thresh=2L,qvalueThreshFunct=.2,
                       pvalueThreshCluster=.2,fileTag=paste0("spearman500F_commmunity_",c),
                       saveDir="/home/kplaney/ovarian_analysis/DAVID/")
  
  message("Significantly enriched DAVID functional groups for overexpressed genes in meta-cluster ",c," with soft q-value thresholds:")
  DAVID_out[[c]]$sigFunctGroups
  
    message("Significantly enriched DAVID clusters for overexpressed genes in meta-cluster ",c,":")
   DAVID_out[[c]]$sigClusters
   # plot2D(DAVID_out[[c]]$functionalAnnotationChart[
    #  na.omit(match(DAVID_out[[c]]$sigFunctGroups,DAVID_out[[c]]$functionalAnnotationChart[,1]))],
     #    color=c("FALSE"="black", "TRUE"="green"))
}
```

```
## 
##  381  genes for commmunity number  1 (name)  2  have a positive effect size above .2 with a wilcoxon rank q-value below .1
```

```
## Warning in mapGeneIDs(geneticFeatures = testGeneVector, bioMartSource =
## "ensembl", : duplicate entries for a certain BiomaRt source may be
## returned.
```

```
## Warning in mapGeneIDs(geneticFeatures = testGeneVector, bioMartSource =
## "ensembl", : if you're feeding in HGNC gene symbols,make sure these are
## the latest updated version, otherwise Biomart may not recognize them.
```

```
## 
##  index 11  gene  CDKN2A  had multiple records returned.
## 
##  index 27  gene  GLDC  had multiple records returned.
## 
##  index 31  gene  APOA1  had multiple records returned.
## 
##  index 38  gene  COL3A1  had multiple records returned.
## 
##  index 44  gene  MMP11  had multiple records returned.
## 
##  index 51  gene  CFB  had multiple records returned.
## 
##  index 53  gene  SERPINA1  had multiple records returned.
## 
##  index 72  gene  PDGFRA  had multiple records returned.
## 
##  index 82  gene  TDO2  had multiple records returned.
## 
##  index 95  gene  CRIP1  had multiple records returned.
## 
##  index 100  gene  COL6A2  had multiple records returned.
## 
##  index 104  gene  COL6A3  had multiple records returned.
## 
##  index 108  gene  WT1  had multiple records returned.
## 
##  index 109  gene  C1S  had multiple records returned.
## 
##  index 110  gene  C1QB  had multiple records returned.
## 
##  index 116  gene  HLA-DPB1  had multiple records returned.
## 
##  index 117  gene  IFI27  had multiple records returned.
## 
##  index 119  gene  CCL5  had multiple records returned.
## 
##  index 127  gene  HLA-DPA1  had multiple records returned.
## 
##  index 136  gene  PMP22  had multiple records returned.
## 
##  index 137  gene  LY6E  had multiple records returned.
## 
##  index 138  gene  CBS  had multiple records returned.
## 
##  index 153  gene  COL1A2  had multiple records returned.
## 
##  index 158  gene  PHLDA2  had multiple records returned.
## 
##  index 163  gene  LPHN2  had multiple records returned.
## 
##  index 164  gene  DPYD  had multiple records returned.
## 
##  index 181  gene  TAP1  had multiple records returned.
## 
##  index 201  gene  IER3  had multiple records returned.
## 
##  index 204  gene  CXCR4  had multiple records returned.
## 
##  index 205  gene  CFI  had multiple records returned.
## 
##  index 215  gene  COL6A1  had multiple records returned.
## 
##  index 223  gene  FBLN5  had multiple records returned.
## 
##  index 228  gene  SOD2  had multiple records returned.
## 
##  index 236  gene  CD55  had multiple records returned.
## 
##  index 245  gene  ALOX5  had multiple records returned.
## 
##  index 252  gene  STAT1  had multiple records returned.
## 
##  index 254  gene  BTG3  had multiple records returned.
## 
##  index 291  gene  NDRG1  had multiple records returned.
## 
##  index 320  gene  DSC2  had multiple records returned.
## 
##  index 328  gene  SSPN  had multiple records returned.
## 
##  index 332  gene  PTPRC  had multiple records returned.
## 
##  index 334  gene  IL7R  had multiple records returned.
## 
##  index 339  gene  RUNX1  had multiple records returned.
## 
##  index 343  gene  IDH2  had multiple records returned.
## 
##  index 344  gene  EZH2  had multiple records returned.
## 
##  index 369  gene  IRF7  had multiple records returned.
## 
##  index 370  gene  CORO1A  had multiple records returned.
## 
##  index 377  gene  IL10RB  had multiple records returned.
```

```
## Connecting to DAVID web service. Occassionally it times out and you need to re-run this function.
## Significantly enriched DAVID functional groups for overexpressed genes in meta-cluster 1 with soft q-value thresholds:
## Significantly enriched DAVID clusters for overexpressed genes in meta-cluster 1:
```

```
## 
##  119  genes for commmunity number  2 (name)  3  have a positive effect size above .2 with a wilcoxon rank q-value below .1
```

```
## Warning in mapGeneIDs(geneticFeatures = testGeneVector, bioMartSource =
## "ensembl", : duplicate entries for a certain BiomaRt source may be
## returned.
```

```
## Warning in mapGeneIDs(geneticFeatures = testGeneVector, bioMartSource =
## "ensembl", : if you're feeding in HGNC gene symbols,make sure these are
## the latest updated version, otherwise Biomart may not recognize them.
```

```
## 
##  index 1  gene  C7  had multiple records returned.
## 
##  index 5  gene  SERPINA1  had multiple records returned.
## 
##  index 11  gene  PDGFRA  had multiple records returned.
## 
##  index 30  gene  PMP22  had multiple records returned.
## 
##  index 48  gene  IER3  had multiple records returned.
## 
##  index 61  gene  MYH11  had multiple records returned.
## 
##  index 94  gene  GATM  had multiple records returned.
## 
##  index 96  gene  RUNX1  had multiple records returned.
```

```
## Connecting to DAVID web service. Occassionally it times out and you need to re-run this function.
## Significantly enriched DAVID functional groups for overexpressed genes in meta-cluster 2 with soft q-value thresholds:
## Significantly enriched DAVID clusters for overexpressed genes in meta-cluster 2:
```

```
## 
##  66  genes for commmunity number  3 (name)  4  have a positive effect size above .2 with a wilcoxon rank q-value below .1
```

```
## Warning in mapGeneIDs(geneticFeatures = testGeneVector, bioMartSource =
## "ensembl", : duplicate entries for a certain BiomaRt source may be
## returned.
```

```
## Warning in mapGeneIDs(geneticFeatures = testGeneVector, bioMartSource =
## "ensembl", : if you're feeding in HGNC gene symbols,make sure these are
## the latest updated version, otherwise Biomart may not recognize them.
```

```
## 
##  index 3  gene  APOA1  had multiple records returned.
## 
##  index 5  gene  SCG5  had multiple records returned.
## 
##  index 23  gene  LPHN2  had multiple records returned.
## 
##  index 25  gene  CDKN1C  had multiple records returned.
## 
##  index 39  gene  SLC16A1  had multiple records returned.
## 
##  index 41  gene  GATM  had multiple records returned.
## 
##  index 45  gene  EZH2  had multiple records returned.
## 
##  index 66  gene  MPZ  had multiple records returned.
```

```
## Warning in mapGeneIDs(geneticFeatures = testGeneVector, bioMartSource =
## "ensembl", : Found more than 1 gene symbol for some input features -
## printed them out above. The first symbol was chosen/kept each time -
## usually the best bet.
```

```
## Connecting to DAVID web service. Occassionally it times out and you need to re-run this function.
## Significantly enriched DAVID functional groups for overexpressed genes in meta-cluster 3 with soft q-value thresholds:
## Significantly enriched DAVID clusters for overexpressed genes in meta-cluster 3:
```
**Survival Analyses**
I'm still working on determing which long-term and binary variables have the least amount of NAs across the samples in these meta-clusters, but it does look like the binary vital_status variable (alive or dead) and continous days to death variables provide survival curves that significantly stratify patients. It was less significant when I used a 5-year cutoff.


```r
library("survival")
```

```
## 
## Attaching package: 'survival'
## 
## The following object is masked from 'package:RDAVIDWebService':
## 
##     cluster
```

```r
#already have esets loaded:
 #load("/home/kplaney/ovarian_analysis/esets_proc_TCGAcombat.RData.gzip")
phenoMasterDF <- createPhenoMasterTableFromMatrixList(esetList=esets)
```

```
## This function assumes samples/patient clinical data rows are not duplicated
```

```r
#save(phenoMasterDF,file="/home/kplaney/ovarian_analysis/curatedOvarian_phenoMasterDF.RData.gzip",compress="gzip")
#load("/home/kplaney/ovarian_analysis/curatedOvarian_phenoMasterDF.RData.gzip")
#study numbers won't align here because some filtered out (had no robust clusters); need to "translate"

origToNewIndexMap <- data.frame(origToNewIndexMap,stringsAsFactors=FALSE)
colnames(origToNewIndexMap) <- c("studyNum","origStudyNum")
sampleClustCommKey <- join(sampleClustCommKey,origToNewIndexMap,by="studyNum",type="full",match="all")
sampleClustCommPhenoData <- addClinicalVarToNodeAttributes(sampleClustCommKey,phenoMasterDF=phenoMasterDF)

#survival analysis
outcomesVarBinary="vital_status"
outcomesVarCont = "days_to_death"
CutoffPointYears=5
uniquePatientID="unique_patient_ID"
groupingTerm="community"
   
     
  #only take samples with the groupingTerm you're looking at.
  sampleClustCommPhenoData <- sampleClustCommPhenoData[which(!is.na(sampleClustCommPhenoData[, groupingTerm])), ]
  #remove samples with NA values.
  groupings <- sampleClustCommPhenoData[, groupingTerm]
  
  #groupings <- groupings[which(!is.na(sampleClustCommPhenoData[,outcomesVarBinary]))];
  #outcomesData <- sampleClustCommPhenoData[which(!is.na(sampleClustCommPhenoData[,outcomesVarBinary])),];
  
  #keep samples with NA days to event for now?
#hmm...one meta-cluster is left out if use "days_to_death"...
  #groupings <- as.numeric(as.factor(groupings[which(!is.na(outcomesData[,outcomesVarCont]))]))
  #outcomesData <- outcomesData[which(!is.na(outcomesData[,outcomesVarCont])),];
  
  #if binary is character string categories: make it a factor first, then numeric,
  #otherwise coxph function will throw errors.
#  nonCensoredTerm=1
sampleClustCommPhenoData[which(sampleClustCommPhenoData[ ,outcomesVarBinary]=="deceased"),outcomesVarBinary] <- 1
sampleClustCommPhenoData[which(sampleClustCommPhenoData[ ,outcomesVarBinary]=="living"),outcomesVarBinary] <- 0
  outcomesDataShort <- data.frame(as.numeric(sampleClustCommPhenoData[,outcomesVarBinary]),as.numeric(sampleClustCommPhenoData[,outcomesVarCont])
  );
  
  #sometimes the names are duplicated across studies - remove this line
  #rownames(outcomesDataShort ) <- outcomesData[,uniquePatientID];
  colnames(outcomesDataShort) <- c("Censoring","TimeToLastContactOrEvent")
  
  nonCensoredTerm=1
  censoredTerm=0
  Survival <- outcomesDataShort
  #creating the survival objects with the time and censoring variables
  OverallSurvival <- Surv(Survival$TimeToLastContactOrEvent,Survival$Censoring==nonCensoredTerm);
  #creating a survival object cutoff at a certain point
   CutoffPoint <- CutoffPointYears*365;
   CutoffSamples=Survival$TimeToLastContactOrEvent>CutoffPoint & !is.na(Survival$TimeToLastContactOrEvent)
   SurvivalCutoff=Survival
   SurvivalCutoff$TimeToLastContactOrEvent[CutoffSamples]=CutoffPoint
   SurvivalCutoff$Censoring[CutoffSamples]=censoredTerm
#   #"Surv" creates a survival object. really for binary outcomes data.
   OverallSurvivalCutoff=Surv(SurvivalCutoff$TimeToLastContactOrEvent,SurvivalCutoff$Censoring==nonCensoredTerm)

    coxfit=coxph(OverallSurvival~groupings, data=Survival)
  message("coxfit summary for overall survival.")
```

```
## coxfit summary for overall survival.
```

```r
  summary(coxfit)
```

```
## Call:
## coxph(formula = OverallSurvival ~ groupings, data = Survival)
## 
##   n= 1551, number of events= 846 
##    (557 observations deleted due to missingness)
## 
##               coef exp(coef) se(coef)      z Pr(>|z|)
## groupings -0.01275   0.98733  0.04540 -0.281    0.779
## 
##           exp(coef) exp(-coef) lower .95 upper .95
## groupings    0.9873      1.013    0.9033     1.079
## 
## Concordance= 0.507  (se = 0.008 )
## Rsquare= 0   (max possible= 0.999 )
## Likelihood ratio test= 0.08  on 1 df,   p=0.7782
## Wald test            = 0.08  on 1 df,   p=0.7788
## Score (logrank) test = 0.08  on 1 df,   p=0.7788
```

```r
  #plot(cox.zph(coxfit))
 kmfit=survdiff(OverallSurvival ~ groupings)
message("kaplan meier p-value for overall survival:")
```

```
## kaplan meier p-value for overall survival:
```

```r
 1 - pchisq(kmfit$chisq, length(kmfit$n) - 1)
```

```
## [1] 0.1361693
```

```r
message("calculating the sign of the survival relationship")
```

```
## calculating the sign of the survival relationship
```

```r
  mfit=survfit(OverallSurvival ~ groupings)
  plot(mfit,main="overall survival")
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png) 

```r
mfitCut <- survfit(OverallSurvivalCutoff ~ groupings)
  plot(mfitCut,main=paste0("overall survival with cutoff at ",
                           CutoffPointYears," years"))
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-2.png) 

```r
coxfit=coxph(OverallSurvivalCutoff~groupings, data= SurvivalCutoff)
  message("coxfit summary for survival cutoff at ",CutoffPointYears,"years:")
```

```
## coxfit summary for survival cutoff at 5years:
```

```r
  summary(coxfit)
```

```
## Call:
## coxph(formula = OverallSurvivalCutoff ~ groupings, data = SurvivalCutoff)
## 
##   n= 1551, number of events= 746 
##    (557 observations deleted due to missingness)
## 
##               coef exp(coef) se(coef)      z Pr(>|z|)
## groupings -0.02433   0.97596  0.04886 -0.498    0.619
## 
##           exp(coef) exp(-coef) lower .95 upper .95
## groupings     0.976      1.025    0.8868     1.074
## 
## Concordance= 0.507  (se = 0.008 )
## Rsquare= 0   (max possible= 0.998 )
## Likelihood ratio test= 0.25  on 1 df,   p=0.6168
## Wald test            = 0.25  on 1 df,   p=0.6185
## Score (logrank) test = 0.25  on 1 df,   p=0.6185
```

```r
 #plot(cox.zph(coxfit))
  
  kmfit=survdiff(OverallSurvivalCutoff ~ groupings)
message("kaplan meier p-value for survival cutoff:")
```

```
## kaplan meier p-value for survival cutoff:
```

```r
  1 - pchisq(kmfit$chisq, length(kmfit$n) - 1)  
```

```
## [1] 0.6408672
```

```r
#message("a chi-square test looking at the binary recurrence status variable, as this data was recorded at least in some patients in all 4 meta-clusters:")

#message: "now with recurrence:"
#chisq.test(sampleClustCommPhenoData[,"recurrence_status"],groupings)
#message("We see a trend by just tabling the recurrence status variable too, but it's still pretty 50-50 (for the samples that didn't have NA values; it looks like a lot of samples were still missing this variable from certain meta-clusters")
#table(sampleClustCommPhenoData[,"recurrence_status"],groupings)


#interesting: two of the clusters are predominantly serous
#but careful: how counts NAs
table(sampleClustCommPhenoData[,groupingTerm],sampleClustCommPhenoData[,"histological_type"],exclude=NULL)
```

```
##       
##        clearcell endo  mix mucinous other  ser undifferentiated <NA>
##   2           12   39    0        3    30 1435                2   47
##   3           31   65    1       41    42   38                1   51
##   4            0   11    0        0     0  259                0    0
##   <NA>         0    0    0        0     0    0                0    0
```

```r
#mixed stages
table(sampleClustCommPhenoData[,groupingTerm],sampleClustCommPhenoData[,"tumorstage"],
      exclude=NULL)
```

```
##       
##           1    2    3    4 <NA>
##   2      63   59 1027  190  229
##   3      76   27   55   19   93
##   4      11    9  215   32    3
##   <NA>    0    0    0    0    0
```

```r
#eh: one group has all zeros
table(sampleClustCommPhenoData[,groupingTerm],sampleClustCommPhenoData[,"recurrence_status"],exclude=NULL)
```

```
##       
##        norecurrence recurrence <NA>
##   2             336        565  667
##   3              46         32  192
##   4              98        151   21
##   <NA>            0          0    0
```

```r
#not grouped by grade
table(sampleClustCommPhenoData[,groupingTerm],sampleClustCommPhenoData[,"grade"],exclude=NULL)
```

```
##       
##          1   2   3   4 <NA>
##   2     30 352 919   9  258
##   3     35  47  59   2  127
##   4      3  49 156   2   60
##   <NA>   0   0   0   0    0
```

```r
#is this only for TCGA?
table(sampleClustCommPhenoData[,groupingTerm],sampleClustCommPhenoData[,"primary_therapy_outcome_success"],exclude=NULL)
```

```
##       
##        completeresponse partialresponse progressivedisease stabledisease
##   2                 293              56                 46            24
##   3                   0               0                  0             0
##   4                  82               9                 15             6
##   <NA>                0               0                  0             0
##       
##        <NA>
##   2    1149
##   3     270
##   4     158
##   <NA>    0
```

```r
#communities spread fairly well through clusters
table(sampleClustCommPhenoData[,groupingTerm],sampleClustCommPhenoData[,"studyNum"],exclude=NULL)
```

```
##       
##          2   3   4   6   8   9  10  11  12  14  15  16  18  19  20  21  22
##   2     15 122  49  18  42  89  99  73  93 112  26  29  54 203  29  57   9
##   3      0   0  14   0   0  51  95  34   0   0   0  31  45   0   0   0   0
##   4      0   0   0   0   0   0   0   0  59   0   0   0   0  64   0  21   6
##   <NA>   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
##       
##         23 <NA>
##   2    449    0
##   3      0    0
##   4    120    0
##   <NA>   0    0
```

```r
#by only TCGA: it does appear that the cluster with the smaller TCGA samples may have worse response, but there's a lot of NAs:
table(sampleClustCommPhenoData[which(sampleClustCommPhenoData[,"studyNum"]==23),"primary_therapy_outcome_success"],sampleClustCommPhenoData[which(sampleClustCommPhenoData[,"studyNum"]==23),"community"],exclude=NULL)
```

```
##                     
##                        2   4 <NA>
##   completeresponse   250  67    0
##   partialresponse     56   9    0
##   progressivedisease  32   9    0
##   stabledisease       24   6    0
##   <NA>                87  29    0
```

```r
#all NAs for binary relapse variable in TCGA data:
table(sampleClustCommPhenoData[which(sampleClustCommPhenoData[,"studyNum"]==23),"relapse_binary"],sampleClustCommPhenoData[which(sampleClustCommPhenoData[,"studyNum"]==23),"community"],exclude=NULL)
```

```
##       
##          2   4 <NA>
##   <NA> 449 120    0
```

```r
table(sampleClustCommPhenoData[which(sampleClustCommPhenoData[,"studyNum"]==23),"relapse_binary"],sampleClustCommPhenoData[which(sampleClustCommPhenoData[,"studyNum"]==23),"community"],exclude=NULL)
```

```
##       
##          2   4 <NA>
##   <NA> 449 120    0
```

As for the days to recurrence available, it turns out there's only one grouping, so we can't run the analysis:


```r
# #re-compute this variable, as subsetted  before:
# sampleClustCommPhenoData <- addClinicalVarToNodeAttributes(sampleClustCommKey,phenoMasterDF=phenoMasterDF)
# 
# #survival analysis
# #booo- only one group has the days to recurrence.
# outcomesVarBinary="primary_therapy_outcome_success"
# outcomesVarCont = "days_to_tumor_recurrence"
# CutoffPointYears=5
# uniquePatientID="unique_patient_ID"
# groupingTerm="community"
#    
#      
#   #only take samples with the groupingTerm you're looking at.
#   sampleClustCommPhenoData <- sampleClustCommPhenoData[which(!is.na(sampleClustCommPhenoData[, groupingTerm])), ]
#   #remove samples with NA values.
#   groupings <- sampleClustCommPhenoData[, groupingTerm]
#   
#   #groupings <- groupings[which(!is.na(sampleClustCommPhenoData[,outcomesVarBinary]))];
#   #outcomesData <- sampleClustCommPhenoData[which(!is.na(sampleClustCommPhenoData[,outcomesVarBinary])),];
#   
#   #keep samples with NA days to event for now?
# #hmm...one meta-cluster is left out if use "days_to_death"...
#   #groupings <- as.numeric(as.factor(groupings[which(!is.na(outcomesData[,outcomesVarCont]))]))
#   #outcomesData <- outcomesData[which(!is.na(outcomesData[,outcomesVarCont])),];
#   
#   #if binary is character string categories: make it a factor first, then numeric,
#   #otherwise coxph function will throw errors.
# #  nonCensoredTerm=1
# sampleClustCommPhenoData[which(sampleClustCommPhenoData[ ,outcomesVarBinary]=="deceased"),outcomesVarBinary] <- 1
# sampleClustCommPhenoData[which(sampleClustCommPhenoData[ ,outcomesVarBinary]=="living"),outcomesVarBinary] <- 0
#   outcomesDataShort <- data.frame(as.numeric(sampleClustCommPhenoData[,outcomesVarBinary]),as.numeric(sampleClustCommPhenoData[,outcomesVarCont])
#   );
#   
#   #sometimes the names are duplicated across studies - remove this line
#   #rownames(outcomesDataShort ) <- outcomesData[,uniquePatientID];
#   colnames(outcomesDataShort) <- c("Censoring","TimeToLastContactOrEvent")
#   
#   nonCensoredTerm=1
#   censoredTerm=0
#   Survival <- outcomesDataShort
#   #creating the survival objects with the time and censoring variables
#   OverallSurvival <- Surv(Survival$TimeToLastContactOrEvent,Survival$Censoring==nonCensoredTerm);
#   #creating a survival object cutoff at a certain point
#    CutoffPoint <- CutoffPointYears*365;
#    CutoffSamples=Survival$TimeToLastContactOrEvent>CutoffPoint & !is.na(Survival$TimeToLastContactOrEvent)
#    SurvivalCutoff=Survival
#    SurvivalCutoff$TimeToLastContactOrEvent[CutoffSamples]=CutoffPoint
#    SurvivalCutoff$Censoring[CutoffSamples]=censoredTerm
# #   #"Surv" creates a survival object. really for binary outcomes data.
#    OverallSurvivalCutoff=Surv(SurvivalCutoff$TimeToLastContactOrEvent,SurvivalCutoff$Censoring==nonCensoredTerm)
# 
#     coxfit=coxph(OverallSurvival~groupings, data=Survival)
#   message("coxfit summary for overall survival.")
#   summary(coxfit)
#   #plot(cox.zph(coxfit))
#  kmfit=survdiff(OverallSurvival ~ groupings)
# message("kaplan meier p-value for overall survival:")
#  1 - pchisq(kmfit$chisq, length(kmfit$n) - 1)
# 
# message("calculating the sign of the survival relationship")
#   mfit=survfit(OverallSurvival ~ groupings)
#   plot(mfit,main="overall survival")
#   
# mfitCut <- survfit(OverallSurvivalCutoff ~ groupings)
#   plot(mfitCut,main=paste0("overall survival with cutoff at ",
#                            CutoffPointYears," years"))
# 
# coxfit=coxph(OverallSurvivalCutoff~groupings, data= SurvivalCutoff)
#   message("coxfit summary for survival cutoff at ",CutoffPointYears,"years:")
#   summary(coxfit)
#  #plot(cox.zph(coxfit))
#   
#   kmfit=survdiff(OverallSurvivalCutoff ~ groupings)
# message("kaplan meier p-value for survival cutoff:")
#   1 - pchisq(kmfit$chisq, length(kmfit$n) - 1)  
```
