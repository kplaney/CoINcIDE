curatedOvarianData meta-clustering with 200 (meta-rank) features using kmeans consensus
========================================================
**Data loading**


This is a report for the curatedOvarianData meta-clustering. Initial clustering was done using k-means and selecting k via the maximum mean cluster consensus (I found this gave stable results, but I am still exploring 1 other metric.) First, I loaded up my data just past the CoINcIDE" "getAdjMatrices() function (that's the function that takes  a really long time to run and is hogging up server space, as this is the part that computes the edge p-values):


```r
#grab data matrix list, clust features list
CoINcIDE_rankFeatures <- "/home/kplaney/ovarian_analysis/metaFeatures_200.RData.gzip"
load(CoINcIDE_rankFeatures)
metaFeatures <- metaFeatures
CoINcIDE_clusterOutput <-  "/home/kplaney/ovarian_analysis/curatedOvarianData_kmeansConsensus_200Features_2015-04-21.RData.gzip"
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

clustSampleIndexList <-  kmeansConsensus$clustSampleIndexList_meanConsensusCluster
clustFeatureIndexList <- kmeansConsensus$clustFeatureIndexList_meanConsensusCluster


CoINcIDE_computeEdgesObject <- "/home/kplaney/ovarian_analysis/kmeansConsensus_200F_meanMatrix_distCor.RData.gzip"
load(CoINcIDE_computeEdgesObject)
output <- kmeansConsensus_200F_meanMatrix_distCor
inputVariablesDF <- output$inputVariablesDF
computeTrueSimilOutput <- output$computeTrueSimilOutput
pvalueMatrix <- output$pvalueMatrix
clustIndexMatrix <- output$clustIndexMatrix

###inputs for edge detection
meanEdgePairPvalueThresh <- .05
indEdgePvalueThresh <- .1
minTrueSimilThresh <- .8
maxTrueSimilThresh <- Inf
clustSizeFractThresh <- inputVariablesDF$clustSizeFractThresh
clustSizeThresh <- inputVariablesDF$clustSizeThresh
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
## 1 2015-04-21 14:47:54    distCor                8               0.25
##   maxTrueSimilThresh  sigMethod maxNullFractSize numSims
## 1                Inf meanMatrix              0.2     500
##   includeRefClustInNull fractFeatIntersectThresh numFeatIntersectThresh
## 1                  TRUE                      0.8                    150
##   clustSizeThresh clustSizeFractThresh
## 1               5                 0.05
```

```r
message("There were ",nrow(clustIndexMatrix), " total input clusters from ",length(unique(clustIndexMatrix[,2])), " studies")
```

```
## There were 53 total input clusters from 24 studies
```

```r
message("The total number of input features was ",length(metaFeatures$finalFeatures))
```

```
## The total number of input features was 240
```

```r
message("Across the entire square (nonsymmetric) p-value matrix, there are ",length(which(pvalueMatrix<=.1))," pvalues less than or equal to .1")
```

```
## Across the entire square (nonsymmetric) p-value matrix, there are 765 pvalues less than or equal to .1
```

```r
message("Across the entire square (nonsymmetric) p-value matrix, there are ",length(which(pvalueMatrix<=.05))," pvalues less than or equal to .05")
```

```
## Across the entire square (nonsymmetric) p-value matrix, there are 678 pvalues less than or equal to .05
```

```r
message("Across the entire square (symmetric) similarity matrix, there are ",length(which(computeTrueSimilOutput$similValueMatrix>=inputVariablesDF$minTrueSimilThresh))," similarities greater than or equal to ",inputVariablesDF$minTrueSimilThresh)
```

```
## Across the entire square (symmetric) similarity matrix, there are 2588 similarities greater than or equal to 0.25
```

```r
message("Across the entire square (symmetric) similarity matrix, there are ",length(which(computeTrueSimilOutput$similValueMatrix>=.8))," similarities greater than or equal to ",.8)
```

```
## Across the entire square (symmetric) similarity matrix, there are 442 similarities greater than or equal to 0.8
```

**Network Analysis**

As we can see in these plots, 4 meta-clusters remained after edge filtering. 

```r
finalEdgeInfo <- assignFinalEdges(computeTrueSimilOutput=computeTrueSimilOutput,pvalueMatrix=pvalueMatrix,indEdgePvalueThresh=indEdgePvalueThresh,
                              meanEdgePairPvalueThresh=meanEdgePairPvalueThresh,
                              minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                              fractFeatIntersectThresh=fractFeatIntersectThresh,numFeatIntersectThresh=numFeatIntersectThresh ,
                              clustSizeThresh=clustSizeThresh, clustSizeFractThresh= clustSizeFractThresh,saveDir=saveDir,fileTag="CoINcIDE_edges_"
)
```

```
## 1 clusters dropped because they were below the clust size thresh clustSizeThresh
##           threshold of 5
```

```
## Warning in dir.create(saveDir): '/home/kplaney/ovarian_analysis' already
## exists
```

```
## A total of 21 clusters removed because they have no significant edges.
## A total of 32 clusters have significant edges.
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
# advancedNetworkPlots(communityMembership=commInfo,
#                                  brewPal = c("Set3"),
#                                  saveDir=saveDir,saveName="network",
#                      plotToScreen=TRUE)$network_stats
```

**Gene meta-rank Analysis**


I ranked genes within each meta-cluster for all samples, and then ran a Kruskal test to see which genes significantly stratified/differentiated patients across the 4 meta-clusters. I still need to implement GSEA; it turns out there's a base GSEA package in Biocondcutor so I've decided to just adapt my code and use their baseline functions.


In the heatmap: red means that gene was ranked high in terms of expression level for patients in that meta-cluster (I took the median rank across all samples in a meta-cluster to create the heatmap. Only significant genes are shown in the heatmap but at over 80 significant genes, of course the gene names are illegible...I print out the top 20 genes below.) The actual meta-cluster numbers are not 1:4 because meta-clusters with only 2 studies were thresholded out earlier.


```r
aggregateData <- returnSampleMemberMatrix(clustSampleIndexList,dataMatrixList,communityInfo=commInfo)
binInfo <- binarizeMetaclustStudyStatus(aggregateData$sampleClustCommKey)
rankInfo <- computeRankMatrix(metaClustSampleNames=binInfo$metaClustSampleNames,
                              featureNames=metaFeatures$finalFeatures,dataMatrixList,
                              sampleClustCommKey=aggregateData$sampleClustCommKey,onlyIntersectingFeat=TRUE)
pvalueInfo <- computeFeaturePvalues(rankMatrix=rankInfo$rankMatrix,featureNames=rankInfo$filteredFeatures,groupings=rankInfo$groupings)
message("There are ",length(which(pvalueInfo$fdr.qvalue<=.05)), " genes with an FDR corrected p-value below .05 whose ranks significantly differed among all of the meta-clusters (using Kruskal's test)")
```

```
## There are 88 genes with an FDR corrected p-value below .05 whose ranks significantly differed among all of the meta-clusters (using Kruskal's test)
```

```r
message("Top 20 significant genes:\n")
```

```
## Top 20 significant genes:
```

```r
rownames(pvalueInfo[which(pvalueInfo$fdr.qvalue<=.05)[1:20],])
```

```
##  [1] "COL11A1" "MMP7"    "DEFB1"   "C7"      "MAL"     "LUM"     "SST"    
##  [8] "NNMT"    "VCAN"    "MFAP5"   "INHBA"   "CDKN2A"  "CXCL10"  "FOS"    
## [15] "KLK10"   "CHI3L1"  "TFAP2A"  "GPX3"    "RARRES1" "TAGLN"
```

```r
cat("\n")
```

```r
metaMatrix <- commMedianRank(rankInfo$rankMatrix[which(pvalueInfo$fdr.qvalue<=.05),],rankInfo$groupings)
message("Red in heatmap means genes were ranked higher across all samples in that meta-cluster.")
```

```
## Red in heatmap means genes were ranked higher across all samples in that meta-cluster.
```

```r
plotMetaAnalysis(metaMatrix,saveFile=FALSE,plotToScreen=TRUE,
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
#save variable
sampleClustCommKey<-aggregateData$sampleClustCommKey
```

**GSEA Analysis**

It turns out "GSEABase" isn't all that great..so just using some Broad gene sets I already downloaded and the standard hypergeometric tests:


```r
GSEA_out <- GSEA(testGeneVector=rownames(metaMatrix),method=c("hypergeometric","fisher"),genomeSize=20000,
                  refGeneListDir="/home/data/MSigDB/GSEA_base_MSigDB_lists_merged.RData.gzip")
```

```
## Warning in GSEA(testGeneVector = rownames(metaMatrix), method = c("hypergeometric", : 
## This code assumes that all of your genes in your test gene list and ref gene list are in the genome.
```

```
## 
## Using default MSigDB lists: MSigDB_onco_symbols, MSigDB_CanPath_symbols,MSigDB_TFT_symbols,MSigDB_immun_symbols,
##        and MSigDB_cancerNeigh_symbols.
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (method == "hypergeometric") {: the condition has length > 1
## and only the first element will be used
```

```r
message("Gene lists enriched with q-value below .05:")
```

```
## Gene lists enriched with q-value below .05:
```

```r
load("/home/data/MSigDB/GSEA_base_MSigDB_lists_merged.RData.gzip")
names(GSEA_base_MSigDB_lists_merged)[which(GSEA_out$qvalues<=.05)]
```

```
##   [1] "E2F1_UP.V1_DN"                                                           
##   [2] "EGFR_UP.V1_UP"                                                           
##   [3] "HINATA_NFKB_MATRIX"                                                      
##   [4] "VEGF_A_UP.V1_UP"                                                         
##   [5] "ATF2_S_UP.V1_DN"                                                         
##   [6] "WNT_UP.V1_UP"                                                            
##   [7] "RELA_DN.V1_DN"                                                           
##   [8] "P53_DN.V1_DN"                                                            
##   [9] "P53_DN.V1_UP"                                                            
##  [10] "SNF5_DN.V1_DN"                                                           
##  [11] "LTE2_UP.V1_DN"                                                           
##  [12] "MEK_UP.V1_DN"                                                            
##  [13] "RAF_UP.V1_DN"                                                            
##  [14] "PRC1_BMI_UP.V1_DN"                                                       
##  [15] "ESC_J1_UP_LATE.V1_UP"                                                    
##  [16] "ESC_V6.5_UP_EARLY.V1_DN"                                                 
##  [17] "ESC_V6.5_UP_LATE.V1_DN"                                                  
##  [18] "BMI1_DN_MEL18_DN.V1_DN"                                                  
##  [19] "BMI1_DN_MEL18_DN.V1_UP"                                                  
##  [20] "BMI1_DN.V1_DN"                                                           
##  [21] "BMI1_DN.V1_UP"                                                           
##  [22] "MEL18_DN.V1_UP"                                                          
##  [23] "PTEN_DN.V1_DN"                                                           
##  [24] "RB_P107_DN.V1_UP"                                                        
##  [25] "CAHOY_ASTROGLIAL"                                                        
##  [26] "PDGF_UP.V1_UP"                                                           
##  [27] "HOXA9_DN.V1_UP"                                                          
##  [28] "SINGH_KRAS_DEPENDENCY_SIGNATURE_"                                        
##  [29] "KRAS.DF.V1_UP"                                                           
##  [30] "KRAS.LUNG_UP.V1_DN"                                                      
##  [31] "KRAS.PROSTATE_UP.V1_DN"                                                  
##  [32] "LEF1_UP.V1_UP"                                                           
##  [33] "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION"                             
##  [34] "KEGG_P53_SIGNALING_PATHWAY"                                              
##  [35] "KEGG_FOCAL_ADHESION"                                                     
##  [36] "KEGG_ECM_RECEPTOR_INTERACTION"                                           
##  [37] "KEGG_CELL_ADHESION_MOLECULES_CAMS"                                       
##  [38] "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES"                                
##  [39] "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY"                               
##  [40] "KEGG_PRION_DISEASES"                                                     
##  [41] "KEGG_PATHWAYS_IN_CANCER"                                                 
##  [42] "BIOCARTA_LAIR_PATHWAY"                                                   
##  [43] "BIOCARTA_CLASSIC_PATHWAY"                                                
##  [44] "BIOCARTA_COMP_PATHWAY"                                                   
##  [45] "BIOCARTA_FIBRINOLYSIS_PATHWAY"                                           
##  [46] "BIOCARTA_IL6_PATHWAY"                                                    
##  [47] "BIOCARTA_IL10_PATHWAY"                                                   
##  [48] "BIOCARTA_PLATELETAPP_PATHWAY"                                            
##  [49] "BIOCARTA_HER2_PATHWAY"                                                   
##  [50] "SA_REG_CASCADE_OF_CYCLIN_EXPR"                                           
##  [51] "PID_SMAD2_3NUCLEARPATHWAY"                                               
##  [52] "PID_INTEGRIN1_PATHWAY"                                                   
##  [53] "PID_E2F_PATHWAY"                                                         
##  [54] "PID_INTEGRIN3_PATHWAY"                                                   
##  [55] "PID_FRA_PATHWAY"                                                         
##  [56] "PID_P53DOWNSTREAMPATHWAY"                                                
##  [57] "PID_AVB3_INTEGRIN_PATHWAY"                                               
##  [58] "PID_ATF2_PATHWAY"                                                        
##  [59] "PID_AP1_PATHWAY"                                                         
##  [60] "PID_UPA_UPAR_PATHWAY"                                                    
##  [61] "PID_FOXM1PATHWAY"                                                        
##  [62] "PID_SYNDECAN_1_PATHWAY"                                                  
##  [63] "PID_ERBB_NETWORK_PATHWAY"                                                
##  [64] "PID_INTEGRIN5_PATHWAY"                                                   
##  [65] "PID_HIF1_TFPATHWAY"                                                      
##  [66] "PID_TOLL_ENDOGENOUS_PATHWAY"                                             
##  [67] "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION"                              
##  [68] "REACTOME_COLLAGEN_FORMATION"                                             
##  [69] "REACTOME_CS_DS_DEGRADATION"                                              
##  [70] "REACTOME_CHONDROITIN_SULFATE_BIOSYNTHESIS"                               
##  [71] "REACTOME_CHONDROITIN_SULFATE_DERMATAN_SULFATE_METABOLISM"                
##  [72] "REACTOME_HEPARAN_SULFATE_HEPARIN_HS_GAG_METABOLISM"                      
##  [73] "REACTOME_GLYCOSAMINOGLYCAN_METABOLISM"                                   
##  [74] "REACTOME_A_TETRASACCHARIDE_LINKER_SEQUENCE_IS_REQUIRED_FOR_GAG_SYNTHESIS"
##  [75] "REACTOME_PEPTIDE_LIGAND_BINDING_RECEPTORS"                               
##  [76] "REACTOME_CLASS_A1_RHODOPSIN_LIKE_RECEPTORS"                              
##  [77] "REACTOME_CHEMOKINE_RECEPTORS_BIND_CHEMOKINES"                            
##  [78] "REACTOME_G_ALPHA_I_SIGNALLING_EVENTS"                                    
##  [79] "REACTOME_INTERFERON_GAMMA_SIGNALING"                                     
##  [80] "AAANWWTGC_UNKNOWN"                                                       
##  [81] "V$IRF1_01"                                                               
##  [82] "V$SRF_Q6"                                                                
##  [83] "V$GR_Q6"                                                                 
##  [84] "V$SRF_C"                                                                 
##  [85] "CATRRAGC_UNKNOWN"                                                        
##  [86] "V$SPZ1_01"                                                               
##  [87] "V$STAT5A_02"                                                             
##  [88] "V$PAX2_02"                                                               
##  [89] "V$CEBPDELTA_Q6"                                                          
##  [90] "V$TFIIA_Q6"                                                              
##  [91] "V$CACCCBINDINGFACTOR_Q6"                                                 
##  [92] "V$SRF_Q4"                                                                
##  [93] "V$GR_Q6_01"                                                              
##  [94] "V$SRF_Q5_01"                                                             
##  [95] "TAATTA_V$CHX10_01"                                                       
##  [96] "TTGTTT_V$FOXO4_01"                                                       
##  [97] "TGATTTRY_V$GFI1_01"                                                      
##  [98] "STTTCRNTTT_V$IRF_Q6"                                                     
##  [99] "YTATTTTNR_V$MEF2_02"                                                     
## [100] "TGGAAA_V$NFAT_Q4_01"                                                     
## [101] "TGASTMAGC_V$NFE2_01"                                                     
## [102] "TATAAA_V$TATA_01"                                                        
## [103] "GSE10325_MYELOID_VS_LUPUS_MYELOID_DN"                                    
## [104] "GSE13306_RA_VS_UNTREATED_MEM_CD4_TCELL_UP"                               
## [105] "GSE13484_UNSTIM_VS_YF17D_VACCINE_STIM_PBMC_DN"                           
## [106] "GSE13485_CTRL_VS_DAY3_YF17D_VACCINE_PBMC_UP"                             
## [107] "GSE13485_CTRL_VS_DAY3_YF17D_VACCINE_PBMC_DN"                             
## [108] "GSE13485_CTRL_VS_DAY7_YF17D_VACCINE_PBMC_DN"                             
## [109] "GSE13485_DAY3_VS_DAY7_YF17D_VACCINE_PBMC_DN"                             
## [110] "GSE13485_PRE_VS_POST_YF17D_VACCINATION_PBMC_DN"                          
## [111] "GSE14000_UNSTIM_VS_4H_LPS_DC_TRANSLATED_RNA_DN"                          
## [112] "GSE14000_UNSTIM_VS_16H_LPS_DC_TRANSLATED_RNA_DN"                         
## [113] "GSE14000_UNSTIM_VS_4H_LPS_DC_DN"                                         
## [114] "GSE1432_CTRL_VS_IFNG_6H_MICROGLIA_DN"                                    
## [115] "GSE1432_CTRL_VS_IFNG_24H_MICROGLIA_DN"                                   
## [116] "GSE14769_UNSTIM_VS_80MIN_LPS_BMDM_DN"                                    
## [117] "GSE15930_NAIVE_VS_48H_IN_VITRO_STIM_CD8_TCELL_UP"                        
## [118] "GSE16755_CTRL_VS_IFNA_TREATED_MAC_DN"                                    
## [119] "GSE17721_CTRL_VS_GARDIQUIMOD_12H_BMDM_DN"                                
## [120] "GSE17721_CPG_VS_GARDIQUIMOD_6H_BMDM_UP"                                  
## [121] "GSE17721_12H_VS_24H_LPS_BMDM_UP"                                         
## [122] "GSE17974_IL4_AND_ANTI_IL12_VS_UNTREATED_48H_ACT_CD4_TCELL_UP"            
## [123] "GSE18791_CTRL_VS_NEWCASTLE_VIRUS_DC_6H_DN"                               
## [124] "GSE18791_CTRL_VS_NEWCASTLE_VIRUS_DC_8H_DN"                               
## [125] "GSE18791_UNSTIM_VS_NEWCATSLE_VIRUS_DC_6H_DN"                             
## [126] "GSE18791_UNSTIM_VS_NEWCATSLE_VIRUS_DC_10H_DN"                            
## [127] "GSE20715_0H_VS_6H_OZONE_LUNG_DN"                                         
## [128] "GSE22886_DC_VS_MONOCYTE_UP"                                              
## [129] "GSE22886_CTRL_VS_LPS_24H_DC_DN"                                          
## [130] "GSE24634_TREG_VS_TCONV_POST_DAY7_IL4_CONVERSION_DN"                      
## [131] "GSE24634_TREG_VS_TCONV_POST_DAY10_IL4_CONVERSION_DN"                     
## [132] "GSE24634_IL4_VS_CTRL_TREATED_NAIVE_CD4_TCELL_DAY3_DN"                    
## [133] "GSE24634_IL4_VS_CTRL_TREATED_NAIVE_CD4_TCELL_DAY5_DN"                    
## [134] "GSE2706_UNSTIM_VS_2H_R848_DC_DN"                                         
## [135] "GSE2706_UNSTIM_VS_8H_R848_DC_DN"                                         
## [136] "GSE2706_UNSTIM_VS_2H_LPS_DC_DN"                                          
## [137] "GSE2706_UNSTIM_VS_2H_LPS_AND_R848_DC_DN"                                 
## [138] "GSE2706_UNSTIM_VS_8H_LPS_AND_R848_DC_DN"                                 
## [139] "GSE2706_R848_VS_LPS_2H_STIM_DC_DN"                                       
## [140] "GSE2706_R848_VS_R848_AND_LPS_2H_STIM_DC_DN"                              
## [141] "GSE2706_2H_VS_8H_R848_AND_LPS_STIM_DC_DN"                                
## [142] "GSE2826_XID_VS_BTK_KO_BCELL_DN"                                          
## [143] "GSE29614_DAY3_VS_DAY7_TIV_FLU_VACCINE_PBMC_UP"                           
## [144] "GSE29615_CTRL_VS_LAIV_FLU_VACCINE_PBMC_UP"                               
## [145] "GSE29617_DAY3_VS_DAY7_TIV_FLU_VACCINE_PBMC_2008_UP"                      
## [146] "GSE29618_PRE_VS_DAY7_FLU_VACCINE_MONOCYTE_UP"                            
## [147] "GSE29618_LAIV_VS_TIV_FLU_VACCINE_DAY7_PDC_DN"                            
## [148] "GSE3337_4H_VS_16H_IFNG_IN_CD8POS_DC_DN"                                  
## [149] "GSE360_CTRL_VS_L_MAJOR_DC_DN"                                            
## [150] "GSE360_CTRL_VS_T_GONDII_DC_DN"                                           
## [151] "GSE360_CTRL_VS_B_MALAYI_HIGH_DOSE_DC_DN"                                 
## [152] "GSE360_CTRL_VS_M_TUBERCULOSIS_DC_DN"                                     
## [153] "GSE360_DC_VS_MAC_T_GONDII_DN"                                            
## [154] "GSE360_DC_VS_MAC_B_MALAYI_HIGH_DOSE_DN"                                  
## [155] "GSE360_DC_VS_MAC_M_TUBERCULOSIS_UP"                                      
## [156] "GSE360_DC_VS_MAC_M_TUBERCULOSIS_DN"                                      
## [157] "GSE360_L_DONOVANI_VS_B_MALAYI_HIGH_DOSE_DC_UP"                           
## [158] "GSE360_L_DONOVANI_VS_B_MALAYI_LOW_DOSE_DC_UP"                            
## [159] "GSE360_L_DONOVANI_VS_M_TUBERCULOSIS_DC_DN"                               
## [160] "GSE360_L_MAJOR_VS_T_GONDII_DC_DN"                                        
## [161] "GSE360_L_MAJOR_VS_B_MALAYI_LOW_DOSE_DC_UP"                               
## [162] "GSE360_T_GONDII_VS_B_MALAYI_LOW_DOSE_DC_DN"                              
## [163] "GSE360_T_GONDII_VS_M_TUBERCULOSIS_DC_DN"                                 
## [164] "GSE360_HIGH_VS_LOW_DOSE_B_MALAYI_DC_DN"                                  
## [165] "GSE360_HIGH_DOSE_B_MALAYI_VS_M_TUBERCULOSIS_DC_DN"                       
## [166] "GSE360_LOW_DOSE_B_MALAYI_VS_M_TUBERCULOSIS_DC_DN"                        
## [167] "GSE360_L_DONOVANI_VS_T_GONDII_MAC_UP"                                    
## [168] "GSE360_L_DONOVANI_VS_B_MALAYI_HIGH_DOSE_MAC_UP"                          
## [169] "GSE360_L_DONOVANI_VS_B_MALAYI_HIGH_DOSE_MAC_DN"                          
## [170] "GSE360_L_DONOVANI_VS_M_TUBERCULOSIS_MAC_DN"                              
## [171] "GSE360_T_GONDII_VS_M_TUBERCULOSIS_MAC_DN"                                
## [172] "GSE360_HIGH_DOSE_B_MALAYI_VS_M_TUBERCULOSIS_MAC_DN"                      
## [173] "GSE36476_YOUNG_VS_OLD_DONOR_MEMORY_CD4_TCELL_40H_TSST_ACT_UP"            
## [174] "GSE3982_CTRL_VS_LPS_48H_DC_DN"                                           
## [175] "GSE3982_BCELL_VS_CENT_MEMORY_CD4_TCELL_DN"                               
## [176] "GSE6269_HEALTHY_VS_STREP_PNEUMO_INF_PBMC_DN"                             
## [177] "GSE6269_FLU_VS_STREP_PNEUMO_INF_PBMC_UP"                                 
## [178] "GSE7764_IL15_TREATED_VS_CTRL_NK_CELL_24H_DN"                             
## [179] "GSE9006_HEALTHY_VS_TYPE_1_DIABETES_PBMC_AT_DX_DN"                        
## [180] "GSE9006_1MONTH_VS_4MONTH_POST_TYPE_1_DIABETES_DX_PBMC_UP"                
## [181] "GSE9988_ANTI_TREM1_VS_LPS_MONOCYTE_DN"                                   
## [182] "GSE9988_ANTI_TREM1_VS_LOW_LPS_MONOCYTE_DN"                               
## [183] "GSE9988_LPS_VS_CTRL_TREATED_MONOCYTE_UP"                                 
## [184] "GSE9988_LPS_VS_VEHICLE_TREATED_MONOCYTE_UP"                              
## [185] "GSE9988_LOW_LPS_VS_CTRL_TREATED_MONOCYTE_UP"                             
## [186] "GNF2_CD33"                                                               
## [187] "GNF2_CDH11"                                                              
## [188] "GNF2_FOS"                                                                
## [189] "GNF2_MMP1"                                                               
## [190] "GNF2_PTX3"
```
**Survival Analyses**
I'm still working on determing which long-term and binary variables have the least amount of NAs across the samples in these meta-clusters, but it does look like the binary vital_status variable (alive or dead) and continous days to death variables provide survival curves that significantly stratify patients. It was less significant when I used a 5-year cutoff.

 
 ```r
 library("survival")
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
 ##   n= 1324, number of events= 717 
 ##    (517 observations deleted due to missingness)
 ## 
 ##               coef exp(coef) se(coef)      z Pr(>|z|)  
 ## groupings -0.14661   0.86363  0.05788 -2.533   0.0113 *
 ## ---
 ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
 ## 
 ##           exp(coef) exp(-coef) lower .95 upper .95
 ## groupings    0.8636      1.158     0.771    0.9674
 ## 
 ## Concordance= 0.527  (se = 0.011 )
 ## Rsquare= 0.005   (max possible= 0.999 )
 ## Likelihood ratio test= 6.66  on 1 df,   p=0.009844
 ## Wald test            = 6.42  on 1 df,   p=0.01131
 ## Score (logrank) test = 6.41  on 1 df,   p=0.01133
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
 ## [1] 0.0006151971
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
 ##   n= 1324, number of events= 622 
 ##    (517 observations deleted due to missingness)
 ## 
 ##               coef exp(coef) se(coef)      z Pr(>|z|)
 ## groupings -0.09497   0.90940  0.06268 -1.515     0.13
 ## 
 ##           exp(coef) exp(-coef) lower .95 upper .95
 ## groupings    0.9094        1.1    0.8043     1.028
 ## 
 ## Concordance= 0.526  (se = 0.011 )
 ## Rsquare= 0.002   (max possible= 0.998 )
 ## Likelihood ratio test= 2.35  on 1 df,   p=0.1253
 ## Wald test            = 2.3  on 1 df,   p=0.1297
 ## Score (logrank) test = 2.29  on 1 df,   p=0.1298
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
 ## [1] 0.0004601185
 ```
 
 ```r
 message("a chi-square test looking at the binary recurrence status variable, as this data was recorded at least in some patients in all 4 meta-clusters:")
 ```
 
 ```
 ## a chi-square test looking at the binary recurrence status variable, as this data was recorded at least in some patients in all 4 meta-clusters:
 ```
 
 ```r
 chisq.test(sampleClustCommPhenoData[,"recurrence_status"],groupings)
 ```
 
 ```
 ## 
 ## 	Pearson's Chi-squared test
 ## 
 ## data:  sampleClustCommPhenoData[, "recurrence_status"] and groupings
 ## X-squared = 19.6412, df = 3, p-value = 0.0002014
 ```
 
 ```r
 message("We see a trend by just tabling the recurrence status variable too, but it's still pretty 50-50 (for the samples that didn't have NA values; it looks like a lot of samples were still missing this variable from certain meta-clusters")
 ```
 
 ```
 ## We see a trend by just tabling the recurrence status variable too, but it's still pretty 50-50 (for the samples that didn't have NA values; it looks like a lot of samples were still missing this variable from certain meta-clusters
 ```
 
 ```r
 table(sampleClustCommPhenoData[,"recurrence_status"],groupings)
 ```
 
 ```
 ##               groupings
 ##                  3   4   5   6
 ##   norecurrence 172 225  16  38
 ##   recurrence   341 365  16  25
 ```
