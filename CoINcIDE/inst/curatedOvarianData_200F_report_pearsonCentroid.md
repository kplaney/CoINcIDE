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
#source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_metaClusterAnalysis.R")

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


CoINcIDE_computeEdgesObject <- "/home/kplaney/ovarian_analysis/kmeansConsensus_200F_centroid_pearson.RData.gzip"
load(CoINcIDE_computeEdgesObject)
output <- kmeansConsensus_200F_centroid_pearson
inputVariablesDF <- output$inputVariablesDF
computeTrueSimilOutput <- output$computeTrueSimilOutput
pvalueMatrix <- output$pvalueMatrix
clustIndexMatrix <- output$clustIndexMatrix

###inputs for edge detection
meanEdgePairPvalueThresh <- .2
indEdgePvalueThresh <- .3
minTrueSimilThresh <- .25
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
## 1 2015-04-21 14:42:30    pearson                8               0.25
##   maxTrueSimilThresh sigMethod maxNullFractSize numSims
## 1                Inf  centroid              0.2     500
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
## Across the entire square (nonsymmetric) p-value matrix, there are 240 pvalues less than or equal to .1
```

```r
message("Across the entire square (nonsymmetric) p-value matrix, there are ",length(which(pvalueMatrix<=.05))," pvalues less than or equal to .05")
```

```
## Across the entire square (nonsymmetric) p-value matrix, there are 189 pvalues less than or equal to .05
```

```r
message("Across the entire square (symmetric) similarity matrix, there are ",length(which(computeTrueSimilOutput$similValueMatrix>=inputVariablesDF$minTrueSimilThresh))," similarities greater than or equal to ",inputVariablesDF$minTrueSimilThresh)
```

```
## Across the entire square (symmetric) similarity matrix, there are 2142 similarities greater than or equal to 0.25
```

```r
message("Across the entire square (symmetric) similarity matrix, there are ",length(which(computeTrueSimilOutput$similValueMatrix>=.8))," similarities greater than or equal to ",.8)
```

```
## Across the entire square (symmetric) similarity matrix, there are 0 similarities greater than or equal to 0.8
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
## A total of 33 clusters removed because they have no significant edges.
## A total of 20 clusters have significant edges.
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
```

```
## Error in eval(expr, envir, enclos): could not find function "binarizeMetaclustStudyStatus"
```

```r
rankInfo <- computeRankMatrix(metaClustSampleNames=binInfo$metaClustSampleNames,
                              featureNames=metaFeatures$finalFeatures,dataMatrixList,
                              sampleClustCommKey=aggregateData$sampleClustCommKey,onlyIntersectingFeat=TRUE)
```

```
## Error in eval(expr, envir, enclos): could not find function "computeRankMatrix"
```

```r
pvalueInfo <- computeFeaturePvalues(rankMatrix=rankInfo$rankMatrix,featureNames=rankInfo$filteredFeatures,groupings=rankInfo$groupings)
```

```
## Error in eval(expr, envir, enclos): could not find function "computeFeaturePvalues"
```

```r
message("There are ",length(which(pvalueInfo$fdr.qvalue<=.05)), " genes with an FDR corrected p-value below .05 whose ranks significantly differed among all of the meta-clusters (using Kruskal's test)")
```

```
## Error in which(pvalueInfo$fdr.qvalue <= 0.05): object 'pvalueInfo' not found
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
## Error in rownames(pvalueInfo[which(pvalueInfo$fdr.qvalue <= 0.05)[1:20], : error in evaluating the argument 'x' in selecting a method for function 'rownames': Error: object 'pvalueInfo' not found
```

```r
cat("\n")
```

```r
metaMatrix <- commMedianRank(rankInfo$rankMatrix[which(pvalueInfo$fdr.qvalue<=.05),],rankInfo$groupings)
```

```
## Error in eval(expr, envir, enclos): could not find function "commMedianRank"
```

```r
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
## Error in eval(expr, envir, enclos): could not find function "plotMetaAnalysis"
```

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
## Error in eval(expr, envir, enclos): could not find function "GSEA"
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
## Error in which(GSEA_out$qvalues <= 0.05): object 'GSEA_out' not found
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
 ## Error in eval(expr, envir, enclos): could not find function "createPhenoMasterTableFromMatrixList"
 ```
 
 ```r
 #save(phenoMasterDF,file="/home/kplaney/ovarian_analysis/curatedOvarian_phenoMasterDF.RData.gzip",compress="gzip")
 #load("/home/kplaney/ovarian_analysis/curatedOvarian_phenoMasterDF.RData.gzip")
 #study numbers won't align here because some filtered out (had no robust clusters); need to "translate"
 
 origToNewIndexMap <- data.frame(origToNewIndexMap,stringsAsFactors=FALSE)
 colnames(origToNewIndexMap) <- c("studyNum","origStudyNum")
 sampleClustCommKey <- join(sampleClustCommKey,origToNewIndexMap,by="studyNum",type="full",match="all")
 ```
 
 ```
 ## Error in eval(expr, envir, enclos): could not find function "join"
 ```
 
 ```r
 sampleClustCommPhenoData <- addClinicalVarToNodeAttributes(sampleClustCommKey,phenoMasterDF=phenoMasterDF)
 ```
 
 ```
 ## Error in eval(expr, envir, enclos): could not find function "addClinicalVarToNodeAttributes"
 ```
 
 ```r
 #survival analysis
 outcomesVarBinary="vital_status"
 outcomesVarCont = "days_to_death"
 CutoffPointYears=5
 uniquePatientID="unique_patient_ID"
 groupingTerm="community"
   
     
  #only take samples with the groupingTerm you're looking at.
  sampleClustCommPhenoData <- sampleClustCommPhenoData[which(!is.na(sampleClustCommPhenoData[, groupingTerm])), ]
 ```
 
 ```
 ## Error in eval(expr, envir, enclos): object 'sampleClustCommPhenoData' not found
 ```
 
 ```r
  #remove samples with NA values.
  groupings <- sampleClustCommPhenoData[, groupingTerm]
 ```
 
 ```
 ## Error in eval(expr, envir, enclos): object 'sampleClustCommPhenoData' not found
 ```
 
 ```r
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
 ```
 
 ```
 ## Error in sampleClustCommPhenoData[which(sampleClustCommPhenoData[, outcomesVarBinary] == : object 'sampleClustCommPhenoData' not found
 ```
 
 ```r
 sampleClustCommPhenoData[which(sampleClustCommPhenoData[ ,outcomesVarBinary]=="living"),outcomesVarBinary] <- 0
 ```
 
 ```
 ## Error in sampleClustCommPhenoData[which(sampleClustCommPhenoData[, outcomesVarBinary] == : object 'sampleClustCommPhenoData' not found
 ```
 
 ```r
  outcomesDataShort <- data.frame(as.numeric(sampleClustCommPhenoData[,outcomesVarBinary]),as.numeric(sampleClustCommPhenoData[,outcomesVarCont])
  );
 ```
 
 ```
 ## Error in data.frame(as.numeric(sampleClustCommPhenoData[, outcomesVarBinary]), : object 'sampleClustCommPhenoData' not found
 ```
 
 ```r
  #sometimes the names are duplicated across studies - remove this line
  #rownames(outcomesDataShort ) <- outcomesData[,uniquePatientID];
  colnames(outcomesDataShort) <- c("Censoring","TimeToLastContactOrEvent")
 ```
 
 ```
 ## Error in colnames(outcomesDataShort) <- c("Censoring", "TimeToLastContactOrEvent"): object 'outcomesDataShort' not found
 ```
 
 ```r
  nonCensoredTerm=1
  censoredTerm=0
  Survival <- outcomesDataShort
 ```
 
 ```
 ## Error in eval(expr, envir, enclos): object 'outcomesDataShort' not found
 ```
 
 ```r
  #creating the survival objects with the time and censoring variables
  OverallSurvival <- Surv(Survival$TimeToLastContactOrEvent,Survival$Censoring==nonCensoredTerm);
 ```
 
 ```
 ## Error in Surv(Survival$TimeToLastContactOrEvent, Survival$Censoring == : object 'Survival' not found
 ```
 
 ```r
  #creating a survival object cutoff at a certain point
   CutoffPoint <- CutoffPointYears*365;
   CutoffSamples=Survival$TimeToLastContactOrEvent>CutoffPoint & !is.na(Survival$TimeToLastContactOrEvent)
 ```
 
 ```
 ## Error in eval(expr, envir, enclos): object 'Survival' not found
 ```
 
 ```r
   SurvivalCutoff=Survival
 ```
 
 ```
 ## Error in eval(expr, envir, enclos): object 'Survival' not found
 ```
 
 ```r
   SurvivalCutoff$TimeToLastContactOrEvent[CutoffSamples]=CutoffPoint
 ```
 
 ```
 ## Error in SurvivalCutoff$TimeToLastContactOrEvent[CutoffSamples] = CutoffPoint: object 'SurvivalCutoff' not found
 ```
 
 ```r
   SurvivalCutoff$Censoring[CutoffSamples]=censoredTerm
 ```
 
 ```
 ## Error in SurvivalCutoff$Censoring[CutoffSamples] = censoredTerm: object 'SurvivalCutoff' not found
 ```
 
 ```r
 #   #"Surv" creates a survival object. really for binary outcomes data.
   OverallSurvivalCutoff=Surv(SurvivalCutoff$TimeToLastContactOrEvent,SurvivalCutoff$Censoring==nonCensoredTerm)
 ```
 
 ```
 ## Error in Surv(SurvivalCutoff$TimeToLastContactOrEvent, SurvivalCutoff$Censoring == : object 'SurvivalCutoff' not found
 ```
 
 ```r
    coxfit=coxph(OverallSurvival~groupings, data=Survival)
 ```
 
 ```
 ## Error in terms.formula(formula, special, data = data): object 'Survival' not found
 ```
 
 ```r
  message("coxfit summary for overall survival.")
 ```
 
 ```
 ## coxfit summary for overall survival.
 ```
 
 ```r
  summary(coxfit)
 ```
 
 ```
 ## Error in summary(coxfit): object 'coxfit' not found
 ```
 
 ```r
  #plot(cox.zph(coxfit))
 kmfit=survdiff(OverallSurvival ~ groupings)
 ```
 
 ```
 ## Error in eval(expr, envir, enclos): object 'OverallSurvival' not found
 ```
 
 ```r
 message("kaplan meier p-value for overall survival:")
 ```
 
 ```
 ## kaplan meier p-value for overall survival:
 ```
 
 ```r
 1 - pchisq(kmfit$chisq, length(kmfit$n) - 1)
 ```
 
 ```
 ## Error in pchisq(kmfit$chisq, length(kmfit$n) - 1): object 'kmfit' not found
 ```
 
 ```r
 message("calculating the sign of the survival relationship")
 ```
 
 ```
 ## calculating the sign of the survival relationship
 ```
 
 ```r
  mfit=survfit(OverallSurvival ~ groupings)
 ```
 
 ```
 ## Error in eval(expr, envir, enclos): object 'OverallSurvival' not found
 ```
 
 ```r
  plot(mfit,main="overall survival")
 ```
 
 ```
 ## Error in plot(mfit, main = "overall survival"): object 'mfit' not found
 ```
 
 ```r
  coxfit=coxph(OverallSurvivalCutoff~groupings, data= SurvivalCutoff)
 ```
 
 ```
 ## Error in terms.formula(formula, special, data = data): object 'SurvivalCutoff' not found
 ```
 
 ```r
  message("coxfit summary for survival cutoff at ",CutoffPointYears,"years:")
 ```
 
 ```
 ## coxfit summary for survival cutoff at 5years:
 ```
 
 ```r
  summary(coxfit)
 ```
 
 ```
 ## Error in summary(coxfit): object 'coxfit' not found
 ```
 
 ```r
 #plot(cox.zph(coxfit))
  
  kmfit=survdiff(OverallSurvivalCutoff ~ groupings)
 ```
 
 ```
 ## Error in eval(expr, envir, enclos): object 'OverallSurvivalCutoff' not found
 ```
 
 ```r
 message("kaplan meier p-value for survival cutoff:")
 ```
 
 ```
 ## kaplan meier p-value for survival cutoff:
 ```
 
 ```r
  1 - pchisq(kmfit$chisq, length(kmfit$n) - 1)  
 ```
 
 ```
 ## Error in pchisq(kmfit$chisq, length(kmfit$n) - 1): object 'kmfit' not found
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
 ## Error in is.data.frame(x): object 'sampleClustCommPhenoData' not found
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
 ## Error in eval(expr, envir, enclos): object 'sampleClustCommPhenoData' not found
 ```
