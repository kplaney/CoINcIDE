pkgname <- "Coincide"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('Coincide')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("CDF_CoINcIDE")
### * CDF_CoINcIDE

flush(stderr()); flush(stdout())

### Name: CDF_CoINcIDE
### Title: TEST
### Aliases: CDF_CoINcIDE
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.




cleanEx()
nameEx("CoINcIDE-package")
### * CoINcIDE-package

flush(stderr()); flush(stdout())

### Name: Coincide-package
### Title: \packageTitleCoincide
### Aliases: Coincide-package Coincide
### Keywords: package

### ** Examples

cat("this is an example")



cleanEx()
nameEx("ConsensusClusterPlus_CoINcIDE")
### * ConsensusClusterPlus_CoINcIDE

flush(stderr()); flush(stdout())

### Name: ConsensusClusterPlus_CoINcIDE
### Title: TEST
### Aliases: ConsensusClusterPlus_CoINcIDE
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.




cleanEx()
nameEx("GSEA")
### * GSEA

flush(stderr()); flush(stdout())

### Name: GSEA
### Title: TEST
### Aliases: GSEA
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.




cleanEx()
nameEx("addClinicalVarToNodeAttributes")
### * addClinicalVarToNodeAttributes

flush(stderr()); flush(stdout())

### Name: addClinicalVarToNodeAttributes
### Title: TEST
### Aliases: addClinicalVarToNodeAttributes
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.





cleanEx()
nameEx("advancedNetworkPlots")
### * advancedNetworkPlots

flush(stderr()); flush(stdout())

### Name: advancedNetworkPlots
### Title: TEST
### Aliases: advancedNetworkPlots
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.



cleanEx()
nameEx("assignCentroidSubtype")
### * assignCentroidSubtype

flush(stderr()); flush(stdout())

### Name: assignCentroidSubtype
### Title: TEST
### Aliases: assignCentroidSubtype
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.




cleanEx()
nameEx("assignFinalEdges")
### * assignFinalEdges

flush(stderr()); flush(stdout())

### Name: assignFinalEdges
### Title: TEST
### Aliases: assignFinalEdges
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.




cleanEx()
nameEx("batchNormalization")
### * batchNormalization

flush(stderr()); flush(stdout())

### Name: batchNormalization
### Title: TEST
### Aliases: batchNormalization
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.




cleanEx()
nameEx("binarizeMetaclustStudyStatus")
### * binarizeMetaclustStudyStatus

flush(stderr()); flush(stdout())

### Name: binarizeMetaclustStudyStatus
### Title: TEST
### Aliases: binarizeMetaclustStudyStatus
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.




cleanEx()
nameEx("clustMatrixListWrapper")
### * clustMatrixListWrapper

flush(stderr()); flush(stdout())

### Name: clustMatrixListWrapper
### Title: TEST
### Aliases: clustMatrixListWrapper
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.





cleanEx()
nameEx("collapseDupProbes")
### * collapseDupProbes

flush(stderr()); flush(stdout())

### Name: collapseDupProbes
### Title: Collapse/handle duplicated probes (genes) in a dataset
### Aliases: collapseDupProbes

### ** Examples

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




cleanEx()
nameEx("computeAdjMatrices")
### * computeAdjMatrices

flush(stderr()); flush(stdout())

### Name: computeAdjMatrices
### Title: TEST
### Aliases: computeAdjMatrices
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.




cleanEx()
nameEx("computeAdjMatricesNullMatrixList")
### * computeAdjMatricesNullMatrixList

flush(stderr()); flush(stdout())

### Name: computeAdjMatricesNullMatrixList
### Title: TEST
### Aliases: computeAdjMatricesNullMatrixList
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.





cleanEx()
nameEx("computeMetaclustEffectSizes")
### * computeMetaclustEffectSizes

flush(stderr()); flush(stdout())

### Name: computeMetaclustEffectSizes
### Title: TEST
### Aliases: computeMetaclustEffectSizes
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.




cleanEx()
nameEx("createExpressionSetList")
### * createExpressionSetList

flush(stderr()); flush(stdout())

### Name: createExpressionSetList
### Title: TEST
### Aliases: createExpressionSetList
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.




cleanEx()
nameEx("createNullDataMatrixList")
### * createNullDataMatrixList

flush(stderr()); flush(stdout())

### Name: createNullDataMatrixList
### Title: TEST
### Aliases: createNullDataMatrixList
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.





cleanEx()
nameEx("createPhenoMasterTableFromMatrixList")
### * createPhenoMasterTableFromMatrixList

flush(stderr()); flush(stdout())

### Name: createPhenoMasterTableFromMatrixList
### Title: TEST
### Aliases: createPhenoMasterTableFromMatrixList
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.




cleanEx()
nameEx("createS4exprSet")
### * createS4exprSet

flush(stderr()); flush(stdout())

### Name: createS4exprSet
### Title: TEST
### Aliases: createS4exprSet
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.




cleanEx()
nameEx("createTissueSimDatasets")
### * createTissueSimDatasets

flush(stderr()); flush(stdout())

### Name: createTissueSimDatasets
### Title: TEST
### Aliases: createTissueSimDatasets
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.





cleanEx()
nameEx("exprSetListToMatrixList")
### * exprSetListToMatrixList

flush(stderr()); flush(stdout())

### Name: exprSetListToMatrixList
### Title: TEST
### Aliases: exprSetListToMatrixList
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.




cleanEx()
nameEx("filterAndImputeSamples")
### * filterAndImputeSamples

flush(stderr()); flush(stdout())

### Name: filterAndImputeSamples
### Title: Filter and Impute Samples
### Aliases: filterAndImputeSamples

### ** Examples

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




cleanEx()
nameEx("filterGenesByVariance")
### * filterGenesByVariance

flush(stderr()); flush(stdout())

### Name: filterGenesByVariance
### Title: Filter genes by variance
### Aliases: filterGenesByVariance

### ** Examples

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



cleanEx()
nameEx("findCommunities")
### * findCommunities

flush(stderr()); flush(stdout())

### Name: findCommunities
### Title: TEST
### Aliases: findCommunities
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.




cleanEx()
nameEx("globalFDR")
### * globalFDR

flush(stderr()); flush(stdout())

### Name: globalFDR
### Title: TEST
### Aliases: globalFDR
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.




cleanEx()
nameEx("global_FDR")
### * global_FDR

flush(stderr()); flush(stdout())

### Name: global_FDR
### Title: TEST
### Aliases: global_FDR
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.




cleanEx()
nameEx("meanMetricDensityPlot")
### * meanMetricDensityPlot

flush(stderr()); flush(stdout())

### Name: meanMetricDensityPlot
### Title: TEST
### Aliases: meanMetricDensityPlot
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.




cleanEx()
nameEx("mergeDatasetList")
### * mergeDatasetList

flush(stderr()); flush(stdout())

### Name: mergeDatasetList
### Title: TEST
### Aliases: mergeDatasetList
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.




cleanEx()
nameEx("merge_datasetList")
### * merge_datasetList

flush(stderr()); flush(stdout())

### Name: merge_datasetList
### Title: TEST
### Aliases: merge_datasetList
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.




cleanEx()
nameEx("procExprSet")
### * procExprSet

flush(stderr()); flush(stdout())

### Name: procExprSet
### Title: Post-process a normalized assayData in an ExpressionSet object
### Aliases: procExprSet

### ** Examples

library("curatedBreastData")
#load up our datasets
data(curatedBreastDataExprSetList);

#just perform on one dataset as an example, GSE9893. 
#This dataset does have NA values, so
#you'll see the impute.knn progress printed to the screen.
#also take only genes that fall in 
#the variance percentiles between .75 and 1 
#(i.e. top 75th percentile genes by variance.)

post_procExprSet <- processExpressionSet(exprSet=
curatedBreastDataExprSetList[[5]], 
outputFileDirectory = "./",
minVarPercentile=.75, maxVarPercentile = 1)




cleanEx()
nameEx("procExprSetList")
### * procExprSetList

flush(stderr()); flush(stdout())

### Name: procExprSetList
### Title: Process a list of S4 expressionSet objects.
### Aliases: procExprSetList procExpressionSetList

### ** Examples

## Not run: 
##D library("curatedBreastData")
##D #warning: takes a while to run! you're processing all datasets in the package!
##D #load up our datasets
##D data(curatedBreastDataExprSetList);
##D 
##D #just take top 5000 genes by variance
##D #this will post-process every dataset in the package
##D #to make them ready for downstream analyses.
##D proc_curatedBreastDataExprSetList <- processExpressionSetList(
##D exprSetList=curatedBreastDataExprSetList, 
##D outputFileDirectory = "./", numTopVarGenes=5000)
## End(Not run)



cleanEx()
nameEx("removeDupPatients")
### * removeDupPatients

flush(stderr()); flush(stdout())

### Name: removeDupPatients
### Title: Remove duplicated patient samples (samples from the same
###   patient/column ID)
### Aliases: removeDupPatients

### ** Examples

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



cleanEx()
nameEx("returnSampleMemberMatrix")
### * returnSampleMemberMatrix

flush(stderr()); flush(stdout())

### Name: returnSampleMemberMatrix
### Title: TEST
### Aliases: returnSampleMemberMatrix
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.




cleanEx()
nameEx("selectFeaturesMetaVariance")
### * selectFeaturesMetaVariance

flush(stderr()); flush(stdout())

### Name: selectFeaturesMetaVariance
### Title: TEST
### Aliases: selectFeaturesMetaVariance
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.




cleanEx()
nameEx("selectFeaturesMetaVariance_wrapper")
### * selectFeaturesMetaVariance_wrapper

flush(stderr()); flush(stdout())

### Name: selectFeaturesMetaVariance_wrapper
### Title: TEST
### Aliases: selectFeaturesMetaVariance_wrapper
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.




cleanEx()
nameEx("selectMetaclustSigGenes")
### * selectMetaclustSigGenes

flush(stderr()); flush(stdout())

### Name: selectMetaclustSigGenes
### Title: TEST
### Aliases: selectMetaclustSigGenes
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.




cleanEx()
nameEx("simulateClusterData")
### * simulateClusterData

flush(stderr()); flush(stdout())

### Name: simulateClusterData
### Title: TEST
### Aliases: simulateClusterData
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.



cleanEx()
nameEx("summarizePosESMetaclustGenes")
### * summarizePosESMetaclustGenes

flush(stderr()); flush(stdout())

### Name: summarizePosESMetaclustGenes
### Title: TEST
### Aliases: summarizePosESMetaclustGenes
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.




### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
