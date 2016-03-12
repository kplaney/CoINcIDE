
#Script that pre-processes the curatedBreastData, filters out small datasets,
#and then picks meta-rank gene sets of varying length.

#CHANGE THESE PATHS to your own user directories
saveRDSDirGlobal <- "/home/ywrfc09/breast_analysis"
saveRDSDir_20 <- "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes"
saveRDSDir_no20 <-  "/home/ywrfc09/breast_analysis/metaRankNoTop20Genes"
outputFile <- "/home/kplaney/breast_analysis/breast_proc_outMessages.txt"
###esets
library("curatedBreastData")

#LOAD data
#file is also saveRDSd here:
#"/home/kplaney/R/x86_64-redhat-linux-gnu-library/curatedBreastData/data/curatedBreastDataExprSetList.rda"
data("curatedBreastDataExprSetList")

#same parameters as for ovarian
#use this code from Coincide and NOT curatedBreastData - so override some of
#curatedBreasetData's functions (future to-do is to move all curatedBreastData
#proc functions to Coincide)
library("Coincide")
esets <- procExprSetList(exprSetList=curatedBreastDataExprSetList,outputFileDirectory=saveRDSDirGlobal,
                                  minVar=.001,featureDataFieldName="gene_symbol",uniquePDataID="patient_ID")

#saveRDS(esets,file=paste0(saveRDSDirGlobal,"/curatedBreastData_esets_proc.RData.gzip"),compress=TRUE)

#remove curatedBreastData package, curatedBreastDataExprSetList - done with them now.
rm(list="curatedBreastData")
rm(list="curatedBreastDataExprSetList")

###now merge datasets for concatenated clustering in merge_analysis script

#with no norm to see which datasets needed to be dropped

#now format just as a list of data matrices.
dataMatrixList <- exprSetListToMatrixList(esets,featureDataFieldName="gene_symbol")

names(dataMatrixList) <- names(esets)

## merge matrices right now/up front to help decide which studies to keep overall - i.e. greater than
#40 patients and at least 10,000 genes.  will help with merged analyses.

#also merge this one
output <- mergeDatasetList(datasetList=dataMatrixList,minNumGenes = 10000, minNumPatients = 40,batchNormalize = c('none'));
saveRDS(output,file=paste0(saveRDSDirGlobal,"/mergedExprMatrix_minVar001_17_studies_no_norm.RData.gzip"),compress=TRUE)


############get final eset list


##ALSO: merge matrices first to help decide which studies to keep overall

#also merge this one
load(paste0(saveRDSDirGlobal,"/mergedExprMatrix_minVar001_17_studies_no_norm.RData.gzip"),compress=TRUE)

#NOW: also remove these smaller datasets from the processed esets list before saveRDS

if(length(output$removeDatasetIndices>0)){
  
  dataMatrixList <- dataMatrixList[-output$removeDatasetIndices]
  
  esets_minVar001_17_studies <- esets[-output$removeDatasetIndices]
}


saveRDS(esets_minVar001_17_studies,file=paste0(saveRDSDirGlobal,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.rds"),compress=TRUE)

#this is the "core" data matrix list we will use for our analyses:
saveRDS(dataMatrixList,file=paste0(saveRDSDirGlobal,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.rds"),compress=TRUE)




###meta-features

#ran meta-feature analysis for 1000,500,1000,2000 features.
metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=200,
                                                   numTopFeatFromEachDataset=20,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile=outputFile)

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

#(saveRDSd each one after ran)
saveRDS(metaFeatures,file=paste0(saveRDSDir_20,"/metaFeatures_200.rds"),compress=TRUE)



metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=500,
                                                   numTopFeatFromEachDataset=20,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile=outputFile)

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

saveRDS(metaFeatures,file=paste0(saveRDSDir_20,"/metaFeatures_500.rds"),compress=TRUE)

metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=1000,
                                                   numTopFeatFromEachDataset=20,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile=outputFile)

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

saveRDS(metaFeatures,file=paste0(saveRDSDir_20,"/metaFeatures_1000.rds"),compress=TRUE)

metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=2000,
                                                   numTopFeatFromEachDataset=20,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile=outputFile)

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

saveRDS(metaFeatures,file=paste0(saveRDSDir_20,"/metaFeatures_2000.rds"),compress=TRUE)
###########
#now try without top-varying genes from each dataset

#ran meta-feature analysis for 1000,500,1000,2000 features.
metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=200,
                                                   numTopFeatFromEachDataset=0,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile=outputFile)

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

#(saveRDSd each one after ran)
saveRDS(metaFeatures,file=paste0(saveRDSDir_no20,"/metaFeatures_200.rds"),compress=TRUE)



metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=500,
                                                   numTopFeatFromEachDataset=0,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile=outputFile)

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

saveRDS(metaFeatures,file=paste0(saveRDSDir_no20,"/metaFeatures_500.rds"),compress=TRUE)

load(paste0(saveRDSDir,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.rds"))

metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=1000,
                                                   numTopFeatFromEachDataset=0,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile=outputFile)

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

saveRDS(metaFeatures,file=paste0(saveRDSDir_no20,"/metaFeatures_1000.rds",compress=TRUE))

load(paste0(saveRDSDir,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.rds"))

metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=2000,
                                                   numTopFeatFromEachDataset=0,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile=outputFile)

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

saveRDS(metaFeatures,file=paste0(saveRDSDir_no20,"/metaFeatures_2000.rds"),compress=TRUE)

###also try: 250, 300 genes

metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=250,
                                                   numTopFeatFromEachDataset=0,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile=outputFile)

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

saveRDS(metaFeatures,file=paste0(saveRDSDir_no20,"/metaFeatures_250.rds",compress=TRUE))


metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=300,
                                                   numTopFeatFromEachDataset=0,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile=outputFile)

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

saveRDS(metaFeatures,file=paste0(saveRDSDir_no20,"/metaFeatures_300.rds",compress=TRUE))

