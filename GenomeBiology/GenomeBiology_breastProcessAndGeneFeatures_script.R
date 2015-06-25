

saveDirGlobal <- "/home/kplaney/breast_analysis"
saveDirMetaWithTop20Genes <- "/home/kplaney/breast_analysis/metaRankWithTop20Genes"
saveDirMetaNoTop20Genes <-  "/home/kplaney/breast_analysis/metaRankNoTop20Genes"
###esets
library("curatedBreastData")
#LOAD data
#file is also saved here:
#"/home/kplaney/R/x86_64-redhat-linux-gnu-library/curatedBreastData/data/curatedBreastDataExprSetList.rda"
data("curatedBreastDataExprSetList")

#same parameters as for ovarian
#use this code from CoINcIDE
library("CoINcIDE")
esets <- processExprSetList(exprSetList=curatedBreastDataExprSetList,outputFileDirectory=saveDirGlobal,
                                  minVar=.001,featureDataFieldName="gene",uniquePDataID="unique_patient_ID")

save(esets,file=paste0(saveDirGlobal,"/curatedBreastData_esets_proc.RData.gzip"),compress="gzip")

#remove curatedBreastData package so functions don't conflict.
rm(list="curatedBreastData")
###merge datasets - with no norm to see which datasets needed to be dropped

load(paste0(saveDirGlobal,"/curatedBreastData_esets_proc.RData.gzip"))
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_featureSelection.R")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
#combat code:
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_batchCorrection.R")
#now format just as a list of data matrices.
dataMatrixList <- exprSetListToMatrixList(esets,featureDataFieldName="gene_symbol")

names(dataMatrixList) <- names(esets)

## merge matrices right now/up front to help decide which studies to keep overall - i.e. greater than
#40 patients and at least 10,000 genes.  will help with merged analyses.

#also merge this one
output <- merge_datasetList(datasetList=dataMatrixList,minNumGenes = 10000, minNumPatients = 40,batchNormalize = c('none'));
save(output,file=paste0(saveDirGlobal,"/mergedExprMatrix_minVar001_17_studies_no_norm.RData.gzip"),compress="gzip")


############get final eset list
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_featureSelection.R")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")

load(paste0(saveDirGlobal,"/curatedBreastData_esets_proc.RData.gzip"))
#now format just as a list of data matrices.
dataMatrixList <- exprSetListToMatrixList(esets,featureDataFieldName="gene_symbol")

names(dataMatrixList) <- names(esets)


##ALSO: merge matrices first to help decide which studies to keep overall

#also merge this one
load(paste0(saveDirGlobal,"/mergedExprMatrix_minVar001_17_studies_no_norm.RData.gzip"),compress="gzip")

#NOW: also remove these smaller datasets from the esets list before save

if(length(output$removeDatasetIndices>0)){
  
  dataMatrixList <- dataMatrixList[-output$removeDatasetIndices]
  
  esets_minVar001_17_studies <- esets[-output$removeDatasetIndices]
}


save(esets_minVar001_17_studies,file=paste0(saveDirGlobal,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"),compress="gzip")

#this is the "core" data matrix list we will use for our analyses:
save(dataMatrixList,file=paste0(saveDirGlobal,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip"),compress="gzip")




###meta-features

load(paste0(saveDir,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_featureSelection.R")
#ran meta-feature analysis for 1000,500,1000,2000 features.
metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=200,
                                                   numTopFeatFromEachDataset=20,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile="/home/kplaney/breast_analysis_withTop20Genes//selectFeaturesMetaVarianceOut.txt")

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

#(saved each one after ran)
save(metaFeatures,file=paste0(saveDirMetaWithTop20Genes,"/metaFeatures_200.RData.gzip"),compress="gzip")



metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=500,
                                                   numTopFeatFromEachDataset=20,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile="/home/kplaney/breast_analysis_withTop20Genes//selectFeaturesMetaVarianceOut.txt")

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

save(metaFeatures,file=paste0(saveDirMetaWithTop20Genes,"/metaFeatures_500.RData.gzip"),compress="gzip")

metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=1000,
                                                   numTopFeatFromEachDataset=20,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile="/home/kplaney/breast_analysis_withTop20Genes//selectFeaturesMetaVarianceOut.txt")

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

save(metaFeatures,file=paste0(saveDirMetaWithTop20Genes,"/metaFeatures_1000.RData.gzip"),compress="gzip")

metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=2000,
                                                   numTopFeatFromEachDataset=20,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile="/home/kplaney/breast_analysis_withTop20Genes//selectFeaturesMetaVarianceOut.txt")

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

save(metaFeatures,file=paste0(saveDirMetaWithTop20Genes,"/metaFeatures_2000.RData.gzip"),compress="gzip")
###########
#now try without top-varying genes from each dataset
saveDir <- 
  
load(paste0(saveDir,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_featureSelection.R")
#ran meta-feature analysis for 1000,500,1000,2000 features.
metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=200,
                                                   numTopFeatFromEachDataset=0,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile="/home/kplaney/breast_analysis//selectFeaturesMetaVarianceOut.txt")

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

#(saved each one after ran)
save(metaFeatures,file=paste0(saveDirMetaNoTop20Genes,"/metaFeatures_200.RData.gzip"),compress="gzip")



load(paste0(saveDir,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_featureSelection.R")

metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=500,
                                                   numTopFeatFromEachDataset=0,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile="/home/kplaney/breast_analysis//selectFeaturesMetaVarianceOut.txt")

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

save(metaFeatures,file=paste0(saveDirMetaNoTop20Genes,"/metaFeatures_500.RData.gzip"),compress="gzip")

load(paste0(saveDir,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_featureSelection.R")

metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=1000,
                                                   numTopFeatFromEachDataset=0,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile="/home/kplaney/breast_analysis//selectFeaturesMetaVarianceOut.txt")

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

save(metaFeatures,file=paste0(saveDirMetaNoTop20Genes,"/metaFeatures_1000.RData.gzip",compress="gzip"))

load(paste0(saveDir,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_featureSelection.R")
metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=2000,
                                                   numTopFeatFromEachDataset=0,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile="/home/kplaney/breast_analysis//selectFeaturesMetaVarianceOut.txt")

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

save(metaFeatures,file=paste0(saveDirMetaNoTop20Genes,"/metaFeatures_2000.RData.gzip"),compress="gzip")

###also try: 250, 300 genes
load(paste0(saveDir,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_featureSelection.R")
metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=250,
                                                   numTopFeatFromEachDataset=0,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile="/home/kplaney/breast_analysis//selectFeaturesMetaVarianceOut.txt")

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

save(metaFeatures,file=paste0(saveDirMetaNoTop20Genes,"/metaFeatures_250.RData.gzip",compress="gzip"))


load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_featureSelection.R")
metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=300,
                                                   numTopFeatFromEachDataset=0,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile="/home/kplaney/breast_analysis//selectFeaturesMetaVarianceOut.txt")

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

save(metaFeatures,file=paste0(saveDirMetaNoTop20Genes,"/metaFeatures_300.RData.gzip",compress="gzip"))

