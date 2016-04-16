
#Script that pre-processes the curatedBreastData, filters out small datasets,
#and then picks meta-rank gene sets of varying length.

#CHANGE THESE PATHS to your own user directories
#these are the breast-specific saveDir variables and you will use them 
#through the rest of the breast cancer scripts.
saveDirGlobal <- "/home/ywrfc09/breast_analysis"
saveDir_20 <- "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes"
saveDir_no20 <-  "/home/ywrfc09/breast_analysis/metaRankNoTop20Genes"
outputFile <- "/home/kplaney/breast_analysis/breast_proc_outMessages.txt"

###assumes this data is in your path. These are provided on Github:
#load("pam50_centroids_updatedSymbols.RData")
#load("pam50Short_genes.RData")

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
esets <- procExprSetList(exprSetList=curatedBreastDataExprSetList,outputFileDirectory=saveDirGlobal,
                                  minVar=.001,featureDataFieldName="gene_symbol",uniquePDataID="patient_ID")

#saveRDS(esets,file=paste0(saveDirGlobal,"/curatedBreastData_esets_proc.RData.gzip"),compress=TRUE)

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
saveRDS(output,file=paste0(saveDirGlobal,"/mergedExprMatrix_minVar001_17_studies_no_norm.RData.rds"),compress=TRUE)


############get final eset list


##ALSO: merge matrices first to help decide which studies to keep overall

#also merge this one
esets <- readRDS(paste0(saveDirGlobal,"/mergedExprMatrix_minVar001_17_studies_no_norm.RData.rds"))

#NOW: also remove these smaller datasets from the processed esets list before saveRDS

if(length(output$removeDatasetIndices>0)){
  
  dataMatrixList <- dataMatrixList[-output$removeDatasetIndices]
  
  esets_minVar001_17_studies <- esets[-output$removeDatasetIndices]
}


saveRDS(esets_minVar001_17_studies,file=paste0(saveDirGlobal,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.rds"),compress=TRUE)

#this is the "core" data matrix list we will use for our analyses:
saveRDS(dataMatrixList,file=paste0(saveDirGlobal,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.rds"),compress=TRUE)




###meta-features

#ran meta-feature analysis for 1000,500,1000,2000 features.
metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=200,
                                                   numTopFeatFromEachDataset=20,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile=outputFile)

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

#(saveRDSd each one after ran)
saveRDS(metaFeatures,file=paste0(saveDir_20,"/metaFeatures_200.rds"),compress=TRUE)



metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=500,
                                                   numTopFeatFromEachDataset=20,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile=outputFile)

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

saveRDS(metaFeatures,file=paste0(saveDir_20,"/metaFeatures_500.rds"),compress=TRUE)

metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=1000,
                                                   numTopFeatFromEachDataset=20,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile=outputFile)

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

saveRDS(metaFeatures,file=paste0(saveDir_20,"/metaFeatures_1000.rds"),compress=TRUE)

metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=2000,
                                                   numTopFeatFromEachDataset=20,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile=outputFile)

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

saveRDS(metaFeatures,file=paste0(saveDir_20,"/metaFeatures_2000.rds"),compress=TRUE)
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
saveRDS(metaFeatures,file=paste0(saveDir_no20,"/metaFeatures_200.rds"),compress=TRUE)



metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=500,
                                                   numTopFeatFromEachDataset=0,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile=outputFile)

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

saveRDS(metaFeatures,file=paste0(saveDir_no20,"/metaFeatures_500.rds"),compress=TRUE)

load(paste0(saveDir,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.rds"))

metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=1000,
                                                   numTopFeatFromEachDataset=0,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile=outputFile)

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

saveRDS(metaFeatures,file=paste0(saveDir_no20,"/metaFeatures_1000.rds",compress=TRUE))

load(paste0(saveDir,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.rds"))

metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=2000,
                                                   numTopFeatFromEachDataset=0,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile=outputFile)

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

saveRDS(metaFeatures,file=paste0(saveDir_no20,"/metaFeatures_2000.rds"),compress=TRUE)

###also try: 250, 300 genes

metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=250,
                                                   numTopFeatFromEachDataset=0,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile=outputFile)

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

saveRDS(metaFeatures,file=paste0(saveDir_no20,"/metaFeatures_250.rds",compress=TRUE))


metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=300,
                                                   numTopFeatFromEachDataset=0,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile=outputFile)

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

saveRDS(metaFeatures,file=paste0(saveDir_no20,"/metaFeatures_300.rds",compress=TRUE))


###assign PAM50 subtypes to all samples:
###look at Pam50 breakdown 

for(d in 1:length(dataMatrixList)){
  
  tmp <- assignCentroidSubtype(t(dataMatrixList[[d]]),minNumGenes=30,centroidRData="pam50_centroids_updatedSymbols.RData");
 
  tmp <- data.frame(colnames(dataMatrixList[[d]]),rep.int(d,times=ncol(dataMatrixList[[d]])),as.numeric(as.character(tmp$subtypes[,1])),tmp$subtypes[,2],stringsAsFactors=FALSE)
  
  if(d>1){
    
   subtypeDF <- rbind(subtypeDF,tmp)
   
  }else{
    
    subtypeDF <- tmp
    
  }
}

colnames(subtypeDF) <- c("sampleName","studyNum","subtypeNum","subtype")
save(subtypeDF,file=paste0(saveDirGlobal,"pam50Full_subtypeDF.RData.gzip"),compress="gzip")

load("pam50_centroids_updatedSymbols.RData")
load("pam50Short_genes.RData")
centroidMatrix <- centroidMatrix[na.omit(match(pam50Short,centroidMatrix[,1])), ]
save(centroidMatrix,file=paste0(saveDirGlobal, "pam50Short_updatedSymbols_centroidMatrix.RData"))
for(d in 1:length(dataMatrixList)){
  
  tmp <- assignCentroidSubtype(t(dataMatrixList[[d]]),minNumGenes=30,centroidRData=paste0(saveDirGlobal, "/pam50Short_updatedSymbols_centroidMatrix.RData"));
  
  tmp <- data.frame(colnames(dataMatrixList[[d]]),rep.int(d,times=ncol(dataMatrixList[[d]])),as.numeric(as.character(tmp$subtypes[,1])),tmp$subtypes[,2],stringsAsFactors=FALSE)
  
  if(d>1){
    
    subtypeDF_short <- rbind(subtypeDF_short,tmp)
    
  }else{
    
    subtypeDF_short <- tmp
    
  }
}

colnames(subtypeDF_short) <- c("sampleName","studyNum","subtypeNum_short","subtype_short")
#all rows ordered correctly?
all(subtypeDF$sampleName==subtypeDF_short$sampleName)
subtypeDF_master <- cbind(subtypeDF,subtypeDF_short[,c(3:4)])

save(subtypeDF_master,file=paste0(saveDirGlobal, "/pam50FullAndShort_subtypeDF.RData.gzip"),compress="gzip")


