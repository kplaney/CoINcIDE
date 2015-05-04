
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_featureSelection.R")
#ran meta-feature analysis for 1000,500,1000,2000 features.
metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=200,
                                                   numTopFeatFromEachDataset=20,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile="/home/kplaney/breast_analysis//selectFeaturesMetaVarianceOut.txt")

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

#(saved each one after ran)
save(metaFeatures,file="/home/kplaney/breast_analysis/metaFeatures_200.RData.gzip",compress="gzip")



metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=500,
                                                   numTopFeatFromEachDataset=20,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile="/home/kplaney/breast_analysis//selectFeaturesMetaVarianceOut.txt")

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

save(metaFeatures,file="/home/kplaney/breast_analysis/metaFeatures_500.RData.gzip",compress="gzip")

metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=1000,
                                                   numTopFeatFromEachDataset=20,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile="/home/kplaney/breast_analysis//selectFeaturesMetaVarianceOut.txt")

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

save(metaFeatures,file="/home/kplaney/breast_analysis/metaFeatures_1000.RData.gzip",compress="gzip")

metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=2000,
                                                   numTopFeatFromEachDataset=20,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile="/home/kplaney/breast_analysis//selectFeaturesMetaVarianceOut.txt")

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

save(metaFeatures,file="/home/kplaney/breast_analysis/metaFeatures_2000.RData.gzip",compress="gzip")

