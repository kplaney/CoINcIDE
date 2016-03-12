
library("Coincide")

##CHANGE these paths to your user directory
saveDirGlobal <- "/home/kplaney/ovarian_analysis/"
saveDir_20 <- "/home/ywrfc09/ovarian_analysis/metaRankWithTop20Genes"

#load up our processed data from ovarianProcessAndGeneFeatures script.
dataMatrixList <- readRDS(paste0(saveDirGlobal,"/esets_proc_TCGAcombat.rds"))

###########do for each number of features (just change numFeatures variable)

numFeatures <- 2000
metaFeatures <- readRDS(paste0(saveDir_no20,"/metaFeatures_",numFeatures,".rds"))

#remove datasets with too many missing top gene features

dataMatrixList <- exprSetListToMatrixList(esets,featureDataFieldName="gene")

names(dataMatrixList) <- names(esets)

if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- metaFeatures$finalFeatures
  
}

##we ended up using nstart=1 for analyses

#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansConsensus<- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                         pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                         numSims=500,maxNumClusters=10,
                                         outputFile=paste0(saveDirGlobal,"/",numFeatures,"F_clust.txt"),distMethod=c("euclidean"),
                                         hclustAlgorithm=c("average"),
                                         consensusHclustAlgorithm=c("average"),
                                         minClustConsensus=.7, minMeanClustConsensus=.85,
                                         corUse="everything",pItem=.9,maxPAC=.15)

saveRDS(kmeansConsensus,file=paste0(saveDirGlobal,"/curatedOvarianData_kmeansConsensus_nstart1_",numFeatures,"Features_",Sys.Date(),".rds"),compress=TRUE)



###########do for each number of features (just change numFeatures variable)
numFeatures <- 1000
metaFeatures <- readRDS(paste0(saveDir_no20,"/metaFeatures_",numFeatures,".rds"))

##do for each 200,500,1000,2000 (load different metaFeatures RData object each time.)
#remove datasets with too many missing top gene features

dataMatrixList <- exprSetListToMatrixList(esets,featureDataFieldName="gene")

names(dataMatrixList) <- names(esets)

if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- metaFeatures$finalFeatures
  
}

##we ended up using nstart=1 for analyses

#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansConsensus<- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                         pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                         numSims=500,maxNumClusters=10,
                                         outputFile=paste0(saveDirGlobal,"/",numFeatures,"F_clust.txt"),distMethod=c("euclidean"),
                                         hclustAlgorithm=c("average"),
                                         consensusHclustAlgorithm=c("average"),
                                         minClustConsensus=.7, minMeanClustConsensus=.85,
                                         corUse="everything",pItem=.9,maxPAC=.15)

saveRDS(kmeansConsensus,file=paste0(saveDirGlobal,"/curatedOvarianData_kmeansConsensus_nstart1_",numFeatures,"Features_",Sys.Date(),".rds"),compress=TRUE)

###########do for each number of features (just change numFeatures variable)
numFeatures <- 500
metaFeatures <- readRDS(paste0(saveDir_no20,"/metaFeatures_",numFeatures,".rds"))
##do for each 200,500,1000,2000 (load different metaFeatures RData object each time.)
#remove datasets with too many missing top gene features

dataMatrixList <- exprSetListToMatrixList(esets,featureDataFieldName="gene")

names(dataMatrixList) <- names(esets)

if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- metaFeatures$finalFeatures
  
}

##we ended up using nstart=1 for analyses

#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansConsensus<- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                         pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                         numSims=500,maxNumClusters=10,
                                         outputFile=paste0(saveDirGlobal,"/",numFeatures,"F_clust.txt"),distMethod=c("euclidean"),
                                         hclustAlgorithm=c("average"),
                                         consensusHclustAlgorithm=c("average"),
                                         minClustConsensus=.7, minMeanClustConsensus=.85,
                                         corUse="everything",pItem=.9,maxPAC=.15)

saveRDS(kmeansConsensus,file=paste0(saveDirGlobal,"/curatedOvarianData_kmeansConsensus_nstart1_",numFeatures,"Features_",Sys.Date(),".rds"),compress=TRUE)



###########do for each number of features (just change numFeatures variable)
numFeatures <- 200
metaFeatures <- readRDS(paste0(saveDir_no20,"/metaFeatures_",numFeatures,".rds"))
##do for each 200,500,1000,2000 (load different metaFeatures RData object each time.)
#remove datasets with too many missing top gene features

dataMatrixList <- exprSetListToMatrixList(esets,featureDataFieldName="gene")

names(dataMatrixList) <- names(esets)

if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- metaFeatures$finalFeatures
  
}

##we ended up using nstart=1 for analyses

#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansConsensus<- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                         pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                         numSims=500,maxNumClusters=10,
                                         outputFile=paste0(saveDirGlobal,"/",numFeatures,"F_clust.txt"),distMethod=c("euclidean"),
                                         hclustAlgorithm=c("average"),
                                         consensusHclustAlgorithm=c("average"),
                                         minClustConsensus=.7, minMeanClustConsensus=.85,
                                         corUse="everything",pItem=.9,maxPAC=.15)

saveRDS(kmeansConsensus,file=paste0(saveDirGlobal,"/curatedOvarianData_kmeansConsensus_nstart1_",numFeatures,"Features_",Sys.Date(),".rds"),compress=TRUE)


###250
###########do for each number of features (just change numFeatures variable)
numFeatures <- 250
metaFeatures <- readRDS(paste0(saveDir_no20,"/metaFeatures_",numFeatures,".rds"))
##do for each 200,500,1000,2000 (load different metaFeatures RData object each time.)
#remove datasets with too many missing top gene features

dataMatrixList <- exprSetListToMatrixList(esets,featureDataFieldName="gene")

names(dataMatrixList) <- names(esets)

if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- metaFeatures$finalFeatures
  
}

##we ended up using nstart=1 for analyses

#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansConsensus<- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                         pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                         numSims=500,maxNumClusters=10,
                                         outputFile=paste0(saveDirGlobal,"/",numFeatures,"F_clust.txt"),distMethod=c("euclidean"),
                                         hclustAlgorithm=c("average"),
                                         consensusHclustAlgorithm=c("average"),
                                         minClustConsensus=.7, minMeanClustConsensus=.85,
                                         corUse="everything",pItem=.9,maxPAC=.15)

saveRDS(kmeansConsensus,file=paste0(saveDirGlobal,"/curatedOvarianData_kmeansConsensus_nstart1_",numFeatures,"Features_",Sys.Date(),".rds"),compress=TRUE)


###300
numFeatures <- 300
metaFeatures <- readRDS(paste0(saveDir_no20,"/metaFeatures_",numFeatures,".rds"))

#remove datasets with too many missing top gene features

dataMatrixList <- exprSetListToMatrixList(esets,featureDataFieldName="gene")

names(dataMatrixList) <- names(esets)

if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- metaFeatures$finalFeatures
  
}

##we ended up using nstart=1 for analyses

#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansConsensus<- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                         pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                         numSims=500,maxNumClusters=10,
                                         outputFile=paste0(saveDirGlobal,"/",numFeatures,"F_clust.txt"),distMethod=c("euclidean"),
                                         hclustAlgorithm=c("average"),
                                         consensusHclustAlgorithm=c("average"),
                                         minClustConsensus=.7, minMeanClustConsensus=.85,
                                         corUse="everything",pItem=.9,maxPAC=.15)

saveRDS(kmeansConsensus,file=paste0(saveDirGlobal,"/curatedOvarianData_kmeansConsensus_nstart1_",numFeatures,"Features_",Sys.Date(),".rds"),compress=TRUE)










