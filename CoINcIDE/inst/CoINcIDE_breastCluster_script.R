
###note: pam50 consensus was run in the "select K" script.
###########do for each number of features (just change numFeatures variable)
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
numFeatures <- 2000
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load(paste0("/home/kplaney/breast_analysis/metaFeatures_",numFeatures,".RData.gzip"))
saveDir <- "/home/kplaney/breast_analysis/"
##do for each 200,500,1000,2000 (load different metaFeatures RData object each time.)
#remove datasets with too many missing top gene features


if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- metaFeatures$finalFeatures
  
}

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")

#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansConsensus<- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                                                    pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                                                    numSims=500,maxNumClusters=10,
                                                                    outputFile=paste0(saveDir,"/",numFeatures,"F_clust.txt"),distMethod=c("euclidean"),
                                                                    hclustAlgorithm=c("average"),
                                                                    consensusHclustAlgorithm=c("average"),
                                                                    minClustConsensus=.7, minMeanClustConsensus=.85,
                                                                    corUse="everything",pItem=.9,maxPAC=.15)

save(kmeansConsensus,file=paste0(saveDir,"/curatedbreastData_kmeansConsensus_nstart1pItem9",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")

#######1000
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
numFeatures <- 1000
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load(paste0("/home/kplaney/breast_analysis/metaFeatures_",numFeatures,".RData.gzip"))
saveDir <- "/home/kplaney/breast_analysis/"
##do for each 200,500,1000,2000 (load different metaFeatures RData object each time.)
#remove datasets with too many missing top gene features


if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- metaFeatures$finalFeatures
  
}

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")

#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansConsensus<- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                         pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                         numSims=500,maxNumClusters=10,
                                         outputFile=paste0(saveDir,"/",numFeatures,"F_clust.txt"),distMethod=c("euclidean"),
                                         hclustAlgorithm=c("average"),
                                         consensusHclustAlgorithm=c("average"),
                                         minClustConsensus=.7, minMeanClustConsensus=.85,
                                         corUse="everything",pItem=.9,maxPAC=.15)

save(kmeansConsensus,file=paste0(saveDir,"/curatedbreastData_kmeansConsensus_nstart1pItem9",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")


####500
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
numFeatures <- 500
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load(paste0("/home/kplaney/breast_analysis/metaFeatures_",numFeatures,".RData.gzip"))
saveDir <- "/home/kplaney/breast_analysis/"
##do for each 200,500,1000,2000 (load different metaFeatures RData object each time.)
#remove datasets with too many missing top gene features


if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- metaFeatures$finalFeatures
  
}

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")

#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansConsensus<- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                         pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                         numSims=500,maxNumClusters=10,
                                         outputFile=paste0(saveDir,"/",numFeatures,"F_clust.txt"),distMethod=c("euclidean"),
                                         hclustAlgorithm=c("average"),
                                         consensusHclustAlgorithm=c("average"),
                                         minClustConsensus=.7, minMeanClustConsensus=.85,
                                         corUse="everything",pItem=.9,maxPAC=.15)

save(kmeansConsensus,file=paste0(saveDir,"/curatedbreastData_kmeansConsensus_nstart1pItem9",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")


####200
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
numFeatures <- 200
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load(paste0("/home/kplaney/breast_analysis/metaFeatures_",numFeatures,".RData.gzip"))
saveDir <- "/home/kplaney/breast_analysis/"
##do for each 200,500,1000,2000 (load different metaFeatures RData object each time.)
#remove datasets with too many missing top gene features


if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- metaFeatures$finalFeatures
  
}

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")

#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansConsensus<- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                         pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                         numSims=500,maxNumClusters=10,
                                         outputFile=paste0(saveDir,"/",numFeatures,"F_clust.txt"),distMethod=c("euclidean"),
                                         hclustAlgorithm=c("average"),
                                         consensusHclustAlgorithm=c("average"),
                                         minClustConsensus=.7, minMeanClustConsensus=.85,
                                         corUse="everything",pItem=.9,maxPAC=.15)

save(kmeansConsensus,file=paste0(saveDir,"/curatedbreastData_kmeansConsensus_nstart1pItem9",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")


###250,300
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
numFeatures <- 250
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load(paste0("/home/kplaney/breast_analysis/metaFeatures_",numFeatures,".RData.gzip"))
saveDir <- "/home/kplaney/breast_analysis/"
##do for each 200,500,1000,2000 (load different metaFeatures RData object each time.)
#remove datasets with too many missing top gene features


if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- metaFeatures$finalFeatures
  
}

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")

#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansConsensus<- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                         pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                         numSims=500,maxNumClusters=10,
                                         outputFile=paste0(saveDir,"/",numFeatures,"F_clust.txt"),distMethod=c("euclidean"),
                                         hclustAlgorithm=c("average"),
                                         consensusHclustAlgorithm=c("average"),
                                         minClustConsensus=.7, minMeanClustConsensus=.85,
                                         corUse="everything",pItem=.9,maxPAC=.15)

save(kmeansConsensus,file=paste0(saveDir,"/curatedbreastData_kmeansConsensus_nstart1pItem9",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")

##300
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
numFeatures <- 300
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load(paste0("/home/kplaney/breast_analysis/metaFeatures_",numFeatures,".RData.gzip"))
saveDir <- "/home/kplaney/breast_analysis/"
##do for each 200,500,1000,2000 (load different metaFeatures RData object each time.)
#remove datasets with too many missing top gene features


if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- metaFeatures$finalFeatures
  
}

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")

#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansConsensus<- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                         pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                         numSims=500,maxNumClusters=10,
                                         outputFile=paste0(saveDir,"/",numFeatures,"F_clust.txt"),distMethod=c("euclidean"),
                                         hclustAlgorithm=c("average"),
                                         consensusHclustAlgorithm=c("average"),
                                         minClustConsensus=.7, minMeanClustConsensus=.85,
                                         corUse="everything",pItem=.9,maxPAC=.15)

save(kmeansConsensus,file=paste0(saveDir,"/curatedbreastData_kmeansConsensus_nstart1pItem9",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")


#####k-means clustering, with gap test. use larger number of starts.
###note: pam50 consensus was run in the "select K" script.
###########do for each number of features (just change numFeatures variable)
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
numFeatures <- 2000
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load(paste0("/home/kplaney/breast_analysis/metaFeatures_",numFeatures,".RData.gzip"))
saveDir <- "/home/kplaney/breast_analysis/"
##do for each 200,500,1000,2000 (load different metaFeatures RData object each time.)
#remove datasets with too many missing top gene features


if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- metaFeatures$finalFeatures
  
}

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")

#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansGap<- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                         pickKMethod=c("gap"),iter.max=20,nstart=25,
                                         numSims=500,maxNumClusters=10,
                                         outputFile=paste0(saveDir,"/",numFeatures,"F_clust.txt"),distMethod=c("euclidean"),
                                         hclustAlgorithm=c("average"),
                                         consensusHclustAlgorithm=c("average"),
                                         minClustConsensus=.7, minMeanClustConsensus=.85,
                                         corUse="everything",pItem=.9,maxPAC=.15)

save(kmeansGap,file=paste0(saveDir,"/curatedbreastData_kmeansGap_nstart25_",numFeatures,"_features_",Sys.Date(),".RData.gzip"),compress="gzip")

#######1000
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
numFeatures <- 1000
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load(paste0("/home/kplaney/breast_analysis/metaFeatures_",numFeatures,".RData.gzip"))
saveDir <- "/home/kplaney/breast_analysis/"
##do for each 200,500,1000,2000 (load different metaFeatures RData object each time.)
#remove datasets with too many missing top gene features


if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- metaFeatures$finalFeatures
  
}

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")

#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansGap<- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                         pickKMethod=c("gap"),iter.max=20,nstart=25,
                                         numSims=500,maxNumClusters=10,
                                         outputFile=paste0(saveDir,"/",numFeatures,"F_clust.txt"),distMethod=c("euclidean"),
                                         hclustAlgorithm=c("average"),
                                         consensusHclustAlgorithm=c("average"),
                                         minClustConsensus=.7, minMeanClustConsensus=.85,
                                         corUse="everything",pItem=.9,maxPAC=.15)

save(kmeansGap,file=paste0(saveDir,"/curatedbreastData_kmeansGap_nstart25_",numFeatures,"_features_",Sys.Date(),".RData.gzip"),compress="gzip")


####500
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
numFeatures <- 500
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load(paste0("/home/kplaney/breast_analysis/metaFeatures_",numFeatures,".RData.gzip"))
saveDir <- "/home/kplaney/breast_analysis/"
##do for each 200,500,1000,2000 (load different metaFeatures RData object each time.)
#remove datasets with too many missing top gene features


if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- metaFeatures$finalFeatures
  
}

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")

#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansGap<- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                         pickKMethod=c("gap"),iter.max=20,nstart=25,
                                         numSims=500,maxNumClusters=10,
                                         outputFile=paste0(saveDir,"/",numFeatures,"F_clust.txt"),distMethod=c("euclidean"),
                                         hclustAlgorithm=c("average"),
                                         consensusHclustAlgorithm=c("average"),
                                         minClustConsensus=.7, minMeanClustConsensus=.85,
                                         corUse="everything",pItem=.9,maxPAC=.15)

save(kmeansGap,file=paste0(saveDir,"/curatedbreastData_kmeansGap_nstart25_",numFeatures,"_features_",Sys.Date(),".RData.gzip"),compress="gzip")


####200
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
numFeatures <- 200
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load(paste0("/home/kplaney/breast_analysis/metaFeatures_",numFeatures,".RData.gzip"))
saveDir <- "/home/kplaney/breast_analysis/"
##do for each 200,500,1000,2000 (load different metaFeatures RData object each time.)
#remove datasets with too many missing top gene features


if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- metaFeatures$finalFeatures
  
}

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")

#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansGap<- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                         pickKMethod=c("gap"),iter.max=20,nstart=25,
                                         numSims=500,maxNumClusters=10,
                                         outputFile=paste0(saveDir,"/",numFeatures,"F_clust.txt"),distMethod=c("euclidean"),
                                         hclustAlgorithm=c("average"),
                                         consensusHclustAlgorithm=c("average"),
                                         minClustConsensus=.7, minMeanClustConsensus=.85,
                                         corUse="everything",pItem=.9,maxPAC=.15)

save(kmeansGap,file=paste0(saveDir,"/curatedbreastData_kmeansGap_nstart25_",numFeatures,"_features_",Sys.Date(),".RData.gzip"),compress="gzip")


#####hierarchical clustering, with gap test.
###note: pam50 consensus was run in the "select K" script.
###########do for each number of features (just change numFeatures variable)
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
numFeatures <- 2000
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load(paste0("/home/kplaney/breast_analysis/metaFeatures_",numFeatures,".RData.gzip"))
saveDir <- "/home/kplaney/breast_analysis/"
##do for each 200,500,1000,2000 (load different metaFeatures RData object each time.)
#remove datasets with too many missing top gene features


if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- metaFeatures$finalFeatures
  
}

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")

#we know these are strong clusters. have  minMeanClustConsensus around .85
hclustOut <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("hc"),
                       pickKMethod=c("gap"),iter.max=20,nstart=25,distMethod=c("euclidean"),hclustAlgorithm=c("average"),
                       numSims=500,maxNumClusters=15,
                       outputFile="/home/kplaney/breast_analysis/test.txt"
)


save(hclustOut,file=paste0(saveDir,"/curatedbreastData_hclust_",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")

#######1000
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
numFeatures <- 1000
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load(paste0("/home/kplaney/breast_analysis/metaFeatures_",numFeatures,".RData.gzip"))
saveDir <- "/home/kplaney/breast_analysis/"
##do for each 200,500,1000,2000 (load different metaFeatures RData object each time.)
#remove datasets with too many missing top gene features


if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- metaFeatures$finalFeatures
  
}

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")

#we know these are strong clusters. have  minMeanClustConsensus around .85
hclustOut <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("hc"),
                                    pickKMethod=c("gap"),iter.max=20,nstart=25,distMethod=c("euclidean"),hclustAlgorithm=c("average"),
                                    numSims=500,maxNumClusters=15,
                                    outputFile="/home/kplaney/breast_analysis/test.txt"
)


save(hclustOut,file=paste0(saveDir,"/curatedbreastData_hclust_",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")


####500
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
numFeatures <- 500
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load(paste0("/home/kplaney/breast_analysis/metaFeatures_",numFeatures,".RData.gzip"))
saveDir <- "/home/kplaney/breast_analysis/"
##do for each 200,500,1000,2000 (load different metaFeatures RData object each time.)
#remove datasets with too many missing top gene features


if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- metaFeatures$finalFeatures
  
}

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")

#we know these are strong clusters. have  minMeanClustConsensus around .85
hclustOut <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("hc"),
                                    pickKMethod=c("gap"),iter.max=20,nstart=25,distMethod=c("euclidean"),hclustAlgorithm=c("average"),
                                    numSims=500,maxNumClusters=15,
                                    outputFile="/home/kplaney/breast_analysis/test.txt"
)


save(hclustOut,file=paste0(saveDir,"/curatedbreastData_hclust_",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")



####200
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
numFeatures <- 200
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load(paste0("/home/kplaney/breast_analysis/metaFeatures_",numFeatures,".RData.gzip"))
saveDir <- "/home/kplaney/breast_analysis/"
##do for each 200,500,1000,2000 (load different metaFeatures RData object each time.)
#remove datasets with too many missing top gene features


if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- metaFeatures$finalFeatures
  
}

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")


#we know these are strong clusters. have  minMeanClustConsensus around .85
hclustOut <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("hc"),
                                    pickKMethod=c("gap"),iter.max=20,nstart=25,distMethod=c("euclidean"),hclustAlgorithm=c("average"),
                                    numSims=500,maxNumClusters=15,
                                    outputFile="/home/kplaney/breast_analysis/test.txt"
)


save(hclustOut,file=paste0(saveDir,"/curatedbreastData_hclust_",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")


