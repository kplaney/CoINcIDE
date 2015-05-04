
###########do for each number of features (just change numFeatures variable)
load("/home/kplaney/ovarian_analysis/esets_proc_TCGAcombat.RData.gzip")
numFeatures <- 2000
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load(paste0("/home/kplaney/ovarian_analysis/metaFeatures_",numFeatures,".RData.gzip"))
saveDir <- "/home/kplaney/ovarian_analysis/"
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

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")



#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansConsensus<- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                                                    pickKMethod=c("consensus"),iter.max=20,nstart=15,
                                                                    numSims=500,maxNumClusters=10,
                                                                    outputFile=paste0(saveDir,"/",numFeatures,"F_clust.txt"),distMethod=c("euclidean"),
                                                                    hclustAlgorithm=c("average"),
                                                                    consensusHclustAlgorithm=c("average"),
                                                                    minClustConsensus=.7, minMeanClustConsensus=.85,
                                                                    corUse="everything",pItem=.9,maxPAC=.15)

save(kmeansConsensus,file=paste0(saveDir,"/curatedOvarianData_kmeansConsensus_",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")




##also try nstart=1

#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansConsensus<- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                         pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                         numSims=500,maxNumClusters=10,
                                         outputFile=paste0(saveDir,"/",numFeatures,"F_clust.txt"),distMethod=c("euclidean"),
                                         hclustAlgorithm=c("average"),
                                         consensusHclustAlgorithm=c("average"),
                                         minClustConsensus=.7, minMeanClustConsensus=.85,
                                         corUse="everything",pItem=.9,maxPAC=.15)

save(kmeansConsensus,file=paste0(saveDir,"/curatedOvarianData_kmeansConsensus_nstart1_",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")







hclustConsensus <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList=clustFeaturesList,clustMethod=c("hc"),
                                          pickKMethod=c("consensus"),
                                          numSims=500,maxNumClusters=10,
                                          outputFile="/home/kplaney/ovarian_analysis/test.txt",iter.max=30,nstart=25,distMethod=c("euclidean"),
                                          hclustAlgorithm=c("ward.D2"),
                                          consensusHclustAlgorithm=c("average"),
                                          minClustConsensus=.7, minMeanClustConsensus=.7,corUse="everything",
                                          pItem=.9)

numFeatures <- 1000
save(hclustConsensus,file=paste0("/home/kplaney/ovarian_analysis/curatedOvarianData_hclustConsensus_",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")




######now to gap test ovarian