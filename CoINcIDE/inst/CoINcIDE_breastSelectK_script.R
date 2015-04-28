

load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
load("/home/kplaney/breast_analysis/pam50Short_genes.RData")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")
saveDir <- "/home/kplaney/breast_analysis/"
numFeatures <- "pam50Short"

clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- pam50Short
  
}

#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansConsensuspam50_short_Nstart15pItem9 <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                          pickKMethod=c("consensus"),iter.max=20,nstart=15,
                                          numSims=500,maxNumClusters=10,
                                          outputFile="/home/kplaney/breast_analysis/test.txt",distMethod=c("euclidean"),
                                          hclustAlgorithm=c("average"),
                                          consensusHclustAlgorithm=c("average"),
                                          minClustConsensus=.7, minMeanClustConsensus=.85,
                                          corUse="everything",pItem=.9,maxPAC=.15)


save(kmeansConsensuspam50_short_Nstart15pItem9,
     file=paste0(saveDir,"/kmeansConsensuspam50_short_Nstart15pItem9_",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")


kmeansConsensuspam50_short_Nstart1pItem9 <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                                                    pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                                                    numSims=500,maxNumClusters=10,
                                                                    outputFile="/home/kplaney/breast_analysis/test.txt",distMethod=c("euclidean"),
                                                                    hclustAlgorithm=c("average"),
                                                                    consensusHclustAlgorithm=c("average"),
                                                                    minClustConsensus=.7, minMeanClustConsensus=.85,
                                                                    corUse="everything",pItem=.9,maxPAC=.15)


save(kmeansConsensuspam50_short_Nstart1pItem9,
     file=paste0(saveDir,"/kmeansConsensuspam50_short_Nstart1pItem9_",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")


kmeansConsensuspam50_short_Nstart15pItem8 <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                                                pickKMethod=c("consensus"),iter.max=20,nstart=15,
                                                                numSims=500,maxNumClusters=10,
                                                                outputFile="/home/kplaney/breast_analysis/test.txt",distMethod=c("euclidean"),
                                                                hclustAlgorithm=c("average"),
                                                                consensusHclustAlgorithm=c("average"),
                                                                minClustConsensus=.7, minMeanClustConsensus=.85,
                                                                corUse="everything",pItem=.8,maxPAC=.15)


save(kmeansConsensuspam50_short_Nstart15pItem8 ,
     file=paste0(saveDir,"/kmeansConsensuspam50_short_Nstart15pItem8_",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")


##hierarchical
hclustConsensuspam50_short_pItem9 <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("hc"),
                                                            pickKMethod=c("consensus"),iter.max=20,nstart=15,
                                                            numSims=500,maxNumClusters=10,
                                                            outputFile="/home/kplaney/breast_analysis/test.txt",distMethod=c("euclidean"),
                                                            hclustAlgorithm=c("ward.D"),
                                                            consensusHclustAlgorithm=c("average"),
                                                            minClustConsensus=.7, minMeanClustConsensus=.85,
                                                            corUse="everything",pItem=.9,maxPAC=.15)


save(hclustConsensuspam50_short_pItem9,file=paste0(saveDir,"/hclustConsensuspam50_short_pItem9_",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")

hclustConsensuspam50_short_pItem8 <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("hc"),
                                                            pickKMethod=c("consensus"),iter.max=20,nstart=15,
                                                            numSims=500,maxNumClusters=10,
                                                            outputFile="/home/kplaney/breast_analysis/test.txt",distMethod=c("euclidean"),
                                                            hclustAlgorithm=c("ward.D"),
                                                            consensusHclustAlgorithm=c("average"),
                                                            minClustConsensus=.7, minMeanClustConsensus=.85,
                                                            corUse="everything",pItem=.8,maxPAC=.15)


save(hclustConsensuspam50_short_pItem8,file=paste0(saveDir,"/hclustConsensuspam50_short_pItem8_",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")

#############
####pam50 full
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
load("/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData")
pam50GenesFull <- centroidMatrix[,1]
saveDir <- "/home/kplaney/breast_analysis/"
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")
numFeatures <- "pam50Full"


#COME BACK: remove clusters with too few pam50 genes?

clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- pam50GenesFull
  
}

kmeansConsensuspam50_full_Nstart15pItem9 <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                          pickKMethod=c("consensus"),iter.max=20,nstart=15,
                          numSims=500,maxNumClusters=10,
                          outputFile="/home/kplaney/breast_analysis/test.txt",distMethod=c("euclidean"),
                          hclustAlgorithm=c("average"),
                          consensusHclustAlgorithm=c("average"),
                          minClustConsensus=.7, minMeanClustConsensus=.85,
                          corUse="everything",pItem=.9,maxPAC=.15)


save(kmeansConsensuspam50_full_Nstart15pItem9,
     file=paste0(saveDir,"/kmeansConsensuspam50_full_Nstart15pItem9_",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")


hclustConsensuspam50_full_pItem9 <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("hc"),
                                          pickKMethod=c("consensus"),iter.max=20,nstart=15,
                                          numSims=500,maxNumClusters=10,
                                          outputFile="/home/kplaney/breast_analysis/test.txt",distMethod=c("euclidean"),
                                          hclustAlgorithm=c("ward.D"),
                                          consensusHclustAlgorithm=c("average"),
                                          minClustConsensus=.7, minMeanClustConsensus=.85,
                                          corUse="everything",pItem=.9,maxPAC=.15)


save(hclustConsensuspam50_full_pItem9 ,file=paste0(saveDir,"/hclustConsensuspam50_full_pItem9_",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")

