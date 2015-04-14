
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_simulation.R")
#create two simulated high quality datasets with a small amount of noise
lungSimData <- createLungSimDatasets(numSimDatasets=1, stddevNoise=.1)
dim(lungSimData$dataMatrixList[[1]])

dataMatrixList <- lungSimData$dataMatrixList
clustFeaturesList <- lungSimData$clustFeaturesList

#let's try the four ways I've coded up clustering, selecting K.
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")

clusterOut_hclust <- clusterMatrixListHclustGap(dataMatrixList=dataMatrixList,clustFeaturesList=clustFeaturesList,maxNumClusters=30,algorithm="ward.D",
                                                distMethod=c("euclidean"),outputFile="/home/kplaney/lungSims/cluster_hclustGap_output.txt",
                                                corUse=c("everything"),numSims=100)


message(paste0("Best K as determined by hclust and gap test is: ", unlist(clusterOut_hclust$bestK)))

clusterOut_kmeans <- clustMatrixListKmeansGap(dataMatrixList,clustFeaturesList,maxNumClusters=30,iter.max=30,nstart=25,numSims=100,algorithm="Hartigan-Wong",
                                              outputFile="/home/kplaney//home/kplaney/lungSims//cluster_kmeansGap_output.txt"
                                           )

message(paste0("Best K as determined by kmeans and gap test is: ", unlist(clusterOut_kmeans$bestK)))

clusterOut_kmeansConsensus <- consensusClusterMatrixList(dataMatrixList,clustFeaturesList,maxNumClusters = 30, numSims=100, pItem=0.8, pFeature=1, clusterAlg=c("km"),
                                                         innerLinkage="average", finalLinkage_hclust="average",
                                                         corUse=c("everything"),
                                                         distMethod=c("euclidean"),
                                                         minClustConsensus=.7,bestKmethod=c("highestMinConsensus"),
                                                         outputFile="/home/kplaney/lungSims/consensus_kmeans.txt")

message(paste0("Best K as determined by kmeans and consensus clustering is: ", unlist(clusterOut_kmeansConsensus$bestK)))

#ward.D seems to work well. I tried average consensus at first and it didn't work so well.
clusterOut_hclustConsensus <- consensusClusterMatrixList(dataMatrixList,clustFeaturesList,maxNumClusters = 30, numSims=100, pItem=0.8, pFeature=1, clusterAlg=c("hc"),
                                                         hclustAlgorithm="ward.D", consensusHclustAlgorithm="ward.D",
                                                         corUse=c("everything"),
                                                         distMethod=c("euclidean"),
                                                         minClustConsensus=.8,bestKmethod=c("highestMinConsensus"),
                                                         outputFile="/home/kplaney/lungSims/consensus_hclust.txt")

message(paste0("Best K as determined by hclust and consensus clustering is: ", unlist(clusterOut_hclustConsensus$bestK)))


