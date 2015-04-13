
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_simulation.R")
#create two simulated high quality datasets with a small amount of noise
lungSimData <- createLungSimDatasets(numSimDatasets=1, stddevNoise=.1)
dim(lungSimData$dataMatrixList[[1]])

dataMatrixList <- lungSimData$dataMatrixList
clustFeaturesList <- lungSimData$clustFeaturesList

#let's try the four ways I've coded up clustering, selecting K.
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")
clusterOut_hclust <- clusterMatrixListHclustGap(dataMatrixList=dataMatrixList,clustFeaturesList=clustFeaturesList,maxNumClusters=30,algorithm="ward.D",
                                                distMethod=c("euclidean"),outputFile="/home/kplaney/ovarian_analysis//cluster_hclustGap_output.txt",
                                                corUse=c("everything"),numSims=100)


message(paste0("Best K as determined by hclust and gap test is: ", unlist(clusterOut_hclust$bestK)))
clusterOut_kmeans <- clustMatrixListKmeansGap(dataMatrixList,clustFeaturesList,maxNumClusters=30,iter.max=30,nstart=25,numSims=100,algorithm="Hartigan-Wong",outputFile="/home/kplaney/ovarian_analysis//cluster_kmeansGap_output.txt",
                                              bestK=c("firstMatch"))

message(paste0("Best K as determined by kmeans and gap test is: ", unlist(clusterOut_kmeans$bestK)))

CoINcIDE_getAdjMatrices(dataMatrixList,clustSampleIndexList,clustFeatureIndexList,
                                    edgeMethod=c("distCor","spearman","pearson","kendall","Euclidean","cosine",
                                                 "Manhattan","Minkowski","Mahalanobis"),numParallelCores=1,minTrueSimilThresh=-Inf,maxTrueSimilThresh=Inf,
                                    sigMethod=c("meanMatrix","centroid"),maxNullFractSize=.1,numSims=100,includeRefClustInNull=TRUE,                        
                                    outputFile="./CoINcIDE_messages.txt",fractFeatIntersectThresh=0,numFeatIntersectThresh=0 ,clustSizeThresh=0, clustSizeFractThresh=0)



clusterOut_kmeansConsensus <- consensusClusterMatrixList(dataMatrixList,clustFeaturesList,maxNumClusters = 30, numSims=100, pItem=0.8, pFeature=1, clusterAlg=c("km"),
                                                         innerLinkage="average", finalLinkage_hclust="average",
                                                         corUse=c("everything"),
                                                         distMethod=c("euclidean"),
                                                         minClustConsensus=.7,bestKmethod=c("highestMinConsensus"),
                                                         outputFile="/home/kplaney/ovarian_analysis/consensus_kmeans.txt")

message(paste0("Best K as determined by kmeans and consensus clustering is: ", unlist(clusterOut_kmeansConsensus$bestK)))

clusterOut_hclustConsensus <- consensusClusterMatrixList(dataMatrixList,clustFeaturesList,maxNumClusters = 30, numSims=100, pItem=0.8, pFeature=1, clusterAlg=c("hc"),
                                                         innerLinkage="average", finalLinkage_hclust="average",
                                                         corUse=c("everything"),
                                                         distMethod=c("euclidean"),
                                                         minClustConsensus=.7,bestKmethod=c("highestMinConsensus"),
                                                         outputFile="/home/kplaney/ovarian_analysis/consensus_hclust.txt")


