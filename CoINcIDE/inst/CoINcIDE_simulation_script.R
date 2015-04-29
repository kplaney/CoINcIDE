
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_simulation.R")
#create two simulated high quality datasets with a small amount of noise
lungSimData <- createLungSimDatasets(numSimDatasets=2, stddevNoise=.1)
dim(lungSimData$dataMatrixList[[1]])

dataMatrixList <- lungSimData$dataMatrixList
clustSampleIndexList <- lungSimData$clustSampleIndexList
clustFeatureIndexList <- lungSimData$clustFeatureIndexList
clustFeaturesList <- lungSimData$clustFeaturesList

#let's try the four ways I've coded up clustering, selecting K.
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")

clusterOut_hclust <- clusterMatrixListHclustGap(dataMatrixList=dataMatrixList,clustFeaturesList=clustFeaturesList,maxNumClusters=8,algorithm="ward.D",
                                                distMethod=c("euclidean"),outputFile="/home/kplaney/lungSims/cluster_hclustGap_output.txt",
                                                corUse=c("everything"),numSims=20)


message(paste0("Best K as determined by hclust and gap test is: ", unlist(clusterOut_hclust$bestK)))

#compute edges
fractFeatIntersectThresh=.7
numFeatIntersectThresh=100
clustSizeThresh=3
clustSizeFractThresh=.05
numParallelCores=8
minTrueSimilThresh=.2
maxNullFractSize=.2
maxTrueSimilThresh=Inf
includeRefClustInNull=TRUE
numSims=10
sigMethod=c("centroid")
edgeMethod=c("pearson")


indEdgePvalueThresh=.3
meanEdgePairPvalueThresh=.2
maxTrueSimilThresh=Inf

clustSampleIndexList <- clusterOut_hclust$clustSampleIndexList
clustFeatureIndexList <- clusterOut_hclust$clustFeatureIndexList

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")
test <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                sigMethod=sigMethod,maxNullFractSize=maxNullFractSize,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                
                                outputFile="/home/kplaney/lungSims/CoINcIDE_hclustConsensus_messages.txt",fractFeatIntersectThresh=fractFeatIntersectThresh,
                                numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, clustSizeFractThresh=clustSizeFractThresh)
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_communityDetection.R")

testEdges <- assignFinalEdges(computeTrueSimilOutput=test$computeTrueSimilOutput,pvalueMatrix=test$pvalueMatrix,indEdgePvalueThresh=indEdgePvalueThresh,
                              meanEdgePairPvalueThresh=meanEdgePairPvalueThresh,
                              minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                              minFractFeatureIntersect=minFractFeatureIntersect,fractFeatIntersectThresh=fractFeatIntersectThresh,numFeatIntersectThresh=numFeatIntersectThresh ,
                              clustSizeThresh=clustSizeThresh, clustSizeFractThresh= clustSizeFractThresh,saveDir="/home/kplaney/lungSims//",fileTag="CoINcIDE_edges_hclustConsensus"
)

communityInfo <- findCommunities(testEdges$filterEdgeOutput$edgeMatrix,testEdges$filterEdgeOutput$edgeWeightMatrix,test$clustIndexMatrix,fileTag="lung_communityNodeAttributes_",
                                             saveDir="/home/kplaney/lungSims/",minNumUniqueStudiesPerCommunity=2,clustMethodName="hclust",
                                             commMethod=c("edgeBetween"),
                                             makePlots=TRUE,saveGraphData=TRUE)
#get meta-communities.



#now: analyze genes





#####
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


