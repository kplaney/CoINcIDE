
#create two simulated high quality datasets with a small amount of noise
lungSimData <- createLungSimDatasets(numSimDatasets=2, stddevNoise=.1)
dim(lungSimData$dataMatrixList[[1]])

clustFeatures <- lungSimData$clustFeaturesList[[1]]

dataset <- dataMatrix[rownames(lungSimData$dataMatrixList[[1]]) %in% clustFeatures, , drop=FALSE]
##test my gap statistic
 source('~/Dropbox/Gevaert_lab/code_RNAseq_pipeline/RP_gapTest.R')
 test <- gapTest(dataset,k=c(1:40),nstart=25,iter.max=10,clusterMethod=c("kmeans"),
                     hclustAlgorithm="complete",distMatrix=NULL,stopAtBestK=TRUE,
                    kmeansAlgorithm="Hartigan-Wong",clusterColumns=TRUE,numSimulations=10,refDistMethod=c("svd"))

 clustF <- function(x,k){
    #k-means clusters the ROWs
    #it turns out that R will recognize the upper-level function input value in this function.
    #transpose BEFORE feed in here; I found it to return odd results if
    #I feed in t(x) in the kmeans function that is feed into clusGap
    kmeans(x, centers=k,iter.max=10,nstart=25,algorithm="Hartigan-Wong");
  
  }
  
  #k-means clusters the ROWs - so transpose dataset
  gapTestOutput <- clusGap(t(dataset), FUNcluster=clustF,K.max=maxNumClusters, B = numSims, verbose = interactive());

#best K?  4 still - i.e. the methods match up. might as well use 'cluster' package and not your own code.
maxSE(f=gapTestOut$Tab[,"gap"], SE.f=gapTestOut$Tab[,"SE.sim"],
       method = c("globalSEmax"))


#let's try hclust gap
lung_hclust <- clusterMatrixListHclustGap(dataMatrixList,clustFeaturesList,maxNumClusters=10,algorithm="ward.D",
                              distMethod=c("euclidean"),outputFile="~/Desktop/cluster_hclustGap_output.txt",numSims=100)

lung_hclust$bestK
