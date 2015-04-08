#use default Ward's method
library("cluster")
library("ConsensusClusterPlus")

clustMatrixListKmeansGap <- function(dataMatrixList,clustFeaturesList,maxNumClusters=30,iter.max=30,nstart=25,numSims=1000,algorithm="Hartigan-Wong",outputFile="./cluster_kmeansGap_output.txt"){
  
  
  #lapply doens't work because also need to loop over clust features list.
  #mapply didn't quite work either.
  clustOutput <- list()
  for(d in 1:length(dataMatrixList)){
    
    outputList[[d]] <- clusterMatrixKmeansGap(dataMatrixList[[d]],clustFeaturesList[[d]],
                                          maxNumClusters=maxNumClusters,iter.max=iter.max,
                                          nstart=nstart,algorithm=algorithm,numSims=numSims,outputFile=outputFile)
    
  }

  clustSampleIndexList <- list()
  clustFeatureIndexList <- list()
  bestK <- list()
  gapTest <- list()
  
  for(d in 1:length(outputList)){
    
    clustSampleIndexList[[d]] <- outputList[[d]]$clustSampleIndexList
    clustFeatureIndexList[[d]] <- outputList[[d]]$clustSampleIndexList
    bestK[[d]] <- outputList[[d]]$bestK
    gapTest[[d]] <- outputList[[d]]$gapTest
    
  }
  
  output <- list(clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                 bestK=bestK,gapTest=gapTest)
  
  return(output)

}

clusterMatrixKmeansGap <- function(dataMatrix,clustFeatures,maxNumClusters=30,iter.max=30,nstart=25,numSims=1000,algorithm="Hartigan-Wong",outputFile="./cluster_kmeansGap_output.txt"){
  
  dataset <- dataMatrix[rownames(dataMatrix) %in% clustFeatures, , drop=FALSE]
  
  if(nrow(dataset)==0){
    
    stop("\nIn clusterMatrixKmeansGap function: no clustFeatures were found in the data matrix inputted.")
    
  }

  clustF <- function(x,k){
    #k-means clusters the ROWs
    #it turns out that R will recognize the upper-level function input value in this function.
    #transpose BEFORE feed in here; I found it to return odd results if
    #I feed in t(x) in the kmeans function that is feed into clusGap
    kmeans(x, centers=k,iter.max=iter.max,nstart=nstart,algorithm=algorithm);
  
  }
  
  #k-means clusters the ROWs - so transpose dataset
  gapTest <- clusGap(t(dataset), FUNcluster=clustF,K.max=maxNumClusters, B = numSims, verbose = interactive());

  #f(k) is the gap statistic.
  # method = c("firstSEmax", "Tibs2001SEmax", "globalSEmax", "firstmax", "globalmax")

  bestK <- maxSE(f=gapTest$Tab[,"gap"], SE.f=gapTest$Tab[,"SE.sim"],
      method = c("Tibs2001SEmax"),
      SE.factor = 1)
  #way to confirm bestK is actually working specifically for Tibs2001SEmax - compute by hand:
  select_metric <- gapTest$Tab[c(1:(maxNumClusters-1)),"gap"] > gapTest$Tab[c(2:maxNumClusters),"gap"] - gapTest$Tab[c(2:maxNumClusters),"SE.sim"]
  #take first one that is true.
  best_k <- which(select_metric)[1]
  
  
  if(bestK != best_k){
    
    stop("\nError when selecting best K in clusGap function in cluster package.")
  }

  cat(paste0("\nBest K as determined by k-means and gap test is: ", bestK ,"\n"),
      append=TRUE,file=outputFile);

  clustFeatureIndexList <- list()
  clustSampleIndexList <- list()

   clusterAssignments <- clustF(dataset,bestK)$cluster
  
 for(k in 1:bestK){
   
   clustFeatureIndexList[[k]] <- clustFeatures
   clustSampleIndexList[[k]] <- which(clusterAssignments==k)
   
   
 }

  output <- list(bestK=bestK,clustFeatureIndexList=clustFeatureIndexList,
                 clustSampleIndexList=clustSampleIndexList, gapTest= gapTest)
 return(output)
  
}

clusterMatrixListHclustGap <- function(dataMatrixList,clustFeaturesList,maxNumClusters=30,algorithm="complete",
                              distMethod=c("pearson","spearman","euclidean", "binary", "maximum", "canberra", "minkowski"),outputFile="./cluster_kmeansGap_output.txt",
                              corUse=c("everything","pairwise.complete.obs", "complete.obs"),numSims=1000){
  
  clustOutput <- list()
  for(d in 1:length(dataMatrixList)){
    
   outputList[[d]] <- clusterMatrixHclustGap(dataMatrixList[[d]],clustFeaturesList[[d]],maxNumClusters=maxNumClusters,algorithm=algorithm,numSims=numSims,distMethod=distMethod,outputFile=outputFile,corUse=corUse)
    
  }

  clustSampleIndexList <- list()
  clustFeatureIndexList <- list()
  bestK <- list()
  gapTest <- list()

  for(d in 1:length(outputList)){
    
    clustSampleIndexList[[d]] <- outputList[[d]]$clustSampleIndexList
    clustFeatureIndexList[[d]] <- outputList[[d]]$clustSampleIndexList
    bestK[[d]] <- outputList[[d]]$bestK
    gapTest[[d]] <- outputList[[d]]$gapTest
    
  }
   
  output <- list(clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                 bestK=bestK,gapTest=gapTest)

return(output)

}

clusterMatrixHclustGap <- function(dataMatrix,clustFeatures,maxNumClusters=30,
                                   algorithm=c("complete","ward.D", "ward.D2", "single", "average","mcquitty","median","centroid"),                     
                              distMethod=c("pearson","spearman","euclidean", "binary", "maximum", "canberra", "minkowski"),outputFile="./cluster_kmeansGap_output.txt",
                              corUse=c("everything","pairwise.complete.obs", "complete.obs"),numSims=1000){
  
  dataset <- dataMatrix[rownames(dataMatrix) %in% clustFeatures, , drop=FALSE]

  if(nrow(dataset)==0){
    
    stop("\nIn clusterMatrixHclustGap function: no clustFeatures were found in the data matrix inputted.")
    
  } 
  
  if(distMethod==("pearson") || distMethod=="spearman"){
    
  clustF <- function(x,k){
   #it turns out that R will recognize the upper-level function input value in this function (algorith, corUse, etc.)
    clustObject <- hclust(as.dist((1-cor(x,use=corUse,method=distMethod))), method=algorithm);
    #clustGap needs a list output
    output <- list()
    output$cluster <- cutree(clustObject,k=k)
    return(output)
    
  }

     gapTest <- clusGap(dataset, FUNcluster=clustF,K.max=maxNumClusters, B = numSims, verbose = interactive());

  }else{
    
        #dist (but not cor) computes across the rows, not columns.
 #dist (but not cor) computes across the rows, not columns.
  clustF <- function(x,k){
       #it turns out that R will recognize the upper-level function input value in this function (algorith, corUse, etc.)
    #tried passing in parent.frame() to make a more elegant solution but didn't work.

    clustObject <- hclust(dist(x,method=distMethod), method=algorithm);
    output <- list()
    output$cluster <- cutree(clustObject,k=k)
    return(output)
  
  }
         #dist (but not cor) computes across the rows, not columns. works better when do t() in outside clustGap function.
    gapTest <- clusGap(t(dataset), FUNcluster=clustF,K.max=maxNumClusters, B = numSims, verbose = interactive());

  }

  #f(k) is the gap statistic.
  # method = c("firstSEmax", "Tibs2001SEmax", "globalSEmax", "firstmax", "globalmax")
  bestK <- maxSE(f=gapTest$Tab[,"gap"], SE.f=gapTest$Tab[,"SE.sim"],
      method = c("Tibs2001SEmax"),
      SE.factor = 1)
  #way to confirm bestK is actually working specifically for Tibs2001SEmax - compute by hand:
  select_metric <- gapTest$Tab[c(1:(maxNumClusters-1)),"gap"] > gapTest$Tab[c(2:maxNumClusters),"gap"] - gapTest$Tab[c(2:maxNumClusters),"SE.sim"];
  #take first one that is true.
  best_k <- which(select_metric)[1];

  if(bestK != best_k){
    
    stop("\nError when selecting best K in clusGap function in cluster package.")
  }

  cat(paste0("\nBest K as determined by k-means and gap test is: ", bestK ,"\n"),
      append=TRUE,file=outputFile);

  clustFeatureIndexList <- list()
  clustSampleIndexList <- list()

  clusterAssignments <- clustF(dataset,bestK)$cluster
  
 for(k in 1:bestK){
   
   clustFeatureIndexList[[k]] <- clustFeatures
   clustSampleIndexList[[k]] <- which(clusterAssignments==k)
   
   
 }

  output <- list(bestK=bestK,clustFeatureIndexList=clustFeatureIndexList,
                 clustSampleIndexList=clustSampleIndexList, gapTest=gapTest)
 
 return(output)
  
}

consensusClusterMatrixList <- function(dataMatrixList,clustFeatures,maxNumClusters = 30, numSims=10, pItem=0.8, pFeature=1, clusterAlg=c("hc","km","pam","kmdist"),
innerLinkage="average", finalLinkage_hclust="average",
corUse=c("everything","pairwise.complete.obs", "complete.obs"),
distMethod=c("pearson","spearman","euclidean", "binary", "maximum", "canberra", "minkowski"),
minClustConsensus=.7){
  
  outputList <- lapply(dataMatrixList,FUN=function(dataMatrix,clustFeatures,maxNumClusters, numSims, pItem, pFeature, clusterAlg,
innerLinkage, finalLinkage_hclust,
corUse,
distMethod,
minClustConsensus){
    
    clusterOutput <- consensusClusterMatrix(dataMatrix,clustFeatures,maxNumClusters, numSims, pItem, pFeature, clusterAlg,
innerLinkage, finalLinkage_hclust,
corUse,
distMethod,
minClustConsensus)
    
    },clustFeatures,maxNumClusters, numSims, pItem, pFeature, clusterAlg,
innerLinkage, finalLinkage_hclust,
corUse,
distMethod,
minClustConsensus)
  
clustSampleIndexList <- list()
clustFeatureIndexList <- list()
bestK <- list()
consensusInfo <- list()

for(d in 1:length(outputList)){
  
  clustSampleIndexList[[d]] <- outputList[[d]]$clustSampleIndexList
  clustFeatureIndexList[[d]] <- outputList[[d]]$clustSampleIndexList
  bestK[[d]] <- outputList[[d]]$bestK
  consensusInfo[[d]] <- outputList[[d]]$consensusCalc
  
}
 
output <- list(clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
               bestK=bestK,consensusInfo=consensusInfo)
return(output)
  
}
#biobase:
#hclustMethod ==finalLinkage
#looks like uses the default Hartigan-Wong for kmeans clustering.
#but that has a pretty low iter.max (10) and nstarts...
#but your own "homemade" cluster functions can only take in a distance matrix,
#so couldn't write a separate kmeans function where I could tweak this.
consensusClusterMatrix <- function(dataMatrix, clustFeatures,maxNumClusters = 30, numSims=10, pItem=0.8, pFeature=1, clusterAlg=c("hc","km","pam","kmdist"),
innerLinkage="average", finalLinkage_hclust="average",
corUse=c("everything","pairwise.complete.obs", "complete.obs"),
distMethod=c("pearson","spearman","euclidean", "binary", "maximum", "canberra", "minkowski"),
minClustConsensus=.7
){

  
  dataset <- dataMatrix[rownames(dataMatrix) %in% clustFeatures, , drop=FALSE]
  
    if(nrow(dataset)==0){
    
    stop("\nIn clusterMatrixHclustGap function: no clustFeatures were found in the data matrix inputted.")
    
  } 
 
  
    consensusClustOutput <- ConsensusClustPlus(d=dataset,
                                               maxK=maxNumClusters,reps=numSims,
                                               pItem=pItem, pFeature=pFeature, clusterAlg=clustF,
  innerLinkage=innerLinkage, finalLinkage=finalLinkage_hclust, ml=NULL,
  tmyPal=NULL,seed=NULL,plot=NULL,writeTable=FALSE,weightsItem=NULL,weightsFeature=NULL,verbose=F,corUse=corUse)
  

  #get consensus score for each cluster for each k: 
  #hmm..doesn't work for k=1?? not so great...
  icl <- calcICL(consensusClustOutput,title=NULL)

  consensusByK <- split(icl[["clusterConsensus"]],f=icl[["clusterConsensus"]][,"k"])
  #if one cluster's consensus is NA: will return NA. but want this.
  #we don't want to pick a K that results in an NaN consensus value; 
  #this probably means should stick with a lower k.
  minConsensusClusterByK <- lapply(consensusByK,FUN=function(consensusUnit){min(consensusUnit)})
  #starts at k=2
  if(length(which(unlist(minConsensusClusterByK)>=minClustConsensus))>0){
  
    #want to pick highest K we can (lowest K would basically always return 2.)
    bestK <- max(which(unlist(minConsensusClusterByK)>=minClustConsensus)) +1

  }else{
    #none passed threshold...let k=1?
    bestK <- 1
  }

    clustFeatureIndexList <- list()
  clustSampleIndexList <- list()

  #COME BACK: is this the correct list index name?
  clusterAssignments <- consensusClustOutput$cluster
  
 for(k in 1:bestK){
   
   clustFeatureIndexList[[k]] <- clustFeatureIndex
   clustSampleIndexList[[k]] <- which(clusterAssignments==k)
   
   
 }
  output <- list(bestK=bestK,clustFeatureIndexList=clustFeatureIndexList,
                 clustSampleIndexList=clustSampleIndexList, consensusCalc=icl,
                 consensusClustOutput=consensusClustOutput)

 return(output)
 
}