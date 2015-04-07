#use default Ward's method
library("cluster")
library("ConsensusClusterPlus")

cluster_kmeansGap <- function(dataMatrix,clustFeatureIndex,maxNumClusters=30,iter.max=30,nstart=25,numSims=1000,algorithm="Hartigan-Wong",outputFile="./cluster_kmeansGap_output.txt"){
  
  if(any(clustFeatureIndex)>nrow(dataMatrix)){
    
    stop("\nIn cluster_kmeansGap: your clustFeatureIndex list contains at least one indice that is beyond the number of rows in your input dataMatrix.\n")
  }
  
  dataset <- dataMatrix[clustFeatureIndex, , drop=FALSE]
  
  clustF <- function(dataset,k,env=parent.frame()){
  
    kmeans(dataset, centers=k,iter.max=env$iter.max,nstart=env$nstart,algorithm=env$algorithm);
  
  }

  gapTest <- clusGap(dataset, FUNcluster=clustF,K.max=maxNumClusters, B = numSims, verbose = interactive());

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
   
   clustFeatureIndexList[[k]] <- clustFeatureIndex
   clustSampleIndexList[[k]] <- which(clusterAssignments==k)
   
   
 }

  output <- list(bestK=bestK,clustFeatureIndexList=clustFeatureIndexList,
                 clustSampleIndexList=clustSampleIndexList, gapTest= gapTest)
 return(output)
  
}


cluster_hclustGap <- function(dataMatrix,clustFeatureIndex,maxNumClusters=30,algorithm="complete",
                              distMethod=c("pearson","spearman","euclidean", "binary", "maximum", "canberra", "minkowski"),outputFile="./cluster_kmeansGap_output.txt",
                              corUse=c("everything","pairwise.complete.obs", "complete.obs")){
  
  if(any(clustFeatureIndex)>nrow(dataMatrix)){
    
    stop("\nIn cluster_kmeansGap: your clustFeatureIndex list contains at least one indice that is beyond the number of rows in your input dataMatrix.\n")
  }
  
  dataset <- dataMatrix[clustFeatureIndex, , drop=FALSE]
  
  if(distMethod==("pearson") || distMethod=="spearman"){
    
  clustF <- function(dataset,k,env=parent.frame()){
  
    clustObject <- hclust(as.dist((1-cor(dataset,use=corUse,method=distMethod))), method=env$algorithm);
    clustObject <- cutree(clustObject,k=k)
    #must be in list format for cluster package to recognize.
    output <- list()
    output$cluster <- clustObject
    names(output$cluster) <- colnames(dataset)
    return(output)
    
  }

  }else{
    
  clustF <- function(dataset,k,env=parent.frame()){
  
    clustObject <- hclust(dist(dataset,method=env$distMethod), method=env$algorithm);
    #ConsensusClustPlus works with cutree output (see vignette.)
    output <- cutree(clustObject,k=k)
    #must be in list format for cluster package to recognize.
    output <- list()
    output$cluster <- clustObject
    names(output$cluster) <- colnames(dataset)
    return(output)
  }
    
  }
  gapTest <- clusGap(dataset, FUNcluster=clustF(distMethod=distMethod,algorithm=algorithm),K.max=maxNumClusters, B = numSims, verbose = interactive());

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
   
   clustFeatureIndexList[[k]] <- clustFeatureIndex
   clustSampleIndexList[[k]] <- which(clusterAssignments==k)
   
   
 }

   output <- list(bestK=bestK,clustFeatureIndexList=clustFeatureIndexList,
                 clustSampleIndexList=clustSampleIndexList, gapTest= gapTest)
 return(output)
  
}


#biobase:
#hclustMethod ==finalLinkage
#looks like uses the default Hartigan-Wong for kmeans clustering.
#but that has a pretty low iter.max (10) and nstarts...
#but your own "homemade" cluster functions can only take in a distance matrix,
#so couldn't write a separate kmeans function where I could tweak this.
consensusClusterPlus_wrapper(dataMatrix, clustFeatureIndex,maxNumClusters = 30, numSims=10, pItem=0.8, pFeature=1, clusterAlg=c("hc","km","pam","kmdist"),
innerLinkage="average", finalLinkage_hclust="average",
corUse=c("everything","pairwise.complete.obs", "complete.obs"),
distMethod=c("pearson","spearman","euclidean", "binary", "maximum", "canberra", "minkowski"),
minClustConsensus=.7
){

  if(any(clustFeatureIndex)>nrow(dataMatrix)){
    
    stop("\nIn cluster_kmeansGap: your clustFeatureIndex list contains at least one indice that is beyond the number of rows in your input dataMatrix.\n")
  }
  
  dataset <- dataMatrix[clustFeatureIndex, , drop=FALSE]
 
  
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

  output <- list(bestK=bestK,clustFeatureIndexList=clustFeatureIndexList,
                 clustSampleIndexList=clustSampleIndexList, consensusCalc=icl,
                 consensusClustOutput=consensusClustOutput)

 return(output)
 
}