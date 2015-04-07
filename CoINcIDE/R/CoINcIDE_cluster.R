#use default Ward's method
library("cluster")
library("ConsensusClusterPlus")

cluster_kmeansGap <- function(dataMatrix,clustFeatureIndex,maxNumClusters=30,iter.max=30,nstart=25,numSims=1000,algorithm="Hartigan-Wong",outputFile="./cluster_kmeansGap_output.txt"){
  
  if(any(clustFeatureIndex)>nrow(dataset)){
    
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


  
}


cluster_hclustGap <- function(dataMatrix,clustFeatureIndex,maxNumClusters=30,algorithm="complete",
                              distMethod=c("cor","euclidean", "maximum", "manhattan", "canberra", "binary" , "minkowski"),outputFile="./cluster_kmeansGap_output.txt"){
  
  if(any(clustFeatureIndex)>nrow(dataset)){
    
    stop("\nIn cluster_kmeansGap: your clustFeatureIndex list contains at least one indice that is beyond the number of rows in your input dataMatrix.\n")
  }
  
  dataset <- dataMatrix[clustFeatureIndex, , drop=FALSE]
  
  if(distMethod=="cor"){
    
  clustF <- function(dataset,k,env=parent.frame()){
  
    clustObject <- hclust(as.dist((1-cor(dataset))), method=env$algorithm);
    clustObject <- cutree(clustObject,k=k)
    #must be in list format for cluster package to recognize.
    output <- list()
    output$cluster <- clustObject
    return(output)
    
  }

  }else{
    
  clustF <- function(dataset,k,env=parent.frame()){
  
    clustObject <- hclust(dist(dataset,method=env$distMethod), method=env$algorithm);
    clustObject <- cutree(clustObject,k=k)
    #must be in list format for cluster package to recognize.
    output <- list()
    output$cluster <- clustObject
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

  
}


#biobase:
consensusClusterPlus