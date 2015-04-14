#use default Ward's method
library("cluster")
library("ConsensusClusterPlus")

#for more specific options: use the individual clustMatrixList functions
clustMatrixListWrapper <- function(dataMatrixList,clustFeaturesList,clustMethod=c("km","hc"),pickKMethod=c("gap","consensus"),numSims=1000,maxNumClusters=30,
                                   outputFile="./cluster_output.txt",iter.max=30,nstart=25,distMethod=c("pearson","spearman","euclidean", "binary", "maximum", "canberra", "minkowski"),
                                   hclustAlgorithm=c("average","complete","ward.D", "ward.D2", "single", "mcquitty","median","centroid"), 
                                   consensusHclustAlgorithm=c("average","complete","ward.D", "ward.D2", "single", "mcquitty","median","centroid"),
                                   bestKconsensusMethod=c("highestMinConsensus","bestKOverThreshBeforeNAclust","maxBestKOverThresh"),minClustConsensus=.9){
  
  if(pickKMethod=="gap"){
    
    if(clustMethod=="kmeans"){
      
      clusterOutputList <- clustMatrixListKmeansGap(dataMatrixList=dataMatrixList,clustFeaturesList=clustFeaturesList,iter.max=iter.max,nstart=nstart,maxNumClusters=maxNumClusters,numSims=numSims,outputFile=outputFile)
        
    }else if(clustMethod=="hc"){
      
      clusterOutputList <- clusterMatrixListHclustGap(dataMatrixList=dataMatrixList,clustFeaturesList=clustFeaturesList,maxNumClusters=maxNumClusters,algorithm=hclustAlgorithm,
                                                              distMethod=distMethod,outputFile=outputFile,
                                                              corUse=corUse,numSims=numSims)      
    }else{
      
      stop("In clustMatrixListWrapper: did not pick hc or kmeans as clustMethod input")
    }
    
  }else if(pickKMethod=="consensus"){
    
    if(clustMethod=="kmeans"){
      
      clusterAlg <- "km"
      
    }else{
      
      clusterAlg <- clustMethod
      
    }
    
    if(length( bestKconsensusMethod)>1){
      
      warning("In clustMatrixListWrapper: bestKconsensusMethod longer than 1; picking default \"highestMinConsensus\"")
      bestKconsensusMethod <- "highestMinConsensus"
      
    }
    
    clusterOutputList <- consensusClusterMatrixList(dataMatrixList=dataMatrixList,clustFeaturesList=clustFeaturesList,maxNumClusters = maxNumClusters, numSims=numSims, pItem=0.8, pFeature=1, clusterAlg=clusterAlg,
                                           hclustAlgorithm=hclustAlgorithm, consensusHclustAlgorithm=consensusHclustAlgorithm,
                                           corUse=corUse,
                                           distMethod=distMethod,
                                           minClustConsensus=minClustConsensus,bestKmethod=bestKconsensusMethod,
                                           outputFile=outputFile)
      
  }else if(is.numeric(pickKMethod)){
    
      
      clusterOutputList <- list()
      
      for(c in 1:length(dataMatrixList)){
        
        clusterOutputList[[c]] <- list()
        dataset <- dataMatrixList[[c]][rownames(dataMatrixList[[c]]) %in% clustFeaturesList[[c]], , drop=FALSE]
        
        if(nrow(dataset)==0){
          
          stop("\nI function clustMatrixListWrapper: no clustFeatures were found in the data matrix inputted.")
          
        }
        
        
        if(ncol(dataset)<pickKMethod){
          #hclust usually returns NA gap test if K.max = ncol(dataset)
          #will also mest up maxSE calculations
          K <- ncol(dataset)-1
          
        }else{
          
          K <- pickKMethod
          
        }

        if(clustMethod=="kmeans"){
          
        clusterAssignments <-   kmeans(t(dataset), centers=K,iter.max=iter.max,nstart=nstart,algorithm=algorithm)$cluster;
        
        clustFeatureIndexList <- list()
        clustSampleIndexList <- list()

        
        }else if(clustMethod=="hc"){
          

          
          if(distMethod==("pearson") || distMethod=="spearman"){
            
              #for cor: don't transpose dataset.
              #it turns out that R will recognize the upper-level function input value in this function (algorith, corUse, etc.)
              clustObject <- hclust(as.dist((1-cor(dataset,use=corUse,method=distMethod))), method=algorithm);
              #clustGap needs a list output
              clusterAssignments <- cutree(clustObject,k=K)

          }else{

              #it turns out that R will recognize the upper-level function input value in this function (algorith, corUse, etc.)
              #tried passing in parent.frame() to make a more elegant solution but didn't work.
              
              clustObject <- hclust(dist(t(dataset),method=distMethod), method=algorithm);
              clusterAssignments <- cutree(clustObject,k=K)

        
          }
        }else{
          
          stop("In clustMatrixListWrapper: did not pick hc or kmeans as clustMethod input")
          
        }
        
        
        clustFeatureIndexList <- list()
        clustSampleIndexList <- list()
        
        for(k in 1:K){
          
          clustFeatureIndexList[[k]] <- clustFeaturesList[[c]]
          clustSampleIndexList[[k]] <- which(clusterAssignments==k)
          
          
        }
        clusterOutputList[[c]]$clustFeatureIndexList <- clustFeatureIndexList
        clusterOutputList[[c]]$clustSampleIndexList <- clustSampleIndexList

        
      }                                            
    
    
  }else{
    
    stop("In clustMatrixListWrapper: did not pick one of the allowed input character strings, or a number, for pickKMethod")
  }
  
  return(clusterOutputList)
  
}

clustMatrixListKmeansGap <- function(dataMatrixList,clustFeaturesList,maxNumClusters=30,iter.max=30,nstart=25,numSims=1000,algorithm="Hartigan-Wong",outputFile="./cluster_kmeansGap_output.txt"
                                     ){
  #lapply doens't work because also need to loop over clust features list.
  #mapply didn't quite work either.
  outputList <- list()
  clustSampleIndexList <- list()
  clustFeatureIndexList <- list()
  bestK <- list()
  gapTest <- list()
  for(d in 1:length(dataMatrixList)){
    
    message(paste0("\nClustering dataset ",d, ": ",names(dataMatrixList)[d]))
    outputList[[d]] <- clusterMatrixKmeansGap(dataMatrix=dataMatrixList[[d]],clustFeatures=clustFeaturesList[[d]],
                                          maxNumClusters=maxNumClusters,iter.max=iter.max,
                                          nstart=nstart,algorithm=algorithm,numSims=numSims,outputFile=outputFile)
    
    clustSampleIndexList[[d]] <- outputList[[d]]$clustSampleIndexList
    clustFeatureIndexList[[d]] <- outputList[[d]]$clustFeatureIndexList
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
  
  if(ncol(dataset)<maxNumClusters){
    #clustGap may return NA gap test if K.max = ncol(dataset)
    #this will also mest up maxSE calculations
    K.max <- ncol(dataset) - 1
    
  }else{
    
    K.max <- maxNumClusters
    
  }
  
  gapTest <- clusGap(t(dataset), FUNcluster=clustF,K.max=K.max, B = numSims, verbose = interactive());

  #f(k) is the gap statistic.
  # method = c("firstSEmax", "Tibs2001SEmax", "globalSEmax", "firstmax", "globalmax")

  bestK <- maxSE(f=gapTest$Tab[,"gap"], SE.f=gapTest$Tab[,"SE.sim"],
      method = c("Tibs2001SEmax"),
      SE.factor = 1)
  #way to confirm bestK is actually working specifically for Tibs2001SEmax - compute by hand:
  select_metric <- gapTest$Tab[c(1:(K.max-1)),"gap"] > gapTest$Tab[c(2:K.max),"gap"] - gapTest$Tab[c(2:K.max),"SE.sim"]
  #take first one that is true.
  best_k <- which(select_metric)[1]
  
  
  if(is.na(best_k)){
    
    warning(paste0("\nDid not find an optimal k below ",maxNumClusters," setting K=1\n"))
    cat(paste0("\nDid not find an optimal k below ",maxNumClusters," setting K=1\n"),
        append=TRUE,file=outputFile)
    best_k <- 1
  
  }
 
  
  if(bestK != best_k){
    
    stop("\nError when selecting best K in clusGap function in cluster package.")
  }

  cat(paste0("\nBest K as determined by k-means and gap test is: ", bestK ,"\n"),
      append=TRUE,file=outputFile);

  clustFeatureIndexList <- list()
  clustSampleIndexList <- list()

  #remember: transpose!
   clusterAssignments <- clustF(t(dataset),bestK)$cluster
  
 for(k in 1:bestK){
   
   clustFeatureIndexList[[k]] <-  na.omit(match(clustFeatures,rownames(dataMatrix)))
   clustSampleIndexList[[k]] <- which(clusterAssignments==k)
   
   
 }

  output <- list(bestK=bestK,clustFeatureIndexList=clustFeatureIndexList,
                 clustSampleIndexList=clustSampleIndexList, gapTest= gapTest)
 return(output)
  
}

clusterMatrixListHclustGap <- function(dataMatrixList,clustFeaturesList,maxNumClusters=30,algorithm="ward.D",
                              distMethod=c("pearson","spearman","euclidean", "binary", "maximum", "canberra", "minkowski"),outputFile="./cluster_hclustGap_output.txt",
                              corUse=c("everything","pairwise.complete.obs", "complete.obs"),numSims=1000){
  
  outputList <- list()
  clustSampleIndexList <- list()
  clustFeatureIndexList <- list()
  bestK <- list()
  gapTest <- list()

  for(d in 1:length(dataMatrixList)){
      message(paste0("\nClustering dataset ",d, ": ",names(dataMatrixList)[d]))
   outputList[[d]] <- clusterMatrixHclustGap(dataMatrix=dataMatrixList[[d]],clustFeatures=clustFeaturesList[[d]],maxNumClusters=maxNumClusters,algorithm=algorithm,numSims=numSims,distMethod=distMethod,outputFile=outputFile,corUse=corUse)
        clustSampleIndexList[[d]] <- outputList[[d]]$clustSampleIndexList
    clustFeatureIndexList[[d]] <- outputList[[d]]$clustFeatureIndexList
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
  
  if(ncol(dataset)<maxNumClusters){
    #hclust usually returns NA gap test if K.max = ncol(dataset)
    #will also mest up maxSE calculations
    K.max <- ncol(dataset)-1
    
  }else{
    
    K.max <- maxNumClusters
    
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
  
     gapTest <- clusGap(dataset, FUNcluster=clustF,K.max=K.max, B = numSims, verbose = interactive());

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
    gapTest <- clusGap(t(dataset), FUNcluster=clustF,K.max=K.max, B = numSims, verbose = interactive());

  }

  #f(k) is the gap statistic.
  # method = c("firstSEmax", "Tibs2001SEmax", "globalSEmax", "firstmax", "globalmax")
  bestK <- maxSE(f=gapTest$Tab[,"gap"], SE.f=gapTest$Tab[,"SE.sim"],
      method = c("Tibs2001SEmax"),
      SE.factor = 1)
  #way to confirm bestK is actually working specifically for Tibs2001SEmax - compute by hand:
  select_metric <- gapTest$Tab[c(1:(K.max-1)),"gap"] > gapTest$Tab[c(2:K.max),"gap"] - gapTest$Tab[c(2:K.max),"SE.sim"];
  #take first one that is true.
  best_k <- which(select_metric)[1];

  if(is.na(best_k)){
    
    warning(paste0("\nDid not find an optimal k below ",maxNumClusters," setting K=1\n"))
    cat(paste0("\nDid not find an optimal k below ",maxNumClusters," setting K=1\n"),
        append=TRUE,file=outputFile)
    best_k <- 1
  }
 
  if(bestK != best_k){
    
    stop("\nError when selecting best K in clusGap function in cluster package.")
  }

  cat(paste0("\nBest K as determined by hclust and gap test is: ", bestK ,"\n"),
      append=TRUE,file=outputFile);

  clustFeatureIndexList <- list()
  clustSampleIndexList <- list()

 if(distMethod==("pearson") || distMethod=="spearman"){
   
  clusterAssignments <- clustF(dataset,bestK)$cluster
  
  }else{
    #mus tranpose dataset for traditional distance matrix hclustering
    clusterAssignments <- clustF(t(dataset),bestK)$cluster
  }
  
 for(k in 1:bestK){
   
   clustFeatureIndexList[[k]] <-  na.omit(match(clustFeatures,rownames(dataMatrix)))
   clustSampleIndexList[[k]] <- which(clusterAssignments==k)
   
   
 }

  output <- list(bestK=bestK,clustFeatureIndexList=clustFeatureIndexList,
                 clustSampleIndexList=clustSampleIndexList, gapTest=gapTest)
 
 return(output)
  
}

consensusClusterMatrixList <- function(dataMatrixList,clustFeaturesList,maxNumClusters = 30, numSims=10, pItem=0.8, pFeature=1, clusterAlg=c("hc","km","pam","kmdist"),
                                       hclustAlgorithm=c("average","complete","ward.D", "ward.D2", "single", "mcquitty","median","centroid"), consensusHclustAlgorithm=c("average","complete","ward.D", "ward.D2", "single", "mcquitty","median","centroid"),
corUse=c("everything","pairwise.complete.obs", "complete.obs"),
distMethod=c("pearson","spearman","euclidean", "binary", "maximum", "canberra", "minkowski"),
minClustConsensus=.7,bestKmethod=c("highestMinConsensus","bestKOverThreshBeforeNAclust","maxBestKOverThresh"),
outputFile="./consensusOut.txt"){
  
 
  outputList <- list()
  
  for(d in 1:length(dataMatrixList)){
    
    message(paste0("\nClustering dataset ",d, " ",names(dataMatrixList)[d]))
    cat(paste0("\nClustering dataset ",d, " ",names(dataMatrixList)[d]),
        append=TRUE,file=outputFile);
    outputList[[d]] <- consensusClusterMatrix(dataMatrixList[[d]],clustFeaturesList[[d]],
                                                     maxNumClusters=maxNumClusters, numSims=numSims, 
                                                     pItem=pItem, pFeature=pFeature, clusterAlg=clusterAlg,
                                              hclustAlgorithm=hclustAlgorithm, consensusHclustAlgorithm=consensusHclustAlgorithm,corUse=corUse,distMethod=distMethod,
  minClustConsensus=minClustConsensus,bestKmethod=bestKmethod,outputFile=outputFile)
    
  }
  
  clustSampleIndexList <- list()
  clustFeatureIndexList <- list()
  bestK <- list()
  consensusInfo <- list()
  minConsensusClusterByK <- list()
  
  for(d in 1:length(outputList)){
    
    clustSampleIndexList[[d]] <- outputList[[d]]$clustSampleIndexList
    clustFeatureIndexList[[d]] <- outputList[[d]]$clustFeatureIndexList
    bestK[[d]] <- outputList[[d]]$bestK
    consensusInfo[[d]] <- outputList[[d]]$consensusCalc
    minConsensusClusterByK[[d]] <- outputList[[d]]$minConsensusClusterByK
    
  }
   
  output <- list(clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                 bestK=bestK,consensusInfo=consensusInfo,minConsensusClusterByK=minConsensusClusterByK)
  return(output)
  
}
#biobase:
#hclustMethod ==finalLinkage
#looks like uses the default Hartigan-Wong for kmeans clustering.
#but that has a pretty low iter.max (10) and nstarts...
#but your own "homemade" cluster functions can only take in a distance matrix,
#so couldn't write a separate kmeans function where I could tweak this.

#bestKOverThreshBeforeNAclust: what the is largest k value whose min cluster consensus passes the minClustConsensus threshold BEFORE the first instance of a cluster with NAs (presumably "too many" clusters then.)
#maxBestKOverThresh: what the is largest k value whose min cluster consensus passes the minClustConsensus threshold?
#highestMinConsensus: what is the k value with the maximum min cluster consensus value?

#warning: can't really silence the ConsensusClusterPlus plotting functions,
#so testing a large maxNumClusters gets slow.
consensusClusterMatrix <- function(dataMatrix, clustFeatures,maxNumClusters = 30, numSims=10, pItem=0.8, pFeature=1, clusterAlg=c("km","hc","pam","kmdist"),
                                   hclustAlgorithm=c("average","complete","ward.D", "ward.D2", "single", "mcquitty","median","centroid"), consensusHclustAlgorithm=c("average","complete","ward.D", "ward.D2", "single", "mcquitty","median","centroid"),
corUse=c("everything","pairwise.complete.obs", "complete.obs"),
distMethod=c("pearson","spearman","euclidean", "binary", "maximum", "canberra", "minkowski"),
minClustConsensus=.8,bestKmethod=c("highestMinConsensus","bestKOverThreshBeforeNAclust","maxBestKOverThresh"),
outputFile="./consensusOut.txt"
){
  
 innerLinkage <- hclustAlgorithm
 finalLinkage <- consensusHclustAlgorithm
 
 if(length(innerLinkage)>1){
   
   warning("Multiple inputs for innerLinkage used so setting to defaule \'average\'")
   innerLinkage <- "average"
 }
 
 if(length(finalLinkage)>1){
   
   warning("Multiple inputs for innerLinkage used so setting to defaule \'average\'")
   finalLinkage <- "average"
   
 }
  #come back: add defaults to make more user-friendly.
  #if(length(clusterAlg)>1 || length(bestKmethod)>1)
   #ConsensusClustPlus: assumes clustering the columns.
  
  dataset <- dataMatrix[rownames(dataMatrix) %in% clustFeatures, , drop=FALSE]
  
    if(nrow(dataset)==0){
    
    stop("\nIn clusterMatrixHclustGap function: no clustFeatures were found in the data matrix inputted.")
    
  } 
  
  #remember: we're clustering pItem*ncol(dataset) each time...so
  #our max K is NOT ncol(dataset), but rather pItem*ncol(dataset)
  if(floor(ncol(dataset)*pItem)<maxNumClusters){
    #hclust usually returns NA gap test if K.max = pItem*ncol(dataset) as opposed to pItem*ncol(dataset)-1
    #will also mest up maxSE calculations
    K.max <- floor(pItem*ncol(dataset))-1
    
  }else{
    
    K.max <- maxNumClusters
    
  }

    consensusClustOutput <- ConsensusClusterPlus(d=dataset,
                                               maxK=K.max,reps=numSims,
                                               pItem=pItem, pFeature=pFeature, clusterAlg=clusterAlg,
  innerLinkage=innerLinkage, finalLinkage=finalLinkage, ml=NULL,
  tmyPal=NULL,seed=NULL,plot=FALSE,writeTable=FALSE,weightsItem=NULL,weightsFeature=NULL,verbose=F,corUse=corUse)
  

  #get consensus score for each cluster for each k: 
  #hmm..doesn't work for k=1?? not so great...
  #weird...can't NOT have the plots outputted...tried title=NULL
  icl <- calcICL(consensusClustOutput,plot=FALSE)

  consensusByK <- split(icl[["clusterConsensus"]],f=icl[["clusterConsensus"]][,"k"])
  #if one cluster's consensus is NA: will return NA. but want this.
  #we don't want to pick a K that results in an NaN consensus value; 
  #this probably means should stick with a lower k.
  minConsensusClusterByK <- lapply(consensusByK,FUN=function(consensusUnit){min(consensusUnit)})
  #starts at k=2
  if(length(which(unlist(minConsensusClusterByK)>=minClustConsensus))>0){

    if(bestKmethod=="bestKOverThreshBeforeNAclust"){  
    #want to pick highest K we can before an NA? (lowest K would basically always return 2.)
    #take ONE step back then from this NA_limit (-1 part) 
    NA_limit <- which(is.na(minConsensusClusterByK))[1] -1
    bestK <- max(which(unlist(minConsensusClusterByK)[c(1:NA_limit)]>=minClustConsensus)) +1
    
    }else if(bestKmethod=="maxBestKOverThresh"){
      
         bestK <- max(which(unlist(minConsensusClusterByK)[c(1:length(minConsensusClusterByK))]>=minClustConsensus)) +1 
      
    }else if(bestKmethod=="highestMinConsensus"){
      
       bestK <- which.max(minConsensusClusterByK) + 1
      
    }

  }else{
    #none passed threshold...let k=1?
    bestK <- 1
  }

  message(paste0("Best K is: ",bestK))
  cat(paste0("\nBest K as determined by ",clusterAlg, " and consensus metric  is: ", bestK ,"\n"),
      append=TRUE,file=outputFile);
    clustFeatureIndexList <- list()
  clustSampleIndexList <- list()

  if(bestK >1){

  clusterAssignments <- consensusClustOutput[[bestK]]$consensusClass
  
  }else{
    
    clusterAssignments <- rep.int(1,times=ncol(dataset))
  }
  
 for(k in 1:bestK){
   
   clustFeatureIndexList[[k]] <- na.omit(match(clustFeatures,rownames(dataMatrix)))
   clustSampleIndexList[[k]] <- which(clusterAssignments==k)
   
   
 }
  output <- list(bestK=bestK,clustFeatureIndexList=clustFeatureIndexList,minConsensusClusterByK=minConsensusClusterByK,
                 clustSampleIndexList=clustSampleIndexList, consensusCalc=icl,
                 consensusClustOutput=consensusClustOutput)

 return(output)
 
}