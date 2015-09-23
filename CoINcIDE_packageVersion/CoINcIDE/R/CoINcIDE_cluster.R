#use default Ward's method
#library("cluster")
#proxy not working...
#library("proxy")

#for more specific options: use the individual clustMatrixList functions
#note: for consenus clustering: need a high # sims otherwise will get NaN values in your consensus matrix (two samples were next chosen in the same resampling run.)
clustMatrixListWrapper <- function(dataMatrixList,clustFeaturesList,clustMethod=c("km","hc"),pickKMethod=c("gap","consensus"),numSims=1000,maxNumClusters=30,
                                   outputFile="./cluster_output.txt",iter.max=10,nstart=1,distMethod=c("euclidean","pearson","spearman", "binary", "maximum", "canberra", "minkowski"),
                                   hclustAlgorithm=c("average","complete","ward.D", "ward.D2", "single", "mcquitty","median","centroid"), 
                                   consensusHclustAlgorithm=c("average","complete","ward.D", "ward.D2", "single", "mcquitty","median","centroid"),
                                   minClustConsensus=.7, minMeanClustConsensus=.7,maxPAC=.1,corUse="everything",pItem=.9){
  
  
  message("This code assumes your patients are in the columns of the data matrix.")
  
  if(pickKMethod=="gap"){
    
    if(clustMethod=="km"){
      

      clusterOutputList <- clustMatrixListKmeansGap(dataMatrixList=dataMatrixList,clustFeaturesList=clustFeaturesList,iter.max=iter.max,nstart=nstart,maxNumClusters=maxNumClusters,numSims=numSims,outputFile=outputFile)
        
    }else if(clustMethod=="hc"){
      
      clusterOutputList <- clusterMatrixListHclustGap(dataMatrixList=dataMatrixList,clustFeaturesList=clustFeaturesList,maxNumClusters=maxNumClusters,algorithm=hclustAlgorithm,
                                                              distMethod=distMethod,outputFile=outputFile,
                                                              corUse=corUse,numSims=numSims)      
    }else{
      
      stop("In clustMatrixListWrapper: did not pick hc or km as clustMethod input")
    }
    
  }else if(pickKMethod=="consensus"){
      
      clusterAlg <- clustMethod
    
    clusterOutputList <- consensusClusterMatrixList(dataMatrixList=dataMatrixList,clustFeaturesList=clustFeaturesList,maxNumClusters = maxNumClusters, 
                                                    numSims=numSims, clusterAlg=clusterAlg,
                                           hclustAlgorithm=hclustAlgorithm, consensusHclustAlgorithm=consensusHclustAlgorithm,
                                           corUse=corUse,pItem=pItem,pFeature=1,iter.max=iter.max,nstart=nstart,
                                           distMethod=distMethod,minMeanClustConsensus=minMeanClustConsensus,
                                           minClustConsensus=minClustConsensus,
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

            if(distMethod=="euclidean"){
              #capitalize
              distMethod <- "Euclidean"
              
            }
                          #convert similarities: does 1- similarity.
              distMatrix <- dist(dataset,method=distMethod,by_rows=FALSE,
                       convert_similarities = TRUE)
              #it turns out that R will recognize the upper-level function input value in this function (algorith, corUse, etc.)
              #tried passing in parent.frame() to make a more elegant solution but didn't work.
              
              clustObject <- hclust(distMatrix, method=algorithm);
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
    bestK <- 1
  }else{
 
  if(bestK != best_k){
    
    stop("\nError when selecting best K in clusGap function in cluster package.")
  }

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
                              distMethod=c("pearson","spearman","Euclidean", "Pearson", "cosine", "Manhattan", "Minkowski"),outputFile="./cluster_hclustGap_output.txt",
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
                              distMethod=c("pearson","spearman","Euclidean", "cosine", "Manhattan", "Minkowski"),outputFile="./cluster_kmeansGap_output.txt",
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
  
    if(distMethod=="euclidean"){
      #capitalize
      distMethod <- "Euclidean"
              
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
    
    #convert similarities: does 1- similarity.
    distMatrix <- dist(dataset,method=distMethod,by_rows=FALSE,
                       convert_similarities = TRUE)
       
    clustObject <- hclust(distMatrix, method=algorithm);
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
    bestK <- 1
    
  }else{
 
  if(bestK != best_k){
    
    stop("\nError when selecting best K in clusGap function in cluster package.")
  }
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
distMethod=c("pearson","spearman","euclidean", "binary", "maximum", "canberra", "minkowski"),minMeanClustConsensus=.7,minClustConsensus=.7,
outputFile="./consensusOut.txt",iter.max=30,nstart=25,maxPAC=.15){
  
 
  outputList <- list()
  
  for(d in 1:length(dataMatrixList)){
    
    message(paste0("\nClustering dataset ",d, " ",names(dataMatrixList)[d]))
    cat(paste0("\nClustering dataset ",d, " ",names(dataMatrixList)[d]),
        append=TRUE,file=outputFile);
    

    outputList[[d]] <- consensusClusterMatrix(dataMatrix=dataMatrixList[[d]],clustFeatures=clustFeaturesList[[d]],
                                                     maxNumClusters=maxNumClusters, numSims=numSims, nstart=nstart,iter.max=iter.max,
                                                     pItem=pItem, pFeature=pFeature, clusterAlg=clusterAlg,
                                              hclustAlgorithm=hclustAlgorithm, consensusHclustAlgorithm=consensusHclustAlgorithm,
                                              corUse=corUse,distMethod=distMethod,
                                              minMeanClustConsensus=minMeanClustConsensus,
                                              minClustConsensus=minClustConsensus,outputFile=outputFile,
                                              maxPAC=maxPAC)
    

  
}
  clustSampleIndexList_consensusFrac <- list()
  clustFeatureIndexList_consensusFrac <- list()
  bestK_consensusFrac <- list()
  clustSampleIndexList_meanConsensusCluster <- list()
  clustFeatureIndexList_meanConsensusCluster <- list()
  bestK_meanConsensusCluster <- list()
  clustSampleIndexList_minConsensusCluster <- list()
  clustFeatureIndexList_minConsensusCluster <- list()
  bestK_minConsensusCluster <- list()
  bestK_PAC <- list()
  clustSampleIndexList_PAC <- list()
  clustFeatureIndexList_PAC <- list()
bestK_PACR <- list()
clustSampleIndexList_PACR <- list()
clustFeatureIndexList_PACR <- list()
  consensusInfo <- list()
  minConsensusClusterByK <- list()
  meanConsensusClusterByK <- list()
  PAC <- list()
PACR <- list()
consensusClustOutput <- list()
consensusFracByK <- list()

for(d in 1:length(outputList)){
  
  clustSampleIndexList_consensusFrac[[d]] <- outputList[[d]]$selectK_consensusFrac$clustSampleIndexList
  clustFeatureIndexList_consensusFrac[[d]] <- outputList[[d]]$selectK_consensusFrac$clustFeatureIndexList
  bestK_consensusFrac[[d]] <- outputList[[d]]$selectK_consensusFrac$bestK
  clustSampleIndexList_meanConsensusCluster[[d]] <- outputList[[d]]$selectK_meanConsensusClust$clustSampleIndexList
  clustFeatureIndexList_meanConsensusCluster[[d]] <- outputList[[d]]$selectK_meanConsensusClust$clustFeatureIndexList
  bestK_meanConsensusCluster[[d]] <- outputList[[d]]$selectK_meanConsensusClust$bestK
  clustSampleIndexList_minConsensusCluster[[d]] <- outputList[[d]]$selectK_minConsensusClust$clustSampleIndexList
  clustFeatureIndexList_minConsensusCluster[[d]] <- outputList[[d]]$selectK_minConsensusClust$clustFeatureIndexList
  bestK_minConsensusCluster[[d]] <- outputList[[d]]$selectK_minConsensusClust$bestK
  bestK_PAC[[d]] <- outputList[[d]]$selectK_PAC$bestK
  clustSampleIndexList_PAC[[d]] <- outputList[[d]]$selectK_PAC$clustSampleIndexList
  clustFeatureIndexList_PAC[[d]] <- outputList[[d]]$selectK_PAC$clustFeatureIndexList
  PAC[[d]] <- outputList[[d]]$PAC
  clustSampleIndexList_PACR[[d]] <- outputList[[d]]$selectK_PACR$clustSampleIndexList
  clustFeatureIndexList_PACR[[d]] <- outputList[[d]]$selectK_PACR$clustFeatureIndexList
  PACR[[d]] <- outputList[[d]]$PACR
  bestK_PACR[[d]] <- outputList[[d]]$selectK_PACR$bestK
  consensusInfo[[d]] <- outputList[[d]]$consensusCalc
  minConsensusClusterByK[[d]] <- outputList[[d]]$minConsensusClusterByK
  meanConsensusClusterByK[[d]] <- outputList[[d]]$meanConsensusClusterByK
  consensusClustOutput[[d]] <- outputList[[d]]$consensusClustOutput
  consensusFracByK[[d]] <- outputList[[d]]$consensusFrac
  
}
   
names(clustSampleIndexList_consensusFrac) <- names(dataMatrixList)
names(clustFeatureIndexList_consensusFrac) <- names(dataMatrixList)
names(bestK_consensusFrac) <- names(dataMatrixList)
names(clustSampleIndexList_meanConsensusCluster) <- names(dataMatrixList)
names(clustFeatureIndexList_meanConsensusCluster) <- names(dataMatrixList)
names(bestK_meanConsensusCluster) <- names(dataMatrixList)
names(clustSampleIndexList_minConsensusCluster) <- names(dataMatrixList)
names(clustFeatureIndexList_minConsensusCluster) <- names(dataMatrixList)
names(bestK_minConsensusCluster) <- names(dataMatrixList)
names(clustSampleIndexList_PAC) <- names(dataMatrixList)
names(clustFeatureIndexList_PAC) <- names(dataMatrixList)
names(bestK_PAC) <- names(dataMatrixList)
names(PAC) <- names(dataMatrixList)
names(clustSampleIndexList_PACR) <- names(dataMatrixList)
names(clustFeatureIndexList_PACR) <- names(dataMatrixList)
names(bestK_PACR) <- names(dataMatrixList)
names(PACR) <- names(dataMatrixList)
names(consensusClustOutput) <- names(dataMatrixList)
names(minConsensusClusterByK)<-  names(dataMatrixList)
names(meanConsensusClusterByK) <-  names(dataMatrixList)
names(consensusFracByK) <- names(dataMatrixList)
      
  output <- list(  clustSampleIndexList_consensusFrac=clustSampleIndexList_consensusFrac,clustFeatureIndexList_consensusFrac =clustFeatureIndexList_consensusFrac,
                   bestK_consensusFrac=bestK_consensusFrac,
                   clustSampleIndexList_meanConsensusCluster=clustSampleIndexList_meanConsensusCluster,
                   clustFeatureIndexList_meanConsensusCluster=clustFeatureIndexList_meanConsensusCluster,
                   bestK_meanConsensusCluster=bestK_meanConsensusCluster,clustSampleIndexList_minConsensusCluster=clustSampleIndexList_minConsensusCluster,
                   clustFeatureIndexList_minConsensusCluster=clustFeatureIndexList_minConsensusCluster,
                   bestK_minConsensusCluster=bestK_minConsensusCluster,consensusInfo=consensusInfo,minConsensusClusterByK=minConsensusClusterByK,
                   meanConsensusClusterByK=meanConsensusClusterByK,    bestK_PAC= bestK_PAC,
                   clustSampleIndexList_PAC= clustSampleIndexList_PAC,
                   clustFeatureIndexList_PAC=clustFeatureIndexList_PAC,PAC=PAC,consensusClustOutput=consensusClustOutput,
                   bestK_PACR= bestK_PACR,
                   clustSampleIndexList_PACR= clustSampleIndexList_PACR,
                   clustFeatureIndexList_PACR=clustFeatureIndexList_PACR,PACR=PACR,consensusFracByK=consensusFracByK
                   )
  
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
#rep of around 100 recommended.
consensusClusterMatrix <- function(dataMatrix, clustFeatures,maxNumClusters = 30, numSims=10, pItem=0.8, pFeature=1, clusterAlg=c("km","hc","pam","kmdist"),
                                   hclustAlgorithm=c("average","complete","ward.D", "ward.D2", "single", "mcquitty","median","centroid"), consensusHclustAlgorithm=c("average","complete","ward.D", "ward.D2", "single", "mcquitty","median","centroid"),
corUse=c("everything","pairwise.complete.obs", "complete.obs"),iter.max=30,nstart=25,
distMethod=c("euclidean","spearman","euclidean", "pearson","binary", "maximum", "canberra", "minkowski"),
minMeanClustConsensus=.7,
outputFile="./consensusOut.txt",minClustConsensus=.7,maxPAC=.15,studyName="test"){
  
 innerLinkage <- hclustAlgorithm
 finalLinkage <- consensusHclustAlgorithm
 
 
 if(length(distMethod)>1){
   warning("Multiple inputs for innerLinkage  so setting to default \'euclidean\'")
   distMethod <- "euclidean"
   
 }
 if(length(innerLinkage)>1){
   
   warning("Multiple inputs for innerLinkage (hclust) method so setting to default \'average\'")
   innerLinkage <- "average"
 }
 
 if(length(finalLinkage)>1){
   
   warning("Multiple inputs for finalLinkage so setting to default \'average\'")
   finalLinkage <- "average"
   
 }
  #come back: add defaults to make more user-friendly.
  #if(length(clusterAlg)>1 || length(bestKmethod)>1)
   #ConsensusClustPlus: assumes clustering the columns.
  
  dataset <- dataMatrix[rownames(dataMatrix) %in% clustFeatures, , drop=FALSE]
  
    if(nrow(dataset)==0){
    
    stop("\nIn clusterMatrixHclustGap function: no clustFeatures were found in the data matrix inputted.")
    
  }
 
     if(clusterAlg==c("hc")){
      #make distance matrix ahead of time and feed in as object; have more distance options with proxy package.
            if(distMethod==("pearson") || distMethod=="spearman"){
            
              #for cor: don't transpose dataset.
              #it turns out that R will recognize the upper-level function input value in this function (algorith, corUse, etc.)
              clustDataset <- as.dist((1-cor(dataset,use=corUse,method=distMethod)), method=algorithm)
              #clustGap needs a list output
       #       clusterAssignments <- cutree(clustObject,k=K)
        #      distMethod=distMethod

         }else{

            if(distMethod=="euclidean"){
              #capitalize
              distMethod <- "Euclidean"
              
            }
             clustDataset <- dist(dataset,method=distMethod,by_rows = FALSE,
                       convert_similarities = TRUE)
 
            
          }
          
    }else{
      
      clustDataset <- dataset
    }


  #remember: we're clustering pItem*ncol(dataset) each time...so
  #our max K is NOT ncol(dataset), but rather pItem*ncol(dataset)
  #add the -1 for rounding...consensus cluster plus appears to do that?
  if((floor(ncol(dataset)*pItem)-1)<maxNumClusters){
    #hclust usually returns NA gap test if K.max = pItem*ncol(dataset) as opposed to pItem*ncol(dataset)-1
    #will also mest up maxSE calculations
    K.max <- floor(ncol(dataset)*pItem)-1
    
  }else{
    
    K.max <- maxNumClusters
    
  }

 #   if(nstart==1 && iter.max==10){
#  library("ConsensusClusterPlus")  
  #using default k-means settings; can use ConsensusClusterPlus package.
 # consensusClustOutput <- ConsensusClusterPlus_CoINcIDE (d=clustDataset,
  #                                             maxK=K.max,reps=numSims,
   #                                            pItem=pItem, pFeature=pFeature, clusterAlg=clusterAlg,
  #innerLinkage=innerLinkage, finalLinkage=finalLinkage, ml=NULL,distance=distMethod,
  #tmyPal=NULL,seed=NULL,plot='pngBMP',writeTable=FALSE,weightsItem=NULL,weightsFeature=NULL,verbose=F,corUse=corUse)
  
 
      
      
   # }else{
  #using Bioconductor's kmeans consensus code, but I altered it so that the user can specify iter.max
 #and nstart for kmeans algorithms.  otherwise you must always run it with nstart=1...this will
 #naturally not product optimal clusterings for each resampling.
  consensusClustOutput <- ConsensusClusterPlus_CoINcIDE(d=clustDataset,iter.max=iter.max,nstart=nstart,
                                               maxK=K.max,reps=numSims,
                                               pItem=pItem, pFeature=pFeature, clusterAlg=clusterAlg,
  innerLinkage=innerLinkage, finalLinkage=finalLinkage, ml=NULL,distance=distMethod,
  tmyPal=NULL,seed=NULL,plot='pngBMP',writeTable=FALSE,weightsItem=NULL,weightsFeature=NULL,verbose=F,corUse=corUse)
  
  #  }
  #get consensus score for each cluster for each k: 
  #hmm..doesn't work for k=1?? not so great...
  #weird...can't NOT have the plots outputted...tried title=NULL
  icl <- calcICL(consensusClustOutput,plot='pngBMP')
  
  consensusByK <- split(icl[["clusterConsensus"]][,"clusterConsensus"],f=icl[["clusterConsensus"]][,"k"])
  #if one cluster's consensus is NA: will return NA. but want this.
  #we don't want to pick a K that results in an NaN consensus value; 
  #this probably means should stick with a lower k.
  #don't remove NAs here.
  meanConsensusClusterByK <- lapply(consensusByK,FUN=function(consensusUnit){mean(consensusUnit)})
  minConsensusClusterByK <- lapply(consensusByK,FUN=function(consensusUnit){min(consensusUnit)})
  
  #is there at least one clustering that passes our minimum meanClust threshold?
   
  #consensusFract calculation
  #are they all NAs?
 
 
 selectK_consensusFrac <- list()
 selectK_minConsensusClust <- list()
 selectK_meanConsensusClust <- list()
 selectK_PAC <- list()
 selectK_PACR <- list()
 
  if(length(meanConsensusClusterByK)> 0 && length(minConsensusClusterByK)> 0 && !(all(unlist(meanConsensusClusterByK)=="NaN")) && any(unlist(meanConsensusClusterByK)>=minMeanClustConsensus,na.rm=TRUE) && any(unlist(minConsensusClusterByK)>=minClustConsensus,na.rm=TRUE)){
    
    meanBetweenConsensusClusterByK <- list()
    consensusFrac <- list()
    consensusMetric <- list()
    PAC <- list()

   
    for(i in 2:K.max){
      
      N <- nrow(consensusClustOutput[[i]]$consensusMatrix)
      #technically could just take mean - is a symmetric matrix. diag=FALSE bc is a sample compared against itself (should always be 1.)
      upperTri <- consensusClustOutput[[i]]$consensusMatrix[upper.tri(consensusClustOutput[[i]]$consensusMatrix,diag=FALSE)]
      #http://www.nature.com/srep/2014/140827/srep06207/full/srep06207.html
      #basically looking at the numbers between 0 and 1
      #paper uses .9, ,.1 indices for gene expression.
      
      PAC[[i]] <- length(upperTri[which(upperTri<=.9)])/(N*(N-1)/2) - length(upperTri[which(upperTri<=.1)])/(N*(N-1)/2) 
      #return NaN if there are NaN values.
      
      #calculate between sum of squares
      sampleAssignments <- consensusClustOutput[[i]]$consensusClass
      
      tmp <- consensusClustOutput[[i]]$consensusMatrix
      for(c in 1:i){
        #zero out samples that are in the same cluster.
        tmp[which(sampleAssignments==c),which(sampleAssignments==c)] <- 0
        
      }
      #technically could just take mean - is a symmetric matrix. diag=FALSE bc it is zero here.
      upperTri <- tmp[upper.tri(tmp,diag=FALSE)]
    #return NaN if there are NaN values.
      #the # of upperTri that is zero will always be the same; samples can only belong to one cluster.
      meanBetweenConsensusClusterByK[[i]] <- mean(upperTri)
      #meanConsensusClusterByK: indices start at 1, even though that's for k=2
      #the meanBetween is almost always pretty low, but it's not monotonic - i.e. a larger K is not automatically favored.
      #not sure this logic works for hclust; in this case, the meanBetweenConsensusClusterByK appears
      consensusFrac[[i]] <- meanConsensusClusterByK[[(i-1)]]/meanBetweenConsensusClusterByK[[i]]

    }
  
    #add NA in index 1, otherwise when unlist, will remove first index and mess up which.max/min
    consensusFrac[[1]] <- NA
    PAC[[1]] <- NA
    ##consensusFrac calculations
    #which.max(): Missing and NaN values are discarded.
    selectK_consensusFrac$bestK <- which.max(unlist(consensusFrac))
    selectK_minConsensusClust$bestK <- which.max(unlist(minConsensusClusterByK))+1
    selectK_meanConsensusClust$bestK <- which.max(unlist(meanConsensusClusterByK))+1
    
    #round PAC - so difference at hundreth level rounded and then pick highest K at the tenth decimal level
    
    PAC_noNAs <- unlist(PAC)[which(!is.na(PAC))]
    
    if(length(PAC_noNAs)>0 && any(PAC_noNAs<maxPAC)){
 
      #this will make .00x 0 (e.g. .002 becomes zero.)
      #anything with zero in tenth, hundreth place counted as zero.
      PACR <- round(unlist(PAC),digits=2)
      possibleKs <- which(PACR==min(PACR,na.rm=TRUE))
      maxPossibleK <- possibleKs[length(possibleKs)]
      selectK_PACR$bestK  <-  maxPossibleK
      #highly unlikely two unrounded values will match, but still keep this:
      possibleKs <- which(unlist(PAC)==min(unlist(PAC),na.rm=TRUE))
      maxPossibleK <- possibleKs[length(possibleKs)]
      selectK_PAC$bestK <- maxPossibleK
    
    }else{
      
      selectK_PAC$bestK <- 1
      selectK_PACR$bestK <- 1
      PACR <- NA
      
    }
    
    if(length(selectK_consensusFrac$bestK)==0){
      #all NAs returned; set K=1
      selectK_consensusFrac$bestK <- 1
      
    }
    
    if(length(selectK_minConsensusClust$bestK)==0){
      #all NAs returned; set K=1
      selectK_minConsensusClust$bestK <- 1
      
    }
    
    if(length(selectK_meanConsensusClust$bestK)==0){
      #all  NAs returned; set K=1
      selectK_meanConsensusClust$bestK <- 1
      
    }
    
    if(length(selectK_PAC$bestK)==0){
      #all NAs returned; set K=1
      selectK_PAC$bestK <- 1
      
    }
    
    
    if(length(selectK_PACR$bestK)==0){
      #all NAs returned; set K=1
      selectK_PACR$bestK <- 1
      PACR <- NA
      
    }

  
  }else{
    #no clusterings passed the minMeanConsensus threshold
    selectK_consensusFrac$bestK <- 1
    consensusFrac <- NA
    meanBetweenConsensusClusterByK <- NA
    PAC <- NA
    PACR <- NA
    selectK_meanConsensusClust$bestK <- 1
    selectK_minConsensusClust$bestK  <- 1
    selectK_PAC$bestK <- 1
    selectK_PACR$bestK <- 1
    
  }
 
 message(paste0("Best K as determined by consensusFract is: ",selectK_consensusFrac$bestK))
 cat(paste0("\nBest K as determined by ",clusterAlg, " and consensusFract  is: ", selectK_consensusFrac$bestKbestK ,"\n"),
     append=TRUE,file=outputFile);
 
 message(paste0("Best K as determined by meanConsensus is: ",selectK_meanConsensusClust$bestK))
 cat(paste0("\nBest K as determined by ",clusterAlg, " and meanConsensus  is: ", selectK_meanConsensusClust$bestK ,"\n"),
     append=TRUE,file=outputFile);
 
 message(paste0("Best K as determined by minConsensus is: ",selectK_minConsensusClust$bestK))
 cat(paste0("\nBest K as determined by ",clusterAlg, " and minConsensus  is: ", selectK_minConsensusClust$bestK ,"\n"),
     append=TRUE,file=outputFile);
 
 
 message(paste0("Best K as determined by unrounded PAC is: ",selectK_PAC$bestK))
 cat(paste0("\nBest K as determined by ",clusterAlg, " and rounded PAC  is: ", selectK_PAC$bestK ,"\n"),
     append=TRUE,file=outputFile);
 
 
 message(paste0("Best K as determined by rounded PACR is: ",selectK_PACR$bestK))
 cat(paste0("\nBest K as determined by ",clusterAlg, " and rounded PAC  is: ", selectK_PACR$bestK ,"\n"),
     append=TRUE,file=outputFile);
 
 
 
packageClustOutput <-  function(listOutputObject,clustFeatures,dataMatrix){
  
 listOutputObject$clustFeatureIndexList <- list()
 listOutputObject$clustSampleIndexList <- list()
  
  if(listOutputObject$bestK >1){

    listOutputObject$clusterAssignments <- consensusClustOutput[[listOutputObject$bestK]]$consensusClass
  
  }else{
    
    listOutputObject$clusterAssignments <- rep.int(1,times=ncol(dataset))
  }
  
 for(k in 1:listOutputObject$bestK){
   
   listOutputObject$clustFeatureIndexList[[k]] <- na.omit(match(clustFeatures,rownames(dataMatrix)))
   listOutputObject$clustSampleIndexList[[k]] <- which(listOutputObject$clusterAssignments==k)
   
   
 }
 
 return(listOutputObject)
 
}

selectK_consensusFrac <- packageClustOutput(selectK_consensusFrac,clustFeatures,dataMatrix)
selectK_meanConsensusClust <- packageClustOutput(selectK_meanConsensusClust,clustFeatures,dataMatrix)
selectK_minConsensusClust <- packageClustOutput(selectK_minConsensusClust,clustFeatures,dataMatrix)
selectK_PAC <- packageClustOutput(selectK_PAC,clustFeatures,dataMatrix)
selectK_PACR <- packageClustOutput(selectK_PACR,clustFeatures,dataMatrix)

  output <- list(selectK_consensusFrac=selectK_consensusFrac, selectK_PAC=selectK_PAC ,
                 selectK_meanConsensusClust=selectK_meanConsensusClust,
                 meanConsensusClusterByK=meanConsensusClusterByK,selectK_minConsensusClust=selectK_minConsensusClust,
                 consensusCalc=icl,consensusByK=consensusByK,
                 consensusClustOutput=consensusClustOutput,consensusFrac=consensusFrac,
                 meanBetweenConsensusClusterByK=meanBetweenConsensusClusterByK,minConsensusClusterByK=minConsensusClusterByK,
                 PAC=PAC,PACR=PACR,selectK_PACR=selectK_PACR)


 return(output)
 
}

#Katie: added iter.max, nstart as options for kmeans (changed nothing else.)
ConsensusClusterPlus_CoINcIDE <- function( d=NULL,iter.max=20,nstart=15,
                                  maxK = 3,
                                  reps=10,
                                  pItem=0.8,
                                  pFeature=1,
                                  clusterAlg="hc",
                                  title="untitled_consensus_cluster",
                                  innerLinkage="average",
                                  finalLinkage="average",
                                  distance="pearson",
                                  ml=NULL,
                                  tmyPal=NULL,
                                  seed=NULL,
                                  plot=NULL,
                                  writeTable=FALSE,
                                  weightsItem=NULL,
                                  weightsFeature=NULL,
                                  verbose=F,
                                  corUse="everything" ) {
  ##description: runs consensus subsamples 
  if(is.null(seed)==TRUE){
    seed=timeSeed = as.numeric(Sys.time())
  }
  set.seed(seed)
  
  #distance=ifelse( inherits(d,"dist"), attr( d, "method" ), "pearson" )
  
  
  if(is.null(ml)==TRUE){
    
    if ( ! class( d ) %in% c( "dist", "matrix", "ExpressionSet" ) ) {
      stop("d must be a matrix, distance object or ExpressionSet (eset object)")
    }
    
    if ( inherits( d, "dist" ) ) {
      ## if d is a distance matrix, fix a few things so that they don't cause problems with the analysis
      ##  Note, assumption is that if d is a distance matrix, the user doesn't want to sample over the row features
      if ( is.null( attr( d, "method" ) ) ) {
        attr( d, "method" ) <- distance <- "unknown - user-specified"
      }
      if ( is.null( distance ) || ( distance != attr( d, "method" ) ) ) {
        distance <- attr( d, "method" ) 
      }
      
      if ( ( ! is.null( pFeature ) ) && ( pFeature < 1 ) ) {
        message( "Cannot use the pFeatures parameter when specifying a distance matrix as the data object\n" )
        pFeature <- 1
      }
      if ( ! is.null( weightsFeature ) ) {
        message( "Cannot use the weightsFeature parameter when specifying a distance matrix as the data object\n" )
        weightsFeature <- NULL
      }
      if ( clusterAlg == "km" ) {
        message( "Note: k-means will cluster the distance matrix you provided.  This is similar to kmdist option when suppling a data matrix")
        ##d <- as.matrix( d )  #this is now done w/in ccRun
      }
    } else {
      if ( is.null( distance ) ) {
        ## we should never get here, but just in case
        distance <- "pearson"
      }
    }
    
    if ( ( clusterAlg == "km" ) && inherits( distance, "character" ) && ( distance != "euclidean" ) ) {
      message( "Note: The km (kmeans) option only supports a euclidean distance metric when supplying a data matrix.  If you want to cluster a distance matrix using k-means use the 'kmdist' option, or use a different algorithm such as 'hc' or 'pam'.  Changing distance to euclidean")
      distance <- 'euclidean'
    }
    
    
    if ( inherits( d,"ExpressionSet" ) ) {
      d <- exprs(d)
    }
    
    ml <- ccRun(d=d,
                 maxK=maxK,
                 repCount=reps,
                 diss=inherits(d,"dist"),
                 pItem=pItem,
                 pFeature=pFeature,
                 innerLinkage=innerLinkage,
                 clusterAlg=clusterAlg,
                 weightsFeature=weightsFeature,
                 weightsItem=weightsItem,
                 distance=distance,
                 verbose=verbose,
                 corUse=corUse)
  }
  res=list();
  
  ##make results directory
  #KP: comment out
  #if((is.null(plot)==FALSE | writeTable) & !file.exists(paste(title,sep=""))){
  #  dir.create(paste(title,sep=""))
  #}
  
  ##write log file
  log <- matrix( ncol=2,
                 byrow=T,
                 c("title",title,
                   "maxK",maxK,
                   "input matrix rows",ifelse ( inherits( d, "matrix" ), nrow(d), "dist-mat" ), 
                   "input matrix columns",ifelse ( inherits( d, "matrix" ), ncol(d), ncol( as.matrix(d) ) ), 
                   "number of bootstraps",reps,
                   "item subsampling proportion",pItem,
                   "feature subsampling proportion",ifelse( is.null(pFeature), 1, pFeature ),
                   "cluster algorithm",clusterAlg,
                   "inner linkage type",innerLinkage,
                   "final linkage type",finalLinkage,
                   "correlation method",distance,
                   "plot",if(is.null(plot)) NA else plot,
                   "seed",if(is.null(seed)) NA else seed))
  colnames(log) = c("argument","value")
  #Katie: commented out.
#   if(writeTable){
#     write.csv(file=paste(title,"/",title,".log.csv",sep=""), log,row.names=F)
#   }
#   if(is.null(plot)){
#     ##nothing
#   }else if(plot=="pngBMP"){
#     bitmap(paste(title,"/","consensus%03d.png",sep=""))
#   }else if(plot=="png"){
#     png(paste(title,"/","consensus%03d.png",sep=""))
#     
#   }else if (plot=="pdf"){
#     pdf(onefile=TRUE, paste(title,"/","consensus.pdf",sep=""))
#   }else if (plot=="ps"){
#     postscript(onefile=TRUE, paste(title,"/","consensus.ps",sep=""))
#   }	
# #   
#   colorList=list()
#   colorM = rbind() #matrix of colors.
#   
#   #18 colors for marking different clusters
#   thisPal <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928",
#                "#bd18ea", #magenta
#                "#2ef4ca", #aqua
#                "#f4cced", #pink,
#                "#f4cc03", #lightorange
#                "#05188a", #navy,
#                "#e5a25a", #light brown
#                "#06f106", #bright green
#                "#85848f", #med gray
#                "#000000", #black
#                "#076f25", #dark green
#                "#93cd7f",#lime green
#                "#4d0776", #dark purple
#                "#ffffff" #white
#   )
#   
#   ##plot scale
#   colBreaks=NA
#   if(is.null(tmyPal)==TRUE){
#     colBreaks=10
#     tmyPal = myPal(colBreaks)
#   }else{
#     colBreaks=length(tmyPal)
#   }
#   sc = cbind(seq(0,1,by=1/( colBreaks) )); rownames(sc) = sc[,1]
#   sc = cbind(sc,sc)
#KP: no heatmap
 # heatmap(sc, Colv=NA, Rowv=NA, symm=FALSE, scale='none', col=tmyPal, na.rm=TRUE,labRow=rownames(sc),labCol=F,main="consensus matrix legend")
  
  for (tk in 2:maxK){
    if(verbose){
      message(paste("consensus ",tk))
    }
    fm = ml[[tk]]
    hc=hclust( as.dist( 1 - fm ), method=finalLinkage);
    message("clustered")	
    ct = cutree(hc,tk)
    names(ct) = colnames(d)
    if(class(d)=="dist"){
      names(ct) = colnames(as.matrix(d))
    }
    c = fm
    #Katie
    #colorList = setClusterColors(res[[tk-1]][[3]],ct,thisPal,colorList)
    pc = c
    pc=pc[hc$order,] #pc is matrix for plotting, same as c but is row-ordered and has names and extra row of zeros.
 
    #Katie:
#     
#     
#     if(!is.null(plot) && plot=="pngBMP"){
#       pc = pc[,hc$order ] #mod for no tree
#       pc = rbind(pc,0)
#       #no dendrogram if pngBMP
#       oc = colorList[[1]][hc$order] #mod for no tree
#       heatmap(pc, Colv = NA, Rowv = NA, symm = FALSE, scale = "none", col = tmyPal, na.rm = TRUE, labRow = F, labCol = F, mar = c(5, 5), main = paste("consensus matrix k=", 
#                                                                                                                                                       tk, sep = ""), ColSideCol = oc)
#     }else{
#       pc = rbind(pc,0)
#       #former with tree:
#       heatmap(pc, Colv=as.dendrogram(hc), Rowv=NA, symm=FALSE, scale='none', col=tmyPal, na.rm=TRUE,labRow=F,labCol=F,mar=c(5,5),main=paste("consensus matrix k=",tk,sep="") , ColSideCol=colorList[[1]])
#     }
#     
#     legend("topright",legend=unique(ct),fill=unique(colorList[[1]]),horiz=FALSE )
    #KP: removed clrs=colorList
#Katie: also returned fullmlList for CDF plotting (yes, a bit redundant to save at each step.)
    res[[tk]] = list(consensusMatrix=c,consensusTree=hc,consensusClass=ct,ml=ml[[tk]])
    #colorM = rbind(colorM,colorList[[1]]) 
  }
  #CDF(ml)
  #clusterTrackingPlot(colorM[,res[[length(res)]]$consensusTree$order])
  #if(is.null(plot)==FALSE){
  #  dev.off();
  #}
  #res[[1]] = colorM
#KP:
res[[1]] = NA

  #if(writeTable){
  #  for(i in 2:length(res)){
  #    write.csv(file=paste(title,"/",title,".k=",i,".consensusMatrix.csv",sep=""), res[[i]]$consensusMatrix)
  #    write.table(file=paste(title,"/",title,".k=",i,".consensusClass.csv",sep=""), res[[i]]$consensusClass,col.names = F,sep=",")
  #  }
  #}
  return(res)
}


calcICL = function(res,title="untitled_consensus_cluster",plot=NULL,writeTable=FALSE){
  #calculates and plots cluster consensus and item consensus
  cc=rbind()
  cci = rbind()
  sumRes=list()
  colorsArr=c()
  
#   #make results directory
#KP: commented out.
#   if((is.null(plot)==FALSE | writeTable) & !file.exists(paste(title,sep=""))){
#     dir.create(paste(title,sep=""))
#   }
#   if(is.null(plot)){
#     #to screen
#   }else if(plot=="pdf"){
#     pdf(onefile=TRUE, paste(title,"/","icl.pdf",sep=""))
#   }else if(plot=="ps"){
#     postscript(onefile=TRUE, paste(title,"/","icl.ps",sep=""))
#   }else if (plot=="png"){
#     png(paste(title,"/","icl%03d.png",sep=""))
#   }else if (plot=="pngBMP"){
#     bitmap(paste(title,"/","icl%03d.png",sep=""))
#   }
#   
  par(mfrow=c(3,1),mar=c(4,3,2,0))
  
  for (k in 2:length(res)){ #each k
    eiCols = c();
    o = res[[k]]
    m = o$consensusMatrix
    m = triangle(m,mode=2)
    for (ci in sort(unique(o$consensusClass))){ #each cluster in k
      items = which(o$consensusClass==ci)
      nk = length(items)
      mk = sum( m[items,items], na.rm=T)/((nk*(nk-1))/2)
      cc=rbind(cc,c(k,ci,mk)) #cluster-consensus
      
      for (ei in rev(res[[2]]$consensusTree$order) ){
        denom = if (ei %in% items) { nk - 1} else { nk }
        mei = sum( c(m[ei,items],m[items,ei]), na.rm=T)/denom  # mean item consensus to a cluster.
        cci = rbind(cci,c(k,ci,ei,mei)) #cluster, cluster index, item index, item-consensus
      }
      eiCols = c(eiCols, rep(ci,length(o$consensusClass)) )
    }
    
    cck = cci[which(cci[,1]==k),] #only plot the new k data.
    
    #group by item, order by cluster i
    w=lapply(split(cck,cck[,3]), function(x) { y=matrix(unlist(x),ncol=4); y[order(y[,2]),4] }) 
    q = matrix(as.numeric(unlist(w)),ncol=length(w),byrow=F)
    q = q[,res[[2]]$consensusTree$order] #order by leave order of k=2
    #q is a matrix of k rows and sample columns, values are item consensus of sample to the cluster.
    #KP: commented out:
#     thisColors = unique(cbind(res[[k]]$consensusClass,res[[k]]$clrs[[1]]))
#     thisColors=thisColors[order(as.numeric(thisColors[,1])),2]
#     colorsArr=c(colorsArr,thisColors)
#     sumRes[[k]] = rankedBarPlot(q,thisColors,cc=res[[k]]$consensusClass[res[[2]]$consensusTree$order],paste("k=",k,sep="") )
  }
  
#   ys=cs=lab=c()
#   lastk=cc[1,1]
#   for(i in 1:length(colorsArr)){
#     if(lastk != cc[i,1]){
#       ys=c(ys,0,0)
#       cs=c(cs,NA,NA)
#       lastk=cc[i,1]
#       lab=c(lab,NA,NA)
#     }
#     ys=c(ys,cc[i,3])
#     cs=c(cs,colorsArr[i])
#     lab=c(lab,cc[i,1])
#   }
#   names(ys) = lab
#   par(mfrow=c(3,1),mar=c(4,3,2,0))
#   barplot(ys,col=cs,border=cs,main="cluster-consensus",ylim=c(0,1),las=1)
#   if(is.null(plot)==FALSE){
#     dev.off()
#   }
  colnames(cc) = c("k","cluster","clusterConsensus")
  colnames(cci) = c("k","cluster","item","itemConsensus")
  cci[,"item"] = names(res[[2]]$consensusClass)[ cci[,"item"] ]
  #type cci
  cci = data.frame( k=as.numeric(cci[,"k"]), cluster=as.numeric(cci[,"cluster"]), item=cci[,"item"], itemConsensus=as.numeric(cci[,"itemConsensus"])) 
  
#   #write to file.
#   if(writeTable){
#     write.csv(file=paste(title,"/",title,".summary.cluster.consensus.csv",sep=""),row.names=F, cc)
#     write.csv(file=paste(title,"/",title,".summary.item.consensus.csv",sep=""), row.names=F, cc)
#   }
  return(list(clusterConsensus=cc,itemConsensus=cci))
}


ccRun <- function( d=d,nstart=15,iter.max=20,
                   maxK=NULL,
                   repCount=NULL,
                   diss=inherits( d, "dist" ),
                   pItem=NULL,
                   pFeature=NULL,
                   innerLinkage=NULL,
                   distance=NULL, #ifelse( inherits(d,"dist"), attr( d, "method" ), "euclidean" ),@@@@@
                   clusterAlg=NULL,
                   weightsItem=NULL,
                   weightsFeature=NULL,
                   verbose=NULL,
                   corUse=NULL) {
  m = vector(mode='list', repCount)
  ml = vector(mode="list",maxK)
  n <- ifelse( diss, ncol( as.matrix(d) ), ncol(d) )
  mCount = mConsist = matrix(c(0),ncol=n,nrow=n)
  ml[[1]] = c(0);
  
  if (is.null( distance ) ) distance <- 'euclidean'  ## necessary if d is a dist object and attr( d, "method" ) == NULL
  
  acceptable.distance <- c( "euclidean", "maximum", "manhattan", "canberra", "binary","minkowski",
                            "pearson", "spearman" )
  
  main.dist.obj <- NULL
  if ( diss ){
    main.dist.obj <- d
    ## reset the pFeature & weightsFeature params if they've been set (irrelevant if d is a dist matrix)
    if ( ( !is.null(pFeature) ) &&
           ( pFeature < 1 ) ) {
      message( "user-supplied data is a distance matrix; ignoring user-specified pFeature parameter\n" )
      pFeature <- 1 # set it to 1 to avoid problems with sampleCols
    }
    if ( ! is.null( weightsFeature ) ) {
      message( "user-supplied data is a distance matrix; ignoring user-specified weightsFeature parameter\n" )
      weightsFeature <- NULL  # set it to NULL to avoid problems with sampleCols
    }
  } else { ## d is a data matrix
    ## we're not sampling over the features
    if ( ( clusterAlg != "km" ) &&
           ( is.null( pFeature ) ||
               ( ( pFeature == 1 ) && is.null( weightsFeature ) ) ) ) {
      ## only generate a main.dist.object IFF 1) d is a matrix, 2) we're not sampling the features, and 3) the algorithm isn't 'km'
      if ( inherits( distance, "character" ) ) {
        if ( ! distance %in%  acceptable.distance  &  ( class(try(get(distance),silent=T))!="function") ) stop("unsupported distance.")
        
        if(distance=="pearson" | distance=="spearman"){
          main.dist.obj <- as.dist( 1-cor(d,method=distance,use=corUse ))
        }else if( class(try(get(distance),silent=T))=="function"){
          main.dist.obj <- get(distance)( t( d )   )
        }else{
          main.dist.obj <- dist( t(d), method=distance )
        }
        attr( main.dist.obj, "method" ) <- distance  
      } else stop("unsupported distance specified.")
    } else {
      ## pFeature < 1 or a weightsFeature != NULL
      ## since d is a data matrix, the user wants to sample over the gene features, so main.dist.obj is left as NULL
    }
  }
  
  
  for (i in 1:repCount){
    if(verbose){
      message(paste("random subsample",i));
    }
    ## take expression matrix sample, samples and genes
    sample_x = sampleCols( d, pItem, pFeature, weightsItem, weightsFeature )
    
    this_dist = NA
    if ( ! is.null( main.dist.obj ) ) {
      boot.cols <- sample_x$subcols
      this_dist <- as.matrix( main.dist.obj )[ boot.cols, boot.cols ]
      if ( clusterAlg != "km" ) {
        ## if this isn't kmeans, then convert to a distance object
        this_dist <- as.dist( this_dist )
        attr( this_dist, "method" ) <- attr( main.dist.obj, "method" )
      }
    } else {
      ## if main.dist.obj is NULL, then d is a data matrix, and either:
      ##   1) clusterAlg is 'km'
      ##   2) pFeatures < 1 or weightsFeatures have been specified, or
      ##   3) both
      ## so we can't use a main distance object and for every iteration, we will have to re-calculate either
      ##   1) the distance matrix (because we're also sampling the features as well), or
      ##   2) the submat (if using km) 
      
      if ( clusterAlg != "km" )  {
        if ( ! distance %in% acceptable.distance &  ( class(try(get(distance),silent=T))!="function")  ) stop("unsupported distance.")
        if( ( class(try(get(distance),silent=T))=="function") ){
          this_dist <- get(distance)( t( sample_x$submat ) )
        }else{
          if( distance == "pearson" | distance == "spearman"){
            this_dist <- as.dist( 1-cor(sample_x$submat,use=corUse,method=distance) )
          }else{
            this_dist <- dist( t( sample_x$submat ), method= distance  )
          }
        }
        attr( this_dist, "method" ) <- distance  
      } else {
        ## if we're not sampling the features, then grab the colslice
        if ( is.null( pFeature ) ||
               ( ( pFeature == 1 ) && is.null( weightsFeature ) ) ) {
          this_dist <- d[, sample_x$subcols ]
        } else {
          if ( is.na( sample_x$submat ) ) {
            stop( "error submat is NA" )
          }
          
          this_dist <- sample_x$submat
        } 
      }
    }
    
    ## cluster samples for HC.
    this_cluster=NA
    if(clusterAlg=="hc"){
      this_cluster = hclust( this_dist, method=innerLinkage)
    }
    ##mCount is possible number of times that two sample occur in same random sample, independent of k
    ##mCount stores number of times a sample pair was sampled together.
    mCount <- connectivityMatrix( rep( 1,length(sample_x[[3]])),
                                  mCount,
                                  sample_x[[3]] ) 
    
    ##use samples for each k		
    for (k in 2:maxK){
      if(verbose){
        message(paste("  k =",k))
      }
      if (i==1){
        ml[[k]] = mConsist #initialize
      }
      this_assignment=NA
      if(clusterAlg=="hc"){
        ##prune to k for hc
        this_assignment = cutree(this_cluster,k)
        
      }else if(clusterAlg=="kmdist"){
        this_assignment = kmeans(this_dist, k, iter.max = iter.max, nstart = nstart, algorithm = c("Hartigan-Wong") )$cluster
        
      }else if(clusterAlg=="km"){
        ##this_dist should now be a matrix corresponding to the result from sampleCols
        this_assignment <- kmeans( t( this_dist ),
                                   k,
                                   iter.max = iter.max,
                                   nstart = nstart,
                                   algorithm = c("Hartigan-Wong") )$cluster
      }else if ( clusterAlg == "pam" ) {
        this_assignment <- pam( x=this_dist,
                                k,
                                diss=TRUE,
                                metric=distance, 
                                cluster.only=TRUE )
      } else{
        ##optional cluterArg Hook.
        this_assignment <- get(clusterAlg)(this_dist, k)
      }
      ##add to tally				
      ml[[k]] <- connectivityMatrix( this_assignment,
                                     ml[[k]],
                                     sample_x[[3]] )
    }
  }
  
  
  ##consensus fraction
  res = vector(mode="list",maxK)
  for (k in 2:maxK){
    ##fill in other half of matrix for tally and count.
    tmp = triangle(ml[[k]],mode=3)
    tmpCount = triangle(mCount,mode=3)
    res[[k]] = tmp / tmpCount
    res[[k]][which(tmpCount==0)] = 0
  }
  message("end fraction")
  return(res)
}


connectivityMatrix <- function( clusterAssignments, m, sampleKey){
  ##input: named vector of cluster assignments, matrix to add connectivities
  ##output: connectivity matrix
  names( clusterAssignments ) <- sampleKey 
  cls <- lapply( unique( clusterAssignments ), function(i) as.numeric( names( clusterAssignments[ clusterAssignments %in% i ] ) ) )  #list samples by clusterId
  
  for ( i in 1:length( cls ) ) {
    nelts <- 1:ncol( m )
    cl <- as.numeric( nelts %in% cls[[i]] ) ## produces a binary vector
    updt <- outer( cl, cl ) #product of arrays with * function; with above indicator (1/0) statement updates all cells to indicate the sample pair was observed int the same cluster;
    m <- m + updt
  }
  return(m)
}



sampleCols <- function( d,
                        pSamp=NULL,
                        pRow=NULL,
                        weightsItem=NULL,
                        weightsFeature=NULL ){
  ## returns a list with the sample columns, as well as the sub-matrix & sample features (if necessary)
  ##  if no sampling over the features is performed, the submatrix & sample features are returned as NAs
  ##  to reduce memory overhead
  
  
  space <- ifelse( inherits( d, "dist" ), ncol( as.matrix(d) ), ncol(d) )
  sampleN <- floor(space*pSamp)
  sampCols <- sort( sample(space, sampleN, replace = FALSE, prob = weightsItem) )
  
  this_sample <- sampRows <- NA
  if ( inherits( d, "matrix" ) ) {
    if ( (! is.null( pRow ) ) &&
           ( (pRow < 1 ) || (! is.null( weightsFeature ) ) ) ) {
      ## only sample the rows and generate a sub-matrix if we're sampling over the row/gene/features
      space = nrow(d)
      sampleN = floor(space*pRow)
      sampRows = sort( sample(space, sampleN, replace = FALSE, prob = weightsFeature) )
      this_sample <- d[sampRows,sampCols]
      dimnames(this_sample) <- NULL
    } else {
      ## do nothing
    }
  }
  return( list( submat=this_sample,
                subrows=sampRows,
                subcols=sampCols ) )
}

CDF=function(ml,breaks=100){
  #plot CDF distribution
 plot(c(0),xlim=c(0,1),ylim=c(0,1),col="white",bg="white",xlab="consensus index",ylab="CDF",main="consensus CDF", las=2)
  k=length(ml)
  this_colors = rainbow(k-1)
  areaK = c()
  for (i in 2:length(ml)){
    v=triangle(ml[[i]],mode=1)
    
    #empirical CDF distribution. default number of breaks is 100    
    h = hist(v, plot=FALSE, breaks=seq(0,1,by=1/breaks))
    h$counts = cumsum(h$counts)/sum(h$counts)
    
    #calculate area under CDF curve, by histogram method.
    thisArea=0
    for (bi in 1:(length(h$breaks)-1)){
      thisArea = thisArea + h$counts[bi]*(h$breaks[bi+1]-h$breaks[bi]) #increment by height by width
      bi = bi + 1
    }
    areaK = c(areaK,thisArea)
    lines(h$mids,h$counts,col=this_colors[i-1],lwd=2,type='l')
  }
  legend(0.8,0.5,legend=paste(rep("",k-1),seq(2,k,by=1),sep=""),fill=this_colors)
  
  #plot area under CDF change.
  deltaK=areaK[1] #initial auc at k=2
  for(i in 2:(length(areaK))){
    #proportional increase relative to prior K.
    deltaK = c(deltaK,( areaK[i] - areaK[i-1])/areaK[i-1])
  }
  plot(1+(1:length(deltaK)),y=deltaK,xlab="k",ylab="relative change in area under CDF curve",main="Delta area",type="b")
}


myPal = function(n=10){
  #returns n colors
  seq = rev(seq(0,255,by=255/(n)))
  palRGB = cbind(seq,seq,255)
  rgb(palRGB,maxColorValue=255)
}

setClusterColors = function(past_ct,ct,colorU,colorList){
  #description: sets common color of clusters between different K
  newColors = c()
  if(length(colorList)==0){
    #k==2
    newColors = colorU[ct]
    colori=2
  }else{
    newColors = rep(NULL,length(ct))
    colori = colorList[[2]]
    mo=table(past_ct,ct)
    m=mo/apply(mo,1,sum)
    for(tci in 1:ncol(m)){ # for each cluster
      maxC = max(m[,tci])
      pci = which(m[,tci] == maxC)				
      if( sum(m[,tci]==maxC)==1 & max(m[pci,])==maxC & sum(m[pci,]==maxC)==1  )  {
        #if new column maximum is unique, same cell is row maximum and is also unique
        ##Note: the greatest of the prior clusters' members are the greatest in a current cluster's members.
        newColors[which(ct==tci)] = unique(colorList[[1]][which(past_ct==pci)]) # one value
      }else{ #add new color.
        colori=colori+1
        newColors[which(ct==tci)] = colorU[colori]
      }
    }
  }
  return(list(newColors,colori,unique(newColors) ))
}

clusterTrackingPlot = function(m){
  #description: plots cluster tracking plot
  #input: m - matrix where rows are k, columns are samples, and values are cluster assignments.
  plot(NULL,xlim=c(-0.1,1),ylim=c(0,1),axes=FALSE,xlab="samples",ylab="k",main="tracking plot")
  for(i in 1:nrow(m)){
    rect(  xleft=seq(0,1-1/ncol(m),by=1/ncol(m)),  ybottom=rep(1-i/nrow(m),ncol(m)) , xright=seq(1/ncol(m),1,by=1/ncol(m)), ytop=rep(1-(i-1)/nrow(m),ncol(m)), col=m[i,],border=NA)   
  }
  #hatch lines to indicate samples
  xl = seq(0,1-1/ncol(m),by=1/ncol(m))
  segments(  xl, rep(-0.1,ncol(m)) , xl, rep(0,ncol(m)), col="black")    #** alt white and black color?
  ypos = seq(1,0,by=-1/nrow(m))-1/(2*nrow(m))
  text(x=-0.1,y=ypos[-length(ypos)],labels=seq(2,nrow(m)+1,by=1))
}

triangle = function(m,mode=1){
  #mode=1 for CDF, vector of lower triangle.
  #mode==3 for full matrix.
  #mode==2 for calcICL; nonredundant half matrix coun
  #mode!=1 for summary 
  n=dim(m)[1]
  nm = matrix(0,ncol=n,nrow=n)
  fm = m
  
  
  nm[upper.tri(nm)] = m[upper.tri(m)] #only upper half
  
  fm = t(nm)+nm
  diag(fm) = diag(m)
  
  nm=fm
  nm[upper.tri(nm)] = NA
  diag(nm) = NA
  vm = m[lower.tri(nm)]
  
  if(mode==1){
    return(vm) #vector 		
  }else if(mode==3){
    return(fm) #return full matrix
  }else if(mode == 2){
    return(nm) #returns lower triangle and no diagonal. no double counts.
  }
  
}


# rankedBarPlot=function(d,myc,cc,title){
#   colors = rbind() #each row is a barplot series
#   byRank = cbind()
#   
#   spaceh = 0.1 #space between bars
#   for(i in 1:ncol(d)){
#     byRank = cbind(byRank,sort(d[,i],na.last=F))
#     colors = rbind(colors,order(d[,i],na.last=F))
#   }
#   maxH = max(c(1.5,apply(byRank,2,sum)),na.rm=T) #maximum height of graph
#   
#   #barplot largest to smallest so that smallest is in front.
#   barp = barplot( apply(byRank,2,sum) ,  col=myc[colors[,1]] ,space=spaceh,ylim=c(0,maxH),main=paste("item-consensus", title),border=NA,las=1  )
#   for(i in 2:nrow(byRank)){
#     barplot( apply(matrix(byRank[i:nrow(byRank),],ncol=ncol(byRank))  ,2,sum), space=spaceh,col=myc[colors[,i]],ylim=c(0,maxH), add=T,border=NA,las=1  )
#   }
#   xr=seq(spaceh,ncol(d)+ncol(d)*spaceh,(ncol(d)+ncol(d)*spaceh)/ncol(d)  )
#   #class labels as asterisks
#   text("*",x=xr+0.5,y=maxH,col=myc[cc],cex=1.4) #rect(xr,1.4,xr+1,1.5,col=myc[cc] )
# }

plotConsensusHeatmap_CoINcIDE <- function(consensusClustOutput,k,plotSave=TRUE,saveDir="./",fileTag=""){
  
  #   #18 colors for marking different clusters
  thisPal <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928",
               "#bd18ea", #magenta
               "#2ef4ca", #aqua
               "#f4cced", #pink,
               "#f4cc03", #lightorange
               "#05188a", #navy,
               "#e5a25a", #light brown
               "#06f106", #bright green
               "#85848f", #med gray
               "#000000", #black
               "#076f25", #dark green
               "#93cd7f",#lime green
               "#4d0776", #dark purple
               "#ffffff" #white
  )
  
  
  res <- consensusClustOutput
  
  #ml in output list is ml[[tk]] inside the ConsensusClusterPlus code.
  #this is the consenusus matrix
  fm = consensusClustOutput[[k]]$ml
  hc=hclust( as.dist( 1 - fm ), method="average");
  message("clustered")  
  ct <- cutree(hc,k)
  
  
  colorList <- list()
  for(t in 2:k){
    tmp <- cutree(hc,t)
   # names(ct) = colnames(exprMatrix)

    colorList = setClusterColors(res[[t-1]][[3]],tmp,thisPal,colorList)
    
  }
  
  if(plotSave){
  #plot heatmap
  png(paste0(plotSaveDir,"/",fileTag,"_consensus_heatmap.png"),width=1000,height=1000,res=200)
      
  }
  
  pc <- fm#for k; from a list 
  pc=pc[hc$order,] #pc is matrix for plotting, same as c but is row-ordered and has names and extra row of zeros.
  
  colBreaks=10
  tmyPal = myPal(colBreaks)
  
  pc = rbind(pc,0)
  # with tree:
  heatmap(pc, Colv=as.dendrogram(hc), Rowv=NA, symm=FALSE, scale='none', col=tmyPal, na.rm=TRUE,labRow=F,labCol=F,mar=c(5,5),main=paste("consensus matrix k=",k,sep="") , ColSideCol=colorList[[1]] )
  legend("topright",legend=unique(ct),fill=unique(colorList[[1]]),horiz=FALSE )
  
  if(plotSave){
    dev.off()
    
  }
  
} 
  
setClusterColors = function(past_ct,ct,colorU,colorList){
    #description: sets common color of clusters between different K
    newColors = c()
    if(length(colorList)==0){
      #k==2
      newColors = colorU[ct]
      colori=2
    }else{
      newColors = rep(NULL,length(ct))
      colori = colorList[[2]]
      mo=table(past_ct,ct)
      m=mo/apply(mo,1,sum)
      for(tci in 1:ncol(m)){ # for each cluster
        maxC = max(m[,tci])
        pci = which(m[,tci] == maxC)				
        if( sum(m[,tci]==maxC)==1 & max(m[pci,])==maxC & sum(m[pci,]==maxC)==1  )  {
          #if new column maximum is unique, same cell is row maximum and is also unique
          ##Note: the greatest of the prior clusters' members are the greatest in a current cluster's members.
          newColors[which(ct==tci)] = unique(colorList[[1]][which(past_ct==pci)]) # one value
        }else{ #add new color.
          colori=colori+1
          newColors[which(ct==tci)] = colorU[colori]
        }
      }
    }
    return(list(newColors,colori,unique(newColors) ))
}
  
myPal = function(n=10){
    #returns n colors
    seq = rev(seq(0,255,by=255/(n)))
    palRGB = cbind(seq,seq,255)
    rgb(palRGB,maxColorValue=255)
}
  



CDF_CoINcIDE <- function(consensusClustOutput,breaks=100,plotSave=TRUE,plotSaveDir="./",fileTag=""){
   
  ml <- list()
  #need to grab all ml values from each k:
  for(c in 2:length(consensusClustOutput)){
    ml[[c]] <- consensusClustOutput[[c]]$ml
    }
  
  if(plotSave){
  #plot CDF distribution
  png(paste0(plotSaveDir,"/",fileTag,"_consensus_CDF.png"),width=2000,height=1000,res=200)
      
  }
  par(oma = c(1, 1, 1, 4))
  
  plot(c(0),xlim=c(0,1),ylim=c(0,1),col="white",bg="white",xlab="consensus index",ylab="CDF",main="consensus CDF", las=2)
  k=length(ml)
  this_colors = rainbow(k-1)
  #grr...still can't get legend to show up if save! oh well...
  legend(1.1,1,legend=paste(rep("",k-1),seq(2,k,by=1),sep=""),inset=c(0,1),
         xpd = TRUE,fill=this_colors)
  #legend(0.7,0.35, # places a legend at the appropriate place ,
  #      legend= as.character(c(1:length(ml)+1)))
  #col=this_colors
  areaK = c()
  for (i in 2:length(ml)){
    v=triangle(ml[[i]],mode=1)
    
    #empirical CDF distribution. default number of breaks is 100    
    h = hist(v, plot=FALSE, breaks=seq(0,1,by=1/breaks))
    h$counts = cumsum(h$counts)/sum(h$counts)
    
    #calculate area under CDF curve, by histogram method.
    thisArea=0
    for (bi in 1:(length(h$breaks)-1)){
      thisArea = thisArea + h$counts[bi]*(h$breaks[bi+1]-h$breaks[bi]) #increment by height by width
      bi = bi + 1
    }
    areaK = c(areaK,thisArea)
    lines(h$mids,h$counts,col=this_colors[i-1],lwd=2,type='l')
  }

  if(plotSave){
    
    dev.off()
  
  }
  #plot area under CDF change.
 # deltaK=areaK[1] #initial auc at k=2
#  for(i in 2:(length(areaK))){
 #   #proportional increase relative to prior K.
#    deltaK = c(deltaK,( areaK[i] - areaK[i-1])/areaK[i-1])
#  }
# plot(1+(1:length(deltaK)),y=deltaK,xlab="k",ylab="relative change in area under CDF curve",main="Delta area",type="b")
 
}

#for consensus clustering.
summarizeChooseK <- function(studyName="study_25055_GPL96_MDACC_M",gapKmeans_pam50Short,
                             consensusClusterOutput){
  
  if(missing(studyName)){
    
    stop("Please provide a studyName character vector that is the name of one of the dataset indices in your consensusClusterOutputList.")
  }
  if(!missing(gapKmeans_pam50Short)){
    message("gap test K: ",length(gapKmeans_pam50Short[[studyName]]))
    
    gapTest <- length(gapKmeans_pam50Short[[studyName]])
    
  }else{
    
    gapTest <- NA
    
  }
  

  message("bestK using consensus fract: ",consensusClusterOutput$bestK_consensusFrac[[studyName]])
  
  consensusFrac <- consensusClusterOutput$bestK_consensusFrac[[studyName]]
  
  message("consensus fract scores for K=1:10: ")
  cat("\n",unlist(consensusClusterOutput$consensus_fract[[studyName]]),"\n")
  
  message("bestK using mean consensus: ",consensusClusterOutput$bestK_meanConsensusCluster[[studyName]])
  
  meanConsensus <- consensusClusterOutput$bestK_meanConsensusCluster[[studyName]]
  
  message("mean consensus for K=1:10: ")
  index <- which(names(consensusClusterOutput$bestK_consensusFrac)==studyName)
  cat("\n",unlist(consensusClusterOutput$meanConsensusClusterByK[[studyName]]),"\n")
  message("bestK using unrounded PAC : ",consensusClusterOutput$bestK_PAC[[studyName]])
  
  PAC <- consensusClusterOutput$bestK_PAC[[studyName]]
  message("unrounded PAC scores: ")
  cat("\n",unlist(consensusClusterOutput$PAC[[studyName]]))
  message("bestK using rounded PAC : ",consensusClusterOutput$bestK_PACR[[studyName]])
  PACR <- consensusClusterOutput$bestK_PACR[[studyName]]
  #just for comparision
  message("rounded PAC scores: ")
  cat("\n",unlist(consensusClusterOutput$PACR[[studyName]]))
  
  
  output <- data.frame(gapTest,consensusFrac,meanConsensus,PAC,PACR)
  
  return(output)
  
}
