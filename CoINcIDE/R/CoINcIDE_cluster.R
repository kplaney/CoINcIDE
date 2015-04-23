#use default Ward's method
library("cluster")
library("ConsensusClusterPlus")

#for more specific options: use the individual clustMatrixList functions
#note: for consenus clustering: need a high # sims otherwise will get NaN values in your consensus matrix (two samples were next chosen in the same resampling run.)
clustMatrixListWrapper <- function(dataMatrixList,clustFeaturesList,clustMethod=c("km","hc"),pickKMethod=c("gap","consensus"),numSims=1000,maxNumClusters=30,
                                   outputFile="./cluster_output.txt",iter.max=30,nstart=25,distMethod=c("euclidean","pearson","spearman", "binary", "maximum", "canberra", "minkowski"),
                                   hclustAlgorithm=c("average","complete","ward.D", "ward.D2", "single", "mcquitty","median","centroid"), 
                                   consensusHclustAlgorithm=c("average","complete","ward.D", "ward.D2", "single", "mcquitty","median","centroid"),
                                   minClustConsensus=.7, minMeanClustConsensus=.7,corUse="everything",pItem=.9){
  
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
      
      clusterAlg <- clustMethod
    
    clusterOutputList <- consensusClusterMatrixList(dataMatrixList=dataMatrixList,clustFeaturesList=clustFeaturesList,maxNumClusters = maxNumClusters, 
                                                    numSims=numSims, clusterAlg=clusterAlg,
                                           hclustAlgorithm=hclustAlgorithm, consensusHclustAlgorithm=consensusHclustAlgorithm,
                                           corUse=corUse,pItem=pItem,pFeature=1,
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
distMethod=c("pearson","spearman","euclidean", "binary", "maximum", "canberra", "minkowski"),minMeanClustConsensus=.7,minClustConsensus=.7,
outputFile="./consensusOut.txt"){
  
 
  outputList <- list()
  
  for(d in 1:length(dataMatrixList)){
    
    message(paste0("\nClustering dataset ",d, " ",names(dataMatrixList)[d]))
    cat(paste0("\nClustering dataset ",d, " ",names(dataMatrixList)[d]),
        append=TRUE,file=outputFile);
    outputList[[d]] <- consensusClusterMatrix(dataMatrix=dataMatrixList[[d]],clustFeatures=clustFeaturesList[[d]],
                                                     maxNumClusters=maxNumClusters, numSims=numSims, 
                                                     pItem=pItem, pFeature=pFeature, clusterAlg=clusterAlg,
                                              hclustAlgorithm=hclustAlgorithm, consensusHclustAlgorithm=consensusHclustAlgorithm,corUse=corUse,distMethod=distMethod,
                                              minMeanClustConsensus=minMeanClustConsensus,minClustConsensus=minClustConsensus,outputFile=outputFile)
    

  
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
  
  consensusInfo <- list()
  minConsensusClusterByK <- list()
  meanConsensusClusterByK <- list()
  
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
    
    consensusInfo[[d]] <- outputList[[d]]$consensusCalc
    minConsensusClusterByK[[d]] <- outputList[[d]]$minConsensusClusterByK
    meanConsensusClusterByK[[d]] <- outputList[[d]]$meanConsensusClusterByK
    
  }
   
  output <- list(  clustSampleIndexList_consensusFrac=clustSampleIndexList_consensusFrac,clustFeatureIndexList_consensusFrac =clustFeatureIndexList_consensusFrac,
                   bestK_consensusFrac=bestK_consensusFrac,
                   clustSampleIndexList_meanConsensusCluster=clustSampleIndexList_meanConsensusCluster,
                   clustFeatureIndexList_meanConsensusCluster=clustFeatureIndexList_meanConsensusCluster,
                   bestK_meanConsensusCluster=bestK_meanConsensusCluster,clustSampleIndexList_minConsensusCluster=clustSampleIndexList_minConsensusCluster,
                   clustFeatureIndexList_minConsensusCluster=clustFeatureIndexList_minConsensusCluster,
                   bestK_meanConsensusCluster=bestK_meanConsensusCluster,consensusInfo=consensusInfo,minConsensusClusterByK=minConsensusClusterByK,
                   meanConsensusClusterByK=meanConsensusClusterByK)
  
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
corUse=c("everything","pairwise.complete.obs", "complete.obs"),
distMethod=c("euclidean","spearman","euclidean", "pearson","binary", "maximum", "canberra", "minkowski"),
minMeanClustConsensus=.7,
outputFile="./consensusOut.txt",minClustConsensus=.7){
  
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
  
  #remember: we're clustering pItem*ncol(dataset) each time...so
  #our max K is NOT ncol(dataset), but rather pItem*ncol(dataset)
  if(floor(ncol(dataset)*pItem)<maxNumClusters){
    #hclust usually returns NA gap test if K.max = pItem*ncol(dataset) as opposed to pItem*ncol(dataset)-1
    #will also mest up maxSE calculations
    K.max <- floor(pItem*ncol(dataset))-1
    
  }else{
    
    K.max <- maxNumClusters
    
  }


  #too bad we can't increase the number of starts for kmeans. (user-defined functions require a distance matrix as input,
 #so can't write own k-means function.)   
  consensusClustOutput <- ConsensusClusterPlus(d=dataset,
                                               maxK=K.max,reps=numSims,
                                               pItem=pItem, pFeature=pFeature, clusterAlg=clusterAlg,
  innerLinkage=innerLinkage, finalLinkage=finalLinkage, ml=NULL,distance=distMethod,
  tmyPal=NULL,seed=NULL,plot='pngBMP',writeTable=FALSE,weightsItem=NULL,weightsFeature=NULL,verbose=F,corUse=corUse)
  

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
  if(!(all(unlist(meanConsensusClusterByK)=="NaN")) && any(unlist(meanConsensusClusterByK)>=minMeanClustConsensus) && any(unlist(minConsensusClusterByK)>=minClustConsensus)){
    
    meanBetweenConsensusClusterByK <- list()
    consensusFrac <- list()
    consensusMetric <- list()
    PAC <- list()
    
    N <- nrow(consensusClustOutput[[i]]$consensusMatrix)
    for(i in 2:K.max){
      #technically could just take mean - is a symmetric matrix. diag=FALSE bc it is zero here.
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
      
      #if(i >2){
      
      #consensusMetric[[i]] <- consensusFrac[[i]] > consensusFrac[[(i-1)]] 
      
      #}
    }
  
  
    consensusFrac[[1]] <- NA

    ##consensusFrac calculations
    #which.max(): Missing and NaN values are discarded.
    selectK_consensusFrac$bestK <- which.max(unlist(consensusFrac))
    selectK_minConsensusClust$bestK <- which.max(unlist(minConsensusClusterByK))+1
    selectK_meanConsensusClust$bestK <- which.max(unlist(meanConsensusClusterByK))+1
    selectK_PAC$bestK  <- min(unlit(PAC))
    
    if(length(selectK_consensusFrac$bestK)==0){
      #all meanConsensus NAs returned; set K=1
      selectK_consensusFrac$bestK <- 1
      
    }
    
    if(length(selectK_minConsensusClust$bestK)==0){
      #all meanConsensus NAs returned; set K=1
      selectK_minConsensusClust$bestK <- 1
      
    }
    
    if(length(selectK_meanConsensusClust$bestK)==0){
      #all meanConsensus NAs returned; set K=1
      selectK_meanConsensusClust$bestK <- 1
      
    }
    
    if(length( selectK_PAC$bestK)==0){
      #all meanConsensus NAs returned; set K=1
      selectK_PAC$bestK <- 1
      
    }

  
  }else{
    #no clusterings passed the minMeanConsensus threshold
    selectK_consensusFrac$bestK <- 1
    consensusFrac <- NA
    meanBetweenConsensusClusterByK <- NA
    PAC <- NA
    selectK_meanConsensusClust$bestK <- 1
    selectK_minConsensusClust$bestK  <- 1
    selectK_PAC$bestK <- 1
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
 
 
 message(paste0("Best K as determined by PAC is: ",selectK_PAC$bestK))
 cat(paste0("\nBest K as determined by ",clusterAlg, " and PAC  is: ", selectK_PAC$bestK) ,"\n"),
     append=TRUE,file=outputFile);
 
 
packageClustOutput <-  function(listOutputObject){
  
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

selectK_consensusFrac <- packageClustOutput(selectK_consensusFrac)
selectK_meanConsensusClust <- packageClustOutput(selectK_meanConsensusClust)
selectK_minConsensusClust <- packageClustOutput(selectK_minConsensusClust)
selectK_PAC <- packageClustOutput(selectK_PAC)

  output <- list(selectK_consensusFrac=selectK_consensusFrac, selectK_PAC=selectK_PAC 
                 selectK_meanConsensusClust=selectK_meanConsensusClust,
                 meanConsensusClusterByK=meanConsensusClusterByK,selectK_minConsensusClust=selectK_minConsensusClust,
                 consensusCalc=icl,consensusByK=consensusByK,
                 consensusClustOutput=consensusClustOutput,consensusFrac=consensusFrac,
                 meanBetweenConsensusClusterByK=meanBetweenConsensusClusterByK,minConsensusClusterByK=minConsensusClusterByK,PAC=PAC)


 return(output)
 
}


