library("proxy")
library("doParallel")
library("foreach")
#suggests:
library("energy")
library("limma")

CoINcIDE_getAdjMatrices <- function(dataMatrixList,clustSampleIndexList,clustFeatureIndexList,
edgeMethod=c("distCor","spearman","pearson","kendall","Euclidean","cosine",
                      "Manhattan","Minkowski","Mahalanobis"),numParallelCores=1,minTrueSimilThresh=-Inf,maxTrueSimilThresh=Inf,
sigMethod=c("meanMatrix","centroid"),maxNullFractSize=.1,numSims=100,includeRefClustInNull=TRUE,

outputFile="./CoINcIDE_messages.txt",minFractFeatureIntersect=0,fractFeatIntersectThresh=0,numFeatIntersectThresh=0 ,clustSizeThresh=0, clustSizeFractThresh=0){
  

  if(length(na.omit(match(edgeMethod,names(summary(pr_DB)[2]$distance))))!=1 && 
              all(is.na(match(edgeMethod,c("distCor","spearman","pearson","kendall"))))){
    
    stop("\nedgeMethod specified does not unique match one of the allowed methods in proxy
       or c(\"distCor\",\"spearman\",\"pearson\",\"kendall\")
       Type \'names(summary(pr_DB)[2]$distance)\' to see permitted string characters.
       Use the formal list name, not the other possible synonyms.")
    
  }

  
  if(length(sigMethod)>1){
    
    sigMethod <- "centroid"
    message("More than one sigMethod given; defaulting to meanMatrix method.")
    
  }

  if(sigMethod != "meanMatrix" && sigMethod != "centroid"){
    
    stop("\nPlease pick \'meanMatrix\' or \'centroid\' for sigMethod variable.")
  }
  
  if(length(clustSampleIndexList) != length(clustFeatureIndexList) || length(clustSampleIndexList) !=
       length(dataMatrixList) || length(clustFeatureIndexList) !=  length(dataMatrixList)){
    
    stop("The lenghts of your clustSampleIndexList,clustFeatureIndexList and dataMatrixList are not equal")
 
  }
  

  ##check: all clustIndex lists have featureIndexLists, and these
  #exist in the fullDataMatrix?
date <- Sys.time();
inputVariablesDF <- data.frame(date,edgeMethod,numParallelCores,minTrueSimilThresh,maxTrueSimilThresh,
sigMethod,maxNullFractSize,numSims,includeRefClustInNull,fractFeatIntersectThresh,
numFeatIntersectThresh,clustSizeThresh, clustSizeFractThresh);

#capture.output prints the data correctly.
capture.output(paste0("\nRunning find CoINcIDE_getAdjMatrices on ",Sys.time()," with the following inputs:\n"),append=TRUE,file=outputFile);
capture.output(inputVariablesDF,append=TRUE,file=outputFile);


message("This code assumes what you're clustering columns, and features are in the rows.");

  #grab total number of clusters
  clustDetailedNames <- array();
  numClust <- 0;
  
   cat("\nComputing cluster index matrix\n",append=TRUE,file=outputFile)

  for(d in 1:length(clustSampleIndexList)){
    
    if(length(clustSampleIndexList[[d]]) != length(clustFeatureIndexList[[d]])){
      
      stopMessage <- paste0("\nFor dataset ",d, ", the lengths of the clustSampleIndexList and clusterFeatureIndexList are not the same.")
      stop(stopMessage)
      
    }
    
    if(length(clustSampleIndexList[[d]])>0){

      #if empty: will skip over it.
      for(a in 1:length(clustSampleIndexList[[d]])){
          
        numClust <- numClust+1;
        clustDetailedNames[numClust] <- paste0(numClust,"_",d,"_",a);
      }
      
    }
    
  }
  
  clustIndexMatrix <- strsplit2(clustDetailedNames,"_");
  cat("\nThere are ",nrow(clustIndexMatrix), " total clusters across ",length(dataMatrixList), "datasets.\n",append=TRUE,file=outputFile)
  
  #compute baseline similiarites, feature/sample stats for clusters.
  cat("\nComputing cluster-cluster true similarities (or distances).\n",append=TRUE,file=outputFile)
  trueSimilData <- computeTrueSimil(clustIndexMatrix=clustIndexMatrix,edgeMethod=edgeMethod,
                                    dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,
                                  clustFeatureIndexList=clustFeatureIndexList,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                  numFeatIntersectThresh=numFeatIntersectThresh ,clustSizeThresh=clustSizeThresh, clustSizeFractThresh=clustSizeFractThresh)
  
 

  pvalueMatrix <- matrix(data=NA,nrow=numClust,ncol=numClust)
  cat("\nComputing p-values for each cluster-cluster similarity using null cluster distributions.\n",append=TRUE,file=outputFile)
  
  for(n in 1:nrow(clustIndexMatrix)){
    
  nullSimilMatrix <- computeNullSimilVector(refClustRowIndex=n,dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,numSims=numSims,
                               trueSimilMatrix=trueSimilData$similValueMatrix,numParallelCores=numParallelCores,sigMethod=sigMethod,edgeMethod=edgeMethod,includeRefClustInNull=includeRefClustInNull,
                          clustIndexMatrix=clustIndexMatrix,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh)
  
   
   for(r in 1:nrow(clustIndexMatrix)){

     if(sigMethod=="meanMatrix"){
       #CHECK: is trueSimilData$similValueMatrix computing the correct numbers??
       thresh <- trueSimilData$similValueMatrix[n,r]
         
     }else if(sigMethod=="centroid"){
       
       thresh <- maxNullFractSize
     }
     #don't analyze clusters from the same dataset.
     if(clustIndexMatrix[n,2] != clustIndexMatrix[r,2]){
       
      if(!any(is.na(nullSimilMatrix[r,]))){
        
        pvalueMatrix[n,r] <- computeEdgePvalue(thresh=thresh,nullSimilVector = nullSimilMatrix[r,],edgeMethod=edgeMethod)
      
      }
     }
     

   }                 
    
  }

 output <- list(trueSimilData=trueSimilData,pvalueMatrix=pvalueMatrix,clustIndexMatrix=clustIndexMatrix)
 return(output)
#EOF
}


######A
#CHECK: is there a bug with indexing here??
computeTrueSimil <- function(clustIndexMatrix,
                                    edgeMethod="correlation",
                                    dataMatrixList,clustSampleIndexList,clustFeatureIndexList,fractFeatIntersectThresh=0,numFeatIntersectThresh=0 ,clustSizeThresh=0, clustSizeFractThresh=0){
  
    numClust <- nrow(clustIndexMatrix)
    similValueMatrix <- matrix(data=NA,nrow=numClust,ncol=numClust)
    numFeatIntersectMatrix <- matrix(data=NA,nrow=numClust,ncol=numClust)
    fractFeatIntersectMatrix <- matrix(data=NA,nrow=numClust,ncol=numClust)
    clustSizeMatrix <- array(data=NA,dim=numClust)
    clustSizeFractMatrix <- array(data=NA,dim=numClust)
    
  for(i in 1:nrow(clustIndexMatrix)){
    
    sampleIndices <- clustSampleIndexList[[as.numeric(clustIndexMatrix[i,2])]][[as.numeric(clustIndexMatrix[i,3])]]
    featureIndices <- clustFeatureIndexList[[as.numeric(clustIndexMatrix[i,2])]][[as.numeric(clustIndexMatrix[i,3])]]
    clustSizeMatrix[as.numeric(clustIndexMatrix[i,1])] <- ncol(dataMatrixList[[as.numeric(clustIndexMatrix[i,2])]][featureIndices,sampleIndices]) 
    clustSizeFractMatrix[as.numeric(clustIndexMatrix[i,1])] <- ncol(dataMatrixList[[as.numeric(clustIndexMatrix[i,2])]][featureIndices,sampleIndices])/ncol(dataMatrixList[[as.numeric(clustIndexMatrix[i,2])]])
  
  }

  
  for(r in 1:nrow(clustIndexMatrix)){
    
    sampleIndices <- clustSampleIndexList[[as.numeric(clustIndexMatrix[r,2])]][[as.numeric(clustIndexMatrix[r,3])]]
    featureIndices <- clustFeatureIndexList[[as.numeric(clustIndexMatrix[r,2])]][[as.numeric(clustIndexMatrix[r,3])]]
    refClust <- dataMatrixList[[as.numeric(clustIndexMatrix[r,2])]][featureIndices,sampleIndices]
    
    refStudyNum <- as.numeric(clustIndexMatrix[r,2])

    for(c in 1:nrow(clustIndexMatrix)){
    
      compareStudyNum <- as.numeric(clustIndexMatrix[c,2])
      
      #in the same study? (if wanted to include clusters in same study, 
      #just do if(r!=c) to avoid comparing the exact same cluster against itself.)
      if(refStudyNum!=compareStudyNum){
        
        #already looked at this pair? won't be NA then.
        if(is.na(numFeatIntersectMatrix[r,c])){
      
          sampleIndices <- clustSampleIndexList[[as.numeric(clustIndexMatrix[c,2])]][[as.numeric(clustIndexMatrix[c,3])]]
          featureIndices <- clustFeatureIndexList[[as.numeric(clustIndexMatrix[c,2])]][[as.numeric(clustIndexMatrix[c,3])]]
          compareClust <- dataMatrixList[[as.numeric(clustIndexMatrix[c,2])]][featureIndices,sampleIndices]  
          numFeatIntersectMatrix[r,c] <- length(intersect(rownames(refClust),rownames(compareClust)))
          numFeatIntersectMatrix[c,r] <- numFeatIntersectMatrix[r,c]
          fractFeatIntersectMatrix[r,c] <- numFeatIntersectMatrix[r,c]/length(union(rownames(refClust),rownames(compareClust)))
          fractFeatIntersectMatrix[c,r] <-    fractFeatIntersectMatrix[r,c]  
       

          if(fractFeatIntersectMatrix[r,c] >= fractFeatIntersectThresh && numFeatIntersectMatrix[r,c] >= numFeatIntersectThresh && clustSizeMatrix[r] >= clustSizeThresh && clustSizeFractMatrix[r] >= clustSizeFractThresh
              && clustSizeMatrix[c] >= clustSizeThresh && clustSizeFractMatrix[c] >= clustSizeFractThresh){
  
            similValueMatrix[r,c] <- computeClusterPairSimil_mean(refClust,compareClust,edgeMethod=edgeMethod)
            similValueMatrix[c,r] <-   similValueMatrix[r,c]  
      
          }
      
        }

      }
    
    }

  }


  output <- list(similValueMatrix=similValueMatrix,numFeatIntersectMatrix=numFeatIntersectMatrix,fractFeatIntersectMatrix=fractFeatIntersectMatrix,
                 clustSizeMatrix=clustSizeMatrix,clustSizeFractMatrix)
  
  return(output)
#EOF  
}


computeEdgePvalue <- function(thresh,nullSimilVector,edgeMethod){
  
  if(is.na(thresh)){
    
    stop("\nthresh variable must be non-NA.")
    
  }
  
  if(any(is.na(nullSimilVector))){

    stop("\nSome of your p-values have NAs - check this out before continuing. 
         This can occur if your data has extremely low variance and you are using correlation metrics.")
    
  }
  
  if(length(na.omit(match(edgeMethod,names(summary(pr_DB)[2]$distance))))!=1 && 
              all(is.na(match(edgeMethod,c("distCor","spearman","pearson","kendall"))))){
    
    stop("\nedgeMethod specified does not unique match one of the allowed methods in proxy
       or c(\"distCor\",\"spearman\",\"pearson\",\"kendall\")
       Type \'names(summary(pr_DB)[2]$distance)\' to see permitted string characters.
       Use the formal list name, not the other possible synonyms.")
    
  }

    #distance or cor matrix?
      if(!is.na(summary(pr_DB)[2]$distance[edgeMethod]) && summary(pr_DB)[2]$distance[edgeMethod]){
      
          #do count NAs in full length?
      pvalue <- length(which(nullSimilVector <= thresh))/length(nullSimilVector)
      
      }else{
        #want greater than for similarity. rest of metrics are simil metrics.
       pvalue <- length(which(nullSimilVector >= thresh))/length(nullSimilVector)  
        
      }
  
  return(pvalue)
  
}

computeNullSimilVector <-   function(refClustRowIndex,dataMatrixList,clustSampleIndexList,clustFeatureIndexList,numSims=100,
                               trueSimilMatrix,numParallelCores=1,sigMethod,edgeMethod,includeRefClustInNull=TRUE,clustIndexMatrix,minTrueSimilThresh=-Inf,maxTrueSimilThresh=Inf){
  
  if(numParallelCores>1){
    
    registerDoParallel(cores=numParallelCores)
  
  }

        refDataMatrix <- dataMatrixList[[as.numeric(clustIndexMatrix[refClustRowIndex,2])]]
        sampleIndices <- clustSampleIndexList[[as.numeric(clustIndexMatrix[refClustRowIndex,2])]][[as.numeric(clustIndexMatrix[refClustRowIndex,3])]]
        featureIndices <- clustFeatureIndexList[[as.numeric(clustIndexMatrix[refClustRowIndex,2])]][[as.numeric(clustIndexMatrix[refClustRowIndex,3])]]
        refClust <- dataMatrixList[[as.numeric(clustIndexMatrix[refClustRowIndex,2])]][featureIndices,sampleIndices]
        
        #will return a matrix of numClusters x numFeatures
        #see documention on CRAN "nesting foreach loops"
        #"nesting" way of dopar with %:% didn't work too well. foreach can be weird with iterative random numbers, so play it safe.
        nullSimilVector <- foreach(i=1:numSims,.combine='cbind') %dopar% {
          
          
          #want same set of nullClusters across ALL comparison clusters.
          #so hold each null cluster "constant", loop through compare clusters.
          if(includeRefClustInNull){
          
            nullClust <- refDataMatrix[rownames(refClust),sample.int(ncol(refDataMatrix),size=ncol(refClust),replace=TRUE), drop=FALSE];
            
          }else{
            
            #note: if only a few samples NOT in the ref cluster, this will mostly be duplicated samples
            colIndicesNoRefClust <- na.omit(match(setdiff(colnames(refDataMatrix),colnames(refClust)),colnames(refDataMatrix)))
            nullClust <- refDataMatrix[rownames(refClust),sample(colIndicesNoRefClust,size=ncol(refClust),replace=TRUE), drop=FALSE];
            
          }
          
          #compare simil against each cluster
          foreach(c=1:nrow(clustIndexMatrix), .combine='c') %dopar%{
            
            #debug: is doParallel keeping the same nullclust in this loop, and is each nullClust truly random?
#             if(i==1){
#               #these should all be equal
#               cat("\n",mean(nullClust),file="~/Desktop/log.txt",append=TRUE)
# 
#             }
            
            if(!is.na(trueSimilMatrix[refClustRowIndex,c])){ 
                
              if(trueSimilMatrix[refClustRowIndex,c] >= minTrueSimilThresh && trueSimilMatrix[refClustRowIndex,c] <= maxTrueSimilThresh)
              
              sampleIndices <- clustSampleIndexList[[as.numeric(clustIndexMatrix[c,2])]][[as.numeric(clustIndexMatrix[c,3])]]
              featureIndices <- clustFeatureIndexList[[as.numeric(clustIndexMatrix[c,2])]][[as.numeric(clustIndexMatrix[c,3])]]
              compareClust <- dataMatrixList[[as.numeric(clustIndexMatrix[c,2])]][featureIndices,sampleIndices]  
              
              if(sigMethod=="meanMatrix"){
                
              similValue <- computeClusterPairSimil_mean(refClust=nullClust,compareClust=compareClust,edgeMethod=edgeMethod)
              
              }else if(sigMethod=="centroid"){
                #even if only one sample, bc used "drop=FALSE" above, should still be a matrix and can take rowMeans.
              refCentroid <- rowMeans(refClust)
              nullCentroid <- rowMeans(nullClust)
              similValue <- computeClusterPairSimil_centroid(refCentroid=refCentroid,nullCentroid=nullCentroid,compareClust=compareClust,edgeMethod=edgeMethod)
              
                
              }

            
          }else{
            
            similValue <- NA
          
          }  
        #end of inner foreach dopar 
        }
    #end of outer foreach dopar
    }

}#EOF

computeClusterPairSimil_mean <- function(refClust,compareClust,
                                    edgeMethod="correlation"){
  
  
  similMatrix <- computeClusterPairSimil(refClust,compareClust,
                                    edgeMethod=edgeMethod)
  
  
  return(mean(similMatrix))
  
  
  
}

#COME BACK: compute centroids correctly?
computeClusterPairSimil_centroid <- function(compareClust,refCentroid,nullCentroid,
                                    edgeMethod="correlation",numParallelCores=1){
  
  centroids <- cbind(refCentroid,nullCentroid)
  #below code adapted from IGP.clusterRepro function in clusterRepro package
  #added other simil metrics to code..
  sampleClass <- rep(NA, ncol(compareClust))
  names(sampleClass) <- colnames(compareClust)
  
  if(!all(rownames(centroids)==rownames(compareClust))){
    
    stop("\nNot all row (feature) names match up.")
    
  }

  
  similMatrix <- computeClusterPairSimil(refClust=compareClust,compareClust=centroids,
                                    edgeMethod=edgeMethod)
  
  for (i in 1:nrow(similMatrix)) {
    
    sampleClass[i] <- which.max(similMatrix[i,])
    
  }
  
  #clusterRepro code we didn't use:
  #result <- c()
  #IGP.clusterRepro from clusterRepro package returns these metrics.
  #we don't need these.
  #result$class <- sampleClass
  #result$IGP <- diag(classification)/rowSums(classification)
  #result$size <- rowSums(classification)
  #Fract is not in clusterRepro, but similar to Size.
  #result$fract <- rowSums(classification)/ncol(compareClust)
  
  #want fraction assigned to second centroid (false centroid.)
  fract <- length(which(sampleClass==2))/ncol(compareClust)
  
  return(fract)

}
 
#NEXT: update simulated data code (are some tissue samples being mixed?  last cluster looked weird...)
#NEXT: still debug a bit: re-run simulation datasets a few times.
#why some p-value results are 0, some are 0.00? Rounds? perhaps..
#centroid method NOT working with "distCor"
#centroid method looks like it has a bug - giving almost reverse results?

#NEXT: community detection, with stats (studies dropped, their attributes, etc.) with 1 p-value

#THEN: visualization

#THEN: meta-analysis.

#THEN:  bagging (bootstrap aggregating) for a classifier? still need to normalize??
  #PERHAPS: normalize by clusters, and NOT entire dataset? (either cluster 6, or not cluster 6.)
  #"univariate" logistic models for each cluster (for all clusters, in all datasets.)
  #group each univariate cluster logistic model into its meta-cluster group.
  #a sample is classified to the meta-cluster whose models most vote it "yes" and not "no".

#Votes: specific cluster, or none at all (how have none at all? perhaps if the probability outputted is very low.)

#THEN: meta-variance code, kmeans.

computeClusterPairSimil <- function(refClust,compareClust,
                                    edgeMethod="correlation"){

  if(!all(rownames(refClust)==rownames(compareClust))){
    
    stop("\nNot all row (feature) names match up.")
    
  }
  
  if(length(features)==0){
    #no features to analyze!
    return(NA)
    
  }
  
  if(length(na.omit(match(edgeMethod,names(summary(pr_DB)[2]$distance))))!=1 && 
              all(is.na(match(edgeMethod,c("distCor","spearman","pearson","kendall"))))){
    
    stop("\nedgeMethod specified does not unique match one of the allowed methods in proxy
       or c(\"distCor\",\"spearman\",\"pearson\",\"kendall\")
       Type \'names(summary(pr_DB)[2]$distance)\' to see permitted string characters.
       Use the formal list name, not the other possible synonyms.")
    
  }

    
    if(edgeMethod!="distCor"){
      
      if(!is.na(summary(pr_DB)[2]$distance[edgeMethod])){
      #distance or cor matrix?
      if(summary(pr_DB)[2]$distance[edgeMethod]){
        
         clusterPairSimil <-   dist(t(refClust),t(compareClust),method=edgeMethod,by_rows=TRUE)
        
      }else{
        
         clusterPairSimil <- simil(t(refClust),t(compareClust),method=edgeMethod,by_rows=TRUE)
        
      }
  
  }else{
    
     clusterPairSimil <- cor(refClust,compareClust,method=edgeMethod)
    
  }
    
  }else{
    
    #number of rows must agree.
    #I have found this metric does NOT actually work too well in separating
    #gene expression clusters.
    clusterPairSimil <- dcor(refClust,compareClust,index=1.0)
    
  }
  
}

