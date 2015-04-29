library("proxy")
library("doParallel")
library("foreach")
#suggests:
library("energy")
library("limma")

CoINcIDE_getAdjMatrices <- function(dataMatrixList,clustSampleIndexList,clustFeatureIndexList,
edgeMethod=c("distCor","spearman","pearson","kendall","Euclidean","cosine",
                      "Manhattan","Minkowski","Mahalanobis"),numParallelCores=1,minTrueSimilThresh=-Inf,maxTrueSimilThresh=Inf,
sigMethod=c("meanMatrix","centroid"),maxNullFractSize=.2,numSims=500,includeRefClustInNull=TRUE,

outputFile="./CoINcIDE_messages.txt",fractFeatIntersectThresh=0,numFeatIntersectThresh=0 ,clustSizeThresh=0, clustSizeFractThresh=0,
checkNA=FALSE){
  
  
   #   
  #check: if rownames are null.
  #this can take a while to run if a long list of datasets.
   if(checkNA){
    if(any(is.na(unlist(dataMatrixList)))){
    
    stop("\nIn CoINcIDE_getAdjMatrices functions found NAs in a data matrix.
             \nPlease impute all NAs using the function filterAndImputeSamples
             function or an imputation method of your choice.")
  
    }
    
    
   }
    
  if(length(na.omit(match(edgeMethod,names(summary(pr_DB)[2]$distance))))!=1 && 
              all(is.na(match(edgeMethod,c("distCor","spearman","pearson","kendall"))))){
    
    stop("\nedgeMethod specified does not unique match one of the allowed methods in proxy
       or c(\"distCor\",\"spearman\",\"pearson\",\"kendall\")
       Type \'names(summary(pr_DB)[2]$distance)\' to see permitted string characters.
       Use the formal list name, not the other possible synonyms.")
    
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


message("This code assumes that you're clustering columns, and features are in the rows.");

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
#this is a symmetric matrix (with diagonal as NA)...so only count half of it.
  message(paste0(length(which(trueSimilData$similValueMatrix>=minTrueSimilThresh))/2," edges above minimum similarity threshold of ",minTrueSimilThresh))
  cat("\n",paste0(length(which(trueSimilData$similValueMatrix>=minTrueSimilThresh))/2," edges above minimum similarity threshold of ",minTrueSimilThresh),append=TRUE,file=outputFile)
  
  for(n in 1:nrow(clustIndexMatrix)){
    #cluster N is our reference cluster.
  nullSimilMatrix <- computeNullSimilVector(refClustRowIndex=n,dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,numSims=numSims,
                               trueSimilMatrix=trueSimilData$similValueMatrix,numParallelCores=numParallelCores,sigMethod=sigMethod,edgeMethod=edgeMethod,includeRefClustInNull=includeRefClustInNull,
                          clustIndexMatrix=clustIndexMatrix,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh)
  
   
  #cluster R is our comparison cluster.
   for(r in 1:nrow(clustIndexMatrix)){

     threshDir <- "greater"
     
     if(sigMethod=="meanMatrix"){
       #CHECK: is trueSimilData$similValueMatrix computing the correct numbers??
       thresh <- trueSimilData$similValueMatrix[n,r]
         
       if(!is.na(summary(pr_DB)[2]$distance[edgeMethod]) && summary(pr_DB)[2]$distance[edgeMethod]){
         
         threshDir <- "less"
         
       }
       
     }else if(sigMethod=="centroid"){
       
       thresh <- maxNullFractSize
     }
     
     
     #don't analyze clusters from the same dataset.
     if(clustIndexMatrix[n,2] != clustIndexMatrix[r,2] && !is.na(thresh)){
       
      if(!any(is.na(nullSimilMatrix[r,]))){
        
        pvalueMatrix[n,r] <- computeEdgePvalue(thresh=thresh,threshDir=threshDir,nullSimilVector = nullSimilMatrix[r,])
      
      }
     }
     

   }                 
    
  }

 output <- list(computeTrueSimilOutput=trueSimilData,pvalueMatrix=pvalueMatrix,clustIndexMatrix=clustIndexMatrix,inputVariablesDF=inputVariablesDF)
 return(output)
#EOF
}


######
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
                 clustSizeMatrix=clustSizeMatrix,clustSizeFractMatrix=clustSizeFractMatrix)
  
  return(output)
#EOF  
}

#nullSimilVector: for that clusterA-nullClust pair, with cluster A as reference,
#all similarities across all nullClusts.
computeEdgePvalue <- function(thresh,threshDir = c("greater","less"),nullSimilVector){
  
  if(is.na(thresh)){
    
    stop("\nthresh variable must be non-NA.")
    
  }
  
  if(any(is.na(nullSimilVector))){

    stop("\nSome of your p-values have NAs - check this out before continuing. 
         This can occur if your data has extremely low variance and you are using correlation metrics.")
    
  }

    #distance or cor matrix?
      if(threshDir=="less"){
      
          #do count NAs in full length?
      pvalue <- length(which(nullSimilVector <= thresh))/length(nullSimilVector)
      
      }else if(threshDir=="greater"){
        #want greater than for similarity. rest of metrics are simil metrics.
       pvalue <- length(which(nullSimilVector >= thresh))/length(nullSimilVector)  
        
      }else{
        
        stop("\nMust provide a value of \'less\' or \'greater\' for threshDir input.")
      
      }
  
  return(pvalue)
  
}

computeNullSimilVector <-   function(refClustRowIndex,dataMatrixList,clustSampleIndexList,clustFeatureIndexList,numSims=100,
                               trueSimilMatrix,numParallelCores=1,sigMethod,edgeMethod,includeRefClustInNull=TRUE,clustIndexMatrix,minTrueSimilThresh=-Inf,maxTrueSimilThresh=Inf){  
    
    registerDoParallel(cores=numParallelCores)


        refDataMatrix <- dataMatrixList[[as.numeric(clustIndexMatrix[refClustRowIndex,2])]]
        sampleIndices <- clustSampleIndexList[[as.numeric(clustIndexMatrix[refClustRowIndex,2])]][[as.numeric(clustIndexMatrix[refClustRowIndex,3])]]
        featureIndices <- clustFeatureIndexList[[as.numeric(clustIndexMatrix[refClustRowIndex,2])]][[as.numeric(clustIndexMatrix[refClustRowIndex,3])]]
        refClust <- refDataMatrix[featureIndices,sampleIndices]
        
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
          
          if(!all(rownames(nullClust) == rownames(refClust))){
            
            stop("In computeNullSimilVector: rownames of nullClust and refClust do not match up.")
            
          }
          
 
          #compare simil against each cluster
          foreach(c=1:nrow(clustIndexMatrix), .combine='c') %dopar%{
         # for(c in 1:nrow(clustIndexMatrix))  {
            #debug: is doParallel keeping the same nullclust in this loop, and is each nullClust truly random?
#             if(i==1){
#               #these should all be equal
#               cat("\n",mean(nullClust),file="~/Desktop/log.txt",append=TRUE)
# 
#             }
            
            if(!is.na(trueSimilMatrix[refClustRowIndex,c])){ 
                
              if(trueSimilMatrix[refClustRowIndex,c] >= minTrueSimilThresh && trueSimilMatrix[refClustRowIndex,c] <= maxTrueSimilThresh){
              
              sampleIndices <- clustSampleIndexList[[as.numeric(clustIndexMatrix[c,2])]][[as.numeric(clustIndexMatrix[c,3])]]
              featureIndices <- clustFeatureIndexList[[as.numeric(clustIndexMatrix[c,2])]][[as.numeric(clustIndexMatrix[c,3])]]
              compareClust <- dataMatrixList[[as.numeric(clustIndexMatrix[c,2])]][featureIndices,sampleIndices,drop=FALSE]  
              
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

  #centroidCompare=TRUE only matter for distCor, as it automatically returns a point value
  #and not a matrix if we feed in the entire refClust and compareClusts at once.
  similMatrix <- computeClusterPairSimil(refClust=compareClust,compareClust=centroids,
                                    edgeMethod=edgeMethod,centroidCompare=TRUE)
  
  
      #distance or cor matrix?
      #running this "if" statement for all simulations does slow down the code...
     #if only using a distance, or a simlarity metric, could remove the if statements to speed this up.
      if(!is.na(summary(pr_DB)[2]$distance[edgeMethod]) && summary(pr_DB)[2]$distance[edgeMethod]){
      
        #want greater than for similarity. rest of metrics are simil metrics.
         for (i in 1:nrow(similMatrix)) {
    
            sampleClass[i] <- which.min(similMatrix[i,])
    
        }
      
      }else{
        #want greater than for similarity. rest of metrics are simil metrics.
        #distCor is also a similarity metric - want highest one.
         for (i in 1:nrow(similMatrix)) {
    
            sampleClass[i] <- which.max(similMatrix[i,])
    
        }
  
        
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
#REMOVE these if statements/debugging
#make a separate function for each distance metric.
#have a function that returns the correct computeClusterPairSimil function, based off
#the edgeMethod
computeClusterPairSimil <- function(refClust,compareClust,
                                    edgeMethod="correlation",centroidCompare=FALSE){
  
  features <- intersect(rownames(refClust),rownames(compareClust))
  refClust <- refClust[rownames(refClust) %in% features, , drop=FALSE]
  compareClust <- compareClust[rownames(compareClust) %in% features, , drop=FALSE]
  
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
    
  }else if(edgeMethod=="distCor"){
    
    #number of rows must agree.

    if(!centroidCompare){
      
      clusterPairSimil <- dcor(refClust,compareClust,index=1.0)
      
    }else{
      #if the compareClust are actually two columns of centroids: want separate measurements for each one against each
      #sample in refClust
      clusterPairSimil <- matrix(data=NA,ncol=2,nrow=ncol(refClust))
      
      for(s in 1:ncol(refClust)){
        
        clusterPairSimil[s, ] <- cbind(dcor(refClust[,s],compareClust[,1],index=1.0),dcor(refClust[,s],compareClust[,2],index=1.0))
        
      }
      
    }
    
  }
  return(clusterPairSimil)
}

# #try having this function returned instead? decide it ahead of time, and pass it as an input
# #to any function that uses computeClusterPairSimil()
# #the returned function will have these inputs: refClust,compareClust,
# computeClusterPairSimilFunctionReturn <- function(
#                                     edgeMethod,centroidCompare=FALSE){
#   
#   func <- function(refClust,compareClust){
#   features <- intersect(rownames(refClust),rownames(compareClust))
#   refClust <- refClust[rownames(refClust) %in% features, , drop=FALSE]
#   compareClust <- compareClust[rownames(compareClust) %in% features, , drop=FALSE]
#   
#   
#   
#    
#   if(edgeMethod!="distCor"){
#     
#     
#     if(!is.na(summary(pr_DB)[2]$distance[edgeMethod])){
#       #distance or cor matrix?
#       if(summary(pr_DB)[2]$distance[edgeMethod]){
#         
#         
#         func <- function(refClust,compareClust){
#           features <- intersect(rownames(refClust),rownames(compareClust))
#           refClust <- refClust[rownames(refClust) %in% features, , drop=FALSE]
#           compareClust <- compareClust[rownames(compareClust) %in% features, , drop=FALSE]
#         clusterPairSimil <-   dist(t(refClust),t(compareClust),method=edgeMethod,by_rows=TRUE)
#         
#         return(clusterPairSimil)
#         
#         }
#         
#       }else{
#         
#         func <- function(refClust,compareClust){
#           features <- intersect(rownames(refClust),rownames(compareClust))
#           refClust <- refClust[rownames(refClust) %in% features, , drop=FALSE]
#           compareClust <- compareClust[rownames(compareClust) %in% features, , drop=FALSE]
#         clusterPairSimil <- simil(t(refClust),t(compareClust),method=edgeMethod,by_rows=TRUE)
#         
#         return(clusterPairSimil)
#         
#         }
#         
#       }
#       
#     }else{
#       
#       func <- function(refClust,compareClust){
#         features <- intersect(rownames(refClust),rownames(compareClust))
#         refClust <- refClust[rownames(refClust) %in% features, , drop=FALSE]
#         compareClust <- compareClust[rownames(compareClust) %in% features, , drop=FALSE]
#       clusterPairSimil <- cor(refClust,compareClust,method=edgeMethod)
#       
#       return(clusterPairSimil)
#       
#       }
#       
#     }
#     
#   }else if(edgeMethod=="distCor"){
#     
#     #number of rows must agree.
#     
#     if(!centroidCompare){
#       
#       func <- function(refClust,compareClust){
#         features <- intersect(rownames(refClust),rownames(compareClust))
#         refClust <- refClust[rownames(refClust) %in% features, , drop=FALSE]
#         compareClust <- compareClust[rownames(compareClust) %in% features, , drop=FALSE]
#       clusterPairSimil <- dcor(refClust,compareClust,index=1.0)
#       
#       return(clusterPairSimil)
#       
#     }else{
#       #if the compareClust are actually two columns of centroids: want separate measurements for each one against each
#       #sample in refClust
#       func <- function(refClust,compareClust){
#         features <- intersect(rownames(refClust),rownames(compareClust))
#         refClust <- refClust[rownames(refClust) %in% features, , drop=FALSE]
#         compareClust <- compareClust[rownames(compareClust) %in% features, , drop=FALSE]
#         
#       clusterPairSimil <- matrix(data=NA,ncol=2,nrow=ncol(refClust))
#       
#       for(s in 1:ncol(refClust)){
#         
#         clusterPairSimil[s, ] <- cbind(dcor(refClust[,s],compareClust[,1],index=1.0),dcor(refClust[,s],compareClust[,2],index=1.0))
#         
#       }
#       
#       return(clusterPairSimil)
#     }
#     
#     }
#   }
#   
#   }
#   return(function)
# }
# 


