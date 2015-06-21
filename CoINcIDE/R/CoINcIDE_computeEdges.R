library("proxy")
library("doParallel")
library("foreach")
#suggests:
library("energy")
library("limma")
library("matrixStats")
#Rclusterpp does not return the distance matrix.
#library("Rclusterpp")

CoINcIDE_getAdjMatrices <- function(dataMatrixList,clustSampleIndexList,clustFeatureIndexList,
edgeMethod=c("distCor","spearman","pearson","kendall","Euclidean","cosine",
                      "Manhattan","Minkowski","Mahalanobis"),numParallelCores=1,minTrueSimilThresh=-Inf,maxTrueSimilThresh=Inf,
sigMethod=c("meanMatrix","centroid"),maxNullFractSize=.2,numSims=500,includeRefClustInNull=TRUE,
outputFile="./CoINcIDE_messages.txt",fractFeatIntersectThresh=0,numFeatIntersectThresh=0 ,clustSizeThresh=0, clustSizeFractThresh=0,
checkNA=FALSE,centroidMethod=c("mean","median")){
  
    if(length(centroidMethod)>1){
      
      centroidMethod<-c("mean")
      
    }
    
  if(edgeMethod=="distCor"){
    
    
    message("Note: because of how it is computed, edgeMethod=distCor may result in a long run time over 24 hours.")
    
  }
  
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
  sigMethod,numSims,includeRefClustInNull,fractFeatIntersectThresh,
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
  message("Computing cluster-cluster true similarities (or distances).")
  trueSimilData <- computeTrueSimil(clustIndexMatrix=clustIndexMatrix,edgeMethod=edgeMethod,
                                    dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,
                                  clustFeatureIndexList=clustFeatureIndexList,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                  numFeatIntersectThresh=numFeatIntersectThresh ,clustSizeThresh=clustSizeThresh, clustSizeFractThresh=clustSizeFractThresh)
  


  message(paste0(length(which(trueSimilData$similValueMatrix>=minTrueSimilThresh))/2," edges above minimum similarity threshold of ",minTrueSimilThresh))
  cat("\n",paste0(length(which(trueSimilData$similValueMatrix>=minTrueSimilThresh))/2," edges above minimum similarity threshold of ",minTrueSimilThresh),append=TRUE,file=outputFile)
  
  cat("\nComputing similarity matrices for null/permutation calculations\n",append=TRUE,file=outputFile)
  message("Computing similarity matrices for null/permutation calculations")

  
  pvalueMatrix <- matrix(data=NA,nrow=numClust,ncol=numClust)

 
  cat("\nComputing p-values for each cluster-cluster similarity using null cluster distributions.\n",append=TRUE,file=outputFile)
  

      #COME BACK: also add median option?
      if(!is.na((summary(pr_DB)[2]$distance[edgeMethod])) && summary(pr_DB)[2]$distance[edgeMethod]){
        #these functions are written later in the script
        
        centroidFunction <-   computeClusterPairAssignFract_matrixStatsDist
 
        
        #this means it's in matrix stats, but is a correlation method
      }else if(!is.na((summary(pr_DB)[2]$distance[edgeMethod])) && !summary(pr_DB)[2]$distance[edgeMethod]){
        
   
          
        centroidFunction <-   computeClusterPairAssignFract_matrixStatsSimil

        
      }else if(edgeMethod=="distCor"){
        
        centroidFunction <- computeClusterPairAssignFract_dcor
        
        
      }else{

        
       centroidFunction <-  computeClusterPairAssignFract_cor
       

      }


  trueFractNNmatrix <- matrix(data=NA,nrow=numClust,ncol=numClust)
 
        #COME BACK: do foreach loop here.
  for(n in 1:nrow(dataMatrixList)){
    
    #for each dataset: make centroid set and null centroid sets
          centroidMatrixOrig <- matrix(data=NA,ncol=length(clustFeatureIndexList[[n]]),nrow=length(clustFeatureIndexList[[1]]))
          
  if(centroidMethod=="mean"){

        for(d in 1:length(clustFeatureIndexList[[n]])){
          
        centroidMatrixOrig[,d] <- rowMeans(dataMatrixList[[n]][ clustFeatureIndexList[[1]], clustSampleIndexList[[n]][[d]] ])
        
        }
        #otherwise: just median. 
      }else{

                for(d in 1:length(clustFeatureIndexList[[n]])){
          
        centroidMatrixOrig[,d] <- rowMedians(dataMatrixList[[n]][ clustFeatureIndexList[[1]], clustSampleIndexList[[n]][[d]] ])
        
        }
  
  }
    nullCentroidList <- createNullCentroidMatrixList(centroidMatrix,numIter=numIter)
        
      #can foreach work here? perhaps if I don't combine.
      # pvalueMatrix[, c] <- foreach(i=1:numSims) %dopar% {
       for(c in 1:nrow(clustIndexMatrix)){

         if(as.numeric(clustIndexMatrix[c,2]) != n){
                    
          message("Running tests for cluster number: ",n," from dataset ",as.numeric(clustIndexMatrix[c,2]))
           ###ALSO: make sure passes feature thresh, etc.
         sampleIndices <- clustSampleIndexList[[as.numeric(clustIndexMatrix[c,2])]][[as.numeric(clustIndexMatrix[c,3])]]
        features <- intersect(dataMatrixList[[as.numeric(clustIndexMatrix[c,2])]][ clustFeatureIndexList[[1]] ],
                                    #just pick first feature index for now.
                                    dataMatrixList[[n]][ clustFeatureIndexList[[1]] ])
        centroidMatrix <- centroidMatrixOrig[features, ]
        
        compareClust <- dataMatrixList[[as.numeric(clustIndexMatrix[c,2])]][features,sampleIndices,drop=FALSE]  
   
         fractNNresults <- centroidFunction(compareMatrix=compareClust,centroidMatrix=centroidMatrix,edgeMethod=edgeMethod)
         bestMatchIndex <- fractNNresults$bestMatch
         globalClustIndex <- intersect(which(clustIndexMatrix[,2,]==n),which(clustIndexMatrix[, ,3]==bestMatchIndex))

         
         trueFractNNmatrix[globalClustIndex,r] <- fractNNresults$bestFract
         #NOW: run through null centroid list.
         
         nullTests <- lapply(nullCentroidList,FUN=function(nullCentroidMatrix,compareClust,edgeMethod,thresh){
           
           passedThresh <- FALSE
           tmp <- centroidFunction(compareMatrix=compareClust,centroidMatrix=centroidMatrix,edgeMethod=edgeMethod)
           
           if(tmp$bestFract>=thresh){
             
             passedThresh <- TRUE
           }    
           
         },compareClust=compareClust,edgeMethod=edgeMethod,thresh=trueFractNNmatrix[globalClustIndex,r])
        
        pvalueMatrix[globalClustIndex,c] <- length(which(unlist(nullTests)))/length(nullCentroidList)
    
  #done computing metrics for this cluster.
  }#else{
    #add else so that foreach loop still counts this column? not needed for now.
  #  pvalueMatrix[globalClustIndex,c] <- NA
  #}

}

#make mean edge p-value? or do this in summary functions?
   
   output <- list(computeTrueSimilOutput=trueSimilData,pvalueMatrix=pvalueMatrix,
                  clustIndexMatrix=clustIndexMatrix,inputVariablesDF=inputVariablesDF,
                  trueFractNNmatrix=trueFractNNmatrix)
   
 }
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
    #drop=FALSE: in case a cluster has 1 sample.
    clustSizeMatrix[as.numeric(clustIndexMatrix[i,1])] <- ncol(dataMatrixList[[as.numeric(clustIndexMatrix[i,2])]][featureIndices,sampleIndices,drop=FALSE]) 
    clustSizeFractMatrix[as.numeric(clustIndexMatrix[i,1])] <- ncol(dataMatrixList[[as.numeric(clustIndexMatrix[i,2])]][featureIndices,sampleIndices,drop=FALSE])/ncol(dataMatrixList[[as.numeric(clustIndexMatrix[i,2])]])
  
  }

  
  for(r in 1:nrow(clustIndexMatrix)){
    
    sampleIndices <- clustSampleIndexList[[as.numeric(clustIndexMatrix[r,2])]][[as.numeric(clustIndexMatrix[r,3])]]
    featureIndices <- clustFeatureIndexList[[as.numeric(clustIndexMatrix[r,2])]][[as.numeric(clustIndexMatrix[r,3])]]
    refClust <- dataMatrixList[[as.numeric(clustIndexMatrix[r,2])]][featureIndices,sampleIndices,drop=FALSE]
    
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
          compareClust <- dataMatrixList[[as.numeric(clustIndexMatrix[c,2])]][featureIndices,sampleIndices,drop=FALSE]  
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

####


computeClusterPairSimil_mean <- function(refClust,compareClust,
                                         edgeMethod="correlation"){
  
  
  similMatrix <- computeClusterPairSimil(refClust,compareClust,
                                         edgeMethod=edgeMethod)
  

  return(mean(similMatrix))
  
  
  
}
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

#warning...this can be slow!
# computeClusterPairAssignFract_dcor <- function(compareMatrix,centroidMatrix,edgeMethod="distCor"){
#   
#       # centroidMatrix is two columns of centroids: want separate measurements for each one against each
#       #sample in refClust
#       clusterPairSimil <- matrix(data=NA,ncol=2,nrow=nrow(compareMatrix))
#       
#       for(s in 1:ncol(compareMatrix)){
#         
#         for(c in 1:ncol(centroidMatrix)){
#         #COME BACK: this is tricky: have to do for each comparison, not just 1 and 2 anymore.  
#         clusterPairSimil[s, c] <- dcor(compareMatrix[,s],centroidMatrix[,1],index=1.0),
#                                        dcor(compareMatrix[,s],centroidMatrix[,2],index=1.0))
#         
#         }
#         
#       }
#       
#       
#      fract <- length(which(clusterPairSimil[,2]>=clusterPairSimil[,1]))/nrow(clusterPairSimil)
# 
#   return(fract)
# 
# }

########
#assumes matching features.
computeClusterPairAssignFract_matrixStatsDist <- function(compareMatrix,centroidMatrix,
                                    edgeMethod="Euclidean"){
  
  Class <- rep(NA, ncol(compareMatrix))
      #assumes is distance (re-run this function in loops,
  #adding an if statement at this level would slow it down, so check in wrapper functions.)       
        sampleCentroidDist <-   dist(t(compareMatrix),t(centroidMatrix),method=edgeMethod,by_rows=TRUE)
            
        for (i in 1:nrow(sampleCentroidDist)) {
          #want minimum distance.
          #which centroid is best?
          Class[i] <- which.min(sampleCentroidDist[i,])
         
        }
  
        fract <- table(Class[i])/ncol(compareMatrix)
        #if zero samples assigned to a cluster: won't be in here.
        bestMatch <- names(fract)[which.max(fract)]
        bestFract <- fract[which.max(fract)]
  
        results <- list(fract=fract,bestMatch=bestMatch,bestFract=bestFract)
    
  return(results)
  
}

computeClusterPairAssignFract_matrixStatsSimil <- function(compareMatrix,centroidMatrix,
                                                          edgeMethod="Euclidean"){
      Class <- rep(NA, ncol(compareMatrix))
    sampleCentroidSimil <- simil(t(compareMatrix),t(centroidMatrix),method=edgeMethod,by_rows=TRUE)

      for (i in 1:nrow(sampleCentroidSimil)) {
          #want minimum distance.
          Class[i] <- which.max(sampleCentroidSimil[i,])
         
        }
  
         fract <- table(Class[i])/ncol(compareMatrix)
        #if zero samples assigned to a cluster: won't be in here.
        bestMatch <- names(fract)[which.max(fract)]
        bestFract <- fract[which.max(fract)]
  
        results <- list(fract=fract,bestMatch=bestMatch,bestFract=bestFract)
  
}

##############
computeClusterPairAssignFract_cor <- function(compareMatrix,centroidMatrix,
                                    edgeMethod="pearson"){
  
  #features <- intersect(rownames(compareMatrix),rownames(centroidMatrix))
  #refClust <- refClust[rownames(refClust) %in% features, , drop=FALSE]
  #compareClust <- compareClust[rownames(compareClust) %in% features, , drop=FALSE]


    Class <- rep(NA, ncol(compareMatrix))
    sampleCentroidSimil <- cor(compareMatrix,centroidMatrix,method=edgeMethod)

      for (i in 1:nrow(sampleCentroidSimil)) {
          #want minimum distance.
          Class[i] <- which.max(sampleCentroidSimil[i,])
         
        }
        fract <- table(Class[i])/ncol(compareMatrix)
        #if zero samples assigned to a cluster: won't be in here.
        bestMatch <- names(fract)[which.max(fract)]
        bestFract <- fract[which.max(fract)]
  
        results <- list(fract=fract,bestMatch=bestMatch,bestFract=bestFract)
    
  return(results)

}


permuteCol <- function(x) {
  dd <- dim(x)
  n <- dd[2]
  p <- dd[1]
  mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
  matrix(x[order(mm)], p, n, byrow = FALSE)
  
}

createNullCentroidMatrixList <- function(centroidMatrix,numIter=100){

  SVD <- svd(centroidMatrix)
  centroidMatrixPrime <- centroidMatrix%*%SVD$v
  nullCentroidMatrixList <- list()
  
  for (i in 1:numIter) {
    
    tmpCentroids <- permuteCol(centroidMatrixPrime)
    nullCentroidMatrixList[[i]] <-  tmpCentroids%*%t(SVD$v)
    rownames(nullCentroidMatrixList[[i]]) <- rownames(centroidMatrix)
    
}
 return(nullCentroidMatrixLst)
}

createNulDataMatrixList <- function(dataMatrixList,numIter=100){

  nullDataMatrixList <- list()
  
  for(d in 1:length(dataMatrixList)){
  
    SVD <- svd(dataMatrixList[[d]])
    nullMatrixPrime <- dataMatrixList[[d]]%*%SVD$v
   
    #just permute once - if want to run again, just re-run this whole function.
    tmp <- permuteCol(nullMatrixPrime )
    nullDataMatrixList[[d]] <-  tmp%*%t(SVD$v)
    rownames(nullDataMatrixList[[d]]) <- rownames(dataMatrixList[[d]])
    colnames(nullDataMatrixList[[d]]) <- colnames(dataMatrixList[[d]])
    
  }

 return(nullDataMatrixList)

}

