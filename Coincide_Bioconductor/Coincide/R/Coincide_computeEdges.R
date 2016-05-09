#library("proxy")
##library("doParallel")
#library("foreach")
#suggests:
#library("energy")
#library("limma")
#library("matrixStats")
#Rclusterpp does not return the distance matrix.
#library("Rclusterpp")

#TO DO: add more filters later so can speed up analyses when clusters don't have any features that intersect, etc.?
#i.e. user has already run FDR function, decided on some of these features.
computeAdjMatrices <- function(dataMatrixList,clustSampleIndexList,clustFeatureIndexList,
edgeMethod=c("spearman","pearson","kendall","Euclidean","cosine",
                      "Manhattan","Minkowski","Mahalanobis"),
numSims=500,outputFile="./CoINcIDE_messages.txt",checkNA=FALSE,centroidMethod=c("mean","median"),seedNum=as.numeric(Sys.time()),
memorySave=FALSE){
  
  
    if(length(centroidMethod)>1){
      
      centroidMethod<-c("mean")
      
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


  
  if(length(clustSampleIndexList) != length(clustFeatureIndexList) || length(clustSampleIndexList) !=
       length(dataMatrixList) || length(clustFeatureIndexList) !=  length(dataMatrixList)){
    
    stop("The lenghts of your clustSampleIndexList,clustFeatureIndexList and dataMatrixList are not equal")
 
  }
  
  ##check: all clustIndex lists have featureIndexLists, and these
  #exist in the fullDataMatrix?
  date <- Sys.time();
  inputVariablesDF <- data.frame(date,edgeMethod,
  numSims,centroidMethod);
  
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



  
  pvalueMatrix <- matrix(data=NA,nrow=numClust,ncol=numClust)
  meanMetricMatrix <-  matrix(data=NA,nrow=numClust,ncol=numClust)
  
 
  cat("\nComputing p-values for each cluster-cluster similarity using null cluster distributions.\n",append=TRUE,file=outputFile)

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
  numDatasets <- length(dataMatrixList)
  
  if(memorySave){
    
    set.seed(seedNum)
    randTag <- rnorm(1)
    save(dataMatrixList,file=paste0(saveDir,"/dataMatrixList.rds"),compress=TRUE)
    #remove the dataMatrixList - not needed anymore and will free up space/memory.
    rm(list="dataMatrixList")
  }
  for(n in 1:numDatasets){
    
    message("On dataset (loop) ",n)
    
  #for each dataset: make centroid set and null centroid sets
  centroidMatrixOrig <- matrix(data=NA,ncol=length(clustFeatureIndexList[[n]]),nrow=length(clustFeatureIndexList[[n]][[1]]))
  
  if(memorySave){
    
    dataMatrixList <- readRDS(file=paste0(saveDir,"/dataMatrixList.rds"))
    
  }
                           
  if(centroidMethod=="mean"){

        for(d in 1:length(clustFeatureIndexList[[n]])){
          
          centroidMatrixOrig[,d] <- rowMeans(dataMatrixList[[n]][ clustFeatureIndexList[[n]][[1]], clustSampleIndexList[[n]][[d]] ,drop=FALSE])
        
        }
        #otherwise: just median. 
      }else{

        
        for(d in 1:length(clustFeatureIndexList[[n]])){
          
          centroidMatrixOrig[,d] <- rowMedians(dataMatrixList[[n]][ clustFeatureIndexList[[n]][[1]], clustSampleIndexList[[n]][[d]] ,drop=FALSE])
        
        }
  
    }

  if(memorySave){
    
    rm(list="dataMatrixList")
    
  }
  rownames( centroidMatrixOrig) <- rownames(dataMatrixList[[n]])[clustFeatureIndexList[[n]][[1]]]
    nullCentroidListOrig <- createNullCentroidMatrixList(centroidMatrixOrig,numIter=numSims)
        
      #can foreach work here? perhaps if I don't combine.
      # pvalueMatrix[, c] <- foreach(i=1:numSims) %dopar% {
       for(c in 1:nrow(clustIndexMatrix)){

         if(as.numeric(clustIndexMatrix[c,2]) != n){
                    
          message("Running tests for cluster number: ",c," or cluster ", as.numeric(clustIndexMatrix[c,3])," from dataset ",as.numeric(clustIndexMatrix[c,2]))
#ADD thresholds here.
        #  if(threshStats$..)
         sampleIndices <- clustSampleIndexList[[as.numeric(clustIndexMatrix[c,2])]][[as.numeric(clustIndexMatrix[c,3])]]
        
         if(all(is.null(rownames(dataMatrixList[[as.numeric(clustIndexMatrix[c,2])]])))){
           
           stop("You have null feature names in at least one of your data matrix list objects.")
           
         }
         features <- intersect(rownames(dataMatrixList[[as.numeric(clustIndexMatrix[c,2])]])[ clustFeatureIndexList[[as.numeric(clustIndexMatrix[c,2])]][[1]] ],
                                    #just pick first feature index for now.
                              
                                    rownames(dataMatrixList[[n]])[ clustFeatureIndexList[[n]][[1]] ])
        centroidMatrix <- centroidMatrixOrig[features, ,drop=FALSE]
        
        nullCentroidList <- lapply(nullCentroidListOrig,FUN=function(nullCentroidSet,featureSet){
          
          result <- nullCentroidSet[features, , drop=FALSE]
          return(result)
          
        },featureSet=features)

        compareClust <- dataMatrixList[[as.numeric(clustIndexMatrix[c,2])]][features,sampleIndices,drop=FALSE]  
   
         #get results for true matrix
         fractNNresults <- centroidFunction(compareMatrix=compareClust,centroidMatrix=centroidMatrix,edgeMethod=edgeMethod)
         bestMatchIndex <- fractNNresults$bestMatch
         globalRefClustIndex <- intersect(which(clustIndexMatrix[,2]==n),which(clustIndexMatrix[,3]==bestMatchIndex))
         allRefDatasetClustIndices <- which(clustIndexMatrix[,2]==n)
         
        if(length(globalRefClustIndex) !=1){
           
           stop("Error in global cluster indexing.")
         
        }

         trueFractNNmatrix[globalRefClustIndex,c] <- fractNNresults$bestFract
         #store alll mean metrics to get a pseudo global distribution
         meanMetricMatrix[allRefDatasetClustIndices,c] <- fractNNresults$allMeanMetrics
         
         #unit test:
         if(meanMetricMatrix[globalRefClustIndex,c] != fractNNresults$meanMetric){
           
           stop("Not indexing/saving mean metrics correctly")
         }
         
        #NOW: run through null centroid list.
         
        #is it a distance metric?
        if(!is.na((summary(pr_DB)[2]$distance[edgeMethod])) && summary(pr_DB)[2]$distance[edgeMethod]){
          
         nullTests <- lapply(nullCentroidList,FUN=function(nullCentroidMatrix,compareClust,edgeMethod,thresh){
           
           passedThresh <- FALSE
           tmp <- centroidFunction(compareMatrix=compareClust,centroidMatrix=nullCentroidMatrix,edgeMethod=edgeMethod)
           
           #does null data meet NN fract thresh and meanMetric (dist) thresh? (ie is null distance smaller than true dist)
           if( (tmp$bestFract>=thresh) && (tmp$meanMetric <= fractNNresults$meanMetric)) {
           #if not looking at correlation: 
         #    if( (tmp$bestFract>=thresh) ) {
             
             passedThresh <- TRUE
           }    
           
               return(passedThresh)
           
         },compareClust=compareClust,edgeMethod=edgeMethod,thresh=trueFractNNmatrix[globalRefClustIndex,c])

        
        
        #similarity, not distance so second boolean is now >=
        }else{
          
           nullTests <- lapply(nullCentroidList,FUN=function(nullCentroidMatrix,compareClust,edgeMethod,thresh){
           
           passedThresh <- FALSE
           tmp <- centroidFunction(compareMatrix=compareClust,centroidMatrix=nullCentroidMatrix,edgeMethod=edgeMethod)
           
           #only compare NN here.  similarity isn't really a fair comparison with this type of data.
           #old code:
          if( (tmp$bestFract>=thresh) && (tmp$meanMetric >= fractNNresults$meanMetric)) {
          #if no similarity comparison: 
        #    if( (tmp$bestFract>=thresh) ) {
             
             passedThresh <- TRUE
             
           } 
           
           return(passedThresh)
           
        } ,compareClust=compareClust,edgeMethod=edgeMethod,thresh=trueFractNNmatrix[globalRefClustIndex,c])

     }
        pvalueMatrix[globalRefClustIndex,c] <- length(which(unlist(nullTests)))/length(nullCentroidList)
    
  #done computing metrics for this cluster.
  }#else{
    #add else so that foreach loop still counts this column? not needed for now.
  #  pvalueMatrix[globalClustIndex,c] <- NA
  #}

#end of looping through all datasets
         #rm(list="nullCentroidList")
}

    rm(list="nullCentroidListOrig")
   
 }

  if(memorySave){
    
    dataMatrixList <- readRDS(file=paste0(saveDir,"/dataMatrixList.rds"))
    
  }

  threshStats <- computeThreshParam(clustIndexMatrix,
                                    dataMatrixList,clustSampleIndexList,clustFeatureIndexList)
  
  #don't need this object anymore- remove to save space
  rm(list="dataMatrixList")
  
  threshStats$meanMetricMatrix <- meanMetricMatrix
  threshStats$trueFractNNmatrix <- trueFractNNmatrix
  
  #to help pick meanMatrix thresh
  meanMatrixQuantiles <- quantile(as.vector(threshStats$meanMetricMatrix),na.rm=TRUE,probs=seq(0,1,.05))
  
  output <- list(pvalueMatrix=pvalueMatrix,meanMetricMatrix=meanMetricMatrix,
                  clustIndexMatrix=clustIndexMatrix,inputVariablesDF=inputVariablesDF,
                  trueFractNNmatrix=trueFractNNmatrix,computeTrueSimilOutput=threshStats,meanMatrixQuantiles=meanMatrixQuantiles)
 return(output)
#EOF
}


######
computeThreshParam <- function(clustIndexMatrix,
                                    dataMatrixList,clustSampleIndexList,clustFeatureIndexList){
  
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
          
        }

      }
    
    }

  }


  output <- list(numFeatIntersectMatrix=numFeatIntersectMatrix,fractFeatIntersectMatrix=fractFeatIntersectMatrix,
                 clustSizeMatrix=clustSizeMatrix,clustSizeFractMatrix=clustSizeFractMatrix)
  
  return(output)
#EOF  
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
  
        fract <- table(Class)/ncol(compareMatrix)
        #if zero samples assigned to a cluster: won't be in here.
        bestMatch <- as.numeric(names(fract)[which.max(fract)])
        bestFract <- fract[which.max(fract)]
  
        #for ALL samples in compareMatrix:  take mean of distance with best centroid.
        meanMetric <- mean(sampleCentroidDist[,bestMatch])
        
        
        allMeanMetrics <- c()
        for(a in 1:ncol(centroidMatrix)){
          
          allMeanMetrics[a] <- mean(sampleCentroidDist[,a])
          
        }
  
        results <- list(fract=fract,bestMatch=bestMatch,bestFract=bestFract,meanMetric=allMeanMetrics[bestMatch],allMeanMetrics=allMeanMetrics)
    
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
  
         fract <- table(Class)/ncol(compareMatrix)
        #if zero samples assigned to a cluster: won't be in here.
         #bestMatch is the index for the best match.
        bestMatch <- as.numeric(names(fract)[which.max(fract)])
        bestFract <- fract[which.max(fract)]

        
        allMeanMetrics <- c()
        
        for(a in 1:ncol(centroidMatrix)){
          #collapse the similarity between this centroid (column a) and all samples (rows) into a single mean value.
          allMeanMetrics[a] <- mean(sampleCentroidSimil[,a])
          
        }
        
  
        results <- list(fract=fract,bestMatch=bestMatch,bestFract=bestFract,meanMetric=allMeanMetrics[bestMatch],allMeanMetrics=allMeanMetrics)
  
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
        fract <- table(Class)/ncol(compareMatrix)
        #if zero samples assigned to a cluster: won't be in here.
        bestMatch <- as.numeric(names(fract)[which.max(fract)])
        bestFract <- fract[which.max(fract)]
  

        allMeanMetrics <- c()
        for(a in 1:ncol(centroidMatrix)){
          
          allMeanMetrics[a] <- mean(sampleCentroidSimil[,a])
          
        }
  
        results <- list(fract=fract,bestMatch=bestMatch,bestFract=bestFract,meanMetric= allMeanMetrics[bestMatch], allMeanMetrics=allMeanMetrics)
    
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
    #needs to be a matrix
    nullCentroidMatrixList[[i]] <-  as.matrix(tmpCentroids%*%t(SVD$v))
    rownames(nullCentroidMatrixList[[i]]) <- rownames(centroidMatrix)
    colnames(nullCentroidMatrixList[[i]]) <- c(1:ncol(centroidMatrix))
    
    
}
 return(nullCentroidMatrixList)

}

createNullDataMatrixList <- function(dataMatrixList){

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

