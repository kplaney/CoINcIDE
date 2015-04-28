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
checkNA=FALSE){
  
    
  if(edgeMethod=="distCor"&& sigMethod=="centroid"){
    
    
    stop("Due to its computational complexity and run time and the actual theory behind distance correlation,
         the joint options of edgeMethod=distCor and sigMethod=centroid are now allowed. Please pick a different combination.")
    
  }
  
  date <- Sys.time();
  inputVariablesDF <- data.frame(date,edgeMethod,numParallelCores,minTrueSimilThresh,maxTrueSimilThresh,sigMethod,maxNullFractSize,numSims,
                                 includeRefClustInNull,outputFile,fractFeatIntersectThresh,numFeatIntersectThresh,clustSizeThresh, clustSizeFractThresh,
                                 checkNA);
  
  #capture.output prints the data correctly.
  capture.output(paste0("\nRunning find similar clusters on ",Sys.time()," with the following inputs:\n"),append=TRUE,file=outputFile);
  capture.output(inputVariablesDF,append=TRUE,file=outputFile);
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
message("nComputing cluster-cluster true similarities (or distances).\n")
#NOTE: we are indeed taking the mean of an entire symmetric matrix and taking the diagonal term.
#(same size similarityMatrix, #diagonal entries regardless of whether it's a null cluster or the true cluster.)
#but this will not affect our final p-value computations, and subsetting to get the upper triangle
#and remove the diagonal will slow down computations
  trueSimilData <- computeTrueSimil(clustIndexMatrix=clustIndexMatrix,edgeMethod=edgeMethod,
                                    dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,
                                  clustFeatureIndexList=clustFeatureIndexList,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                  numFeatIntersectThresh=numFeatIntersectThresh ,clustSizeThresh=clustSizeThresh, clustSizeFractThresh=clustSizeFractThresh)
  


  message(paste0(length(which(trueSimilData$similValueMatrix>=minTrueSimilThresh))/2," edges above minimum similarity threshold of ",minTrueSimilThresh))
  cat("\n",paste0(length(which(trueSimilData$similValueMatrix>=minTrueSimilThresh))/2," edges above minimum similarity threshold of ",minTrueSimilThresh),append=TRUE,file=outputFile)
  
cat("\nComputing similarity matrices for null/permutation calculations\n",append=TRUE,file=outputFile)
message("Computing similarity matrices for null/permutation calculations")
clustSimilMatrixList <- computeSimilMatrices(clustIndexMatrix=clustIndexMatrix,dataMatrixList=dataMatrixList,
                                             edgeMethod=edgeMethod,clustFeatureIndexList=clustFeatureIndexList.
                                             trueSimilMatrix=trueSimilData$similValueMatrix,
                                             minTrueSimilThresh= minTrueSimilThresh,
                                             maxTrueSimilThresh=maxTrueSimilThresh)

pvalueMatrix <- matrix(data=NA,nrow=numClust,ncol=numClust)

if(sigMethod=="centroid"){
  
  pvalueMatrix3 <- pvalueMatrix
  pvalueMatrix2 <- pvalueMatrix
  
}
cat("\nComputing p-values for each cluster-cluster similarity using null cluster distributions.\n",append=TRUE,file=outputFile)

  if(edgeMethod=="centroid"){
    
    if(!is.na(summary(pr_DB)[2]$distance[edgeMethod])){
      #these functions are written later in the script
      centroidFunction <-   computeClusterPairSimilMatrixStats
      
    }else if(edgeMethod=="distCor"){
      
      centroidFunction <- computeClusterPairSimil_dcor
      
    }else{
      
     centroidFunction <-  computeClusterPairSimil_cor
    }
      
  }
  for(n in 1:nrow(clustIndexMatrix)){
    #cluster N is our reference cluster.
    if(edgeMethod=="meanMatrix"){
     

  nullSimilMatrix <- computeNullSimilVector_mean(refClustRowIndex=n,dataMatrixList=dataMatrixList,clustSimilMatrixList=clustSimilMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,numSims=numSims,
                               trueSimilMatrix=trueSimilData$similValueMatrix,numParallelCores=numParallelCores,edgeMethod=edgeMethod,includeRefClustInNull=includeRefClustInNull,
                          clustIndexMatrix=clustIndexMatrix,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh)
  
    }else if(edgeMethod=="centroid"){

      nullSimilMatrix <- computeNullSimilVector_centroid(refClustRowIndex=n,clustSimilMatrixList=clustSimilMatrixList,
                                                         dataMatrixList=dataMatrixList, clustSampleIndexList= clustSampleIndexList,
                                                         clustFeatureIndexList= clustFeatureIndexList,numSims=numSims,
                                                    trueSimilMatrix=trueSimilData$similValueMatrix,numParallelCores=numParallelCores,edgeMethod,includeRefClustInNull=includeRefClustInNull,
                                                    clustIndexMatrix,minTrueSimilThresh=minTrueSimilThresh,
                                                    maxTrueSimilThresh=maxTrueSimilThresh,centroidFunction=centroidFunction)
      
    }
  
 #this is a symmetric matrix (with diagonal as NA)...so only count half of it.
  
  #cluster R is our comparison cluster.
 #COME BACK: test 3 c(.3,.2,.1) as maxNullFractSize
   for(r in 1:nrow(clustIndexMatrix)){

     if(sigMethod=="meanMatrix"){
       #CHECK: is trueSimilData$similValueMatrix computing the correct numbers??
       thresh <- trueSimilData$similValueMatrix[n,r]
       
       if(!is.na(summary(pr_DB)[2]$distance[edgeMethod]) && summary(pr_DB)[2]$distance[edgeMethod]){
         
         threshDir <- "less"
         
       }else{
         
         threshDir <- "greater"
       }
       
     }else if(sigMethod=="centroid"){
       

       threshDir <- "greater"
       
           #don't analyze clusters from the same dataset.
     if(clustIndexMatrix[n,2] != clustIndexMatrix[r,2] && !is.na(thresh)){
       
      if(!any(is.na(nullSimilMatrix[r,]))){
        
        pvalueMatrix[n,r] <- computeEdgePvalue(thresh=.1,threshDir=threshDir,nullSimilVector = nullSimilMatrix[r,])
        pvalueMatrix2[n,r] <- computeEdgePvalue(thresh=.2,threshDir=threshDir,nullSimilVector = nullSimilMatrix[r,])
        pvalueMatrix3[n,r] <- computeEdgePvalue(thresh=.3,threshDir=threshDir,nullSimilVector = nullSimilMatrix[r,])
      
      
      }
     }
       
     }else{
     #don't analyze clusters from the same dataset.
     if(clustIndexMatrix[n,2] != clustIndexMatrix[r,2] && !is.na(thresh)){
       
      if(!any(is.na(nullSimilMatrix[r,]))){
        
        pvalueMatrix[n,r] <- computeEdgePvalue(thresh=thresh,threshDir=threshDir,nullSimilVector = nullSimilMatrix[r,])
      
      
      }
     }
     
}
   }
 

  }

  if(sigMethod=="centroid"){
    
     output <- list(computeTrueSimilOutput=trueSimilData,pvalueMatrix1=pvalueMatrix,
                    pvalueMatrix2=pvalueMatrix2,pvalueMatrix3=pvalueMatrix3,
                    clustIndexMatrix=clustIndexMatrix,inputVariablesDF=inputVariablesDF)
    
 }else{
   
   output <- list(computeTrueSimilOutput=trueSimilData,pvalueMatrix=pvalueMatrix,clustIndexMatrix=clustIndexMatrix,inputVariablesDF=inputVariablesDF)
   
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


####
computeSimilMatrices <- function(clustIndexMatrix,dataMatrixList,edgeMethod,clustFeatureIndexList,
                                 trueSimilMatrix,minTrueSimilThresh,maxTrueSimilThresh){
  #NOTE: assumes features are the same across all clusters in a dataset.
  #if want to use different features within one dataset: can't optimize the code as much
  #could just have a fake "second" dataset for clusters with the different rows.
  clustSimilMatrixList <- list()

  if(edgeMethod!="distCor"){
    
    for(d in 1:length(dataMatrixList)){
      #what clusters are in this dataset?
      row_indices <- as.numeric(clustIndexMatrix[which(clustIndexMatrix[,2]==d) ,1])
      
      for(c in 1:length(dataMatrixList)){
        
        #has this been calculated yet? (just in c_d, not d_c direction)
        if(is.null(clustSimilMatrixList[[paste0(c,"_",d)]])){
          
          
          col_indices <- as.numeric(clustIndexMatrix[which(clustIndexMatrix[,2]==c) ,1])
          
          if(any(!is.na(trueSimilMatrix[row_indices,col_indices]))){
          #diagonal will be NA/not counted here.
          if(any(trueSimilMatrix[row_indices,col_indices]>= minTrueSimilThresh) && any(trueSimilMatrix[row_indices,col_indices]<=maxTrueSimilThresh)){


            #NOTE: assumes features are the same across all clusters in a dataset (why: once get similarity correlation
            #or dist matrices, can't subset further on gene indices.)
            #if want to use different features within one dataset: can't optimize the code as much
            #could just have a fake "second" dataset for clusters with the different rows. 
              features1 <- rownames(dataMatrixList[[d]][clustFeatureIndexList[[d]][[1]], ])
              
              features2 <- rownames(dataMatrixList[[c]][clustFeatureIndexList[[c]][[1]], ])

            featuresIntersect <- intersect(features2,features2)
            
            if(!is.na(summary(pr_DB)[2]$distance[edgeMethod])){
              #distance or cor matrix?
              if(summary(pr_DB)[2]$distance[edgeMethod]){
                
                clustSimilMatrixList[[paste0(d,"_",c)]] <-   as.matrix(dist(t(dataMatrixList[[d]][featuresIntersect,]),
                                                                  t(dataMatrixList[[c]][featuresIntersect,]),method=edgeMethod,by_rows=TRUE))

                
              }else{
                
                clustSimilMatrixList[[paste0(d,"_",c)]] <- simil(t(dataMatrixList[[d]][featuresIntersect,]),
                                                                 t(dataMatrixList[[c]][featuresIntersect,]),method=edgeMethod,by_rows=TRUE)
                
              }
              
            }else{
              
              clustSimilMatrixList[[paste0(d,"_",c)]] <- cor(dataMatrixList[[d]][featuresIntersect,],
                                                             dataMatrixList[[c]][featuresIntersect,],method=edgeMethod)
              
            }
            
            
          }
          
        }
     
      }
      
      }
      
    }
    
    }else{
      
      #just compute dist for each cluster.
      for(d in 1:length(dataMatrixList)){
        
        row_indices <- clustIndexMatrix[which(clustIndexMatrix[,2]==d) ,1]
        
        for(r in 1:length(row_indices)){
        #diagonal will be NA/not counted here (i.e. row_indices,row_indices)
          
        if(any(!is.na(trueSimilMatrix[row_indices, ]))){
            
        if(any(trueSimilMatrix[row_indices[r],]>= minTrueSimilThresh) && any(trueSimilMatrix[row_indices[r],]<=maxTrueSimilThresh)){

          #want gene by gene: a summary of the genetic pattern, across all samples.
          clustSimilMatrixList[[row_indices]] <- as.matrix(dist(dataMatrixList[[d]][features1,],
                                                      method="euclidean"))
          
          
        }
        
      }
      
        }
    }
  
  return(clustSimilMatrixList)
  
}

computeNullSimilVector_mean <-   function(refClustRowIndex,clustSimilMatrixList,dataMatrixList,clustSampleIndexList,clustFeatureIndexList,numSims=100,
                               trueSimilMatrix,numParallelCores=1,edgeMethod,includeRefClustInNull=TRUE,clustIndexMatrix,minTrueSimilThresh=-Inf,maxTrueSimilThresh=Inf){  
    
    registerDoParallel(cores=numParallelCores)

    
    refStudyNum <- as.numeric(clustIndexMatrix[refClustRowIndex,2])
    refDataMatrix <- dataMatrixList[[as.numeric(clustIndexMatrix[refClustRowIndex,2])]]
    sampleIndices <- clustSampleIndexList[[as.numeric(clustIndexMatrix[refClustRowIndex,2])]][[as.numeric(clustIndexMatrix[refClustRowIndex,3])]]
    featureIndices1 <- clustFeatureIndexList[[as.numeric(clustIndexMatrix[refClustRowIndex,2])]][[as.numeric(clustIndexMatrix[refClustRowIndex,3])]]
    refClust <- refDataMatrix[featureIndices1,sampleIndices]

  
  if(edgeMethod!="distCor"){
     
        #will return a matrix of numClusters x numFeatures
        #see documention on CRAN "nesting foreach loops"
        #"nesting" way of dopar with %:% didn't work too well. foreach can be weird with iterative random numbers, so play it safe.
        nullSimilVector <- foreach(i=1:numSims,.combine='cbind') %dopar% {
          
          
          #want same set of nullClusters across ALL comparison clusters.
          #so hold each null cluster "constant", loop through compare clusters.
          if(includeRefClustInNull){
          
            nullClustIndices <-sample.int(ncol(refDataMatrix),size=ncol(refClust),replace=TRUE)
            
          }else{
            
            #note: if only a few samples NOT in the ref cluster, this will mostly be duplicated samples
            colIndicesNoRefClust <- na.omit(match(setdiff(colnames(refDataMatrix),colnames(refClust)),colnames(refDataMatrix)))
            nullClustIndices <- sample(colIndicesNoRefClust,size=ncol(refClust),replace=TRUE)
            
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
              
              compareStudyNum <- as.numeric(clustIndexMatrix[c,2])
              compareSampleIndices <- clustSampleIndexList[[as.numeric(clustIndexMatrix[c,2])]][[as.numeric(clustIndexMatrix[c,3])]]
                
  
            if(!is.null( clustSimilMatrixList[[paste0(refStudyNum,"_",compareStudyNum )]])){
                
                #NOTE: we are indeed taking the mean of an entire symmetric matrix and taking the diagonal term.
                #(same size similarityMatrix, #diagonal entries regardless of whether it's a null cluster or the true cluster.)
                #but this will not affect our final p-value computations, and subsetting to get the upper triangle
                #and remove the diagonal will slow down computations
                similValue <- mean(clustSimilMatrixList[[paste0(refStudyNum,"_",compareStudyNum )]][nullClustIndices,compareSampleIndices])
                
              }else if(!is.null( clustSimilMatrixList[[paste0(compareStudyNum,"_",refStudyNum)]])){
                
                similValue <- mean(clustSimilMatrixList[[paste0(paste0(compareStudyNum,"_",refStudyNum))]][compareSampleIndices,nullClustIndices])
                
              }else{
                
                stop(paste0("\n\n finding an non-null index for studies ",compareStudyNum, " and ",refStudyNum))
             
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
#end of if not distCor (don't want this if statement embedded in doParallel loop - might slow things down.)
}else if(edgeMethod=="distCor"){
  #will return a matrix of numClusters x numFeatures
  #see documention on CRAN "nesting foreach loops"
  #"nesting" way of dopar with %:% didn't work too well. foreach can be weird with iterative random numbers, so play it safe.
  nullSimilVector <- foreach(i=1:numSims,.combine='cbind') %dopar% {
    
    
    #want same set of nullClusters across ALL comparison clusters.
    #so hold each null cluster "constant", loop through compare clusters.
    if(includeRefClustInNull){
      
      nullClustIndices <-sample.int(ncol(refDataMatrix),size=ncol(refClust),replace=TRUE)
      
    }else{
      
      #note: if only a few samples NOT in the ref cluster, this will mostly be duplicated samples
      colIndicesNoRefClust <- na.omit(match(setdiff(colnames(refDataMatrix),colnames(refClust)),colnames(refDataMatrix)))
      nullClustIndices <- sample(colIndicesNoRefClust,size=ncol(refClust),replace=TRUE)
      
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
          
          compareStudyNum <- as.numeric(clustIndexMatrix[c,2])
          featureIndices2 <- clustFeatureIndexList[[as.numeric(clustIndexMatrix[c,2])]][[as.numeric(clustIndexMatrix[c,3])]]
          featureIndices <- intersect(featureIndices1,featureIndices2)

          
          if(!is.null(clustSimilMatrixList[[c]])){
            #need to make clustSimilMatrixList a dist object again.
          similValue <- dcor(dist(refDataMatrix[featureIndices,nullClustIndices],method="euclidean"),
                             as.dist(clustSimilMatrixList[[c]][featureIndices, featureIndices]),index=1.0)
          
          stop(paste0("\n\n finding an non-null index for global cluster ",c))
               
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
}

}#EOF


computeNullSimilVector_centroid <-   function(refClustRowIndex,clustSimilMatrixList,dataMatrixList,
                                              clustSampleIndexList,clustFeatureIndexList,numSims=100,
                                          trueSimilMatrix,numParallelCores=1,edgeMethod,includeRefClustInNull=TRUE,
                                          clustIndexMatrix,minTrueSimilThresh=-Inf,maxTrueSimilThresh=Inf,
                                          centroidFunction){  
  
  registerDoParallel(cores=numParallelCores)
  
  
  refStudyNum <- as.numeric(clustIndexMatrix[refClustRowIndex,2])
  refDataMatrix <- dataMatrixList[[as.numeric(clustIndexMatrix[refClustRowIndex,2])]]
  refSampleIndices <- clustSampleIndexList[[as.numeric(clustIndexMatrix[refClustRowIndex,2])]][[as.numeric(clustIndexMatrix[refClustRowIndex,3])]]
  featureIndices1 <- clustFeatureIndexList[[as.numeric(clustIndexMatrix[refClustRowIndex,2])]][[as.numeric(clustIndexMatrix[refClustRowIndex,3])]]
  refClust <- refDataMatrix[featureIndices1,  refSampleIndices]
  
    
    #will return a matrix of numClusters x numFeatures
    #see documention on CRAN "nesting foreach loops"
    #"nesting" way of dopar with %:% didn't work too well. foreach can be weird with iterative random numbers, so play it safe.
    nullSimilVector <- foreach(i=1:numSims,.combine='cbind') %dopar% {
      
      
      #want same set of nullClusters across ALL comparison clusters.
      #so hold each null cluster "constant", loop through compare clusters.
      if(includeRefClustInNull){
        
        nullClustIndices <-sample.int(ncol(refDataMatrix),size=ncol(refClust),replace=TRUE)
        
      }else{
        
        #note: if only a few samples NOT in the ref cluster, this will mostly be duplicated samples
        colIndicesNoRefClust <- na.omit(match(setdiff(colnames(refDataMatrix),colnames(refClust)),colnames(refDataMatrix)))
        nullClustIndices <- sample(colIndicesNoRefClust,size=ncol(refClust),replace=TRUE)
        
      }
      
      nullClustMatrix <- refDataMatrix[featureIndices1,nullClustIndices]
      
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
            
            compareStudyNum <- as.numeric(clustIndexMatrix[c,2])
            #only take intersecting features
            features <- intersect(features1,clustFeatureIndexList[[as.numeric(clustIndexMatrix[c,2])]][[as.numeric(clustIndexMatrix[c,3])]])
            compareSampleIndices <- clustSampleIndexList[[as.numeric(clustIndexMatrix[c,2])]][[as.numeric(clustIndexMatrix[c,3])]]
            #COME BACK:
            compareClust <- dataMatrixList[[compareStudyNum]][features,  compareSampleIndices]
            centroidMatrix <- cbind(rowMeans(refClust),rowMeans(nullClust))[features, ]
           
            fractNullSimil <- centroidFunction(compareClust=compareClust,centroidMatrix=centroidMatrix,
                                               edgeMethod=edgeMethod)
            
            
          }else{
            
            fractNullSimil  <- NA
            
          }  
          
        }else{
          
          fractNullSimil  <- NA
          
        }
        #end of inner foreach dopar 
      }
      #end of outer foreach dopar
    }

  
}#EOF
# 
# computeNullSimilVector_centroid <-   function(refClustRowIndex,clustSimilMatrixList,dataMatrixList,clustSampleIndexList,clustFeatureIndexList,numSims=100,
#                                           trueSimilMatrix,numParallelCores=1,edgeMethod,includeRefClustInNull=TRUE,clustIndexMatrix,minTrueSimilThresh=-Inf,maxTrueSimilThresh=Inf){  
#   
#   registerDoParallel(cores=numParallelCores)
#   
#   
#   refStudyNum <- as.numeric(clustIndexMatrix[refClustRowIndex,2])
#   refDataMatrix <- dataMatrixList[[as.numeric(clustIndexMatrix[refClustRowIndex,2])]]
#   refSampleIndices <- clustSampleIndexList[[as.numeric(clustIndexMatrix[refClustRowIndex,2])]][[as.numeric(clustIndexMatrix[refClustRowIndex,3])]]
#   featureIndices1 <- clustFeatureIndexList[[as.numeric(clustIndexMatrix[refClustRowIndex,2])]][[as.numeric(clustIndexMatrix[refClustRowIndex,3])]]
#   refClust <- refDataMatrix[featureIndices1,  refSampleIndices]
#   
#   
#   if(edgeMethod=="distCor"){
#     
#     
#     stop("Due to its computational complexity and run time and the actual theory behind distance correlation,
#          the joint options of edgeMethod=distCor and sigMethod=centroid are now allowed. Please pick a different combination.")
#     
#   }
#     
#     #will return a matrix of numClusters x numFeatures
#     #see documention on CRAN "nesting foreach loops"
#     #"nesting" way of dopar with %:% didn't work too well. foreach can be weird with iterative random numbers, so play it safe.
#     nullSimilVector <- foreach(i=1:numSims,.combine='cbind') %dopar% {
#       
#       
#       #want same set of nullClusters across ALL comparison clusters.
#       #so hold each null cluster "constant", loop through compare clusters.
#       if(includeRefClustInNull){
#         
#         nullClustIndices <-sample.int(ncol(refDataMatrix),size=ncol(refClust),replace=TRUE)
#         
#       }else{
#         
#         #note: if only a few samples NOT in the ref cluster, this will mostly be duplicated samples
#         colIndicesNoRefClust <- na.omit(match(setdiff(colnames(refDataMatrix),colnames(refClust)),colnames(refDataMatrix)))
#         nullClustIndices <- sample(colIndicesNoRefClust,size=ncol(refClust),replace=TRUE)
#         
#       }
#       
#       #compare simil against each cluster
#       foreach(c=1:nrow(clustIndexMatrix), .combine='c') %dopar%{
#         # for(c in 1:nrow(clustIndexMatrix))  {
#         #debug: is doParallel keeping the same nullclust in this loop, and is each nullClust truly random?
#         #             if(i==1){
#         #               #these should all be equal
#         #               cat("\n",mean(nullClust),file="~/Desktop/log.txt",append=TRUE)
#         # 
#         #             }
#         
#         if(!is.na(trueSimilMatrix[refClustRowIndex,c])){ 
#           
#           if(trueSimilMatrix[refClustRowIndex,c] >= minTrueSimilThresh && trueSimilMatrix[refClustRowIndex,c] <= maxTrueSimilThresh){
#             
#             compareStudyNum <- as.numeric(clustIndexMatrix[c,2])
#             compareSampleIndices <- clustSampleIndexList[[as.numeric(clustIndexMatrix[c,2])]][[as.numeric(clustIndexMatrix[c,3])]]
#             if(!is.null( clustSimilMatrixList[[paste0(refStudyNum,"_",compareStudyNum )]])){
#               
#               #NOTE: we are indeed taking the mean of an entire symmetric matrix and taking the diagonal term.
#               #(same size similarityMatrix, #diagonal entries regardless of whether it's a null cluster or the true cluster.)
#               #but this will not affect our final p-value computations, and subsetting to get the upper triangle
#               #and remove the diagonal will slow down computations
#               similTrue <- colMeans(clustSimilMatrixList[[paste0(refStudyNum,"_",compareStudyNum )]][            
#                 refSampleIndices,   compareSampleIndices])
#               similNull <- colMeans(clustSimilMatrixList[[paste0(refStudyNum,"_",compareStudyNum )]][nullClustIndices, compareSampleIndices])
#   
#               fractNullSimil <- length(which(similNull>=similTrue))/length(similTrue)
#               
#             }else if(!is.null( clustSimilMatrixList[[paste0(compareStudyNum,"_",refStudyNum)]])){
#               
#               similTrue <- rowMeans(clustSimilMatrixList[[paste0(paste0(compareStudyNum,"_",refStudyNum))]][compareSampleIndices,refSampleIndices])
#               similNull <- rowMeans(clustSimilMatrixList[[paste0(paste0(compareStudyNum,"_",refStudyNum))]][compareSampleIndices,nullClustIndices])
#           
#               fractNullSimil  <- length(which(similNull>=similTrue))/length(similTrue)
#               
#                 }else{
#               
#               stop(paste0("\n\n finding an non-null index for studies ",compareStudyNum, " and ",refStudyNum))
#                    
#             }
#             
#             
#             
#           }else{
#             
#             fractNullSimil  <- NA
#             
#           }  
#           
#         }else{
#           
#           fractNullSimil  <- NA
#           
#         }
#         #end of inner foreach dopar 
#       }
#       #end of outer foreach dopar
#     }
# 
#   
# }#EOF

computeClusterPairSimil_mean <- function(refClust,compareClust,
                                         edgeMethod="correlation"){
  
  
  similMatrix <- computeClusterPairSimil(refClust,compareClust,
                                         edgeMethod=edgeMethod)
  
  
  return(mean(similMatrix))
  
  
  
}

#come back: could this be more efficient/run similarly to nullClust sim code?
computeClusterPairSimil_dcor <- function(compareClust,centroidMatrix,edgeMethod="distCor"){
  
      # centroidMatrix is two columns of centroids: want separate measurements for each one against each
      #sample in refClust
      clusterPairSimil <- matrix(data=NA,ncol=2,nrow=ncol(distCompareClust))
      
      for(s in 1:ncol(refClust)){
        
        clusterPairSimil[s, ] <- cbind(dcor(compareClust[,s],centroidMatrix[,1],index=1.0),
                                       dcor(compareClust[,s],centroidMatrix[,2],index=1.0))
        
      }
      
                   
     fract <- length(which(clusterPairSimil[,2]>=clusterPairSimil[,1]))/nrow(clusterPairSimil)
    
  }

  return(fract)

}

########
#come back: could this be more efficient/run similarly to nullClust sim code?
#assumes matching features.
computeClusterPairSimil_matrixStats <- function(compareMatrix,centroidMatrix,
                                    edgeMethod="Euclidean"){
  
  #features <- intersect(rownames(compareMatrix),rownames(centroidMatrix))
  #refClust <- refClust[rownames(refClust) %in% features, , drop=FALSE]
  #compareClust <- compareClust[rownames(compareClust) %in% features, , drop=FALSE]


      #distance or cor matrix?
      if(summary(pr_DB)[2]$distance[edgeMethod]){
        
        clusterPairSimil <-   dist(t(compareMatrix),t(centroidMatrix),method=edgeMethod,by_rows=TRUE)
              
     
        fract <- length(which(clusterPairSimil[,2]<=clusterPairSimil[,1]))/nrow(clusterPairSimil)
        
      }else{
        
        clusterPairSimil <- simil(t(compareMatrix),t(centroidMatrix),method=edgeMethod,by_rows=TRUE)
              
     fract <- length(which(clusterPairSimil[,2]>=clusterPairSimil[,1]))/nrow(clusterPairSimil)
        
      }


    
  return(fract)
  
}

##############
computeClusterPairSimil_cor <- function(compareMatrix,centroidMatrix,
                                    edgeMethod="pearson"){
  
  #features <- intersect(rownames(compareMatrix),rownames(centroidMatrix))
  #refClust <- refClust[rownames(refClust) %in% features, , drop=FALSE]
  #compareClust <- compareClust[rownames(compareClust) %in% features, , drop=FALSE]


        
    clusterPairSimil <- cor(compareMatrix,centroidMatrix,method=edgeMethod,by_rows=TRUE)
              
    fract <- length(which(clusterPairSimil[,2]>=clusterPairSimil[,1]))/nrow(clusterPairSimil)
    
  return(fract)

}
########

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
