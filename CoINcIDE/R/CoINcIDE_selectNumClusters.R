#careful: are your features/sample in correct dimension for yoru clustermembership function? (in rows or columns?)
#usually assume samples in rows for meta-clustering code.

#EDGE case: what if there's a really small outlier? like 2 samples and the rest belong to one cluster? 
#THen the number of clusters will be correctly deemed 1, but it would be nice to know that we should throw out those outliers...

CoINcIDE_selectK <- function(dataMatrix,clustFeatures,
                                            edgeMethod=c("distCor","spearman","pearson","kendall","Euclidean","cosine",
                                                         "Manhattan","Minkowski","Mahalanobis"),numParallelCores=1,minTrueSimilThresh=-Inf,maxTrueSimilThresh=Inf,
                                            sigMethod=c("meanMatrix","centroid"),maxNullFractSize=.1,numSims=100,includeRefClustInNull=TRUE,
                                            
                                            outputFile="./CoINcIDE_messages.txt",fractFeatIntersectThresh=0,numFeatIntersectThresh=0 ,clustSizeThresh=0, clustSizeFractThresh=0,
  distMethod=c("pearson","spearman","euclidean", "binary", "maximum", "canberra", "minkowski"),
  hclustAlgorithm=c("average","complete","ward.D", "ward.D2", "single", "mcquitty","median","centroid"), 
  clustMethod=c("km","hc"), saveDir="/home/kplaney/ovarian_analysis/", indEdgePvalueThresh=.1,corUse="everything"
  meanEdgePairPvalueThresh=.05, commMethod=c("edgeBetween"),bestKSelect=c("max","min"),numIter=20,maxNumClusters=15,
  iter.max=30,nstart=10,splitDataBeforeClust=FALSE
  
  ){
  
  dataset <- dataMatrix[rownames(dataMatrix) %in% clustFeatures, , drop=FALSE]
  
  if(nrow(dataset)==0){
    
    stop("\nIn clusterMatrixKmeansGap function: no clustFeatures were found in the data matrix inputted.")
    
  }
  
  warning("Assumes samples are in the columns.")
  #must have unique row names!
  if(any(duplicated(colnames(dataMatrix)))){   
    colnames(dataMatrix) <- paste0(colnames(dataMatrix),"_",c(1:ncol(dataMatrix)));
  }
  
  if(ncol(dataset)<maxNumClusters){
    #hclust usually returns NA gap test if kMax = ncol(dataset)
    #will also mest up maxSE calculations
    kMax <- ncol(dataset)-1
    
  }else{
    
    kMax <- maxNumClusters
    
  }
  
  if(clustMethod=="hc"){
    
  if(distMethod==("pearson") || distMethod=="spearman"){
    
    datasetClust <- dataset
    
    clustF <- function(x,k){
      #it turns out that R will recognize the upper-level function input value in this function (algorith, corUse, etc.)
      clustObject <- hclust(as.dist((1-cor(x,use=corUse,method=distMethod))), method=hclustAlgorithm);
      #clustGap needs a list output
      output <- cutree(clustObject,k=k)
      return(output)
      
    }

    
  }else{
    
    #dist (but not cor) computes across the rows, not columns.
    #dist (but not cor) computes across the rows, not columns.
    datasetClust <- t(dataset)
    clustInput <- dist(datasetClust,method=distMethod)
    clustF <- function(clustInput,k){
      #it turns out that R will recognize the upper-level function input value in this function (algorith, corUse, etc.)
      #tried passing in parent.frame() to make a more elegant solution but didn't work.
      
      clustObject <- hclust(clustInput, method=hclustAlgorithm);
      output <- cutree(clustObject,k=k)
      return(output)
      
    }
  
  
  }
  
  }else if(clustMethod=="km"){
    #k-means clusters rows.
    datasetClust <- t(dataset)
    clustF <- function(x,k){
  
      #it turns out that R will recognize the upper-level function input value in this function.
      #transpose BEFORE feed in here; I found it to return odd results if
      #I feed in t(x) in the kmeans function that is feed into clusGap
      output <- kmeans(x, centers=k,iter.max=iter.max,nstart=nstart,algorithm="Hartigan-Wong")$cluster
      return(output)
      
    }
    
  }
  
  
    numComm <- list();
    commTest <- list();
    medianComm <- c();
    membershipMatrixList <- list()

    for(k in 2:kMax){
      

      numComm[[k]] <- array(data=NA,dim=numIter);
      #FIRST: main clustering. need to "arrange" patients in clusters to make sure get a reasonably balanced number of patients.
      #clusterAssignments <- clusterMembershipFunction(dataMatrix,k);
      
      membershipMatrixList[[k]] <- list()
      
      if(!splitDataBeforeClust){
        
      if(clustMethod == "hc"){
        #hc always gives same outputs, so can run just once.
        clusterAssignments <- clustF(datasetClust,k);
        
      }
      
      }
      
      for(i in 1:numIter){
        
        message("Tesing out k = ",k, " iteration ",i)
        #NOTE: if want to resample from your dataset, redefine a new datasetClust each time: do it here. 
        #I don't think it's quite necessary. even if have a few outliers, they'll usually be included in the resampling by random chance.
        if(clustMethod != "hc" && !splitDataBeforeClust){
          #must run each time because k-means has random starts.
        #run on each iteration: in case clustering algorithm has a random aspect.?
        clusterAssignments <- clustF(datasetClust,k);
        
        }
        #rnd_indices <- sample(1:nrow(dataMatrix),size=nrow(dataMatrix),replace=FALSE);
        #rnd1 <- rnd_indices[c(1:ceiling(nrow(dataMatrix)/2))];
        #rnd2 <- rnd_indices[c((ceiling(nrow(dataMatrix)/2)+1):nrow(dataMatrix))];
        
        #clusterAssignments1 <- clusterMembershipFunction(dataMatrix[rnd1, ],k);
        #clusterAssignments2 <- clusterMembershipFunction(dataMatrix[rnd2, ],k);
        
        d1Indices <- c();
        d2Indices <- c();
        clustMatrixList <- list();
        clustSampleIndexList <- list();
        clustFeatureIndexList <- list();
        dataMatrixList <- list()
        
        clustSampleIndexList[[1]] <- list()
        clustSampleIndexList[[2]] <- list()
        clustFeatureIndexList[[1]] <- list()
        clustFeatureIndexList[[2]] <- list()
        #go through each cluster.
 
        if(!splitDataBeforeClust){
          
        for(c in 1:k){
          
          clusterAssignmentsC <- which(clusterAssignments==c);
          rnd_indices <- sample(clusterAssignmentsC,size=length(clusterAssignmentsC),replace=FALSE);
          rnd1 <- rnd_indices[c(1:ceiling(length(clusterAssignmentsC)/2))];
          rnd2 <- rnd_indices[c((ceiling(length(clusterAssignmentsC)/2)+1):length(clusterAssignmentsC))];
          
          #this is too random - we want balanced # samples from a cluster split across the two fake datasets.
          #clusterAssignments1 <- clusterMembershipFunction(dataMatrix[rnd1, ],k);
          #clusterAssignments2 <- clusterMembershipFunction(dataMatrix[rnd2, ],k);
          
          
          clustSampleIndexList[[1]][[c]] <- rnd1
          clustFeatureIndexList[[1]][[c]] <- c(1:nrow(dataset))
          clustSampleIndexList[[2]][[c]] <- rnd2
          clustFeatureIndexList[[2]][[c]] <-  c(1:nrow(dataset))
          d1Indices <- append(d1Indices,rnd1);
          d2Indices <- append(d2Indices,rnd2)
          
        }    
        
        if(is.null(d1Indices) || is.null(d2Indices)){
          
          stop("Not computing d1Indices or d2Indices correctly.")
        }
        
        if(length(intersect(d1Indices,d2Indices))>0){
          
          stop("Not computing split dataset indices correctly.")
        }
      
        dataMatrixList[[1]] <- dataMatrix[ ,d1Indices];
        dataMatrixList[[2]] <- dataMatrix[ ,d2Indices];
        
        #NOW: must update indices w.r.t. new split datasets
        
        for(c in 1:k){
                 
          clustSampleIndexList[[1]][[c]] <- na.omit(match(colnames(dataMatrix[,clustSampleIndexList[[1]][[c]]]),colnames(  dataMatrixList[[1]])))
          clustSampleIndexList[[2]][[c]] <- na.omit(match(colnames(dataMatrix[, clustSampleIndexList[[2]][[c]]]),colnames(  dataMatrixList[[2]])))
        
        }
        
      }else{
        
        rnd_indices <- sample(c(1:ncol(dataset)),size=ncol(dataset),replace=FALSE);
        d1Indices <- rnd_indices[c(1:ceiling(ncol(dataset)/2))];
        d2Indices <- rnd_indices[c((ceiling(ncol(dataset)/2)+1):ncol(dataset))];
        
        
        dataMatrixList[[1]] <- dataMatrix[ ,d1Indices];
        dataMatrixList[[2]] <- dataMatrix[ ,d2Indices];
        
        if(length(intersect(d1Indices,d2Indices))>0){
          
          stop("Not computing split dataset indices correctly.")
        }
        
        
        if(all(colnames(dataset)==colnames(datasetClust))){
          
          clusterAssignments1 <- clustF(datasetClust[,d1Indices],k);
          clusterAssignments2 <- clustF(datasetClust[,d2Indices],k);
          
        }else{
          
          clusterAssignments1 <- clustF(datasetClust[d1Indices,],k);
          clusterAssignments2 <- clustF(datasetClust[d2Indices,],k);
        
        }
          
        for(c in 1:k){

          clustSampleIndexList[[1]][[c]] <- which(clusterAssignments1==c);
          clustSampleIndexList[[2]][[c]] <- which(clusterAssignments2==c);
          clustFeatureIndexList[[1]][[c]] <- c(1:nrow(dataset))
          clustFeatureIndexList[[2]][[c]] <- c(1:nrow(dataset))
        
      }
      
      }
        
        #run analysis.

        adjMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                            edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                            sigMethod=sigMethod,maxNullFractSize=maxNullFractSize,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                            
                                            outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                            numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, clustSizeFractThresh=clustSizeFractThresh)
        
  
        edgeOutput <-  assignFinalEdges(computeTrueSimilOutput= adjMatrix$computeTrueSimilOutput,
                                      pvalueMatrix= adjMatrix$pvalueMatrix,indEdgePvalueThresh=indEdgePvalueThresh,
                                      meanEdgePairPvalueThresh=meanEdgePairPvalueThresh,
                                      minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                      fractFeatIntersectThresh=fractFeatIntersectThresh,numFeatIntersectThresh=numFeatIntersectThresh ,
                                      clustSizeThresh=clustSizeThresh, clustSizeFractThresh= clustSizeFractThresh,saveDir=saveDir,fileTag="selectNumTmp"
        )
        
        if(nrow(edgeOutput$filterEdgeOutput$edgeMatrix)>0){
        
          #if 1 community: should be dropped. means clusters got too noisy to form strong edges.
          communityInfo <- findCommunities(edgeMatrix=edgeOutput$filterEdgeOutput$edgeMatrix,edgeWeightMatrix=edgeOutput$filterEdgeOutput$edgeWeightMatrix,
                                         adjMatrix$clustIndexMatrix,fileTag="selectNumTmp",
                                         saveDir=saveDir,minNumUniqueStudiesPerCommunity=2,clustMethodName=clustMethod,
                                         commMethod=commMethod,
                                         makePlots=FALSE,saveGraphData=FALSE)
        
        
  
          numComm[[k]][i] <- communityInfo$numCommunities
         
        #make a membership matrix to get consensus matrices later on...
        #if any communities are dropped: they become NA (stop and take k-1 if this happens?)
        membershipMatrixList[[k]][[i]] <- returnSampleMemberMatrix(clustSampleIndexList=clustSampleIndexList,
                                                                   dataMatrixList=dataMatrixList,communityAttrDF=communityInfo$attrDF)$fullMemberMatrix
      
       #no final communities found 
      }else{
        
        membershipMatrixList[[k]][[i]] <- NA
        numComm[[k]][i] <- 0
        
      }
      message( numComm[[k]][i] ," final communities")
        #end of loop i
      }
      
      medianComm[k] <- median(unlist(numComm[[k]]));
      
      commTest[k] <- medianComm[k]==k;
      
      if(is.null(commTest[k])){
        #if NULL: change to NA so can keep correct indices when take max outside of loop
        commTest[k] <- NA
      }
      
      #end of looping over all K's
    }
  
  if(bestKSelect="max"){
  bestK <- max(which(unlist(commTest)))
  
  }else if(bestKSelect="min"){
    
    bestK <- min(which(unlist(commTest)))
  }
  
  if(length(bestK)==0 && medianComm[2]==1){
    
    bestK <- 1
    consensusClustAssignments <- rep.int(1, times=ncol(dataMatrix))
    
  }else{
    #not even 1 cluster!
    bestK <- NA
    consensusClustAssignments <- rep.int(NA, times=ncol(dataMatrix))
    
  }
  
  realBestK <- bestK
  
  if(!is.na(bestK)){
    
    if(bestK ==1){
      
      bestK <- 2
      #let's take the consensus matrix from k=2; either all clusters will be globbed together, or perhaps there's
      #a chunk of samples that are complete noise we should throw out
      #this could still happen with bestK >1, but we have no way of really confirming this...
      #can always remove these samples and re-run the analysis to see if that makes for clearer and higher # K clusters.
    }
  #now calculate final assignments
  consensusMatrix <- membershipMatrixList[[bestK]][[1]]
  #set all values to NA?
  consensusMatrix[c(1:ncol(consensusMatrix)),c(1:ncol(consensusMatrix))] <- 0
  numNA <- consensusMatrix
  
    for(n in 1:numIter){
      
      #basically: need to sum up whenever there is NOT an NA
      #BUT to divide later on...must keep the number of NA memberships a sample-sample pair ever had
      #temporarily make NAs in consnsusMatrix 0 so can add if this iteration of membershipMatrixList[[n]] has a non-NA values.
      tmp <- consensusMatrix
      tmp[which(is.na(tmp))] <- 0
      numNA[which(is.na(membershipMatrixList[[bestK]][[n]]))] <- numNA[which(is.na(membershipMatrixList[[bestK]][[n]]))] +1
      consensusMatrix[which(!is.na(membershipMatrixList[[bestK]][[n]]))] <- tmp[which(!is.na(membershipMatrixList[[bestK]][[n]]))] + membershipMatrixList[[bestK]][[n]][which(!is.na(membershipMatrixList[[bestK]][[n]]))]
      
    }
  
  #COME BACK:
  #CAREFUL: what if numIter-numNA =numIter? then you're dividing by zero...but top value should be an NA anways and then will just change to Inf
  consensusMatrix <- consensusMatrix/(numIter-numNA)
  if(realBestK != 1 && any(consensusMatrix==Inf)){
    
    stop("Trying running select K method with a greater numIter input; 
         some of your samples never made it into a resampled dataset with the current numIter so the consensusMatrix has Inf values.")
  }else if(any(consensusMatrix==Inf)){
    #these Inf samples are most likely very noisy samples that fall out as cluster, but when split into subclusters,
    #can't find significance between them.
    samplesRemove <- apply(consensusMatrix,MARGIN=1,FUN=function(consensusRow){
      
      if(all(consensusRow==Inf)){
        
        return(TRUE)
        
      }else{
        
        return(FALSE)
        
      }
      
    })
    
  }else{
    
    samplesRemove <- c()
    
  }
  #how cluster these -Inf relationships remove them if ALL of the time this sample never had one relationship with another cluster?
  #that would be odd to happen by random chance...
  #then just 
  
  #reset K if did not alreayd for bestK=1 scenario
  bestK <- realBestK
  
  if(bestK !=1){
    
    consensusClustAssignments <- hclust(consensusMatrix,method="average",k=bestK)
    
  }else{
    
    #remove these samples
    if(length(samplesRemove)>0){
      
      consensusClustAssignments[samplesRemove] <- NA
      
    }
    
  }
  
  
      #WANT: k-1 where first k is k AFTER hit at least one truth, where medianComm[k] < k (went one too far.)
      #CHECK: what if this misses a certain case??? like randomly, a k that's actually too low passes this test?
      #if(stopAtBestK){
      
      #  if(!is.na(medianComm[k]) && medianComm[k] < k){
      
      #    output <- list(bestK=(k-1),numComm=numComm,medianComm=medianComm,kVector=c(2:k));
      #    return(output);
      
      #  }
      
  
  #now: compute a clustering with k clusters on our consensus matrix derived from all matrices in membershipMatrixList[[k]]
  output <- list(bestK=bestK,numComm=numComm,medianComm=medianComm,consensusMatrix=consensusMatrix,consensusClustAssignments=consensusClustAssignments);
  
  }else{
    
    output <- list(bestK=bestK,numComm=numComm,medianComm=medianComm,consensusClustAssignments=consensusClustAssignments)
  }
  return(output);

#EOF
}

####START HERE

CoINcIDE_selectK_hclust <- function(dataMatrix,clustFeatures,
                                            edgeMethod=c("distCor","spearman","pearson","kendall","Euclidean","cosine",
                                                         "Manhattan","Minkowski","Mahalanobis"),numParallelCores=1,minTrueSimilThresh=-Inf,maxTrueSimilThresh=Inf,
                                            sigMethod=c("meanMatrix","centroid"),maxNullFractSize=.1,numSims=100,includeRefClustInNull=TRUE,
                                            
                                            outputFile="./CoINcIDE_messages.txt",fractFeatIntersectThresh=0,numFeatIntersectThresh=0 ,clustSizeThresh=0, clustSizeFractThresh=0,
  distMethod=c("pearson","spearman","euclidean", "binary", "maximum", "canberra", "minkowski"),
  hclustAlgorithm=c("average","complete","ward.D", "ward.D2", "single", "mcquitty","median","centroid"), 
  clustMethod=c("km","hc"), saveDir="/home/kplaney/ovarian_analysis/", indEdgePvalueThresh=.1,corUse="everything"
  meanEdgePairPvalueThresh=.05, commMethod=c("edgeBetween"),bestKSelect=c("max","min"),numIter=20,maxNumClusters=15,
  iter.max=30,nstart=10,numSubDatasets=2,computeDistMatrixOnce=TRUE
  
  ){
  
  dataset <- dataMatrix[rownames(dataMatrix) %in% clustFeatures, , drop=FALSE]
  
  if(nrow(dataset)==0){
    
    stop("\nIn clusterMatrixKmeansGap function: no clustFeatures were found in the data matrix inputted.")
    
  }

  warning("Assumes samples are in the columns.")
  #must have unique row names!
  if(any(duplicated(colnames(dataMatrix)))){   
    colnames(dataMatrix) <- paste0(colnames(dataMatrix),"_",c(1:ncol(dataMatrix)));
  }
  
  if(ncol(dataset)<maxNumClusters){
    #hclust usually returns NA gap test if kMax = ncol(dataset)
    #will also mest up maxSE calculations
    kMax <- ncol(dataset)-1
    
  }else{
    
    kMax <- maxNumClusters
    
  }
  
  
  if(computeDistMatrixOnce){
    
  if(distMethod==("pearson") || distMethod=="spearman"){
    
    datasetClust <- dataset
    distMatrix <- as.matrix(as.dist((1-cor(datasetClust,use=corUse,method=distMethod))))
    
#     ##let's chunk this up?
#     numChunks <- ceiling(ncol(dataset)/2000)
#     
#     if(numChunks>1){
#       
#       distMatrix <- matrix(data=NA,ncol=ncol(dataset),nrow=nrow(dataset))
#     
#     #do diagonals first
#     count <- 0
#     for(n in 1:numChunks){
#       
#     patIndices <-         
#     distMatrix[]  as.matrix(as.dist((1-cor(datasetClust,use=corUse,method=distMethod))))
#    
#       
#     }
# 
#     }else{
#       
#       distMatrix <- as.matrix(as.dist((1-cor(datasetClust,use=corUse,method=distMethod))))
#  
#       
#     }
    
  }else{
    
    #dist (but not cor) computes across the rows, not columns.
    #dist (but not cor) computes across the rows, not columns.
    datasetClust <- t(dataset)
    distMatrix <- as.matrix(dist(datasetClust,method=distMethod))

     
  }
  }
    clustF <- function(distMatrix,k){
      #it turns out that R will recognize the upper-level function input value in this function (algorith, corUse, etc.)
      clustObject <- hclust(distMatrix, method=hclustAlgorithm);
      #clustGap needs a list output
      output <- cutree(clustObject,k=k)
      return(output)
      
    }
  
  
    numComm <- list();
    commTest <- list();
    medianComm <- c();
    membershipMatrixList <- list()
    ml <- list()
    #ml is the consensus matrix for each K run.
  
    mCount = mConsist = matrix(c(0),ncol=ncol(dataset),nrow=ncol(dataset))

      for(i in 1:numIter){
        
        if (i==1){
          
          ml[[k]] = mConsist #initialize
        
        }
    
    ##mCount is possible number of times that two sample occur in same random sample, independent of k
    ##mCount stores number of times a sample pair was sampled together.  
    #will count up over each numIter rep.
    mCount <- connectivityMatrix( rep( 1,ncol(dataset)),
                                  mCount,
                                  c(1:ncol(dataset)) ) 
    
    
        message("Tesing out k = ",k, " iteration ",i)
        #NOTE: if want to resample from your dataset, redefine a new datasetClust each time: do it here. 
        #I don't think it's quite necessary. even if have a few outliers, they'll usually be included in the resampling by random chance.
      
        
        clusterAssignments <- clustF(datasetClust,k);
        
        clustSampleIndexList <- list()
        clustFeatureIndexList <- list()
        dataMatrixList <- list()
        
        rnd_indices <- sample(c(1:ncol(dataset)),size=ncol(dataset),replace=FALSE);
        names(rnd_indices) <- rep.int(c(1:numSubDatasets), times=ceiling(ncol(dataset)/numSubDatasets))[1:ncol(dataset)]
        
        for(k in 2:kMax){
      

        numComm[[k]] <- array(data=NA,dim=numIter);
        #FIRST: main clustering. need to "arrange" patients in clusters to make sure get a reasonably balanced number of patients.
        #clusterAssignments <- clusterMembershipFunction(dataMatrix,k);
      
        membershipMatrixList[[k]] <- list()
      
        #parallelize this step?
        #run on numSubDatasets
        for(r in 1:numSubDatasets){
          
          d_indices <- rnd_indices[which(names(rnd_indices)==r)]
          dataMatrixList[[r]] <- dataset[ ,d_indices]
            
          if(computeDistMatrixOnce){

          distMatrixClust <- as.dist(distMatrix[d_indices,d_indices])
          
            }else{
                #if k-means: this is where you'd compute it too.
              
          if(distMethod==("pearson") || distMethod=="spearman"){
    
          datasetClust <- dataMatrixList[[r]]
          distMatrixClust <- as.dist((1-cor(datasetClust,use=corUse,method=distMethod)))

    
          }else{
    
          #dist (but not cor) computes across the rows, not columns.
          #dist (but not cor) computes across the rows, not columns.
          datasetClust <- t(dataMatrixList[[r]])
          distMatrixClust <- dist(datasetClust,method=distMethod)
   
        }
            
      }
      
      #now cluster.
          clusterAssign <- clustF(distMatrixClust,k)
          #now populate the clustSampleIndices
          clustSampleIndexList[[r]] <- list()
          clustFeatureIndexList[[r]] <- list()
          
          for(c in 1:k){
            
           clustSampleIndexList[[r]][[c]] <- which(clusterAssign==c)
           clustFeatureIndexList[[r]][[c]] <- c(1:nrow(dataset))
           
          }
      #end of loop R
      }
        
        #run analysis for this iteration for this k.

        adjMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                            edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                            sigMethod=sigMethod,maxNullFractSize=maxNullFractSize,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                            
                                            outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                            numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, clustSizeFractThresh=clustSizeFractThresh)
        
  
        edgeOutput <-  assignFinalEdges(computeTrueSimilOutput= adjMatrix$computeTrueSimilOutput,
                                      pvalueMatrix= adjMatrix$pvalueMatrix,indEdgePvalueThresh=indEdgePvalueThresh,
                                      meanEdgePairPvalueThresh=meanEdgePairPvalueThresh,
                                      minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                      fractFeatIntersectThresh=fractFeatIntersectThresh,numFeatIntersectThresh=numFeatIntersectThresh ,
                                      clustSizeThresh=clustSizeThresh, clustSizeFractThresh= clustSizeFractThresh,saveDir=saveDir,fileTag="selectNumTmp"
        )
        
        if(nrow(edgeOutput$filterEdgeOutput$edgeMatrix)>0){
        
          #if 1 community: should be dropped. means clusters got too noisy to form strong edges.
          communityInfo <- findCommunities(edgeMatrix=edgeOutput$filterEdgeOutput$edgeMatrix,edgeWeightMatrix=edgeOutput$filterEdgeOutput$edgeWeightMatrix,
                                         adjMatrix$clustIndexMatrix,fileTag="selectNumTmp",
                                         saveDir=saveDir,minNumUniqueStudiesPerCommunity=2,clustMethodName=clustMethod,
                                         commMethod=commMethod,
                                         makePlots=FALSE,saveGraphData=FALSE)
        
        
  
          numComm[[k]][i] <- communityInfo$numCommunities
        
 
        #make a membership matrix to get consensus matrices later on...
        #if any communities are dropped: they become NA (stop and take k-1 if this happens?)
        #later: will need to line up all membership matrices so patients are in same order.
        membershipMatrixList[[k]][[i]] <- returnSampleMemberMatrix(clustSampleIndexList=clustSampleIndexList,
                                                                   dataMatrixList=dataMatrixList,communityAttrDF=communityInfo$attrDF$fullMemberMatrix
      
   ##add to tally  			
   ml[[k]] <- connectivityMatrix( this_assignment,
                                     ml[[k]],
                                     NAMESOFSAMPLESINORDER)
       #no final communities found 
      }else{
        #COME BACK: update mcount??
        membershipMatrixList[[k]][[i]] <- NA
        #don't change ml[[k]]?  
        numComm[[k]][i] <- 0
        
      }
      message( numComm[[k]][i] ," final communities")
      
        }#end of loop k

        #end of loop i
      }
   
  
   ##compute final consensus fraction for each K
  res = vector(mode="list",maxK)
  for (k in 2:kMax){
    ##fill in other half of matrix for tally and count.
    tmp = triangle(ml[[k]],mode=3)
    tmpCount = triangle(mCount,mode=3)
    #hmm....COME BACK: does this make sense??
    res[[k]] = tmp / tmpCount
    res[[k]][which(tmpCount==0)] = 0
  }

  #re-assign
  ml <- res
  
  #make output list
  res <- list()
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
    #(plotting code from consensus cluster matrix list here.)
    res[[tk]] = list(consensusMatrix=c,consensusTree=hc,consensusClass=ct,ml=ml[[tk]])
    #colorM = rbind(colorM,colorList[[1]]) 
  }


    for(c in 1:Kmax){
  consensusClustOutput <- list()
    
     consensusClustOutput[[c]]$ml <- res[[c]]$ml
  
  }
#now select best K based off of CoINcIDE.
for(k in 2:maxK){
      medianComm[k] <- median(unlist(numComm[[k]]));
      
      commTest[k] <- medianComm[k]==k;
      
      if(is.null(commTest[k])){
        #if NULL: change to NA so can keep correct indices when take max outside of loop
        commTest[k] <- NA
      }
      
}


  if(bestKSelect="max"){
  bestK <- max(which(unlist(commTest)))
  
  }else if(bestKSelect="min"){
    
    bestK <- min(which(unlist(commTest)))
  }
  
  if(length(bestK)==0 && medianComm[2]==1){
    
    bestK <- 1
    consensusClustAssignments <- rep.int(1, times=ncol(dataMatrix))
    
  }else{
    #not even 1 cluster!
    bestK <- NA
    consensusClustAssignments <- rep.int(NA, times=ncol(dataMatrix))
    
  }
  
  realBestK <- bestK
  
  if(!is.na(bestK)){
    
    if(bestK ==1){
      
      bestK <- 2
      #let's take the consensus matrix from k=2; either all clusters will be globbed together, or perhaps there's
      #a chunk of samples that are complete noise we should throw out
      #this could still happen with bestK >1, but we have no way of really confirming this...
      #can always remove these samples and re-run the analysis to see if that makes for clearer and higher # K clusters.
    }
   
# 
#   for(c in 1:length(Kmax)){
#   #now calculate final assignments
#   consensusMatrix <- membershipMatrixList[[c]][[1]]
#   #set all values to NA?
#   consensusMatrix[c(1:ncol(consensusMatrix)),c(1:ncol(consensusMatrix))] <- 0
#   numNA <- consensusMatrix
#   
#     for(n in 1:numIter){
#       
#       #basically: need to sum up whenever there is NOT an NA
#       #BUT to divide later on...must keep the number of NA memberships a sample-sample pair ever had
#       #temporarily make NAs in consnsusMatrix 0 so can add if this iteration of membershipMatrixList[[n]] has a non-NA values.
#       tmp <- consensusMatrix
#       tmp[which(is.na(tmp))] <- 0
#       numNA[which(is.na(membershipMatrixList[[bestK]][[n]]))] <- numNA[which(is.na(membershipMatrixList[[c]][[n]]))] +1
#       consensusMatrix[which(!is.na(membershipMatrixList[[c]][[n]]))] <- tmp[which(!is.na(membershipMatrixList[[c]][[n]]))] + membershipMatrixList[[c]][[n]][which(!is.na(membershipMatrixList[[c]][[n]]))]
# 
#     }
#   
#   #COME BACK:
#   #CAREFUL: what if numIter-numNA =numIter? then you're dividing by zero...but top value should be an NA anways and then will just change to Inf
#   consensusMatrix <- consensusMatrix/(numIter-numNA)
#   
#     for(c in 1:Kmax){
#   consensusClustOutput <- list()
#     
#      consensusClustOutput[[c]]$ml <- consensusMatrix
#   
#   }
#   
#   }
  if(realBestK != 1 && any(res[[realBestK]]$ml==Inf)){
    
    stop("Trying running select K method with a greater numIter input; 
         some of your samples never made it into a resampled dataset with the current numIter so the consensusMatrix has Inf values.")
  }else if(any(consensusMatrix==Inf)){
    #these Inf samples are most likely very noisy samples that fall out as cluster, but when split into subclusters,
    #can't find significance between them.
    samplesRemove <- apply(consensusMatrix,MARGIN=1,FUN=function(consensusRow){
      
      if(all(consensusRow==Inf)){
        
        return(TRUE)
        
      }else{
        
        return(FALSE)
        
      }
      
    })
    
  }else{
    
    samplesRemove <- c()
    
  }
  #how cluster these -Inf relationships remove them if ALL of the time this sample never had one relationship with another cluster?
  #that would be odd to happen by random chance...
  #then just 
  
  #reset K if did not alreayd for bestK=1 scenario
  bestK <- realBestK
  
  if(bestK !=1){
    
    consensusClustAssignments <- hclust(1-consensusMatrix,method="average",k=bestK)
    
  }else{
    
    #remove these samples
    if(length(samplesRemove)>0){
      
      consensusClustAssignments[samplesRemove] <- NA
      
    }
    
  }
  
  
      #WANT: k-1 where first k is k AFTER hit at least one truth, where medianComm[k] < k (went one too far.)
      #CHECK: what if this misses a certain case??? like randomly, a k that's actually too low passes this test?
      #if(stopAtBestK){
      
      #  if(!is.na(medianComm[k]) && medianComm[k] < k){
      
      #    output <- list(bestK=(k-1),numComm=numComm,medianComm=medianComm,kVector=c(2:k));
      #    return(output);
      
      #  }
      
  
  #now: compute a clustering with k clusters on our consensus matrix derived from all matrices in membershipMatrixList[[k]]
  output <- list(bestK=bestK,numComm=numComm,medianComm=medianComm,consensusMatrix=consensusMatrix,consensusClustAssignments=consensusClustAssignments);
  
  }else{
    
    output <- list(bestK=bestK,numComm=numComm,medianComm=medianComm,consensusClustAssignments=consensusClustAssignments)
  }
  return(output);

#EOF
}

#atken from consensus cluster plus
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

library("doParallel")
library("foreach")
library("matrixStats")
#  dataset <- dataMatrix[rownames(dataMatrix) %in% clustFeatures, , drop=FALSE]
chunkUpSymmDistMatrixCalc <- function(dataset,distMethod="Euclidean",numParallelCores=1){
 
                                  registerDoParallel(cores=numParallelCores)

#distMatrixWhole <- matrix(data=NA,ncol=ncol(dataset),nrow=ncol(dataset),dimnames=list(colnames(dataset),colnames(dataset)))
                          
distMatrixWhole <-  foreach(p=1:ncol(dataset),.combine='rbind') %dopar%{
    
      numNAs <- ncol(dataset) - p
      
      if(distMethod==("pearson") || distMethod=="spearman"){
    
    
    distMatrix <-  append(as.matrix(as.dist((1-cor(dataset[,p],dataset[,c(p:ncol(dataset))], use=corUse,method=distMethod))))[1,],
                           rep.int(NA,numNAs))

  
                          
    }else{
    
      #need a lower triangle matrix.
    #dist (but not cor) computes across the rows, not columns.
    #dist (but not cor) computes across the rows, not columns.

    #last numNAs will be zero, but append can handle this.
    distMatrix <- append(as.matrix(dist(dataset[,p],dataset[,c(1:p)], method=distMethod,by_rows=FALSE))[1,],
                         rep.int(NA,numNAs))

    #distMatrixWhole[p, ] <- distMatrix[1,]
     
  }
    
  }
  
   rownames(distMatrixWhole) <- colnames(dataset)
   colnames(distMatrixWhole) <- colnames(dataset)

#save(distMatrixWhole,file="Harvard_distMatrixWhole_cui_6Months_cosine_2015-06-08.RData.gzip",compress="gzip")

  #no need to copy over to create a full distance matrix - as.dist() works on a lower triangle matrix.
  #distMatrixWhole[upper.tri(distMatrixWhole,diag=FALSE)] <- distMatrixWhole[lower.tri(distMatrixWhole,diag=FALSE)]

  return(distMatrixWhole)
  
  }