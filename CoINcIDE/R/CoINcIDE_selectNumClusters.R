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
  clustMethod=c("km","hc"), saveDir="/home/kplaney/ovarian_analysis/", indEdgePvalueThresh=.1,
  meanEdgePairPvalueThresh=.05, commMethod=c("edgeBetween"),bestKSelect=c("max","min")
  
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
    #hclust usually returns NA gap test if K.max = ncol(dataset)
    #will also mest up maxSE calculations
    K.max <- ncol(dataset)-1
    
  }else{
    
    K.max <- maxNumClusters
    
  }
  
  if(clustMethod=="hc"){
    
  if(distMethod==("pearson") || distMethod=="spearman"){
    
    datasetClust <- dataset
    
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
    datasetClust <- t(dataset)
    clustF <- function(x,k){
      #it turns out that R will recognize the upper-level function input value in this function (algorith, corUse, etc.)
      #tried passing in parent.frame() to make a more elegant solution but didn't work.
      
      clustObject <- hclust(dist(x,method=distMethod), method=hclustAlgorithm);
      output <- list()
      output$cluster <- cutree(clustObject,k=k)
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
    kmeans(x, centers=k,iter.max=iter.max,nstart=nstart,algorithm=algorithm)$cluster;
    
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
    
    for(i in 1:numIter){
      #NOTE: if want to resample from your dataset, redefine a new datasetClust each time: do it here. 
      #I don't think it's quite necessary. even if have a few outliers, they'll usually be included in the resampling by random chance.

      #run on each iteration: in case clustering algorithm has a random aspect.?
      clusterAssignments <- clustF(datasetClust,k);
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
      sampleCount1 <- 1
      sampleCount2 <- 1
      
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
        #THIS IS MESSED UP>>>>
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
      
      #run analysis.
      #COME BACK: if can't find sample or feature name in end dataset.
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
    
      #end of loop i
    }
    
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

if(bestK !=1 || !is.na(bestK)){
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
if(any(consensusMatrix==Inf)){
  
  stop("Trying running select K method with a greater numIter input; 
       some of your samples never made it into a resampled dataset with the current numIter so the consensusMatrix has Inf values.")
}
#how cluster these -Inf relationships remove them if ALL of the time this sample never had one relationship with another cluster?
#that would be odd to happen by random chance...
#then just 

consensusClustAssignments <- hclust(consensusMatrix,method="average",k=bestK)


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