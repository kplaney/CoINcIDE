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
    clustF <- function(x,k){
      #it turns out that R will recognize the upper-level function input value in this function (algorith, corUse, etc.)
      #tried passing in parent.frame() to make a more elegant solution but didn't work.
      
      clustObject <- hclust(dist(x,method=distMethod), method=hclustAlgorithm);
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