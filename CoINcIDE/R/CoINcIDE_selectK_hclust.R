

#sourceDir <- "/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/"
#setwd(sourceDir)
#source("CoINcIDE_computeEdges.R")
#source("CoINcIDE_communityDetection.R")
library("proxy")
CoINcIDE_selectK_hclust <- function(dataMatrix,clustFeatures,
                                    edgeMethod=c("distCor","spearman","pearson","kendall","Euclidean","cosine",
                                                 "Manhattan","Minkowski","Mahalanobis"),numParallelCores=1,minTrueSimilThresh=-Inf,maxTrueSimilThresh=Inf,
                                    sigMethod=c("meanMatrix","centroid"),maxNullFractSize=.1,numSims=100,includeRefClustInNull=TRUE,
                                     outputFile="./CoINcIDE_messages.txt",clustSizeThresh=0, clustSizeFractThresh=0,
                                    hclustDistMethod=c("pearson","spearman","euclidean", "binary", "maximum", "canberra", "minkowski"),
                                    hclustAlgorithm=c("average","complete","ward.D", "ward.D2", "single", "mcquitty","median","centroid"), 
                                    clustMethod=c("km","hc"), saveDir="/home/kplaney/ovarian_analysis/", indEdgePvalueThresh=.1,corUse="everything",
                                    meanEdgePairPvalueThresh=.05, commMethod=c("edgeBetween"),bestKSelect=c("max","min"),numIter=20,maxNumClusters=15,
                                    iter.max=30,nstart=10,numSubDatasets=2,computeDistMatrixOnce=TRUE,distMatrix=NULL,
                                    centroidMethod=c("mean","median"),consensusLinkage="average"
                                    
){
  
  adjMatrixList <- list()
  consensusData <- list()
  
  if(length(centroidMethod)>1){
    
    centroidMethod <- "mean"
    
  }
  
  if(length(bestKSelect)>1){
    
    bestKSelect <- "max"
    
  }
  

    
  message("This function assumes you are clustering the columns of your data matrix.")
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
  
 
  
  if(is.null(distMatrix)){
    
    if(computeDistMatrixOnce){
      
      if(distMethod==("pearson") || distMethod=="spearman"){
        
        datasetClust <- dataset
        distMatrix <- as.matrix(as.dist((1-cor(datasetClust,use=corUse,method=hclustDistMethod))))

      }else{
        
        #dist (but not cor) computes across the rows, not columns.
        #dist (but not cor) computes across the rows, not columns.
        distMatrix <- as.matrix(dist(dataset,method=hclustDistMethod,by_rows=FALSE))
        
        
      }
      
    }
   
  }else{
    
    #also: this matrix must be symmetric, not just a lower triangle.
    #subsetting like this would mix up which (upper or lower) had NAs for all rows.
    #diag=0 avoids NaN issues.
    
    distMatrix[upper.tri(distMatrix)] <- 0
    distMatrix <- distMatrix + t(distMatrix)
    
    diag(distMatrix) <- 0
    
  } #else: do nothing. distMatrix already computed
  
  clustF <- function(distMatrixTmp,k){
    #it turns out that R will recognize the upper-level function input value in this function (algorith, corUse, etc.)
    clustObject <- hclust(distMatrixTmp, method=hclustAlgorithm);
    #clustGap needs a list output
    output <- cutree(clustObject,k=k)
    return(output)
    
  }
  
  
  numComm <- list();
  #membershipMatrixList <- list()
  ml <- list()
  #ml is the consensus matrix for each K run.
  
  mConsist <- matrix(c(0),ncol=ncol(dataset),nrow=ncol(dataset))
  mCount <- mConsist
  
  for(d in 1:kMax){
    
    numComm[[d]] <- array(data=NA,dim=numIter);
    
  }
  
  for(i in 1:numIter){
    
    consensusData[[i]] <- list()
    ##mCount is possible number of times that two sample occur in same random sample, independent of k
    ##mCount stores number of times a sample pair was sampled together.  
    #will count up over each numIter rep.
    mCount <- connectivityMatrix( rep( 1,ncol(dataset)),
                                  mCount,
                                  c(1:ncol(dataset)) ) 
    
    
    
    #NOTE: if want to resample from your dataset, redefine a new datasetClust each time: do it here. 
    #I don't think it's quite necessary. even if have a few outliers, they'll usually be included in the resampling by random chance.
    
    
    #clusterAssignments <- clustF(datasetClust,k);

    dataMatrixList <- list()
    
    rnd_indices <- sample(c(1:ncol(dataset)),size=ncol(dataset),replace=FALSE);
    #randomly assign samples to different sub datasets
    names(rnd_indices) <- rep.int(c(1:numSubDatasets), times=ceiling(ncol(dataset)/numSubDatasets))[1:ncol(dataset)]
    
    #parallelize this step?
    #run on numSubDatasets
    distMatrixClust <- list()
    for(r in 1:numSubDatasets){
      
      d_indices <- rnd_indices[which(names(rnd_indices)==r)]
      dataMatrixList[[r]] <- dataset[ ,d_indices]
      
      if(computeDistMatrixOnce || !is.null(distMatrix)){
        
         #also: this matrix is already symmetric. if was only lower or upper,
        #subsetting like this would mix up which (upper or lower) had NAs for all rows.
        distMatrixClust[[r]] <- as.dist(distMatrix[d_indices,d_indices])
        
        
      }else{
        #if k-means: this is where you'd compute it too.
        
        if(distMethod==("pearson") || distMethod=="spearman"){
          
          datasetClust <- dataMatrixList[[r]]
          distMatrixClust[[r]] <- as.dist((1-cor(datasetClust,use=corUse,method=hclustDistMethod)))
          
          
        }else{
          
          #dist (but not cor) computes across the rows, not columns.
          #dist (but not cor) computes across the rows, not columns.
          datasetClust <- dataMatrixList[[r]]
          distMatrixClust[[r]] <- dist(datasetClust,method=hclustDistMethod,by_rows=FALSE)
          
        }
        
      }
      
    }
    
    for(k in 2:kMax){
      
      consensusData[[i]][[k]] <- list()
      clustSampleIndexList <- list()
      clustFeatureIndexList <- list()
      
      if (i==1){
        
        ml[[k]] = mConsist #initialize
        
      }
      
      message("Tesing out k = ",k, " iteration ",i)
      
      #FIRST: main clustering. need to "arrange" patients in clusters to make sure get a reasonably balanced number of patients.
      #clusterAssignments <- clusterMembershipFunction(dataMatrix,k);
      
      for(r in 1:numSubDatasets){
  
        #now cluster.
       # diag(distMatrixClust) <- 0
        clusterAssign <- clustF(distMatrixClust[[r]],k)
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
      #features will always match up here - they're from the same dataset! (well unless you're clustering genes and not patients...)
      adjMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                           edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                           sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,                                  
                                           outputFile=outputFile,fractFeatIntersectThresh=0,numFeatIntersectThresh=0 ,clustSizeThresh=clustSizeThresh,clustSizeFractThresh=clustSizeFractThresh,
                                           checkNA=FALSE,centroidMethod=centroidMethod)
      
      
      edgeOutput <-  assignFinalEdges(computeTrueSimilOutput= adjMatrix$computeTrueSimilOutput,
                                      pvalueMatrix= adjMatrix$pvalueMatrix,indEdgePvalueThresh=indEdgePvalueThresh,
                                      meanEdgePairPvalueThresh=meanEdgePairPvalueThresh,fractFeatIntersectThresh=0,numFeatIntersectThresh=0,
                                      minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                      clustSizeThresh=clustSizeThresh, clustSizeFractThresh= clustSizeFractThresh,saveDir=saveDir,fileTag="selectNumTmp"
      )
      
      #save last one
      if(i==1){

        adjMatrixList[[k]] <- adjMatrix
        
      }
      
      if(nrow(edgeOutput$filterEdgeOutput$edgeMatrix)>0){
        
        #if 1 community: should be dropped. means clusters got too noisy to form strong edges.
        communityInfo <- findCommunities(edgeMatrix=edgeOutput$filterEdgeOutput$edgeMatrix,edgeWeightMatrix=edgeOutput$filterEdgeOutput$edgeWeightMatrix,
                                         adjMatrix$clustIndexMatrix,fileTag="selectNumTmp",
                                         saveDir=saveDir,minNumUniqueStudiesPerCommunity=numSubDatasets,
                                         commMethod=commMethod,
                                         makePlots=FALSE,saveGraphData=FALSE)
        
        
        
        numComm[[k]][i] <- communityInfo$numCommunities
        
        if(length(numComm[[k]][i])>0){
          
          #if not in a community: just set to zero.
          aggregateData <- returnSampleMemberMatrix(clustSampleIndexList=clustSampleIndexList,
                                                    dataMatrixList=dataMatrixList,communityInfo=communityInfo,
                                                    noClustNA=FALSE)
          
          
          
          #get patients that were actually assigned to a community.
          patientNonNAindices <- which(!is.na(aggregateData$sampleClustCommKey$community))
          
          ##add to tally      	
          #first argument: actual cluster assignments
          #third argument: numeric column indices of these samples.
          #run into issues for harvard data, so just saved instead below.
          
          #for consensus calculations: if NA, just set to zero.
      
          connectivityM <- aggregateData$fullMemberMatrix
          
          if(all(colnames(connectivityM)==colnames(ml[[k]]))){
          #for consensus calculations: don't add NAs.
          ml[[k]] <- ml[[k]] + connectivityM
          
          }else{
            
            connectivityM <- connectivityM[match(colnames(ml[[k]]),colnames(connectivityM)),match(colnames(ml[[k]]),colnames(connectivityM))]
            #for consensus calculations: don't add NAs.
            ml[[k]] <- ml[[k]] + connectivityM 
          }
            
          #  ml[[k]] <-    connectivityMatrix(aggregateData$sampleClustCommKey$community[patientNonNAindices],
           #                             ml[[k]],
            #                            na.omit(match(aggregateData$sampleClustCommKey$sampleName[patientNonNAindices],
             #colnames(dataset))))
          
          consensusData[[i]][[k]]$patientAssignment <- aggregateData$sampleClustCommKey$community[patientNonNAindices]
          consensusData[[i]][[k]]$patientIndex <-  na.omit(match(aggregateData$sampleClustCommKey$sampleName[patientNonNAindices],
                                                                 colnames(dataset)))
          
          
        }
        
        #no final communities found 
      }else{
        
        #COME BACK: update ml[[k]]?? do I just not include this in the consensus matrix?
        #OR do I say they were all zero patient-patient concordances?
        #this code wouldn't work as it stands right now:
        #ml[[k]] <- connectivityMatrix(NULL,
        #                             ml[[k]],
        #                            NULL)
        #membershipMatrixList[[k]][[i]] <- NA
        #don't change ml[[k]]?  
        numComm[[k]][i] <- 0
        
      }
      message( numComm[[k]][i] ," final communities")
      
    }#end of loop k
    
    #end of loop i
  }
  
  
  ##compute final consensus fraction for each K
  res = vector(mode="list",kMax)
  for (k in 2:kMax){
    ##fill in other half of matrix for tally and count.
    #Katie: my membership matrix is already symmetric.
    #tmp = triangle(ml[[k]],mode=3)
    #mcount still from original code:
    #this takes a while to run...come back and figure this out!
    tmpCount = triangle(mCount,mode=3)
    # number of times each patient-patient pair put in a cluster/# iterations.
    res[[k]] = ml[[k]] / tmpCount
    res[[k]][which(tmpCount==0)] = 0
  }
  
  #re-assign
  ml <- res
  
  #make output list
  res <- list()
  for (tk in 2:kMax){
    fm = ml[[tk]]
    hc=hclust( as.dist( 1 - fm ), method=consensusLinkage);
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
  
  
  
  commTest <- c()
  medianComm <- c()
  
  #now select best K based off of CoINcIDE.
  for(k in 2:kMax){
    medianComm[k] <- median(unlist(numComm[[k]]));
    
    commTest[k] <- medianComm[k]==k;
    
    if(is.null(commTest[k])){
      #if NULL: change to NA so can keep correct indices when take max outside of loop
      commTest[k] <- NA
    }
    
  }
  
  
  if(bestKSelect=="max"){
    
    bestK <- max(which(unlist(commTest)))
    
  }else if(bestKSelect=="min"){
    
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
  

  
  #now: compute a clustering with k clusters on our consensus matrix derived from all matrices in membershipMatrixList[[k]]
  output <- list(adjMatrixList=adjMatrixList,bestK=bestK,numComm=numComm,medianComm=medianComm,consensusClustOutput=res,consensusClustAssignments=consensusClustAssignments);
  
  return(output);
  
  #EOF
}

#doesn't work well for sizes like Harvard's matrix, but I've already computed the membership matrix.
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

