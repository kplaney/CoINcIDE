#careful: are your features/sample in correct dimension for yoru clustermembership function? (in rows or columns?)
#usually assume samples in rows for meta-clustering code.

#EDGE case: what if there's a really small outlier? like 2 samples and the rest belong to one cluster? 
#THen the number of clusters will be correctly deemed 1, but it would be nice to know that we should throw out those outliers...


# 
#   }else if(clustMethod=="km"){
#     #k-means clusters rows.
#     datasetClust <- t(dataset)
#     clustF <- function(x,k){
#   
#       #it turns out that R will recognize the upper-level function input value in this function.
#       #transpose BEFORE feed in here; I found it to return odd results if
#       #I feed in t(x) in the kmeans function that is feed into clusGap
#       output <- kmeans(x, centers=k,iter.max=iter.max,nstart=nstart,algorithm="Hartigan-Wong")$cluster
#       return(output)
#       
#     }
#     
#   }
#   



#sourceDir <- "/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/"
#setwd(sourceDir)
#source("CoINcIDE_computeEdges.R")
#source("CoINcIDE_communityDetection.R")
#library("proxy")
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
  centroidMethod=c("mean","median"),consensusLinkage="average",kMax=10
  
  ){
  
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
      distMatrix <- as.matrix(dist(dataset,method=hclustDistMethod,by_rows=FALSE))
  
       
    }
  
  }
  #distMatrix already computed
  }
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
        
 
    ##mCount is possible number of times that two sample occur in same random sample, independent of k
    ##mCount stores number of times a sample pair was sampled together.  
    #will count up over each numIter rep.
    mCount <- connectivityMatrix( rep( 1,ncol(dataset)),
                                  mCount,
                                  c(1:ncol(dataset)) ) 
    
    

        #NOTE: if want to resample from your dataset, redefine a new datasetClust each time: do it here. 
        #I don't think it's quite necessary. even if have a few outliers, they'll usually be included in the resampling by random chance.
      
        
        #clusterAssignments <- clustF(datasetClust,k);
        
        clustSampleIndexList <- list()
        clustFeatureIndexList <- list()
        dataMatrixList <- list()
        
        rnd_indices <- sample(c(1:ncol(dataset)),size=ncol(dataset),replace=FALSE);
        #randomly assign samples to different sub datasets
        names(rnd_indices) <- rep.int(c(1:numSubDatasets), times=ceiling(ncol(dataset)/numSubDatasets))[1:ncol(dataset)]
        
        for(k in 2:kMax){
         
          if (i==1){
            
            ml[[k]] = mConsist #initialize
            
          }
    
        message("Tesing out k = ",k, " iteration ",i)

        #FIRST: main clustering. need to "arrange" patients in clusters to make sure get a reasonably balanced number of patients.
        #clusterAssignments <- clusterMembershipFunction(dataMatrix,k);

        #parallelize this step?
        #run on numSubDatasets
        for(r in 1:numSubDatasets){
          
          d_indices <- rnd_indices[which(names(rnd_indices)==r)]
          dataMatrixList[[r]] <- dataset[ ,d_indices]
            
          if(computeDistMatrixOnce || !is.null(distMatrix)){

          distMatrixClust <- as.dist(distMatrix[d_indices,d_indices])
          
            }else{
                #if k-means: this is where you'd compute it too.
              
          if(distMethod==("pearson") || distMethod=="spearman"){
    
          datasetClust <- dataMatrixList[[r]]
          distMatrixClust <- as.dist((1-cor(datasetClust,use=corUse,method=hclustDistMethod)))

    
          }else{
    
          #dist (but not cor) computes across the rows, not columns.
          #dist (but not cor) computes across the rows, not columns.
          datasetClust <- t(dataMatrixList[[r]])
          distMatrixClust <- dist(datasetClust,method=hclustDistMethod)
   
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
        
        if(nrow(edgeOutput$filterEdgeOutput$edgeMatrix)>0){
        
          #if 1 community: should be dropped. means clusters got too noisy to form strong edges.
          communityInfo <- findCommunities(edgeMatrix=edgeOutput$filterEdgeOutput$edgeMatrix,edgeWeightMatrix=edgeOutput$filterEdgeOutput$edgeWeightMatrix,
                                         adjMatrix$clustIndexMatrix,fileTag="selectNumTmp",
                                         saveDir=saveDir,minNumUniqueStudiesPerCommunity=numSubDatasets,
                                         commMethod=commMethod,
                                         makePlots=FALSE,saveGraphData=FALSE)
        
        
  
          numComm[[k]][i] <- communityInfo$numCommunities
        
         if(length(numComm[[k]][i])>0){
           
          aggregateData <- returnSampleMemberMatrix(clustSampleIndexList=clustSampleIndexList,
                                                             dataMatrixList=dataMatrixList,communityInfo=communityInfo)
          
        
 
   #get patients that were actually assigned to a community.
   patientNonNAindices <- which(!is.na(aggregateData$sampleClustCommKey$community))
   
   ##add to tally    		
   #first argument: actual cluster assignments
   #third argument: numeric column indices of these samples.
   ml[[k]] <- connectivityMatrix(aggregateData$sampleClustCommKey$community[patientNonNAindices],
                                     ml[[k]],
                                 na.omit(match(aggregateData$sampleClustCommKey$sampleName[patientNonNAindices],
                                  colnames(dataset))))
   
   
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
  
  
#   if(!is.na(bestK)){
#     
#     if(bestK ==1){
#       
#       bestK <- 2
#       #let's take the consensus matrix from k=2; either all clusters will be globbed together, or perhaps there's
#       #a chunk of samples that are complete noise we should throw out
#       #this could still happen with bestK >1, but we have no way of really confirming this...
#       #can always remove these samples and re-run the analysis to see if that makes for clearer and higher # K clusters.
#     }
   
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

#COME BACK: note if for best K, some samples have zero concordance with ALL other samples. may mean outliers.
#   if(realBestK != 1 && any(res[[realBestK]]$ml==Inf)){
#     
#     stop("Trying running select K method with a greater numIter input; 
#          some of your samples never made it into a resampled dataset with the current numIter so the consensusMatrix has Inf values.")
#   
#     }else if(any(consensusMatrix==Inf)){
#     #these Inf samples are most likely very noisy samples that fall out as cluster, but when split into subclusters,
#     #can't find significance between them.
#     samplesRemove <- apply(consensusMatrix,MARGIN=1,FUN=function(consensusRow){
#       
#       if(all(consensusRow==Inf)){
#         
#         return(TRUE)
#         
#       }else{
#         
#         return(FALSE)
#         
#       }
#       
#     })
#     
#   }else{
#     
#     samplesRemove <- c()
#     
#   }
       
  
  #now: compute a clustering with k clusters on our consensus matrix derived from all matrices in membershipMatrixList[[k]]
  output <- list(bestK=bestK,numComm=numComm,medianComm=medianComm,consensusClustOutput=res,consensusClustAssignments=consensusClustAssignments);

  return(output);

#EOF
}

#atken from consensus cluster plus
connectivityMatrix <- function(clusterAssignments, m, sampleKey){
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




##can run this in a separate function
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
selectK_PACR <- list()

if(length(meanConsensusClusterByK)> 0 && length(minConsensusClusterByK)> 0 && !(all(unlist(meanConsensusClusterByK)=="NaN")) && any(unlist(meanConsensusClusterByK)>=minMeanClustConsensus,na.rm=TRUE) && any(unlist(minConsensusClusterByK)>=minClustConsensus,na.rm=TRUE)){
  
  meanBetweenConsensusClusterByK <- list()
  consensusFrac <- list()
  consensusMetric <- list()
  PAC <- list()
  
  
  for(i in 2:K.max){
    
    N <- nrow(consensusClustOutput[[i]]$consensusMatrix)
    #technically could just take mean - is a symmetric matrix. diag=FALSE bc is a sample compared against itself (should always be 1.)
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
    
  }
  
  #add NA in index 1, otherwise when unlist, will remove first index and mess up which.max/min
  consensusFrac[[1]] <- NA
  PAC[[1]] <- NA
  ##consensusFrac calculations
  #which.max(): Missing and NaN values are discarded.
  selectK_consensusFrac$bestK <- which.max(unlist(consensusFrac))
  selectK_minConsensusClust$bestK <- which.max(unlist(minConsensusClusterByK))+1
  selectK_meanConsensusClust$bestK <- which.max(unlist(meanConsensusClusterByK))+1
  
  #round PAC - so difference at hundreth level rounded and then pick highest K at the tenth decimal level
  
  PAC_noNAs <- unlist(PAC)[which(!is.na(PAC))]
  
  if(length(PAC_noNAs)>0 && any(PAC_noNAs<maxPAC)){
    
    #this will make .00x 0 (e.g. .002 becomes zero.)
    #anything with zero in tenth, hundreth place counted as zero.
    PACR <- round(unlist(PAC),digits=2)
    possibleKs <- which(PACR==min(PACR,na.rm=TRUE))
    maxPossibleK <- possibleKs[length(possibleKs)]
    selectK_PACR$bestK  <-  maxPossibleK
    #highly unlikely two unrounded values will match, but still keep this:
    possibleKs <- which(unlist(PAC)==min(unlist(PAC),na.rm=TRUE))
    maxPossibleK <- possibleKs[length(possibleKs)]
    selectK_PAC$bestK <- maxPossibleK
    
  }else{
    
    selectK_PAC$bestK <- 1
    selectK_PACR$bestK <- 1
    PACR <- NA
    
  }
  
  if(length(selectK_consensusFrac$bestK)==0){
    #all NAs returned; set K=1
    selectK_consensusFrac$bestK <- 1
    
  }
  
  if(length(selectK_minConsensusClust$bestK)==0){
    #all NAs returned; set K=1
    selectK_minConsensusClust$bestK <- 1
    
  }
  
  if(length(selectK_meanConsensusClust$bestK)==0){
    #all  NAs returned; set K=1
    selectK_meanConsensusClust$bestK <- 1
    
  }
  
  if(length(selectK_PAC$bestK)==0){
    #all NAs returned; set K=1
    selectK_PAC$bestK <- 1
    
  }
  
  
  if(length(selectK_PACR$bestK)==0){
    #all NAs returned; set K=1
    selectK_PACR$bestK <- 1
    PACR <- NA
    
  }
  
  
}else{
  #no clusterings passed the minMeanConsensus threshold
  selectK_consensusFrac$bestK <- 1
  consensusFrac <- NA
  meanBetweenConsensusClusterByK <- NA
  PAC <- NA
  PACR <- NA
  selectK_meanConsensusClust$bestK <- 1
  selectK_minConsensusClust$bestK  <- 1
  selectK_PAC$bestK <- 1
  selectK_PACR$bestK <- 1
  
}