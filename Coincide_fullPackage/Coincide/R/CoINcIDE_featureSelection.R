#library("matrixStats")
createTestVarianceDataMatrixList  <- function(){
  
testDataMatrixList <- list();
  #5 genes each so can easily keep track of highest varying ones.
  #let's say 5 patients.
  testDataMatrixList[[1]] <- matrix(data=NA,nrow=5,ncol=5,
                                    dimnames=list(paste0("feat_",c(1:5)),
                                                  paste0("sample_",c(1:5))))
  #first row: no sd
  testDataMatrixList[[1]][1,] <- rnorm(5,mean=5,sd=0)
  #second row: no sd also, just different mean
  testDataMatrixList[[1]][2,] <- rnorm(5,mean=1,sd=0)
  #now: increasing sd. need very large differences to ensure
  #actual variance from random sampling holds.
  testDataMatrixList[[1]][3,] <- rnorm(5,mean=3,sd=.1)
  testDataMatrixList[[1]][4,] <- rnorm(5,mean=3,sd=1)
  testDataMatrixList[[1]][5,] <- rnorm(5,mean=3,sd=5)
  
  #add a very similar second dataset, but one extra gene.
  testDataMatrixList[[2]] <- matrix(data=NA,nrow=6,ncol=5,
                                    dimnames=list(paste0("feat_",c(1:6)),
                                                  paste0("sample_",c(1:5))))
  #first row: no sd
  testDataMatrixList[[2]][1,] <- rnorm(5,mean=5,sd=0)
  #second row: no sd also, just different mean
  testDataMatrixList[[2]][2,] <- rnorm(5,mean=1,sd=0)
  #now: increasing sd. need very large differences to ensure
  #actual variance from random sampling holds.
  testDataMatrixList[[2]][3,] <- rnorm(5,mean=3,sd=.1)
  testDataMatrixList[[2]][4,] <- rnorm(5,mean=3,sd=1)
  testDataMatrixList[[2]][5,] <- rnorm(5,mean=3,sd=5)
  #extra row with lots of variance
  testDataMatrixList[[2]][6,] <- rnorm(5,mean=3,sd=20)

  return(testDataMatrixList)

}

selectFeaturesMetaVariance_wrapper <- function(dataMatrixList,rankMethod=c("sd,mad"), 
                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),numFeatSelectByGlobalRank=1000,
                                       numTopFeatFromEachDataset=20,fractNATopFeatAllowedPerDataset=.2,selectMethod=c("mean","median"),
                                       outputFile="./selectFeaturesMetaVarianceOut.txt"){
  
  rankList <- computeFeatureRanks(dataMatrixList,method=rankMethod)
  
  rankMatrix <- combineFeatureRanks(rankList)
  
  output <- selectFeaturesMetaVariance(rankMatrix=rankMatrix,rankList=rankList,numNAstudiesAllowedPerFeat=numNAstudiesAllowedPerFeat,numFeatSelectByGlobalRank=numFeatSelectByGlobalRank,
                                       numTopFeatFromEachDataset=numTopFeatFromEachDataset,fractNATopFeatAllowedPerDataset=fractNATopFeatAllowedPerDataset,selectMethod=selectMethod,
                                       outputFile=outputFile)
  
  return(output)
  
  
}

computeFeatureRanks <- function(dataMatrixList,method=c("sd,mad")){
  
    if(any(is.na(unlist(dataMatrixList)))){
    
    warning("\nIn computeFeatureRanks functions found NAs in a data matrix.
            \nsd and mad will be calculated leaving out these NAs.
            \nYou will need to impute these NAs before runnin the main CoINcIDE functions.\n")
  
    }
  rankList <- lapply(dataMatrixList,FUN=function(dataMatrix,method){
    
    if(any(duplicated(rownames(dataMatrix)))){
      
      stop("Error: You have duplicated gene names. Please correct this before running this function.")
      
    }
    if(method=="sd"){
      
      output <- rowSds(dataMatrix,na.rm=TRUE)
      
    }else if(method=="mad"){
      
      output <- rowMads(dataMatrix,na.rm=TRUE)
      
    }else{
      
      stop("\nIn function computeFeatureRanks:
           Please pick either \'sd\' or \'mad\' (not both) as input for method variable.")
    
    }
  
  if(all(!is.null(names(dataMatrix)))){
    
    stop("\nError in computeFeatureRanks function: dataMatrix has no unique feature (row) names.")
  
  }
  
  names(output) <- rownames(dataMatrix)

  if(any(is.na(output))){
    
    warning("\nIn computeFeatureRanks function NAs were detected.\n
            These row indices will not be returned in the sorted/ranked output list.\n")
    
  }
  #we want the highest sd/mad FIRST.
  #it looks like sort randomly breaks ties.
  tmp <- sort(output,decreasing = TRUE,na.last=NA)
  #sort.int$ix doesn't work for vectors.
  output <- c(1:length(tmp))
  names(output) <- names(tmp)
  
  return(output)
  
  },method=method)
  
  names(rankList) <- names(dataMatrixList)
  
  return(rankList)
  
}

combineFeatureRanks <- function(rankList){
  
  uniqueFeat <- c()
  
  for(g in 1:length(rankList)){
    
    uniqueFeat <- union(uniqueFeat,names(rankList[[g]]))
    
  }
  
  rankMatrix <- matrix(data=NA,nrow=length(uniqueFeat),ncol=length(rankList),
                       dimnames=list(uniqueFeat,names(rankList)))
  
  for(c in 1:length(rankList)){
    
    rankMatrix[names(rankList[[c]]),c] <- rankList[[c]]
    
    
  }
  
  return(rankMatrix)

}

selectFeaturesMetaVariance <- function(rankMatrix,rankList, numNAstudiesAllowedPerFeat=ceiling(ncol(rankMatrix)/10),numFeatSelectByGlobalRank=1000,
                                       numTopFeatFromEachDataset=10,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("mean","median"),
                                       outputFile="./selectFeaturesMetaVarianceOut.txt"){
  
  if(!is.matrix(rankMatrix)){
    
    stop("rankMatrix input must be a 2D matrix.")
    
  }

  #what genes are not in enough studies?
  tmp <- rankMatrix
  tmp[which(!is.na(tmp))] <- 0
  tmp[which(is.na(tmp))] <- 1
  #count up #NAs across studies for each feature.
  tmp <- rowSums(tmp,na.rm=TRUE)
  
  cat(paste0(length(which(tmp>numNAstudiesAllowedPerFeat)), " features did not pass the numNAstudiesAllowedPerFeat= ",numNAstudiesAllowedPerFeat," threshold."),
      append=TRUE,file=outputFile)
  
  if(length(which(tmp>numNAstudiesAllowedPerFeat))>0){
    
    rankMatrix <- rankMatrix[-which(tmp>numNAstudiesAllowedPerFeat), ]
    
  }
  
  if(numTopFeatFromEachDataset>0){
    
    topFeatures <- lapply(rankList,FUN=function(rankUnit){
      
      if(length(rankUnit>=numTopFeatFromEachDataset)){
        
        output <- names(rankUnit)[c(1:numTopFeatFromEachDataset)]
        
      }else{
        #not numTopFeatFromEachDataset features in this list - just take entire list.
        output <- names(rankUnit)
        
      }
      
      return(output)
      
    })
    
    
    topFeatures <- unique(unlist(topFeatures))
    #don't take features that were only in a few datasets - ie that were
    #already filtered out of the rankMatrix
    topFeatures <- topFeatures[na.omit(match(rownames(rankMatrix),topFeatures))]
    
    if(length(topFeatures)==0){
      
      m <- paste0("\nIn selectFeaturesMetaVariance: No topFeatures from each dataset will be used because all of the features did not pass the numNAstudiesAllowedPerFeat=",numNAstudiesAllowedPerFeat," threshold.\n")
      message(m)
      cat(m,append=TRUE,file=outputFile)
      
    }
  }else{
    
    topFeatures <- NA
  }

  
  if(selectMethod=="mean"){
    #we are looking for the SMALL means...
    metric <- rowMeans(rankMatrix,na.rm=TRUE)
    
  }else if(selectMethod=="median"){
    
    metric <- rowMedians(rankMatrix,na.rm=TRUE)
    
  }
  
  names(metric) <- rownames(rankMatrix)

  if(length(metric)==0){
    
    stop("selectFeaturesMetaVariance function: all your rowMeans or rowMedians resulted in NA values.")
  
  }
  
  #decreasing=FALSE: want small ones (higher mean or median ranks) to be at the top.
  if(numFeatSelectByGlobalRank>length(metric)){
    
    cat(paste0("\nIn selectFeaturesMetaVariance: not enough total features to get numFeatSelectByGlobalRank=",numFeatSelectByGlobalRank,
               " so just taking the total number of features (",length(metric),")"),append=TRUE,file=outputFile)
  
    metric <- sort(metric,decreasing=FALSE)
  
  }else{
    
      metric <- sort(metric,decreasing=FALSE)[c(1:numFeatSelectByGlobalRank)]
  }

  if(numTopFeatFromEachDataset>0){
    
    uniqueFeat <- union(names(metric),topFeatures)
  
  cat(paste0(length(topFeatures)," features from numNAstudiesAllowedPerFeat= ",numNAstudiesAllowedPerFeat," threshold\n",
             "and ",numFeatSelectByGlobalRank, " from general meta-ranking for a total of ",length(uniqueFeat)," unique features.\n"),append=TRUE,file=outputFile)
  
  }else{
    uniqueFeat <- names(metric)
    
    cat(paste0("\nThere is a total of ",length(uniqueFeat)," unique meta-ranked features.\n"),append=TRUE,file=outputFile)
    
  }

  #now "best" feature will be in the first row of the matrix.
  rankMatrix <- rankMatrix[uniqueFeat, ,drop=FALSE]
  
  tmp <- rankMatrix
  tmp[which(!is.na(tmp))] <- 0
  tmp[which(is.na(tmp))] <- 1
  
  #only want to remove studies whose # NAs (count the 1's) is above this fraction.
  datasetsRemove <- which(colSums(tmp) > (nrow(rankMatrix)*fractNATopFeatAllowedPerDataset) )
  
    cat(paste0("\nThere are ",length(datasetsRemove), " datasets that do not match the fractNATopFeatAllowedPerDataset=",fractNATopFeatAllowedPerDataset," threshold.\n"),append=TRUE,file=outputFile)


  
  output <- list(finalRankMatrix=rankMatrix,finalFeatures=uniqueFeat,topFeaturesAcrossDatasets=topFeatures,
                 datasetListIndicesToRemove=datasetsRemove)
  
  return(output)

}