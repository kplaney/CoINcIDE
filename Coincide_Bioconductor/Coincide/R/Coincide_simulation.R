#library("ggplot2")
#library("s4vd")


######
#NOTE: right now, for counting true negatives: DO count clusters from own dataset
#(otherwise would have to go dataset by dataset and figure out which clusters belong to the
#same dataset.) the true positive rate is not inflated this way.
#NOTE: to make random mixtures of samples: just randomly remove 1-3 # of clusters from each data matrix list.
#NOTE: to adjust cluster sizes: just make all clusters that have max cluster size, say 100.
#then subset random # (say 0-100) for each dataset, updating sampleIndex and featureIndex list, for each random dataset
#rows are numbers of features, such as genes, that you are NOT clustering one.
createTissueSimDatasets <-  function(numSimDatasets=10,
                                   eigenValueMin = -.001,simType=c("highQualityClust","mixedClustQualityClust","unevenSizeClust"),
                                   randNumClust=FALSE,
                                   numPerClust = c(50,50,50,50),
                                   stddevNoise=0,numRows=200,minRandSize=1,maxRandSize=100,minRandNumClust=2){
  
  clustSampleIndexList <- list()
  clustFeatureIndexList <- list()
  clustFeaturesList <- list()
  numPerClustList <- list()
  for(d in 1:numSimDatasets){
    
    clustSampleIndexList[[d]] <- list()
    clustFeatureIndexList[[d]] <- list()

    
    if( simType != "unevenSizeClust"){

      numPerClustList[[d]] <- list()
      count <- 1
      
      for(c in 1:length(numPerClust)){
        
        numPerClustList[[d]][[c]] <-   numPerClust[c]
        clustSampleIndexList[[d]][[c]] <- c(count:(count+numPerClustList[[d]][[c]]-1))
        clustFeatureIndexList[[d]][[c]] <- c(1:numRows)
        count <- count+numPerClustList[[d]][[c]]
        
      }
      
    }else{

        
        count <- 1
        numPerClustList[[d]] <- list()
        
        for(c in 1:length(numPerClust)){
          
          numPerClustList[[d]][[c]] <- sample(c(minRandSize:maxRandSize),size=1)
          clustSampleIndexList[[d]][[c]] <- c(count:(count+numPerClustList[[d]][[c]] -1))
          clustFeatureIndexList[[d]][[c]] <- c(1:numRows)
          count <- count+numPerClustList[[d]][[c]]
          
        }
        

    }
  }
  
  clustFeatureIndexListOrig <- clustFeatureIndexList
  clustSampleIndexListOrig <- clustSampleIndexList
  tissueData <- createLungMatrixList()
  tissueSimData <- simulateClusterData(numSimDatasets=numSimDatasets, stddevNoise= stddevNoise,method="eigen",eigenValueMin=eigenValueMin,clustMatrixListOrig=tissueData$clustMatrixList,numRows=numRows,numPerClustList=numPerClustList)

  dataMatrixList <- list()
  clustMatrixList <- list()

  
  if(length(simType)>1){
    
    message("\nMore than one simType selected; defaulting to highQualityClust option.")
    simType <- "highQualityClust"
    
  }
  
  if(simType=="highQualityClust" || simType=="unevenSizeClust"){

      if(randNumClust){
      #need to trim down/randomly remove clusters
      clustSampleIndexList <- list()
      clustFeatureIndexList <- list()
      clustMatrixList <- list()
      dataMatrixList <- list()
      
      for(s in 1:numSimDatasets){
        
        clustSampleIndexList[[s]] <- list()
        clustMatrixList[[s]] <- list()
        clustFeatureIndexList[[s]] <- list()
        #sample from 1 to original number of clusters
        numClust <- sample(c(minRandNumClust:length(tissueSimData$clustMatrixList[[s]])),size=1)
        #message(numClust)
        clustIDs <- sample(c(1:length(tissueSimData$clustMatrixList[[s]])),size=numClust)
        #message(clustIDs)
        clustSampleIndices <- c()
        dataMatrixList[[s]] <- data.matrix(tissueSimData$dataMatrixList[[s]])
        clustMatrixList[[s]] <- list()
        clustFeaturesList[[s]] <- rownames(dataMatrixList[[s]])
        
        counter <- 1
        
        for(e in 1:length(clustIDs)){
          
          clustMatrixList[[s]][[e]] <- tissueSimData$clustMatrixList[[s]][[ clustIDs[e] ]]
          clustSampleIndexList[[s]][[e]] <- c(counter:(counter+ncol(clustMatrixList[[s]][[e]])-1))
          clustFeatureIndexList[[s]][[e]] <-  clustFeatureIndexListOrig[[s]][[  clustIDs[e] ]]
          
          if(e>1){
            
            exprMatrix <- cbind(exprMatrix,clustMatrixList[[s]][[e]])
            
          }else{
            
            exprMatrix <-  clustMatrixList[[s]][[e]]
            
          }
          
          counter <- counter + ncol( clustMatrixList[[s]][[e]] ) 
        #end of loop e
        }
        
        dataMatrixList[[s]] <- exprMatrix
        rownames(dataMatrixList[[s]]) <- clustFeaturesList[[s]]
        #end of loop s
      }
      
      }else{
      
        clustSampleIndexList <- clustSampleIndexListOrig
        clustFeatureIndexList <-  clustFeatureIndexListOrig
        dataMatrixList <- tissueSimData$dataMatrixList
        clustFeaturesList <- tissueSimData$clustFeaturesList
        clustMatrixList <- tissueSimData$clustMatrixList
        
    }
    
  }else if(simType=="mixedClustQualityClust"){
    
    #don't allow random # clusters here - controlled to know that two are noisy.
    for(s in 1:numSimDatasets){

      dataMatrixList[[s]] <- data.matrix(tissueSimData$dataMatrixList[[s]])
      clustMatrixList[[s]] <- list()
      
      clustMatrixList[[s]][[1]] <- tissueSimData$clustMatrixList[[s]][[1]]
      #mix up the two middle clusters.
      clustMatrixList[[s]][[2]] <- tissueSimData$dataMatrixList[[s]][ ,sample.int(ncol(tissueSimData$dataMatrixList[[s]]),size=numPerClust[2],replace=TRUE)]
      clustMatrixList[[s]][[3]] <- tissueSimData$dataMatrixList[[s]][ ,sample.int(ncol(tissueSimData$dataMatrixList[[s]]),size=numPerClust[3],replace=TRUE)]
      clustMatrixList[[s]][[4]] <- tissueSimData$clustMatrixList[[s]][[2]]
      #replace data matrix list too
      dataMatrixList[[s]][ , c((numPerClust[1]+1):(numPerClust[1]+numPerClust[2]))] <- clustMatrixList[[s]][[2]]
      dataMatrixList[[s]][ , c((numPerClust[1]+numPerClust[2]+1):(numPerClust[1]+numPerClust[2]+numPerClust[3]))] <- clustMatrixList[[s]][[3]]
      rownames(dataMatrixList[[s]]) <- rownames(tissueSimData$dataMatrixList[[s]])
      
    }
    
    clustFeatureIndexList <-  clustFeatureIndexListOrig
    clustFeaturesList <- tissueSimData$clustFeaturesList
    

  }else{
    
    stop("\nYou did not specify one of the correct input options for simType.")
    
  }
  

  output <- list(clustMatrixList=clustMatrixList,dataMatrixList=dataMatrixList,clustFeaturesList=clustFeaturesList,
                   clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList)
  
  return(output)
  
}

##grab real tissue data to base our simulations off of.
createLungMatrixList <- function(){
  #ref from s4vd library:
  #Bhattacharjee, A., Richards, W. G., Staunton, J., Li, C., Monti, S., Vasa, P., Ladd, C.,<br> Beheshti, J., Bueno, R., Gillette, M., Loda, M., Weber, G., Mark, E. J., Lander,<br> E. S., Wong, W., Johnson, B. E., Golub, T. R., Sugarbaker, D. J., and Meyerson,<br> M. (2001). Classification of human tissue carcinomas by mRNA expression profiling<br> reveals distinct adenocarcinoma subclasses. Proceedings of the National Academy<br> of Sciences of the United States of America. 

  data(lung200)
  #before use these in simulations, make correlations within and between each subtype.
  carcinoid_indices <- which(colnames(lung200)=="Carcinoid")
  colon_indices <- which(colnames(lung200)=="Colon")
  normal_indices <- which(colnames(lung200)=="Normal")
  smallcell_indices <- which(colnames(lung200)=="SmallCell")
  
  clustMatrixList <- list()
  clustMatrixList[[1]] <- lung200[ ,carcinoid_indices]
  clustMatrixList[[2]] <- lung200[ ,colon_indices]
  clustMatrixList[[3]] <- lung200[ ,normal_indices]
  clustMatrixList[[4]] <- lung200[ ,smallcell_indices]
  
  names(clustMatrixList) <- c("Carcinoid","Colon","Normal","SmallCell")
  output <- list(clustMatrixList=clustMatrixList,lung200=lung200)
  return(output)
  
}
########
calcCorMeanMatrix <- function(clustMatrixList,corMethod=c("pearson")){

  
  clustCorMeans <- matrix(data=NA,nrow=length(clustMatrixList),ncol=length(clustMatrixList),
                          dimnames=list(names(clustMatrixList),
                                        names(clustMatrixList)))
  
  intraClustCorMeans <- array(data=NA,dim=length(clustMatrixList),dimnames=list(names(clustMatrixList)))
  
  diag(clustCorMeans) <- 1
  
  for(r in 1:nrow(clustCorMeans)){
    
    intraClustCorMeans[r] <- mean(cor(clustMatrixList[[r]],method=corMethod))
    
    for(s in 1:nrow(clustCorMeans)){
      
      if(r != s){
        
        clustCorMeans[r,s] <- mean(cor(clustMatrixList[[r]],clustMatrixList[[s]]))
        
      }  
      
    }
    
    
  }
  
  finalDF <- data.frame(clustCorMeans,intraClustCorMeans)
  colnames(finalDF)[ncol(finalDF)] <- "intraCor"
  
  output <- list(clustCorMeans=clustCorMeans,intraClustCorMeans=intraClustCorMeans,finalDF=finalDF)
  return(output)
  
}

#############
#create a strongly correlated fake dataset:
#examples: http://stevencarlislewalker.wordpress.com/2012/06/05/simulating-random-variables-with-a-particular-correlation-structure/
#http://r-forge.r-project.org/R/?group_id=1171 . that example uses a more complicated way than rnorm to make z. 
#made a cleaned-up version in the rmv.R file in this repo.
#similar blog: http://statistical-research.com/simulating-random-multivariate-correlated-data-continuous-variables/
#why I don't use chol: oftentimes results in not a positive semi-definite matrix.

createSimGeneData <- function(sampleExpr,numRows,method=c("eigen","chol"),eigenValueMin=-.001,numSims=1,stddevNoise=0){
  
  if(stddevNoise<=0 && numSims>1){
    
    warning("\nBecause you inputted zero noise, all simulated datasets will be exactly the same.\n")
    
  }
  simDataList <- list()
  warning("\nAssumes you want to capture the correlation of the columns.\n")
  #also: must be in patient dimensions to capture this correlation structure...
  #(correlation is just a scaled version of covriance. if use correlation, may get smaller magnitudes but similar patterns.)
  #we can't add noise directly the the covariance matrix, otherwise it induces large-ish negative eigenvalues.
  covMatrix <- cov(sampleExpr)
  
  #simData <- matrix(data=NA,nrow= nrow(sampleExpr),ncol=numRows)
  z <- matrix(data=rnorm(n=ncol(covMatrix)*numRows,mean=0,sd=1),nrow=nrow(covMatrix),ncol=numRows)
  
  if(method=="eigen"){
    
    
    eigenFact <- eigen(covMatrix)
    #what if eigenvalues are negative??
    if(is.complex(eigenFact$values)){
      
      stop("Complex numbers induced in eigenvalues.")
    }
    
    #looks like minorly negative eigenvalues aren't a huge deal?
    #http://www.mathworks.com/matlabcentral/newsreader/view_thread/317549
    #sqrt the absolute value, add back the sign??
    #or just remove the sign??? most negative values are near zero anyways.
    #http://www.mathworks.com/matlabcentral/newsreader/view_thread/317549
    #so don't keep as negative? just set these small values to positive.
    #sign(eigenFact$values)*sqrt(abs(eigenFact$values))
    
    if(any(eigenFact$values<=eigenValueMin)){
      
      stop("Some eigenvalues are below the pre-set -.001 threshold.  Variables may be overly correlated?")
    }
    
    diagMatrix <- diag(sqrt(abs(eigenFact$values)), nrow=nrow(covMatrix),ncol=ncol(covMatrix))
    #code at http://r-forge.r-project.org/R/?group_id=1171 helped me confirm the correct order...otherwise 
    for(n in 1:numSims){
      
      simDataList[[n]] <- t(z)%*%diagMatrix%*%t(eigenFact$vectors)
      
      #add some noise.
      if(stddevNoise>0){
        
        #really need a lot of noise to put a dent in the correlation structure.
        simDataList[[n]] <- simDataList[[n]] + rnorm(length(simDataList[[n]]),mean=0,sd=stddevNoise)
        
      }
      
    }
    
  }else if(method=="chol"){
    
    chol_M <- chol(covMatrix)
    
    for(n in 1:numSims){
      
      simDataList[[n]] <- t(z)%*%chol_M
      
      #add some noise.
      if(stddevNoise>0){
        
        simDataList[[n]] <- simDataList[[n]] + rnorm(length(simDataList[[n]]),mean=0,sd=stddevNoise)
        
      }
      
    }
    
  }
  
  output <- list(simDataList=simDataList,realData=sampleExpr,normalMatrixList=z)
  
  return(output)
  
  
}
###############

###simulate this correlation structure. can also add noise.
#numRows=200
#numPerClust <- c(50,50,50,50)

simulateClusterData <- function(numSimDatasets=1,clustMatrixListOrig,numRows,numPerClustList,stddevNoise=0,method=c("eigen","chol"),eigenValueMin=-.001){


  maxNumSamples <- rep.int(0,length(clustMatrixListOrig))
  
  for(s in 1:numSimDatasets){

     if(length(numPerClustList[[s]])!=length(clustMatrixListOrig)){
       
       stop("\nnumPerClust and clustMatrixList must have same length.\n")
       
     }
    
    for(c in 1:length(numPerClustList[[s]])){
      
      maxNumSamples[c] <- max(maxNumSamples[c],numPerClustList[[s]][[c]]) 
      
    }
    
 
  } 
 
    numSamples <- sum(maxNumSamples)
    tempExprMatrix <- matrix(data=NA,nrow=numRows,ncol=numSamples)
    clustMatrixList <- list()
    count <- 1
    #clustMatrixList: the original real tissue cluster data list. so only 1 dataset, not numSims datasets, in here.
    for(m in 1:length(clustMatrixListOrig)){

      exprMatrix <- matrix(data=NA,ncol=maxNumSamples[m],nrow=nrow(clustMatrixListOrig[[m]]))
      
      #NOTE: this will just repeat samples a lot if
      #simulated cluster sizes are much larger than true ones - not create truly distinct ones.
      #But for simulations, we want fairly clear baseline signals.
      exprMatrix[ ,c(1:maxNumSamples[m])] <- clustMatrixListOrig[[m]][, sample(c(1:ncol(clustMatrixListOrig[[m]])),
        size=maxNumSamples[m],replace=TRUE)]

      tempExprMatrix[ , c(count:(count+(maxNumSamples[m]-1)))] <- exprMatrix
      count <- count + maxNumSamples[m]
      #end of loop m.
    }
  
    
    dataMatrixList <- createSimGeneData(numSims=numSimDatasets,sampleExpr=tempExprMatrix,stddevNoise=stddevNoise,numRows=numRows,method=method,eigenValueMin=eigenValueMin)$simDataList
    
    for(s in 1:length(dataMatrixList)){
      cNames <- c()
      #NOW pair down using specific cluster sizes.
      clustMatrixList[[s]] <- list()
      count <- 1
      for(d in 1:length(clustMatrixListOrig)){
      
        #no need for replacement here: have at least the maximum number of samples required for any simulation iteration.
        #also want zero noise scenario to be truly similar across datasets.
        totalClustIndices <- c(count:(count+(maxNumSamples[d]-1)))
        subClustIndices <- sample(totalClustIndices,size=numPerClustList[[s]][[d]],replace=FALSE)
        clustMatrixList[[s]][[d]] <- dataMatrixList[[s]][ ,subClustIndices, drop=FALSE]
        
        cNamesTmp <- paste0(rep.int(colnames(clustMatrixListOrig[[d]])[1],times=numPerClustList[[s]][[d]]),"_",c(1:numPerClustList[[s]][[d]]))
        cNames <- append(cNames,cNamesTmp)
        colnames(clustMatrixList[[s]][[d]]) <- cNamesTmp
        names(clustMatrixList[[s]])[d] <- names(clustMatrixList)[d]
        count <- count + maxNumSamples[d]
    
        if(d>1){
          
          tmp <- cbind(tmp,clustMatrixList[[s]][[d]])
          
        }else{
          
          tmp <- clustMatrixList[[s]][[d]]
          
        }

      }
      
      names(dataMatrixList)[s] <- s

      #re-assign with matrix with specific cluster sizes.
      dataMatrixList[[s]] <- tmp
      colnames(dataMatrixList[[s]]) <- cNames
      #create fake gene names
      rownames(dataMatrixList[[s]]) <- paste0("gene_",c(1:nrow(dataMatrixList[[s]])))
      #end of loop s.
      
    }
    
  output <- list(dataMatrixList=dataMatrixList,clustMatrixList=clustMatrixList)
  
  return(output)
  
}

#####

#assumes clusters are in the same order each time.
createTrueEdges <- function(clustIDVector){
  
  #If you have N nodes, there are N - 1 directed edges than can lead from it (going to every other node). 
  #Therefore, the maximum number of edges is N * (N - 1)
  #we do not allow edges to go from a node to itself.
  
  numTotalEdges <- 0
  
  for(c in 1:length(unique(clustIDVector))){
    
    numWithID <- length(which(clustIDVector==unique(clustIDVector)[c]))
    numTotalEdges <- numTotalEdges + numWithID*(numWithID-1)
    
  }
  #we are doing unidirectional, no bidirectional, edges.
  numTotalEdges <- numTotalEdges/2
  trueEdges <- matrix(data=NA,ncol=2,nrow=numTotalEdges);
  counter <- 1;
  
  if(nrow(trueEdges)>0){
    
  for(r in 1:length(clustIDVector)){
    
    connectNodes <- which(clustIDVector==clustIDVector[r])
    #don't take noes that are themselves
    connectNodes <- connectNodes[which(connectNodes!=r)]
    
    if(length(connectNodes)>0){
      
    PASS <- FALSE
    
    for(n in 1:length(connectNodes)){
      
      #already counted? the nodes can have flipped column order (not directed nodes here so this is equivalent.)
      if(length(which(trueEdges[,1]==r)) == 0 &&  length(which(trueEdges[,2]==r))==0 ){
        
       PASS <- TRUE
        
      }else if( (length(which(trueEdges[,1]==r)) != 0) ){
        
        if(all(trueEdges[which(trueEdges[,1]==r) ,2] != connectNodes[n]) ){
          
          #now check in other direction
          if( (length(which(trueEdges[,2]==r)) != 0) ){
            
            if(all(trueEdges[which(trueEdges[,2]==r) ,1] != connectNodes[n]) ){
              
              PASS <- TRUE
              
            }
            
            #none in the second column.
          }else{
            
            PASS <- TRUE
          }
          
        }
        
      }else if((length(which(trueEdges[,2]==r)) != 0)){
        
        if(all(trueEdges[which(trueEdges[,2]==r) ,1] != connectNodes[n]) ){
          
          PASS <- TRUE
          
        }
        
      }else{
        
        stop("Bug: already checked to make sure at least one column had a hit for this node. Shouldn't reach this else statement.")
      }
        
      
      if(PASS){
        
        #this node has not appeared in the edge matrix before.
        trueEdges[counter, 1] <- r
        trueEdges[counter, 2] <- connectNodes[n]
        counter <- counter + 1
        
      }
        
      }

    }
    

    }
    
    }
    

  if(!all(clustIDVector[trueEdges[,1]]==clustIDVector[trueEdges[,2]])){
    
    stop("Not correctly matching up tissue IDS")
    
  }
  
  if((counter-1) != nrow(trueEdges)){
    
    stop("Bug: not looping through all of true edges correctly.")
    
  }
  
  return(trueEdges);
  
}
###########
compute_edge_ROC_metrics <- function(trueEdges,predEdges,numTotalClusters,numNonExistEdges){
  
  trueEdgesFound <- matrix(data=NA,ncol=2,nrow=1);
  trueEdgesMissed <- matrix(data=NA,ncol=2,nrow=1);
  
  for(r in 1:nrow(trueEdges)){
    
    rowIndices <- which(predEdges[,1]==trueEdges[r,1]);
    
    if(length(rowIndices)>0){
      
      if(any(predEdges[rowIndices,2]==trueEdges[r,2])){
        
        #num_truePos <- num_truePos + 1;
        #remove this edge: will help sum up FP later
        if(length(which(predEdges[rowIndices,2]==trueEdges[r,2]))>1){
          
          stop("indexing error.");
          
        }
        #cat("\n",predEdges[rowIndices[which(predEdges[rowIndices,2]==trueEdges[r,2])],],"\n")
        predEdges <- predEdges[-rowIndices[which(predEdges[rowIndices,2]==trueEdges[r,2])], , drop=FALSE];
        
        trueEdgesFound <- rbind(trueEdgesFound,trueEdges[r,c(1:2),drop=FALSE]);
        
        #found correct edge - move on to next row
        next;
        #no row index found.
      }
      #no row index found. 
    }
    # but flip node order, check again
    
    rowIndices <- which(predEdges[,2]==trueEdges[r,1]);
    
    if(length(rowIndices)>0){
      
      if(any(predEdges[rowIndices,1]==trueEdges[r,1])){
        
        #num_truePos <- num_truePos + 1;
        #remove this edge: will help sum up FP later
        if(length(which(predEdges[rowIndices,2]==trueEdges[r,2]))>1){
          
          stop("indexing error.");
          
        }
        #cat("\n",predEdges[rowIndices[which(predEdges[rowIndices,2]==trueEdges[r,2])],],"\n")
        predEdges <- predEdges[-rowIndices[which(predEdges[rowIndices,2]==trueEdges[r,2])], ,drop=FALSE];
        trueEdgesFound <- rbind(trueEdgesFound,trueEdges[r,c(1:2),drop=FALSE]);
        
        #found a true edge - move on to next loop.
        next;
      }
      
    }
    
    
    #didn't find an edge when there was a true edge.
    trueEdgesMissed <- rbind(trueEdgesMissed,trueEdges[r,c(1:2),drop=FALSE]);
    #move on to next r index in true edges.
  }
  
  trueEdgesFound <- trueEdgesFound[-1, ,drop=FALSE];
  trueEdgesMissed <- trueEdgesMissed[-1, ,drop=FALSE];
  num_truePos <- nrow(trueEdgesFound);
  #any predicted left that weren't true edges?
  num_falsePos <- nrow(predEdges);
  #how many true edges did we miss?
  num_falseNeg <- nrow(trueEdgesMissed);
  
  #total number of possible edges is (N * (N-1)) / 2
  #http://stackoverflow.com/questions/5058406/what-is-the-maximum-number-of-edges-in-a-directed-graph-with-n-nodes
  #/2 becuase edges are not directed. (one edge possible between each cluster pair.)
  totalPossEdges <- (numTotalClusters*(numTotalClusters-1))/2

  #NOTE: it's not REALLY this traditional calculation. Must "minus out" edges that are from
  #the same dataset, as these never could be connected to each other using Coincide.
  #calculated this outside of the function when had the information what clusters came from what datasets.
  totalPossEdges <- totalPossEdges- numNonExistEdges
  
  #must also minus out the false positives, and the false negatives - these should not be counted as true negatives.
  num_trueNeg <- totalPossEdges-num_truePos-num_falsePos-num_falseNeg;
  
  if(num_falseNeg+num_truePos != nrow(trueEdges)){
    
    stop("\nUnit test 1: Calculating FN or TP wrong.\n");
    
  }
  
  #num_falsePos: want to minus this:  true edges found + true edges missed + "rest" of edges are the num_trueNeg
  if( (num_falseNeg+num_truePos+num_trueNeg+num_falsePos) != totalPossEdges){
    
    stop("\nUnit test 2: Calculating FP or TN wrong as FP,TN,FN and TP sum does not equal total possible edges.\n");
    
  }

  #TPR is sensitivity, specificity is FPR
  #we mostly care about PPV and TPR. FPR will almost alwways be low.
  if(num_truePos >0){
    
    PPV <- num_truePos/(num_truePos+num_falsePos);
    
  }else{
    
    PPV <- 0;
  }
  
  FPR <- num_falsePos/(num_falsePos+num_trueNeg);
  TPR <- num_truePos/(num_truePos + num_falseNeg);
  output <- list(PPV=PPV,TPR=TPR,FPR=FPR,FP=num_falsePos,TP=num_truePos,FN=num_falseNeg,TN=num_trueNeg,
                 trueEdgesFound=trueEdgesFound,trueEdgesMissed=trueEdgesMissed);
  
  return(output);
  
}

runTissueClusterSimROC <- function(saveDir="./",numSimDatasets=10,
                       eigenValueMin = -.001,simType=c("highQualityClust","mixedClustQualityClust","unevenSizeClust"),
                       noiseVector = seq(from=0,to=2.5,by=.2),numPerClust = c(50,50,50,50),
                       numWrapperSims=100,numSims=500,fractFeatIntersectThresh=.7,numFeatIntersectThresh=190,
                       clustSizeThresh=0,clustSizeFractThresh=0,minTrueSimilThresh=.4,
                       maxTrueSimilThresh=Inf,edgeMethod=c("pearson"),centroidMethod=c("mean"),
                       indEdgePvalueThresh=.05,meanEdgePairPvalueThresh=.1,
                      minFractNN=.8,minRandNumClust=2,randNumClust=FALSE,minRandSize=1,maxRandSize=100
){
  

  numWrapperSimsOrig <- numWrapperSims
  outputFile <- paste0(saveDir,"/simsMessages.txt")
  
  ROC <- list();
  commMedianMaxTissueType <- list();
  commMeanMaxTissueType <- list();
  
  
  #just run n=1, i.e. zero data, once. do by hand.
  for(n in 1:length(noiseVector)){
    cat("\nLoop ",n,"Noise level: ",noiseVector[n],"\n")
    #only run 1 simulation iteration for noise level o.
    
    if(noiseVector[n]==0){
      
      numWrapperSims <- 1;
      
    }else{
      
      numWrapperSims <- numWrapperSimsOrig;
      
    }
    ROC[[n]] <- list();
    commMedianMaxTissueType[[n]] <- list();
    commMeanMaxTissueType[[n]] <- list();
    message("numWrapperSims ",numWrapperSims)
    #    message("numWrapperSims ",numWrapperSimsOrig)
    for(t in c(1:numWrapperSims)){
      
      #don't save these: just recalc each time. this will save the last iteration.
      tissueSimData <- createTissueSimDatasets(numSimDatasets=numSimDatasets,
                                           eigenValueMin =eigenValueMin,simType=simType,
                                           numPerClust = numPerClust,
                                           stddevNoise=noiseVector[n],numRows=200,
                                           minRandSize=minRandSize,maxRandSize=maxRandSize,
                                           minRandNumClust=minRandNumClust,randNumClust=randNumClust);
      
      
      dataMatrixList <- tissueSimData$dataMatrixList
      clustSampleIndexList <- tissueSimData$clustSampleIndexList
      clustFeatureIndexList <- tissueSimData$clustFeatureIndexList
      
      numTotalClusters <- 0
    
      clustNameVector <- c()
      
      for(s in 1:length(clustSampleIndexList)){
        
      for(c in 1:length(clustSampleIndexList[[s]])){
        
        numTotalClusters <- numTotalClusters+1
        #tissue names in same cluster will all be identical before _ tag.
        #drop=FALSE: in case cluster is a size of 1.
        clustNameVector[numTotalClusters] <- strsplit2(colnames(dataMatrixList[[s]][,clustSampleIndexList[[s]][[c]], drop=FALSE])[1],"_")[,1]
      }
        
      }  
      #STEP: create your true edge list
      trueEdgeMatrix <- createTrueEdges(clustIDVector=clustNameVector);

      message("running CoINcIDE")
      adjMatrixOut <- computeAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                         edgeMethod=edgeMethod,centroidMethod=centroidMethod,
                                         numSims=numSims,
                                         outputFile=outputFile)
      
      
      predEdges <- assignFinalEdges(computeTrueSimilOutput=adjMatrixOut$computeTrueSimilOutput,pvalueMatrix=adjMatrixOut$pvalueMatrix,indEdgePvalueThresh=indEdgePvalueThresh,
                                    meanEdgePairPvalueThresh=meanEdgePairPvalueThresh,
                                    minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                    fractFeatIntersectThresh=fractFeatIntersectThresh,numFeatIntersectThresh=numFeatIntersectThresh ,
                                    clustSizeThresh=clustSizeThresh, clustSizeFractThresh= clustSizeFractThresh,saveDir=saveDir,fileTag="simil_edges_",
                                    minFractNN = minFractNN
      )

        #need to figure out edges that could never be drawn because are from same dataset.
      #can't have an edge drawn to anywhere in these datasets
      #note: this code is different than the null dataset calculations, as here,
      #we want to count true negative edges here. so for ROC calculations, "possible" edges
      #will be a larger number, as it will count true negatives
      #BUT an edge could never be drawn between clusters from the same dataset,
      #so don't want to include these numbers in our true negative denominator.
        numNonExistEdges <- 0
        for(c in 1:length(clustSampleIndexList)){
          #each length(clustSampleIndexList[[c]]: # of clusters in dataset C
            #total number of possible edges (here, within a dataset) for a set of N nodes is (N * (N-1)) / 2
            #http://stackoverflow.com/questions/5058406/what-is-the-maximum-number-of-edges-in-a-directed-graph-with-n-nodes
            #/2 becuase edges are not directed. (one edge possible between each cluster pair.)
            numNonExistEdges <- numNonExistEdges + (length(clustSampleIndexList[[c]])*(length(clustSampleIndexList[[c]])-1))/2
          
        }
      
        ROC[[n]][[t]]   <-    compute_edge_ROC_metrics(trueEdges=trueEdgeMatrix,predEdges=predEdges$filterEdgeOutput$edgeMatrix,numTotalClusters=numTotalClusters,numNonExistEdges=numNonExistEdges);

        commInfo <- findCommunities(edgeMatrix=predEdges$filterEdgeOutput$edgeMatrix,edgeWeightMatrix=predEdges$filterEdgeOutput$edgeWeightMatrix,
                                    clustIndexMatrix=adjMatrixOut$clustIndexMatrix,fileTag="sims",
                                    saveDir=saveDir,minNumUniqueStudiesPerCommunity=2,experimentName="sims",
                                    commMethod="edgeBetween",
                                    makePlots=FALSE,saveGraphData=FALSE,plotToScreen=FALSE, findCommWithWeights=TRUE, plotSimilEdgeWeight = FALSE,
                                    fractEdgesInVsOutComm=0, fractEdgesInVsOutEdge=0)

        if(commInfo$numCommunities>0){        
        
          sampleClustCommKey <- returnSampleMemberMatrix(clustSampleIndexList=clustSampleIndexList,
                                                  dataMatrixList=dataMatrixList,communityInfo=commInfo)$sampleClustCommKey
          
          
        commClustTable <- table(sampleClustCommKey$community,sampleClustCommKey$globalClustNum)
        commTissueMakeUp <- list()
        
        for(c in 1:commInfo$numCommunities){
          
          tmp <- c()
          for(r in 1:ncol(commClustTable)){
            
          tmp <- append(tmp,rep.int(clustNameVector[r],commClustTable[c,r]))
            
          }
          
          #fraction of most prevalent tissue type.
          commTissueMakeUp[[c]] <- max(table(tmp))/length(tmp)
          
        }
        
        #take mean across all communities
        commMeanMaxTissueType[[n]][[t]] <- mean(unlist(commTissueMakeUp))
        commMedianMaxTissueType[[n]][[t]] <- median(unlist(commTissueMakeUp))
        
        }else{
          #non communities - technically all match then? use NA as placeholder to symbolize this special case.
          commMeanMaxTissueType[[n]][[t]] <- NA
          commMedianMaxTissueType[[n]][[t]] <- NA
          
        }

      #loop of t simulations
    }

    #loop n
  }


  #convert to matrices
  
  ROC_matrixFull <- do.call(rbind,lapply(ROC,FUN=function(ROC_unit){
    
    #this will make 1 matrix for 1 noise level, with each row a simulation.                     
    oneNoiseLevelData <- do.call(rbind,lapply(ROC_unit,FUN=function(data){
      
      ROC_stats <- c(data$PPV,data$TPR,data$FPR,data$FP,data$TP,data$FN,data$TN);
      return(ROC_stats);
    }
    
    )
    
    )
    
    #take the mean across all simulations for this noise level.
    #if only one: take data.matrix
    oneNoiseLevelData <- colMeans(data.matrix(oneNoiseLevelData));
    return(oneNoiseLevelData);
    
  }
  
  )
  
  );
  
  ROC_matrixFull <- cbind(noiseVector,ROC_matrixFull)
  colnames(ROC_matrixFull) <- c("sd_noise","PPV","TPR","FPR","FP","TP","FN","TN");
  rownames(ROC_matrixFull) <- noiseVector;
  
  ROC_SDmatrixFull <- do.call(rbind,lapply(ROC,FUN=function(ROC_unit){
    
    #this will make 1 matrix for 1 noise level, with each row a simulation.                     
    oneNoiseLevelData <- do.call(rbind,lapply(ROC_unit,FUN=function(data){
      
      ROC_stats <- c(data$PPV,data$TPR,data$FPR,data$FP,data$TP,data$FN,data$TN);
      return(ROC_stats);
    }
    
    )
    
    )
    
    #take standard deviation across all iterations for this noise level.
    oneNoiseLevelData <- colSds(data.matrix(oneNoiseLevelData));
    return(oneNoiseLevelData);
    
  }
  
  )
  
  );
  
  ROC_SDmatrixFull <- cbind(noiseVector, ROC_SDmatrixFull)
  colnames(ROC_SDmatrixFull) <- c("sd_noise","sd_PPV","sd_TPR","sd_FPR","sd_FP","sd_TP","sd_FN","sd_TN");
  rownames(ROC_SDmatrixFull) <- noiseVector;
  
  #a simple plot
  df <- data.frame(ROC_matrixFull);
  #theme(plot.title = element_text(size = rel(2)))
  TPR_plot <- ggplot(data =df,aes(x=sd_noise,y=TPR))+geom_line()+  labs(title = "TPR for simulations \nwith increasing noise.",
                                                                        y="TPR",x="sd noise")+
    theme(panel.background = element_rect(fill='white', colour='black')) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                 plot.title = element_text(size = rel(2)));
  commMedianMaxTissueTypeList <- commMedianMaxTissueType
  commMeanMaxTissueTypeList <- commMeanMaxTissueType
  #aggregate across all t to return for each noise level.
  commMedianMaxTissueType <- lapply(commMedianMaxTissueType,FUN=function(noiseUnit){
    
    #NA: no communities detected.
    return(median(unlist(noiseUnit),na.rm=TRUE))
    
    
  })
  commMedianMaxTissueType <- unlist(commMedianMaxTissueType)
  
  commMeanMaxTissueType <- lapply(commMeanMaxTissueType,FUN=function(noiseUnit){
    
    return(mean(unlist(noiseUnit),na.rm=TRUE))
    
    
  })
  
  commMeanMaxTissueType <- unlist( commMeanMaxTissueType)
  output <- list(TPR_plot=TPR_plot,ROC=ROC, ROC_matrixFull=ROC_matrixFull,
                 ROC_SDmatrixFull= ROC_SDmatrixFull,  commMedianMaxTissueTypeList=commMedianMaxTissueTypeList,
                 commMeanMaxTissueTypeList=commMeanMaxTissueTypeList,commMedianMaxTissueType=commMedianMaxTissueType,
                 commMeanMaxTissueType=commMeanMaxTissueType
  );
  
  return(output);

}
