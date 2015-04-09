
ISMB_lungSimROC <- function(sourceDir="/home/kplaney/gitRepos//IGP_network/igp_network/",
                            saveDir = "/home/kplaney/ISMB/lung_sims/",numSimDatasets=10,tmpSaveName = paste0("tmp",rnorm(1)),
                            eigenValueMin = -.001,simType=c("highQualityClust","mixedClustQualityClust","unevenSizeClust"),simMetric=c("centroid","correlation"),
                                                            noiseVector = seq(from=0,to=2.5,by=.2),numPerClust = c(50,50,50,50),
                                                            numWrapperSims=10,numCores=2,numSims=1000
                          ){
                              
                              
                              if(length(numPerClust) != 4){
                                
                                stop("\nMust pick a cluster size for all 4 tissue types:\nyour numPerClust variable vector was not a length of 4.")
                              }  
                              
                              #library("doParallel");
                              #registerDoParallel(cores=numCores);
                              
                              setwd(sourceDir);
                              source("clust_robust.R");

                              ROC <- list();
                              
                              #fabricate your true edge list
                              trueEdgeMatrix <- createTrueEdges(numRowClust=1,numColClust=4,numSimDatasets=numSimDatasets);
                              
                              if(simType=="unevenSizeClust"){
                                #remove any edges that containi last 2 nodes -they aren't in our simulation.
                                removeNode1 <- 4*numSimDatasets;
                                removeNode2 <- 4*numSimDatasets-1;
                                rowRemove1_1 <- which(trueEdgeMatrix[,1]==removeNode1);
                                if(length(rowRemove1_1 >1)){
                                  trueEdgeMatrix <- trueEdgeMatrix [-rowRemove1_1, ,drop=FALSE];
                                }
                                rowRemove1_2 <- which(trueEdgeMatrix[,2]==removeNode1);
                                if(length(rowRemove1_2 >1)){
                                trueEdgeMatrix <- trueEdgeMatrix [-rowRemove1_2, ,drop=FALSE];
                                }
                                rowRemove2_1 <- which(trueEdgeMatrix[,1]==removeNode2);
                                if(length(rowRemove2_1 >1)){
                                trueEdgeMatrix <- trueEdgeMatrix [-rowRemove2_1, ,drop=FALSE];
                                }
                                rowRemove2_2 <- which(trueEdgeMatrix[,2]==removeNode2);
                                if(length(rowRemove2_2 >1)){
                                trueEdgeMatrix <- trueEdgeMatrix [-rowRemove2_2, ,drop=FALSE];
                                }
                                
                                numTotalClusters <- 4*numSimDatasets-2;
                                
                              }else{
                                
                                numTotalClusters <- 4*numSimDatasets;
                                
                              }
                              
                              numWrapperSimsOrig <- numWrapperSims;
                              #just run n=1, i.e. zero data, once. do by hand.
                              for(n in 1:length(noiseVector)){
                                cat("\nLoop ",n,"Noise level: ",noiseVector[n],"\n")
                                #only run 1 simulation iteration for noise level o.
                                
                                if(n==1){
                                  
                                  numWrapperSims <- 1;
                                  
                                }else if(n>1){
                                  
                                  numWrapperSims <- numWrapperSimsOrig;
                                  
                                }
                                ROC[[n]] <- list();
                                
                                for(t in c(1:numWrapperSims)){
                                
                                #will automatically return a list.
                                  #bc save datasets here...Rparallel doesn't work!
                                #ROC_tmp <-   foreach(c(1:numWrapperSims)) %dopar% {
                                  
                                  #must transpose for my clust robust methods...
                                  #don't save these: just recalc each time. this will save the last iteration.
                                  lungSimData <- createLungSimDatasets(numSimDatasets=numSimDatasets,
                                                                       eigenValueMin =eigenValueMin,simType=simType,
                                                                       numPerClust = numPerClust,
                                                                       stddevNoise=noiseVector[n],numRows=200);
                                  
                                  if(simMetric=="centroid"){
                                    
                                    clustRobust_IGP <- findSimilarClusters(clustMatrixList=lungSimData$clustMatrixList,dataMatrixList=lungSimData$dataMatrixList,minClustSampleSize=10,resampleWithClustSamples=TRUE,compareMetric=c("correlation"),
                                                                           normFeatureWise=c("none"),method=c("NN_distribution"),numSims=numSims,
                                                                           minFractFeatureIntersect=.05,minNumFeatureIntersect=5,fract_nullMax=.1,
                                                                           centroidType=c("mean"),outputFile="/home/kplaney/ISMB/lung_sims/tmp.txt",
                                                                           nullClustTag =tmpSaveName ,sourceDir="/home/kplaney/gitRepos/IGP_network/igp_network/",
                                                                           trueDiffMin=.3,intraSimilMax=.3,simil_pvalueMax=1,tempSaveDir=saveDir,
                                                                           runIntraStabilityTests=FALSE,compareIntraInterDistributions=FALSE
                                    );
                                    
                                    edge_IGPList <- computeEdgeMatrix(clustRobust_output=clustRobust_IGP,restrictEdgesByStudy=TRUE,method=c("NN_distribution"),
                                                                      compareMetric=c("correlation"),edge_trueCorr_thresh=.3,edge_fractFeatIntersect_thresh=.05,
                                                                      edge_numFeatIntersect_thresh=20,edge_IGP_thresh=.2,sourceDir="/home/kplaney/gitRepos/IGP_network/igp_network/",
                                                                      edge_simil_overNull_pvalue_thresh=.05,clustMethodName="tmp")$edgeResults$edgeMatrix;
                                    
                                    #return this.
                                    #ROC[[n]][[t]]  <-
                                    ROC[[n]][[t]] <-     compute_edge_ROC_metrics(trueEdges=trueEdgeMatrix,predEdges=edge_IGPList,numTotalClusters=numTotalClusters);
                                    
                                  }else if(simMetric=="correlation"){
                                    
                                    clustRobust_corr <- findSimilarClusters(clustMatrixList=lungSimData$clustMatrixList,dataMatrixList=lungSimData$dataMatrixList,minClustSampleSize=10,resampleWithClustSamples=TRUE,compareMetric=c("correlation"),
                                                                            normFeatureWise=c("none"),method=c("diff_distribution"),numSims=numSims,
                                                                            minFractFeatureIntersect=.05,minNumFeatureIntersect=5,fract_nullMax=.1,
                                                                            centroidType=c("mean"),outputFile="/home/kplaney/ISMB/lung_sims/tmp.txt",
                                                                            nullClustTag =tmpSaveName ,sourceDir="/home/kplaney/gitRepos/IGP_network/igp_network/",
                                                                            trueDiffMin=.3,intraSimilMax=.3,simil_pvalueMax=1,tempSaveDir=saveDir,
                                                                            runIntraStabilityTests=FALSE,compareIntraInterDistributions=FALSE
                                    );
                                    
                                    
                                    edge_corrList <- computeEdgeMatrix(clustRobust_output=clustRobust_corr,restrictEdgesByStudy=TRUE,method=c("diff_distribution"),
                                                                       compareMetric=c("correlation"),edge_trueCorr_thresh=.3,edge_fractFeatIntersect_thresh=.05,
                                                                       edge_numFeatIntersect_thresh=20,edge_IGP_thresh=.2,sourceDir="/home/kplaney/gitRepos/IGP_network/igp_network/",
                                                                       edge_simil_overNull_pvalue_thresh=.05,clustMethodName="tmp",saveDir="/home/kplaney")$edgeResults$edgeMatrix;
                                    
                                    
                                    #return this.
                                    #  ROC[[n]][[t]] <-
                                    ROC[[n]][[t]]   <-    compute_edge_ROC_metrics(trueEdges=trueEdgeMatrix,predEdges=edge_corrList,numTotalClusters=numTotalClusters);
                                    
                                  }
                                  
                                  
                                  #loop of simulations/foreach
                                }
                                #loop n
                                #if did R parallel.
                                #ROC[[n]] <- ROC_tmp
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
                              
                              #a simple plot
                              library("ggplot2");
                              #only take first 12...
                              df <- data.frame(ROC_matrixFull);
                              #theme(plot.title = element_text(size = rel(2)))
                              TPR_plot <- ggplot(data =df,aes(x=sd_noise,y=TPR))+geom_line()+  labs(title = "TPR for simulations with increasing noise.",
                                                                                                         y="TPR",x="sd noise")+
                                theme(panel.background = element_rect(fill='white', colour='black')) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                             plot.title = element_text(size = rel(2)));
                              
                              
                              output <- list(TPR_plot=TPR_plot,ROC=ROC, ROC_matrixFull=ROC_matrixFull);
                              
                              return(output);
                              
}
######
createLungSimDatasets <-  function(numSimDatasets=10,
                                   eigenValueMin = -.001,simType=c("highQualityClust","mixedClustQualityClust","unevenSizeClust"),
                                   numPerClust = c(50,50,50,50),
                                   stddevNoise=0,numRows=200){
  
  clustSampleIndexList <- list()
  clustFeatureIndexList <- list()
  clustFeaturesList <- list()
  
  for(d in 1:numSimDatasets){
    
    clustSampleIndexList[[d]] <- list()
    clustFeatureIndexList[[d]] <- list()

    count <- 1
    for(c in 1:length(numPerClust)){
      
      clustSampleIndexList[[d]][[c]] <- c(count:(count+numPerClust[c]-1))
      clustFeatureIndexList[[d]][[c]] <- c(1:numRows)
      count <- count+numPerClust[c]
      
    }
    
  }
  
  
  lungData <- createLungMatrixList()
  lungSimData <- simulateClusterData(numSimDatasets=numSimDatasets, stddevNoise= stddevNoise,method="eigen",eigenValueMin=eigenValueMin,clustMatrixList=lungData$clustMatrixList,numRows=numRows,numPerClust=numPerClust)
  
  dataMatrixList <- list()
  clustMatrixList <- list()
  
  if(length(numPerClust)!= 4){
    
    stop("Must specify the number of samples for all 4 tissue types.")
    
  }
  
  if(length(simType)>1){
    
    message("\nMore than one simType selected; defaulting to highQualityClust option.")
    simType <- "highQualityClust"
    
  }
  
  if(simType=="highQualityClust"){
    
    
    #construct your data matrices.
    for(s in 1:numSimDatasets){
      
      dataMatrixList[[s]] <- data.matrix(lungSimData$simMatrixList[[s]])
      clustMatrixList[[s]] <- list()
      clustFeaturesList[[s]] <- rownames(dataMatrixList[[s]])
    
      #COME BACK: simplify?
      for(e in 1:length(lungSimData$simClustList[[s]])){
        
        clustMatrixList[[s]][[e]] <- lungSimData$simClustList[[s]][[e]]
        
      }
      
    }
    
  }else if(simType=="mixedClustQualityClust"){
    
    
    for(s in 1:numSimDatasets){

      dataMatrixList[[s]] <- data.matrix(lungSimData$simMatrixList[[s]])
      clustMatrixList[[s]] <- list()
      
      clustMatrixList[[s]][[1]] <- lungSimData$simClustList[[s]][[1]]
      #mix up the two middle clusters.
      clustMatrixList[[s]][[2]] <- lungSimData$simMatrixList[[s]][ ,sample.int(ncol(lungSimData$simMatrixList[[s]]),size=numPerClust[2],replace=TRUE)]
      clustMatrixList[[s]][[3]] <- lungSimData$simMatrixList[[s]][ ,sample.int(ncol(lungSimData$simMatrixList[[s]]),size=numPerClust[3],replace=TRUE)]
      clustMatrixList[[s]][[4]] <- lungSimData$simClustList[[s]][[2]]
      #replace data matrix list too
      dataMatrixList[[s]][ , c((numPerClust[1]+1):(numPerClust[1]+numPerClust[2]))] <- clustMatrixList[[s]][[2]]
      dataMatrixList[[s]][ , c((numPerClust[1]+numPerClust[2]+1):(numPerClust[1]+numPerClust[2]+numPerClust[3]))] <- clustMatrixList[[s]][[3]]
      
      
    }
    
  }else if(simType=="unevenSizeClust"){
    
    for(s in 1:numSimDatasets){
      
      if(s==numSimDatasets){
        #must transpose for my clust robust methods..
        clustMatrixList[[s]] <- list()
        #only resample from first cluster - make a large one.
        clustMatrixList[[s]][[1]] <- lungSimData$simClustList[[s]][[1]][ ,sample.int(numPerClust[1],size=100,replace=TRUE)]
        #make a small cluster
        clustMatrixList[[s]][[2]] <- lungSimData$simClustList[[s]][[2]][ ,sample.int(numPerClust[2],size=10,replace=TRUE)]
        #replace data matrix list too
        dataMatrixList[[s]] <- data.matrix(cbind(clustMatrixList[[s]][[1]],clustMatrixList[[s]][[2]]))
        
        #just make the rest normal clusters.
      }else{
        
        dataMatrixList[[s]] <- data.matrix(lungSimData$simMatrixList[[s]])
        clustMatrixList[[s]] <- list()
        
        for(e in 1:length(lungSimData$simClustList[[s]])){
          
          clustMatrixList[[s]][[e]] <- lungSimData$simClustList[[s]][[e]]
          
        }
        
      }   
      
    }
    
  }else{
    
    stop("\nYou did not specify one of the correct input options for simType.")
    
  }
  
  output <- list(clustMatrixList=clustMatrixList,dataMatrixList=dataMatrixList,clustFeaturesList=clustFeaturesList,
                   clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList)
  
  return(output)
  
}

##grab real lung data to base our simulations off of.
createLungMatrixList <- function(){
  #ref:
  #Bhattacharjee, A., Richards, W. G., Staunton, J., Li, C., Monti, S., Vasa, P., Ladd, C.,<br> Beheshti, J., Bueno, R., Gillette, M., Loda, M., Weber, G., Mark, E. J., Lander,<br> E. S., Wong, W., Johnson, B. E., Golub, T. R., Sugarbaker, D. J., and Meyerson,<br> M. (2001). Classification of human lung carcinomas by mRNA expression profiling<br> reveals distinct adenocarcinoma subclasses. Proceedings of the National Academy<br> of Sciences of the United States of America. 
  library("s4vd")
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
calcCorMeanMatrix <- function(clustMatrixList,fileSaveName="/home/kplaney/ISMB/lung_sims/lungCorMatrix.txt"){
  
  warning("Assumes your clustMatrixList is in this order: Carcinoid,Colon,Normal,SmallCell\n")
  
  if(!all(names(clustMatrixList) == c("Carcinoid","Colon","Normal","SmallCell"))){
    
    warning("Your clustMatrixList is NOT in this order: Carcinoid,Colon,Normal,SmallCell\n")
    
  }
  
  clustCorMeans <- matrix(data=NA,nrow=4,ncol=4,
                          dimnames=list(names(clustMatrixList),
                                        names(clustMatrixList)))
  
  intraClustCorMeans <- array(data=NA,dim=4,dimnames=list(names(clustMatrixList)))
  
  diag(clustCorMeans) <- 1
  
  for(r in 1:nrow(clustCorMeans)){
    
    intraClustCorMeans[r] <- mean(cor(clustMatrixList[[r]]))
    
    for(s in 1:nrow(clustCorMeans)){
      
      if(r != s){
        
        clustCorMeans[r,s] <- mean(cor(clustMatrixList[[r]],clustMatrixList[[s]]))
        
      }  
      
    }
    
    
  }
  
  finalDF <- data.frame(clustCorMeans,intraClustCorMeans)
  colnames(finalDF)[ncol(finalDF)] <- "intraCor"
  write.table(finalDF,file=fileSaveName)
  
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

simulateClusterData <- function(numSimDatasets=1,clustMatrixList,numRows,numPerClust,stddevNoise=0,method=c("eigen","chol"),eigenValueMin=-.001){
  
  if(length(numPerClust)!=length(clustMatrixList)){
    
    stop("\nnumPerClust and clustMatrixList must have same length.\n")
    
  }
  simMatrix <- matrix(data=NA,nrow=numRows,ncol=sum(numPerClust))
  tempExprMatrix <- matrix(data=NA,nrow=numRows,ncol=sum(numPerClust))
  simClustList <- list()
  count <- 1
  cNames <- c()
  for(m in 1:length(clustMatrixList)){
    
    #the number of patients we have to "create"
    simPatientNum <- numPerClust[m]-ncol(clustMatrixList[[m]])
    exprMatrix <- matrix(data=NA,ncol=numPerClust[m],nrow=nrow(clustMatrixList[[m]]))
    
    exprMatrix[ ,c(1:ncol(clustMatrixList[[m]]))] <- clustMatrixList[[m]]
    
    
    for(e in 1:simPatientNum){
      
      for(r in 1:nrow(clustMatrixList[[m]])){
        
        #pick a random patient index for each gene?
        exprMatrix[r ,(ncol(clustMatrixList[[m]])+e)] <- clustMatrixList[[m]][r,sample.int(ncol(clustMatrixList[[m]]),size=1)]
        
      }
      
    }
    
    tempExprMatrix[ , c(count:(count+(numPerClust[m]-1)))] <- exprMatrix
    cNames <- append(cNames,paste0(rep.int(colnames(clustMatrixList[[m]])[1],times=numPerClust[m]),"_",c(1:numPerClust[m])))
    count <- count + numPerClust[m]
    
    #end of loop m.
  }
  
  
  simMatrixList <- createSimGeneData(numSims=numSimDatasets,sampleExpr=tempExprMatrix,stddevNoise=stddevNoise,numRows=numRows,method=method,eigenValueMin=eigenValueMin)$simDataList
  
  for(s in 1:length(simMatrixList)){
    
    colnames(simMatrixList[[s]]) <- cNames
    #create fake gene names
    rownames(simMatrixList[[s]]) <- paste0("gene_",c(1:nrow(simMatrixList[[s]])))
    
    
    count <- 1
    simClustList[[s]] <- list()
    
    for(d in 1:length(clustMatrixList)){
      
      simClustList[[s]][[d]] <- simMatrixList[[s]][, c(count:(count+(numPerClust[m]-1)))]
      count <- count + numPerClust[m]
      names(simClustList[[s]])[d] <- names(clustMatrixList)[d]
      
      
    }
    
    #end of loop s.
  }
  output <- list(simMatrixList=simMatrixList,simClustList=simClustList)
  
  return(output)
  
}

#####
#simulation metrics
#####
#assumes clusters are in the same order each time.
createTrueEdges <- function(numRowClust=4,numColClust=1,numSimDatasets=2){
  
  trueEdges <- matrix(data=NA,ncol=2,nrow=(numColClust*numRowClust*numSimDatasets)*(numSimDatasets-1)/2)
  counter <- 0
  
  #have a feeling recursion could do this... 1 is an edge with numRowClust*numColClust*i where i = a loop over num sims. then move on toe
  #1+numRowClust*numColClust*i edge with numRowClust*numColClust*(i+1), etc.
  if(numSimDatasets>1){
    
    for(u in 1:(numSimDatasets-1)){
      
      for(m in 2:(numSimDatasets)){
        
        #want to start when simulation number is larger than this one.
        #previous ones "taken care of" by the edges in u-1
        if(m <= u){
          
          next
          
        }
        indexes <- c( (1+(numColClust*numRowClust)*counter):(1+(numColClust*numRowClust)*counter+ numColClust*numRowClust-1))
        
        startEdges <- c(1:(numColClust*numRowClust))+numColClust*numRowClust*(u-1)
        trueEdges[indexes, 1] <- startEdges
        trueEdges[indexes, 2] <- c(1:(numColClust*numRowClust))+numColClust*numRowClust*(m-1)
        
        counter <- counter +1
      }
      
    }
    
  }
  
  return(trueEdges)
  
}
###########
compute_edge_ROC_metrics <- function(trueEdges,predEdges,numTotalClusters){
  
  trueEdgesFound <- matrix(data=NA,ncol=2,nrow=1)
  trueEdgesMissed <- matrix(data=NA,ncol=2,nrow=1)
  
  for(r in 1:nrow(trueEdges)){
    
    rowIndices <- which(predEdges[,1]==trueEdges[r,1])
    
    if(length(rowIndices)>0){
      
      if(any(predEdges[rowIndices,2]==trueEdges[r,2])){
        
        #num_truePos <- num_truePos + 1
        #remove this edge: will help sum up FP later
        if(length(which(predEdges[rowIndices,2]==trueEdges[r,2]))>1){
          
          stop("indexing error.")
          
        }
        #cat("\n",predEdges[rowIndices[which(predEdges[rowIndices,2]==trueEdges[r,2])],],"\n")
        predEdges <- predEdges[-rowIndices[which(predEdges[rowIndices,2]==trueEdges[r,2])], , drop=FALSE]
        
        trueEdgesFound <- rbind(trueEdgesFound,trueEdges[r,c(1:2),drop=FALSE])
        
        #found correct edge - move on to next row
        next
        #no row index found.
      }
      #no row index found. 
    }
    # but flip node order, check again
    
    rowIndices <- which(predEdges[,2]==trueEdges[r,1])
    
    if(length(rowIndices)>0){
      
      if(any(predEdges[rowIndices,1]==trueEdges[r,1])){
        
        #num_truePos <- num_truePos + 1
        #remove this edge: will help sum up FP later
        if(length(which(predEdges[rowIndices,2]==trueEdges[r,2]))>1){
          
          stop("indexing error.")
          
        }
        #cat("\n",predEdges[rowIndices[which(predEdges[rowIndices,2]==trueEdges[r,2])],],"\n")
        predEdges <- predEdges[-rowIndices[which(predEdges[rowIndices,2]==trueEdges[r,2])], ,drop=FALSE]
        trueEdgesFound <- rbind(trueEdgesFound,trueEdges[r,c(1:2),drop=FALSE])
        
        #found a true edge - move on to next loop.
        next
      }
      
    }
    
    
    #didn't find an edge when there was a true edge.
    trueEdgesMissed <- rbind(trueEdgesMissed,trueEdges[r,c(1:2),drop=FALSE])
    #move on to next r index in true edges.
  }
  
  trueEdgesFound <- trueEdgesFound[-1, ,drop=FALSE]
  trueEdgesMissed <- trueEdgesMissed[-1, ,drop=FALSE]
  num_truePos <- nrow(trueEdgesFound)
  #any predicted left that weren't true edges?
  num_falsePos <- nrow(predEdges)
  #how many true edges? did we miss?
  num_falseNeg <- nrow(trueEdgesMissed)
  #total number of possible edges is number of clusters squared (a lot!)
  # minus numTotalClusters is minusing out of the diagonal
  #/2: we only want the upper (or lower) triangle to remove redundancy.
  num_trueNeg <- numTotalClusters^2/2 -numTotalClusters -nrow(trueEdges)
  
  if(num_falseNeg+num_truePos != nrow(trueEdges)){
    
    stop("\nCalculating FN or TP wrong.\n")
    
  }
  
  if( (num_trueNeg+num_falseNeg+num_truePos+num_falsePos) != (numTotalClusters^2/2-numTotalClusters)){
    
    stop("\nCalculating FP or TN wrong.\n")
    
  }
  
  #TPR is sensitivity, specificity is FPR
  #we mostly care about PPV and TPR. FPR will almost alwways be low.
  if(num_truePos >0){
    
    PPV <- num_truePos/(num_truePos+num_falsePos)
    
  }else{
    
    PPV <- 0
  }
  
  FPR <- num_falsePos/(num_falsePos+num_trueNeg)
  TPR <- num_truePos/(num_truePos + num_falseNeg)
  output <- list(PPV=PPV,TPR=TPR,FPR=FPR,FP=num_falsePos,TP=num_truePos,FN=num_falseNeg,TN=num_trueNeg,
                 trueEdgesFound=trueEdgesFound,trueEdgesMissed=trueEdgesMissed)
  
  return(output)
  
}