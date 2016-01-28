#library("plyr")
#library("GSEABase")
#library("Biobase")
#library("fdrtool")
#library("rmeta")
#library("RCurl")
#library("RDAVIDWebService")
#library("limma")
#library("biomaRt")
#library("biomaRt")
#library("Rgraphviz")
#http://david.abcc.ncifcrf.gov/content.jsp?file=WS.html
#http://bioinformatics.oxfordjournals.org/content/29/21/2810

#to do survival, other ggplot stuff (and advanced network plots.)


createPhenoMasterTableFromMatrixList <- function(esetList,dataMatrixList,sampleKeyColName="unique_patient_ID"){
  
  message("This function assumes samples/patient clinical data rows are not duplicated")
  
  if(!missing(esetList)){
    
    #check first object index: is it actually an expression object 
    #only a high-level check; not checking all objects in list.
    if(validObject(esetList[[1]])){
      
      for(d in 1:length(esetList)){
        
        if(d>1){
          
          tmp <- data.frame(pData(esetList[[d]]),rep.int(d,times=nrow(pData(esetList[[d]]))),rownames(pData(esetList[[d]])))
          if(all(colnames(phenoMasterDF)[1:(ncol(phenoMasterDF)-2)]==colnames(pData(esetList[[d]])))){
            
            #just do rbind - faster than join.
            phenoMasterDF <- rbind(phenoMasterDF,tmp)
            
          }else{
            
            phenoMasterDF <- join(phenoMasterDF,tmp,match="all",type="full",by=sampleKeyColName)
            
          }
          
        }else{
          
          phenoMasterDF <- data.frame(pData(esetList[[d]]),rep.int(d,times=nrow(pData(esetList[[d]]))),rownames(pData(esetList[[d]])))
          
        }
        
      }
      
    }
    
  }else if(!missing(dataMatrixList)){
    
    for(d in 1:length(dataMatrixList)){
      
      
      if(is.null(dataMatrixList[[d]]$phenoData)){
        
        stop(paste0("Index ",d,"dataMatrixList[[d]]$phenoData is null. Clinical (pheno) data should have this list index name."))
        
      }
      
      if(d>1){
        
        tmp <- data.frame(dataMatrixList[[d]]$phenoData,rep.int(d,times=nrow(dataMatrixList[[d]]$phenoData)),rownames(dataMatrixList[[d]]$phenoData))
        if(all(colnames(phenoMasterDF)[1:(ncol(phenoMasterDF)-2)]==colnames(pData(esetList[[d]])))){
          
          #just do rbind - faster than join.
          
          phenoMasterDF <- rbind(phenoMasterDF,tmp)
          
        }else{
          
          phenoMasterDF <- join(phenoMasterDF,tmp,match="all",type="full",by=sampleKeyColName)
          
        }
        
      }else{
        
        phenoMasterDF <- data.frame(dataMatrixList[[d]]$phenoData,rep.int(d,times=nrow(dataMatrixList[[d]]$phenoData)),rownames(dataMatrixList[[d]]$phenoData))
        
      }
      
    }
    
    
  }
  
  colnames(phenoMasterDF)[c((ncol(phenoMasterDF)-1):ncol(phenoMasterDF))] <- c("studyNum","sampleName")
  
  return(phenoMasterDF)
  
}

addClinicalVarToNodeAttributes <- function(sampleClustCommKey,phenoMasterDF){
  
  message("Note: assumes there is a column variable \'sampleName\' in both input data frames.")
  #do by study, in case sample names from matrix are actually duplicated.
  colnames(phenoMasterDF)[which(colnames(phenoMasterDF)=="studyNum")] <- "origStudyNum"
  #need to first cluster by studyNum: some names are duplicated across different studies.
  sampleClustCommPhenoData <- join(sampleClustCommKey,phenoMasterDF,by=c("origStudyNum","sampleName"),type="full",
                                   match="all")

  if(nrow( sampleClustCommPhenoData )>nrow(phenoMasterDF)){
    
    stop("Getting more rows than expected.")
  }
  return(sampleClustCommPhenoData)
                                   
}


binarizeMetaclustStudyStatus <- function(sampleClustCommKey){
  
  samplesInNoComm <- sampleClustCommKey[which(is.na(sampleClustCommKey[,"community"])),"sampleName"]
  
  if(!is.data.frame(sampleClustCommKey)){
    
    stop("sampleClustCommKey is not a data frame. Coercing it to a data.frame can alter numeric vs. character values; this is left up to the user before running the function.")
  }
  
  commSplit <- split(sampleClustCommKey,f=sampleClustCommKey$community)
  
  metaClustSampleStatus <- list()
  metaClustSampleNames <- list()
  
  for(n in 1:length(commSplit)){
    #still keep original community name as name of list.
    #careful: split make numbers into character strings.
    commNum <- as.numeric(as.character(unique(commSplit[[n]]$community)))
    metaClustSampleStatus[[commNum]] <- list()
    metaClustSampleNames[[commNum]] <- list()
    studySplit <- split(commSplit[[n]],f=commSplit[[n]]$studyNum)
    
    lengthStudySplit <- lapply(studySplit,FUN=nrow)
    
    if(any(unlist(lengthStudySplit))==0){
      
      studySplit <- studySplit[-which(lengthStudySplit==0)]
      
    }
    
    
    for(s in 1:length(studySplit)){
      
      #just to be safe...add as.character first if somehow still levels in here.
      studyNum <- as.numeric(as.character(unique(studySplit[[s]]$studyNum)))
      inComm <- as.character(studySplit[[s]]$sampleName)
      metaClustSampleNames[[commNum]][[studyNum]] <- inComm
      #community ,studyNum in original DF still numeric, so don't need to coerce to a charater.
      notInCommIndices <- intersect(which(sampleClustCommKey$community != commNum),which(sampleClustCommKey$studyNum==studyNum))
      
      if(length( notInCommIndices)>0){
        
        notInComm <- as.character(sampleClustCommKey[notInCommIndices, "sampleName"])
        #which ones were actually in ANY meta-cluster?
        notInComm <- setdiff(notInComm,samplesInNoComm)
        
      }else{
        
        notInComm <- NULL
        
      }
      
      if(length( notInComm)>0){
        
        #CHECK: did these all come from the same study?
        #remove levels.
        studySampleNames <- as.character(sampleClustCommKey[which(sampleClustCommKey$studyNum==studyNum),"sampleName"])
        #all the inComm,notInComm names should be in here!
        if(any(is.na(match(append(inComm,notInComm),studySampleNames)))){
          
          stop("\nNot calculating inComm,notInComm indices correctly")
          
        }
        metaClustSampleStatus[[commNum]][[studyNum]] <- append(rep.int(1,length(inComm)),rep.int(0,length(notInComm)))
        names(metaClustSampleStatus[[commNum]][[studyNum]]) <- append(inComm,notInComm)
        
      }else{
        
        metaClustSampleStatus[[commNum]][[studyNum]] <- rep.int(1,length(inComm))
        names(metaClustSampleStatus[[commNum]][[studyNum]]) <- inComm
        
      }
      
      names(metaClustSampleStatus[[commNum]])[studyNum] <- studyNum
      names(metaClustSampleNames[[commNum]])[studyNum] <- studyNum
      #end of loop s
    }
    #still keep original community name as name of list.
    names(metaClustSampleStatus)[commNum] <- commNum
    names(metaClustSampleNames)[commNum] <- commNum
    #end of loop N
  }
  
  
  output <- list(metaClustSampleStatus=metaClustSampleStatus,metaClustSampleNames=metaClustSampleNames)
  return(output)
  
}

#NOTE: if lots of NA genes in a study: this will result in skewed gene rankings.
computeRankMatrix <- function(metaClustSampleNames,featureNames,dataMatrixList,sampleClustCommKey,onlyIntersectingFeat=TRUE){
  
  if(length(names(metaClustSampleNames))==0){
    
    stop("No names on metaClustSampleNames")
  }
  
  if(onlyIntersectingFeat){
    
    studiesUsed <- unique(sampleClustCommKey[which(!is.na(sampleClustCommKey$community)),"studyNum"])
    
    for(s in 1:length(studiesUsed)){
      
      featureNames <- intersect(featureNames,rownames(dataMatrixList[[studiesUsed[s]]]))
      
      
    }
    
  }
  communityNames <- names(metaClustSampleNames)[which(names( metaClustSampleNames)!="")]
  
  
  for(c in 1:length(communityNames)){
    
    sampleNames <- unlist(metaClustSampleNames[[communityNames[c]]])
    
    rankMatrix <- matrix(data=NA,nrow=length(featureNames),ncol=length(sampleNames),dimnames=list(featureNames,sampleNames))
    
    
    if(length(names(metaClustSampleNames[[communityNames[c]]]))==0){
      
      stop("No names on metaClustSampleNames[[c]]")
      
    }
    studyNames <- names(metaClustSampleNames[[communityNames[c]]])[which(names( metaClustSampleNames[[communityNames[c]]])!="")]
    
    for(s in 1:length(studyNames)){
      
      #remember: names of data matrix list may not necessarily be just the study number. but the numeric index is.
      fullStudyMatrix <- dataMatrixList[[as.numeric(studyNames[s])]]
      tmp <- fullStudyMatrix[na.omit(match(featureNames,rownames(fullStudyMatrix))) ,metaClustSampleNames[[communityNames[c]]][[studyNames[s]]] ,drop=FALSE]
      
      if(c==1 && s==1){
        
        groupings <-  data.frame(rep.int(as.numeric(as.character(communityNames[c])),times=length( metaClustSampleNames[[communityNames[c]]][[studyNames[s]]])), 
                                 rep.int(as.numeric(as.character(studyNames[s])),times=length( metaClustSampleNames[[communityNames[c]]][[studyNames[s]]])),
                                 metaClustSampleNames[[communityNames[c]]][[studyNames[s]]], stringsAsFactors=FALSE)
        
      }else{
        groupings <- rbind(groupings,data.frame(rep.int(as.numeric(as.character(communityNames[c])),times=length( metaClustSampleNames[[communityNames[c]]][[studyNames[s]]])), 
                                                rep.int(as.numeric(as.character(studyNames[s])),times=length( metaClustSampleNames[[communityNames[c]]][[studyNames[s]]])),
                                                metaClustSampleNames[[communityNames[c]]][[studyNames[s]]], stringsAsFactors=FALSE))
        
      }
      
      if(nrow(tmp)==0){
        
        stop("After selecting for feature names and samples, no data rows left.")
        
      }
      
      #fill in for ranks of these genes
      #gene indices should match up?
      rankMatrix[na.omit(match(rownames(tmp),rownames(rankMatrix))) ,metaClustSampleNames[[communityNames[c]]][[studyNames[s]]] ] <-  apply(tmp, MARGIN=2,FUN=function(sampleUnit){
        
        sort.int(sampleUnit,index.return=TRUE)$ix
        
      })
      
      #end of study s  
    }
    
    if(c >1){
      
      finalRankMatrix <- cbind(finalRankMatrix,rankMatrix)
    }else{
      
      finalRankMatrix <- rankMatrix
    }
    
    #end of community looping.
  }
  
  colnames(groupings) <- c("community","studyNum","sampleName")
  output <- list(rankMatrix=finalRankMatrix,groupings=groupings,filteredFeatures=featureNames)
  return(output)
  
}

#rankMatrix can be actual ranks or effect size.
computeFeaturePvalues <- function(rankMatrix,featureNames,groupings){
  #now do a kruskal wallace test
  featureTests <- array(data=NA,dim=length(featureNames),dimnames=list(featureNames))
  for(r in 1:nrow(rankMatrix)){
    #REMOVE which features are NA.
    rankArray <- rankMatrix[r,which(!is.na(rankMatrix[r,]))]
    groups <- groupings[which(!is.na(rankMatrix[r,])),"community"]
    #what if a whole community is removed??
    featureTests[r] <- kruskal.test(rankArray,g=as.numeric(groups))$p.value
    
  }
  #NOW: fdr correct
  #note: "hochberg" original 1988 method tends to be overly conservative.
  #"fdr" is hochberg 1995 method.
  featureQvalues <- p.adjust(featureTests,method="fdr")
  
  featurePvalues <- data.frame(featureTests,featureQvalues)
  rownames(featurePvalues) <- featureNames
  colnames(featurePvalues) <- c("kruskal.pvalue","fdr.qvalue")
  
  return(featurePvalues)
  
}

#select only sig genes, take median for each meta-cluster

commMedianRank <- function(rankMatrix,groupings){
  
  commNames <- unique(groupings[,"community"])
  featureMedians <- matrix(data=NA,ncol=length(commNames),nrow=nrow(rankMatrix),dimnames=list(rownames(rankMatrix),colnames=commNames))
  
  for(c in 1:length(commNames)){
    
    featureMedians[ ,c] <- rowMedians(rankMatrix[,which(groupings[,"community"]==commNames[c])],na.rm=TRUE)
    
  }
  
  return(featureMedians)
  
}


######compute effect sizes
###"other" groupings include any other cluster in dataset, even if dropped
#COME BACK: try using samr for microarray data?
#remember: samr's version for rna-seq requires you to normalize it
#using their read depth method. so not as readily applicable after clustering.
computeMetaclustEffectSizes <- function(metaClustSampleNames,dataMatrixList,
                                        featureNames,minOtherClass=5,computeWilcoxon=FALSE){
  
  
  fishersMethod = function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower=FALSE);
  #sample names may not be 1,2,3...in linear order if some were dropped earlier bc were too small.
  commNames <- names(metaClustSampleNames)[intersect(which(!is.na(names(metaClustSampleNames))),
                                                     which(names(metaClustSampleNames)!=""))]
  
  
  summHedgeG_ES <- matrix(data=NA,nrow=length(featureNames),        
                          ncol=length(commNames),
                          dimnames=list(featureNames,
                          commNames)); 
  
  summHedgeG_ES_se <- matrix(data=NA,nrow=length(featureNames),                       
                          ncol=length(commNames),
                          dimnames=list(featureNames,
                          commNames)); 
  
  summWilcoxon_qvalue <- matrix(data=NA,nrow=length(featureNames),             
                             ncol=length(commNames),
                             dimnames=list(featureNames,
                             commNames)); 
  
  
  numDatasetsPerGene <- matrix(data=NA,nrow=length(featureNames),        
                           ncol=length(commNames),
                           dimnames=list(featureNames,
                                         paste0("num_datasets_comm",commNames))); 
  
  wilcoxon_qvalueList <- list()
  hedgeGList <- list()
  hedgeG_seList <- list()
                          
                          
  
  for(c in 1:length(commNames)){
    
    studyNames <- names(metaClustSampleNames[[commNames[c]]])[intersect(which(!is.na( names(metaClustSampleNames[[commNames[c]]]))),
                                                                        which( names(metaClustSampleNames[[commNames[c]]])!=""))]
    

                                                                 

    hedgeG <- matrix(data=NA,nrow=length(featureNames),                     
                     ncol=length(studyNames),dimnames=list(featureNames,
                  studyNames)); 

    hedgeG_se <- matrix(data=NA,nrow=length(featureNames),
                    ncol=length(studyNames), dimnames=list(featureNames,
              studyNames)); 


    wilcoxon_qvalue <- matrix(data=NA,nrow=length(featureNames),                                               
                          ncol=length(studyNames), dimnames=list(featureNames,
              studyNames));       

    for(s in 1:length(studyNames)){

      #study name is actually just a number.
      studyIndex <- as.numeric(studyNames[s])                
      dataMatrix <- dataMatrixList[[studyIndex]]
      
      metaClusterSampleNames <- metaClustSampleNames[[commNames[c]]][[studyNames[s]]]
      otherSampleNames <- setdiff(colnames(dataMatrix),metaClusterSampleNames)
      
      
      if(length(otherSampleNames)>minOtherClass){   
        #drop=FALSE in case user allows othr samples to be of length 1
        matrix1 <- dataMatrix[,  metaClusterSampleNames,drop=FALSE]
        matrix2 <- dataMatrix[,  otherSampleNames,drop=FALSE]
        
      for(g in 1:length(featureNames)){

        if(featureNames[g] %in% rownames(dataMatrix)){

          #this assumes variances between the two groups are the same for this gene...hedge's g sort of helps
          #http://www.polyu.edu.hk/mm/effectsizefaqs/effect_size_equations2.html
          #want hedge's g bc may have pretty different sample sizes
          #adapted from Purvesh's code: getES() function.
          #assumes no NAs.
          mean_diff <- mean(matrix1[featureNames[g], ,drop=FALSE])-mean(matrix2[featureNames[g], ,drop=FALSE]);
          sd1  <- sd(matrix1[featureNames[g], ,drop=FALSE]);  
          sd2 <- sd(matrix2[featureNames[g], ,drop=FALSE]);
          n1 <- length(metaClusterSampleNames);
          n2 <- length(otherSampleNames);
          pooled_sd   <- sqrt( ( (n1-1)*sd1^2 + (n2-1)*sd2^2 )/( n1 + n2 - 2 ) );
          cf   <- 1 - 3/( 4*(n1 + n2) - 9 );
          
          hedgeG[g,s]   <- cf * mean_diff/pooled_sd;
          hedgeG_se[g,s] <- sqrt( (n1+n2)/(n1*n2) + 0.5*(hedgeG[g,s] )^2 /(n1+n2-3.94) );
          
          if(computeWilcoxon){
          #get p-value; compute/update to q-value after loop through all genes
          wilcoxon_qvalue[g,s] <- wilcox.test(t(matrix1[featureNames[g], ]), t(matrix2[featureNames[g],]),alternative = c("two.sided"),
                                                    paired=FALSE)$p.value
          
          }
          
          
          
            }else{
              
              hedgeG[g,s]    <- NA
              hedgeG_se[g,s] <- NA
              wilcoxon_qvalue[g,s] <- NA         
           
            }
          
          
      
      #end of looping through all features
      }
      
      }else{
        
        #not enough "other" samples: return NA for all features.
        hedgeG[ ,s]    <- NA
        hedgeG_se[ ,s] <- NA
        wilcoxon_qvalue[ ,s] <- NA         
        
      }
      
      
      if(computeWilcoxon){
      #for this study: fdr correct all of the genes.
      wilcoxon_qvalue[which(!is.na(wilcoxon_qvalue[,s])), s] <- 
        
        p.adjust(wilcoxon_qvalue[which(!is.na(wilcoxon_qvalue[ ,s])), s],method="fdr")
      
      }
      
      
      #end of looping through all studies for this meta-cluster/community
    }
    #now compute summary stats across all communities
    #ours is random? A “group” effect is random if we can think of the levels we
    #observe in that group to be samples from a larger population.
    # Example: if collecting data from different medical centers,
    #“center” might be thought of as random.
    # Example: if surveying students on different campuses,
    #“campus” may be a random effect.
    #http://statweb.stanford.edu/~jtaylo/courses/stats203/notes/fixed+random.pdf
    
    for(g in 1:length(featureNames)){
      
    d <- as.vector(hedgeG[g, which(!is.na(hedgeG[g, ]))]);
    se <- as.vector(hedgeG_se[g, which(!is.na(hedgeG[g, ]))]);
    numDatasetsPerGene[g,c] <- length(which(!is.na(hedgeG[g, ])))
    
    if(length(d)>0 && length(se)>0){
      #weighted mean; inverse weigthing by standard deviation
      ES <- meta.summaries(d=d, se=se, method=c("random"),
                           logscale=TRUE,
                           conf.level=0.95);
      
      summHedgeG_ES[g,c] <- ES$summary;
      summHedgeG_ES_se[g,c] <- ES$se.summary;
    
    }else{
      
      summHedgeG_ES[g,c] <- NA
      summHedgeG_ES_se[g,c] <- NA
      
    }
    #summarize p-value across all studies (can use Fisher's because studies are independent.)
    
    if(computeWilcoxon){
    summWilcoxon_qvalue[g,c] <-  fishersMethod(wilcoxon_qvalue[g, which(!is.na(wilcoxon_qvalue[g,]))])
    
    }
    
    #end of looping over genes
    }

    
    
    if(computeWilcoxon){
    wilcoxon_qvalueList[[commNames[c]]] <- wilcoxon_qvalue
    }
    
    hedgeGList[[commNames[c]]] <- hedgeG
    hedgeG_seList[[commNames[c]]] <- hedgeG_se
    
    #end of looping over communities
  }
  
  
  if(computeWilcoxon){
  #correct across the number of meta-clusters testing.
  #could do apply I guess here...but just sticking for for loop format throughout the function:
  for(g in 1:length(featureNames)){
    
    summWilcoxon_qvalue[g,which(!is.na(summWilcoxon_qvalue[g, ]))] <- p.adjust(summWilcoxon_qvalue[g,which(!is.na(summWilcoxon_qvalue[g, ]))],method="fdr");
  
  
  }
  
  }
  
  output <- list(summWilcoxon_qvalue=summWilcoxon_qvalue,summHedgeG_ES=summHedgeG_ES,
                 summHedgeG_ES_se=summHedgeG_ES_se,wilcoxon_qvalueList=wilcoxon_qvalueList,
                 hedgeGList=hedgeGList,hedgeG_seList=hedgeG_seList,numDatasetsPerGene=numDatasetsPerGene)
  
  return(output)

  #EOF
}

#thresh is greater than.
selectMetaclustSigGenes <- function(computeMetaclustEffectSizesOutput,includeWilcoxon=FALSE,
                                    ESthresh=0,qvalueThresh=1){
  
  commNames <- colnames(computeMetaclustEffectSizesOutput$summHedgeG_ES)
  
  sigMetaclustGenesWilcox <- list()
  sigMetaclustGenesES_pos <- list()
  sigMetaclustGenes_pos <- list()
  sigMetaclustGenesES_neg <- list()
  sigMetaclustGenes_neg <- list()
  

    
  for(c in 1:length(commNames)){
    
    if(includeWilcoxon){
      
    sigMetaclustGenesWilcox[[commNames[c]]] <- rownames(computeMetaclustEffectSizesOutput$summWilcoxon_qvalue[
      which(computeMetaclustEffectSizesOutput$summWilcoxon_qvalue[,commNames[c]]<=qvalueThresh),commNames[c],drop=FALSE])
    
    }
    
    sigMetaclustGenesES_pos[[commNames[c]]] <- rownames(computeMetaclustEffectSizesOutput$summHedgeG_ES[
      which(computeMetaclustEffectSizesOutput$summHedgeG_ES[,commNames[c]]>=ESthresh),commNames[c],drop=FALSE])
    
    if(includeWilcoxon){
      
    sigMetaclustGenes_pos[[commNames[c]]] <- intersect(sigMetaclustGenesWilcox[[commNames[c]]],
                                                       sigMetaclustGenesES_pos[[commNames[c]]])
    
    }
    
    #negative: reverse ESthresh sign.
    sigMetaclustGenesES_neg[[commNames[c]]] <- rownames(computeMetaclustEffectSizesOutput$summHedgeG_ES[
      which(computeMetaclustEffectSizesOutput$summHedgeG_ES[,commNames[c]] <= -ESthresh),commNames[c],drop=FALSE])
    
    if(includeWilcoxon){
      
    sigMetaclustGenes_neg[[commNames[c]]] <- intersect(sigMetaclustGenesWilcox[[commNames[c]]],
                                                       sigMetaclustGenesES_neg[[commNames[c]]])
    }
      
  }
  
  
  
  output <- list( sigMetaclustGenesWilcox= sigMetaclustGenesWilcox, sigMetaclustGenesES_pos= sigMetaclustGenesES_pos,
                  sigMetaclustGenes_pos=sigMetaclustGenes_pos,sigMetaclustGenesES_neg=sigMetaclustGenesES_neg,
                  sigMetaclustGenes_neg= sigMetaclustGenes_neg)
  
  return(output)
  
}

summarizePosESMetaclustGenes <- function(selectMetaclustSigGenesOut,computeMetaclustEffectSizesOutput){
  
  featureNames <- c()
  commNames <- colnames(computeMetaclustEffectSizesOutput$summHedgeG_ES)
  
  for(c in 1:length(commNames)){
    
    featureNames <- append(featureNames,selectMetaclustSigGenesOut$sigMetaclustGenesES_pos[[commNames[c]]])
    
  }
  
  featureNames <- unique(featureNames)
  ESMatrix <- matrix(data=NA,ncol=length(commNames),nrow=length(featureNames),dimnames=list(featureNames,commNames))
  
  #all featureNames will be in this matrix, even if NA.
  for(c in 1:length(commNames)){
    
    ESMatrix[selectMetaclustSigGenesOut$sigMetaclustGenesES_pos[[commNames[c]]],commNames[c]] <-
      computeMetaclustEffectSizesOutput$summHedgeG_ES[selectMetaclustSigGenesOut$sigMetaclustGenesES_pos[[commNames[c]]], commNames[c]]
    
  }
  
  if(any(is.na(rowMeans(ESMatrix,na.rm=TRUE)))){
    
    stop("Not computing ESMatrix correctly; getting all NAs in some rows.")
  }
  
  #now add # datsets per gene
  ESMatrix <- cbind(ESMatrix,computeMetaclustEffectSizesOutput$numDatasetsPerGene[featureNames, ])
  #now can do a heatmap of these genes by effect size.
  return(ESMatrix)
  
}

summarizeNegESMetaclustGenes <- function(selectMetaclustSigGenesOut,computeMetaclustEffectSizesOutput){
  
  featureNames <- c()
  commNames <- colnames(computeMetaclustEffectSizesOutput$summHedgeG_ES)
  
  for(c in 1:length(commNames)){
    
    featureNames <- append(featureNames,selectMetaclustSigGenesOut$sigMetaclustGenesES_neg[[commNames[c]]])
    
  }
  
  featureNames <- unique(featureNames)
  ESMatrix <- matrix(data=NA,ncol=length(commNames),nrow=length(featureNames),dimnames=list(featureNames,commNames))
  
  #all featureNames will be in this matrix, even if NA.
  for(c in 1:length(commNames)){
    
    ESMatrix[selectMetaclustSigGenesOut$sigMetaclustGenesES_neg[[commNames[c]]],commNames[c]] <-
      computeMetaclustEffectSizesOutput$summHedgeG_ES[selectMetaclustSigGenesOut$sigMetaclustGenesES_neg[[commNames[c]]], commNames[c]]
    
  }
  
  if(any(is.na(rowMeans(ESMatrix,na.rm=TRUE)))){
    
    stop("Not computing ESMatrix correctly; getting all NAs in some rows.")
  }
  
  #now add # datsets per gene
  ESMatrix <- cbind(ESMatrix,computeMetaclustEffectSizesOutput$numDatasetsPerGene[featureNames, ])
  #now can do a heatmap of these genes by effect size.
  return(ESMatrix)
  
}


GSEA_posMetaclustGenes <- function(selectMetaclustSigGenesOut,selectMetaclustSigGenesOutIndexName="sigMetaclustGenesES_pos",
                                     genomeSize=2000,refGeneLists=NULL, 
                                   refGeneListDir="/home/data/MSigDB/GSEA_base_MSigDB_lists_merged.RData.gzip",
                                   qvalueThresh=.1,fileTag="test"){
  
  commNames <- names(selectMetaclustSigGenesOut[[1]])[intersect(which(names(selectMetaclustSigGenesOut[[1]])!=""),
                                                           which(!is.na(names(selectMetaclustSigGenesOut[[1]]))))]
  
  if(is.null(commNames)){
    
    stop("error: in GSEA_posMetaclustGenes, getting NULL community names")
  }
  
  GSEAout <- list()
  
  for(c in 1:length(commNames)){
    
    tmp <- GSEA(testGeneVector=selectMetaclustSigGenesOut[[selectMetaclustSigGenesOutIndexName]][[commNames[c]]],method=c("hypergeometric"),genomeSize=genomeSize,
                refGeneListDir= refGeneListDir,refGeneLists=refGeneLists)
    
    GSEAout[[commNames[c]]] <- tmp$qvalues[which(tmp$qvalues<=qvalueThresh)]
    
  }
  
  return(GSEAout)
  
}
#for allowed geneType inputs: http://david.abcc.ncifcrf.gov/content.jsp?file=DAVID_API.html
#email: must be registered with DAVID knowledge database.
GSEA_DAVID <- function(testGeneVector,idType="GENE_SYMBOL",yourEmail="katie.planey@stanford.edu",
                       EASE_thresh=0.1,count_thresh=2L,qvalueThreshFunct=.05,
                       pvalueThreshCluster=.1,fileTag="commmunity1",
                       saveDir="/home/kplaney/ovarian_analysis/DAVID/"){

  #DAVID APIs (alpha) allow other bioinformatics web sites to directly link to DAVID tools and functions ONLY 
  #for light-duty jobs (i.e. a gene list with no more than 400 genes)
  #also need to wait 10 seconds between queries.
  if(length(testGeneVector)>400){
    
    stop("DAVID can only look at gene lists with fewer than 400 genes.")
    
  }

  if(idType=="GENE_SYMBOL"){
    #first: map gene symbols to entrez IDs. R package doesn't take gene symbols as inputs.
    entrezIDs <- mapGeneIDs(geneticFeatures=testGeneVector,bioMartSource="ensembl",bioMartDataset="hsapiens_gene_ensembl",
                         geneticFeatureSource = 'hgnc_symbol',returnOnlyOneSymbolPerFeature=TRUE)
    
    entrezIDs <- entrezIDs[ ,"entrezgene"]
  
    idType="ENTREZ_GENE_ID"
  }
  
  message("Connecting to DAVID web service. Occassionally it times out and you need to re-run this function.")
  david<-DAVIDWebService$new(email=yourEmail)
  connect(david)
  #idTypes: getIdTypes(david)
  if(!idType %in% getIdTypes(david)){
    
    stop("Must provide one of the accepted DAVID webservice ID types: ",getIdTypes(david))
  
  }
  result<-addList(david,inputIds=entrezIDs,
                      idType=idType,
                      listName="DAVID_genes", listType="Gene")
  #need to update.
  #setAnnotationCategories(david, c("GOTERM_BP_ALL",
   #                                 "GOTERM_MF_ALL", "GOTERM_CC_ALL"))

  termCluster<-getClusterReport(david, type="Term")
  getClusterReportFile(david, type="Term",
                           fileName=paste0(saveDir,"/",fileTag,"_termClusterReport.tab"))
  
  #chart option uses statistics, table option does not.
  functAnnotChart <- getFunctionalAnnotationChart(david)
  if(nrow(functAnnotChart)>0){
  
    getFunctionalAnnotationChartFile(david,
                                   fileName=paste0(saveDir,"/",fileTag,"_functAnnotChartReport.tab"),
                                   threshold=EASE_thresh,count=count_thresh)
  
  
  

  sigFunct <- functAnnotChart[which(functAnnotChart[,"Benjamini"]<=qvalueThreshFunct) ,c("Category","Term","Benjamini")]
  
  }else{
    
    stop("Returning an empty functional annotation chart.")
    
  }
  sigClusters <- which(summary(termCluster)[,"Enrichment"]<=pvalueThreshCluster)
  
  output <- list(termCluster=termCluster,termClusterSummary=summary(termCluster),
            clusterReportFilePath=paste0(saveDir,"/",fileTag,"_termClusterReport.tab"),
            functionalAnnotationChart=functAnnotChart,
            functAnnotCategories=categories(functAnnotChart),
            functionalAnnotFilePath=paste0(saveDir,"/",fileTag,"_functAnnotChartReport.tab"),
            sigFunctGroups=sigFunct,sigClusters=sigClusters
            )
  #make it wait at least 10 seconds, in case running in a loop.
  Sys.sleep(10)
  return(output)
  
}

# DAVID_plotClusters <- function(termClusterSummary,clusterReportFilePath){
#   
#   fileName<-system.file(clusterReportFilePath,
#                             package="RDAVIDWebService")
#   #untar(fileName)
#   # termCluster<-DAVIDTermCluster(untar(fileName, list=TRUE))
#   termCluster<-DAVIDTermCluster(fileName
#   plot2D(termCluster, clustNumber)
#   davidGODag<-DAVIDGODag(members(termCluster)[[clustNumber]], pvalueCutoff=0.1, "CC")
#   #requires Rgraphviz :()
#   #plotGOTermGraph(g=goDag(davidGODag), r=davidGODag, max.nchar=40, node.shape="ellipse")
#   
#   
# }
# 
# DAVID_plotFunctAnnotChart(){
#   
#   fileName<-system.file("files/termClusterReport1.tab.tar.gz",
#                         package="RDAVIDWebService")
#   untar(fileName)
#   functChart<- DAVIDFunctionalAnnotationChart(untar(fileName, list=TRUE))
#   
#   plot2D(DAVIDFunctionalAnnotationChart(functChart,
#          color=c("FALSE"="black", "TRUE"="green"))
# }
# #listAttributes(mart) after setting up "mart" tells you all of the feature sources you can use.
# #ex: "refseq_mrna" for refseq mrna ID.
mapGeneIDs <- function(geneticFeatures,bioMartSource="ensembl",bioMartDataset="hsapiens_gene_ensembl",
                                          geneticFeatureSource = 'hgnc_symbol',returnOnlyOneSymbolPerFeature=TRUE){
  

  mart <- useMart(bioMartSource,dataset=bioMartDataset)
  
  #can't have ' in a genetic features name for biomaRt to recognize it.
  if(length(dim(geneticFeatures))>0){
    
    #take first column. otherwise gsub does weird stuff.
    geneticFeatures <- geneticFeatures[,1];
    
  }
  
  if(length(geneticFeatures)>450){
    
    stop("Your gene list is probably too long for biomaRt in R; try looping over chunks of the list.")
    
  }
  #can't have ' in a genetic features name for biomaRt to recognize it.
  geneticFeatures<- gsub("\'","",geneticFeatures);
  #oftentimes this period is added to ensembl IDs - but it's not in biomart.
  geneticFeatures <- strsplit2(geneticFeatures,"\\.")[,1];
  
  #some code alters geneticFeatures - we'll need the original list for some purposes later on down in the code
  #(returning 1 row per original genetic feature source item.)
  geneticFeaturesFull <- geneticFeatures;
  warning("duplicate entries for a certain BiomaRt source may be returned.")
  warning("if you're feeding in HGNC gene symbols,make sure these are the latest updated version, otherwise Biomart may not recognize them.")
  #listAttributes(mart)
  
  #use BM, as BMList is slower and loops over each value you put in. BM looks at all values more efficiently.
   
      linkedData <- getBM(c(geneticFeatureSource,'hgnc_symbol','ensembl_gene_id','start_position','end_position','band','chromosome_name','gene_biotype','entrezgene'),
                          filters=c(geneticFeatureSource),values=geneticFeatures,
                          mart=mart,uniqueRows=TRUE); 
  

  if(dim(linkedData)[1] < length(geneticFeatures)){
    
    warning("couldn't find a biomaRt record for all features.");
    
  }
  
  if(returnOnlyOneSymbolPerFeature){
    
    if(length(geneticFeaturesFull) != length(unique(geneticFeaturesFull))){
      
      warning("Some of your genetic feature inputs were duplicated. So you'll get duplicated output rows by choosing to have only 1 symbol returned per 
              (unique) feature input.")
      
    }
    
    if(length(dim(linkedData))==0){
      #so can grab column ID if it's a vector (shouldn't occur bc I put the feature source along with the gene symbol - so min 2 columns.)
      linkedData <- as.matrix(linkedData);
      
    }
    
    featureInfo <- linkedData;
    linkedData <- matrix(data=NA,nrow=length(geneticFeaturesFull),ncol=dim(linkedData)[2]);
    colnames(linkedData) <- colnames(featureInfo);
    
    #looking for original column we fed in - want this many rows and to match on this.
    #will only return 1 row per unique geneticFeatureSource
    colIndex <- which(colnames(linkedData)==geneticFeatureSource)[1];
    
    for(g in 1:length(geneticFeaturesFull)){
      
      #may have duplicated entries per one feature.
      #technically, if just did match function, would only return first match.
      #but I liked doing this stuff deliberately so I can remember it better later...
      matchIDs <- which(featureInfo[,colIndex ]==geneticFeaturesFull[g]);
      
      #may have returned NULL.
      if(!is.na(all(matchIDs))){
        
        #just take the first match.
        linkedData[g,] <- unlist(featureInfo[matchIDs[1],]);
        
      }else{
        
        cat("nothing for ", geneticFeaturesFull[g]);
        #data is already NA. no need for this.
        #linkedData[g,] <- NA;
        
      }
      
      #COME BACK: is this working OK?
      if(length(matchIDs)>1){
        warn <- TRUE;
        cat("\n index",g," gene ",geneticFeaturesFull[g]," had multiple records returned.\n");
        
      }else{
        warn <- FALSE;
      }
      
      
    }
    
    if(warn){
      
      warning("Found more than 1 gene symbol for some input features - printed them out above. The first symbol was chosen/kept each time - usually the best bet.")
    }
    
  }
  
  
  return(linkedData);
  
  }

# #COME BACK: need to debug how compute this....
 GSEA <- function(testGeneVector,refGeneLists=NULL,method=c("hypergeometric","fisher"),
                  refGeneListDir="/home/data/MSigDB/GSEA_base_MSigDB_lists_merged.RData.gzip"){
   
  #all reference - take union
   #genes in reference gene lists + length test GeneVector
   
   warning("\nThis code assumes that all of your genes in your test gene list and ref gene list are in the genome.");
   
   if(is.null(refGeneLists)){
     
   load(refGeneListDir)
   #cat("\nUsing default MSigDB lists: MSigDB_onco_symbols, MSigDB_CanPath_symbols,MSigDB_TFT_symbols,MSigDB_immun_symbols,
    #   and MSigDB_cancerNeigh_symbols.\n");
   
   }
   
  refGeneLists <- GSEA_base_MSigDB_lists_merged;
  
  genome <- c()
  
  for(r in 1:length(refGeneLists)){
    
    genome <- union(genome,unlist(refGeneLists[[r]]))
  }
  
  #now add on test gene vector
  genome <- union(genome,testGeneVector)
  message("genome size is ",length(genome))
  genomeSize <- length(genome)
  pvalues <- list();
  
  for(r in 1:length(refGeneLists)){
    
    numIntersect <- length(intersect(testGeneVector,refGeneLists[[r]]));
    
    if(numIntersect>0){
     # message(names(refGeneLists)[r])
      if(method=="hypergeometric"){
        #we are taking the area under the curve to calculate a <= scenario.
        #from 1 to the number of genes we've intersected - this is our vector of quantiles
        #"representing the number of white balls drawn without replacement from an urn which contains both black and white balls"
        #length(refGeneLists[[r]]) is the total number of balls we drew from the urn. so the urn has to be the entire genome here?
        #tricky when filter down gene list first.  a larger genomeSize will make for more conservative p-values.
        
        #phyper(x, m, n, k) gives the probability of getting x or fewer
        #phyper(numIntersect-1, length(geneList), genomeSize, length(refGeneLists[[r]]), lower.tail=FALSE)
        #lower tail =FALSE same as 1-:
        #what's the probability of have a larger tail than this?
        pvalues[[r]] <- 1- phyper(numIntersect-1, length(testGeneVector), genomeSize, length(refGeneLists[[r]]));
        
        #this will also give the same answer:
        #we want to know how probably it is by chance that we'd see MORE (or equal?) intersecting genes than we have.
        #pvalue <- min(1-cumsum(dhyper(0:(numIntersect-1), length(testGeneVector), genomeSize-length(testGeneVector), length(refGeneLists[[r]])) ));
      }else if(method=="fisher"){
        
        counts <-  matrix(data = c(numIntersect, length(refGeneLists[[r]])-numIntersect, length(testGeneVector),genomeSize-length(testGeneVector)),nrow = 2);
        pvalues[[r]] <- fisher.test(counts,alternative="greater")$p.value;
        
      }else{
        
        stop("\nPlease input method=hypergeometric or method=fisher as the statistical method.");
        
      }
      
      
      
    }else{
      
      #set equal to one? no sure how to do this here....
      pvalues[[r]] <- NA;
      
    }
    
  }
  
  nonNA_indices <- which(!is.na(unlist(unlist(pvalues))));
  nonNA_pvalues <- unlist(pvalues)[which(!is.na(unlist(unlist(pvalues))))];
  qvalues <- array(data=NA,dim=length(unlist(pvalues)));
  #want to preserve NA values
  #now do multiple hyp testing
  qvalues[nonNA_indices] <- p.adjust(p=nonNA_pvalues,method="fdr");
  
  names(pvalues) <- names(refGeneLists)
  names(qvalues) <- names(refGeneLists)
  output <- list(qvalues=qvalues,pvalues=unlist(pvalues),refGeneListNames=names(refGeneLists));
  
  return(output);
  
}

computeAdjMatricesNullMatrixList <- function(dataMatrixList,numIter=5,numParallelCores=1,
                                   clustSampleIndexList,clustFeatureIndexList,
edgeMethod=c("spearman","pearson","kendall","Euclidean","cosine",
                      "Manhattan","Minkowski","Mahalanobis"),
numSims=500,
outputFile="./CoINcIDE_messages.txt",
centroidMethod=c("mean","median")){
  
  meanMetric <- c()
  indPvalues <- c()
  indFractNN <- c()
  
  CoINcIDE_NullOutputList <- list()
  #use R parallel here??
  for(i in 1:numIter){
    
    message("Running CoINcIDE on a null list of data matrices for iteration ",i)
    nullMatrixList <- createNullDataMatrixList(dataMatrixList)
    CoINcIDE_NullOutputList[[i]] <- computeAdjMatrices(dataMatrixList=nullMatrixList,clustSampleIndexList=clustSampleIndexList,
                                                            clustFeatureIndexList=clustFeatureIndexList,
    edgeMethod=edgeMethod,numSims=numSims, outputFile=outputFile,checkNA=FALSE,centroidMethod=centroidMethod)
    meanMetric <- append(meanMetric, as.vector(CoINcIDE_NullOutputList[[i]]$computeTrueSimilOutput$meanMetricMatrix))
    indPvalues <- append(indPvalues,as.vector(CoINcIDE_NullOutputList[[i]]$pvalueMatrix))
    indFractNN <- append(indFractNN,as.vector(CoINcIDE_NullOutputList[[i]]$computeTrueSimilOutput$trueFractNNmatrix))
  }
  
  meanMatrixQuantiles <- quantile(meanMetric,na.rm=TRUE,probs=seq(0,1,.05))
  indPvalueQuantiles <- quantile(indPvalues,na.rm=TRUE,probs=seq(0,1,.05))
  indFractNNQuantiles <- quantile(indFractNN,na.rm=TRUE,probs=seq(0,1,.05))
  
  nullOutput <- list(CoINcIDE_NullOutputList=CoINcIDE_NullOutputList,meanMatrixQuantiles=meanMatrixQuantiles,
                     indPvalueQuantiles=indPvalueQuantiles,indFractNNQuantiles=indFractNNQuantiles)
  return(nullOutput)
  
}
#do in loop nullMatrixList <- createNulDataMatrixList(dataMatrixList)
globalFDR <- function(CoINcIDE_outputList,
edgeMethod=c("spearman","pearson","kendall","Euclidean","cosine",
                      "Manhattan","Minkowski","Mahalanobis"),minTrueSimilThresh=-Inf,maxTrueSimilThresh=Inf,
outputFile="./CoINcIDE_messages.txt",fractFeatIntersectThresh=0,numFeatIntersectThresh=0 ,clustSizeThresh=0, clustSizeFractThresh=0,
    meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, 
saveDir = "./",experimentName = "nullTest",
    commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=3,minFractNN =.8,findCommWithWeights=FALSE,
minNumEdgesForCluster=1,fractEdgesInVsOutComm=0,
fractEdgesInVsOutEdge=0
   
    
){
  
  numFalseIndPvaluesBeforeOtherThresh <- list()
  numFalseEdges_ind <- list()
  numFalseEdges_afterThresh <- list()
  numFalseEdges_afterGN <- list()
  numFalseMetaclusters_afterGN <- list()
  numFalseClusters_afterGN <- list()
  numFalseDatasetsInNetwork <- list()
  falseEdgeLists <- list()
  meanMeanMetric <- c()
  meanFractNN <- c()
  
  for(i in 1:length(CoINcIDE_outputList)){
    
  CoINcIDE_output <- CoINcIDE_outputList[[i]]
  #look at if just thresholded by indEdgePvalueThresh
  numFalseIndPvaluesBeforeOtherThresh[[i]] <- length(which(CoINcIDE_output$pvalueMatrix<=indEdgePvalueThresh))
  computeTrueSimilOutput = CoINcIDE_output$computeTrueSimilOutput
  pvalueMatrix = CoINcIDE_output$pvalueMatrix
  clustIndexMatrix = CoINcIDE_output$clustIndexMatrix

  finalEdgeInfo <- assignFinalEdges(computeTrueSimilOutput=computeTrueSimilOutput,pvalueMatrix=pvalueMatrix,indEdgePvalueThresh=indEdgePvalueThresh,
                                    meanEdgePairPvalueThresh=meanEdgePairPvalueThresh,
                                    minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                    fractFeatIntersectThresh=fractFeatIntersectThresh,numFeatIntersectThresh=numFeatIntersectThresh ,
                                    clustSizeThresh=clustSizeThresh, clustSizeFractThresh= clustSizeFractThresh,saveDir=saveDir,fileTag="CoINcIDE_edges_",
                                    clustIndexMatrix,minFractNN =minFractNN,minNumEdgesForCluster=minNumEdgesForCluster
  )
  
  meanMeanMetric <- append(meanMeanMetric, as.vector(finalEdgeInfo$meanMeanMetricMatrix))
  meanFractNN <- append(meanFractNN, as.vector(finalEdgeInfo$meanFractNNmatrix))                            
  #NOW count # false edges
  numFalseEdges_afterThresh[[i]] <- nrow(finalEdgeInfo$filterEdgeOutput$edgeMatrix)
  #can be helpful to see if say always the small clusters get dropped out, etc.
  falseEdgeLists[[i]] <- finalEdgeInfo$filterEdgeOutput$edgeMatrix
  
  if(nrow(finalEdgeInfo$filterEdgeOutput$edgeMatrix)>0){
    
    commInfo <- findCommunities(edgeMatrix=finalEdgeInfo$filterEdgeOutput$edgeMatrix,edgeWeightMatrix=finalEdgeInfo$filterEdgeOutput$edgeWeightMatrix,
                              clustIndexMatrix=CoINcIDE_output$clustIndexMatrix,fileTag=experimentName,
                              saveDir=saveDir,minNumUniqueStudiesPerCommunity=minNumUniqueStudiesPerCommunity,experimentName=experimentName,
                              commMethod=commMethod,plotSimilEdgeWeight = FALSE,
                              makePlots=FALSE,saveGraphData=FALSE,plotToScreen=FALSE,findCommWithWeights=findCommWithWeights,
                              fractEdgesInVsOutComm=fractEdgesInVsOutComm,fractEdgesInVsOutEdge=fractEdgesInVsOutEdge)
  
  
  numFalseEdges_afterGN[[i]] <- nrow(commInfo$edgeDF)
  numFalseClusters_afterGN[[i]] <- nrow(commInfo$attrDF)
  numFalseMetaclusters_afterGN[[i]] <- commInfo$numCommunities
  numFalseDatasetsInNetwork[[i]] <- length(unique(commInfo$attrDF$studyNum))
  
  }else{
    
    numFalseEdges_afterGN[[i]] <- 0
    numFalseMetaclusters_afterGN[[i]] <- 0
    numFalseClusters_afterGN[[i]]  <- 0
    numFalseDatasetsInNetwork[[i]] <- 0
    
  }

  
  }
  
  #just use last clustIndexMatrixList - # datasets will always be the same.
  uniqueDatasets <- unique(clustIndexMatrix[,2])
  uniqueClusters <- unique(clustIndexMatrix[,1])
  #each cluster could be assigned to 1 cluster in all other (N-1) datasets
  #N-1 so don't count an edges to clusters in own dataset.
  #2 clusters per edge (hence the divide by 2)
  #this number may end in 0.5 (i.e. with the arrangment of clusters and datasets-
  #so just round up to be safe, but should only cause minor rounding errors.
  totalPossibleEdges <- ceiling((length(uniqueClusters)*(length(uniqueDatasets)-1))/2)
  meanMeanMetricQuantilesBeforeThresh <-  quantile(meanMeanMetric,na.rm=TRUE,probs=seq(0,1,.05))
  meanFractNNQuantilesBeforeThresh <-  quantile(meanFractNN,na.rm=TRUE,probs=seq(0,1,.05))
  numFalseEdges_afterThresh <- mean(unlist(numFalseEdges_afterThresh))
  numFalseEdges_afterGN <- mean(unlist(numFalseEdges_afterGN))
  numFalseIndPvaluesBeforeOtherThresh <- mean(unlist(numFalseIndPvaluesBeforeOtherThresh))
  numFalseDatasetsInNetwork <- mean(unlist( numFalseDatasetsInNetwork))
  #must NOT count edges that are in the same dataset...how do this?


  #each cluster assigned to only one centroid in another dataset.
  #so # non-NA p-values = total # of possible p-values for this experiment.
  #this is actually the same as the total # of edges
  totalPossibleIndPvalues <- length(which(!is.na(CoINcIDE_output$pvalue)))
  FDR_indPvalue_beforeOtherThresh <- numFalseIndPvaluesBeforeOtherThresh/totalPossibleIndPvalues
  #an "edge" means two p-values in essence were kept.
  FDR_indPvalue_afterOtherThresh <-  (numFalseEdges_afterThresh*2)/totalPossibleIndPvalues
  FDR_indPvalue_afterGN <-(numFalseEdges_afterGN*2)/totalPossibleIndPvalues
  
  FDR_edgesAfterAllThresh <- numFalseEdges_afterThresh/totalPossibleEdges
  FDR_edgesAfterGN <- numFalseEdges_afterGN/totalPossibleEdges

  message("FDR returned is a fraction, not a percentage.")
  
  FDR_results <- list(numFalseIndPvaluesBeforeOtherThresh=numFalseIndPvaluesBeforeOtherThresh,numFalseEdges_ind =numFalseEdges_ind,
                      numFalseEdges_afterThres=numFalseEdges_afterThresh,numFalseEdges_afterGN=numFalseEdges_afterGN,
                      numFalseMetaclusters_afterGN=numFalseMetaclusters_afterGN,numFalseClusters_afterGN=numFalseClusters_afterGN ,
                      falseEdgeLists=falseEdgeLists,meanMeanMetricQuantilesBeforeThresh=meanMeanMetricQuantilesBeforeThresh,
                      meanFractNNQuantilesBeforeThresh=meanFractNNQuantilesBeforeThresh,numFalseEdges_afterThresh=numFalseEdges_afterThresh,
                      numFalseEdges_afterGN=numFalseEdges_afterGN,numFalseIndPvaluesBeforeOtherThresh=numFalseIndPvaluesBeforeOtherThresh,
                      FDR_indPvalue_beforeOtherThresh=FDR_indPvalue_beforeOtherThresh,FDR_indPvalue_afterOtherThresh=FDR_indPvalue_afterOtherThresh,
                      FDR_indPvalue_afterGN=FDR_indPvalue_afterGN,FDR_edgesAfterGN=FDR_edgesAfterGN,FDR_edgesAfterAllThresh=FDR_edgesAfterAllThresh,numFalseDatasetsInNetwork=numFalseDatasetsInNetwork)

  
  return(FDR_results)
  
}

#edge output is edgeDF in findCommunities
networkLeaveOutAnalysis <- function(finalNodeMatrix, origEdgeMatrix,
                                    origEdgeWeightsMatrix,finalEdgeMarix,fractLeaveOut=0.1,
                                    numIter=100,commMethod="edgeBetween",
                                    findCommWithWeights=TRUE,clustIndexMatrix,messageSaveDir="./"){
  
  
  #commMethod=c("edgeBetween","fastGreedy","walktrap","eigenvector","optimal","spinglass","multilevel")
  
  if(nrow(origEdgeWeightsMatrix) != nrow(origEdgeMatrix)){
    
    stop("Indices of your original weight and edge matrix aren't matching up.")
    
  }
  #careful! the append function does something weird with the factors, so just remove the levels here
  finalNodeMatrix$clust <- as.character(finalNodeMatrix$clust)
  finalNodeMatrix$igraph_id <- as.character(finalNodeMatrix$igraph_id)
  finalNodeMatrix$community <- as.character(finalNodeMatrix$community)
  origEdgeMatrix[,1] <- as.character(origEdgeMatrix[,1])
  origEdgeMatrix[,2] <- as.character(origEdgeMatrix[,2])
  finalEdgeMatrix[,1] <- as.character(finalEdgeMatrix[,1])
  finalEdgeMatrix[,2] <- as.character(finalEdgeMatrix[,2])


  #replace igraph clust/node ids with original clust ids
  
  if(nrow(finalEdgeMatrix) != length(na.omit(match(finalEdgeMatrix[,1],finalNodeMatrix$igraph_id))) 
     || nrow(finalEdgeMatrix) != length(na.omit(match(finalEdgeMatrix[,2],finalNodeMatrix$igraph_id))) ){
    
    stop("Indexing error going from igraph clust id to global clust id.")
    
  }
  
  finalEdgeMatrix[,1] <- finalNodeMatrix[na.omit(match(finalEdgeMatrix[,1],finalNodeMatrix$igraph_id)), "clust"]
  #do for both columns of clusters/nodes
  finalEdgeMatrix[,2] <- finalNodeMatrix[na.omit(match(finalEdgeMatrix[,2],finalNodeMatrix$igraph_id)), "clust"]
  
  
  
  uniqueNodes <- unique(append(finalEdgeMatrix[,1],finalEdgeMatrix[,2]))
  uniqueOrigNodes <-  unique(append(origEdgeMatrix[,1],origEdgeMatrix[,2]))
  numNodesRemove <- round(fractLeaveOut*length(uniqueNodes))
  #remove nodes (not edges) that are not in the final edge matrix
  #if Girvan-Newman removes edges connected to nodes that stayed in network:
  #want to keep these edges as input for the leave-X-out analysis.
  prunedNodes <- setdiff(uniqueOrigNodes,uniqueNodes)
  
  if(length(setdiff(uniqueNodes,uniqueOrigNodes)) !=0 ){
    
    stop("Error: somehow you have extra clusters (nodes) that were not in your original edge matrix.")
    
  }
  
  removeNodes <- function(nodesToRemove, edgeMatrix,edgeWeightsMatrix,
                          clustIndexMatrix){
  
  if(length(nodesToRemove)>0){
    
    #can have duplicated nodes, so must loop and not use match which only finds first match
    for(n in 1:length(nodesToRemove)){
      
    #look in first column
      removeRows <- which(edgeMatrix[,1]==nodesToRemove[n])
      if(length(removeRows)>0){
        
        edgeMatrix <- edgeMatrix[-removeRows, ,drop=FALSE]
        edgeWeightsMatrix <- edgeWeightsMatrix[-removeRows, ,drop=FALSE]
       
      }
      
      }
    
      
      if(nrow(edgeMatrix)>0){
      #now look in second column
      #can have duplicated nodes, so must loop and not use match which only finds first match
      for(n in 1:length(nodesToRemove)){
        
        removeRows <- which(edgeMatrix[,2]==nodesToRemove[n])
      
      if(length(removeRows)>0){
        
        edgeMatrix <- edgeMatrix[-removeRows, ,drop=FALSE]
        edgeWeightsMatrix <- edgeWeightsMatrix[-removeRows, ,drop=FALSE]
        
      }
  
      
     }
    }

  }
  return(list(edgeMatrix=edgeMatrix,edgeWeightsMatrix=edgeWeightsMatrix))
  
  }
  
  tmp <-  removeNodes(nodesToRemove=prunedNodes, edgeMatrix=origEdgeMatrix, 
                      edgeWeightsMatrix=origEdgeWeightsMatrix)
  origEdgeMatrix <- tmp$edgeMatrix
  origEdgeWeightsMatrix <- tmp$edgeWeightsMatrix
  rm(list="tmp")
  
  
  #create a membership matrix by community (meta-cluster)
  membershipMatrix <- matrix(data=0,ncol=length(unique(finalNodeMatrix$clust)), nrow=length(unique(finalNodeMatrix$clust)),
                             dimnames = list(unique(finalNodeMatrix$clust),unique(finalNodeMatrix$clust)))
  #want to keep diagonals zero - don't count this as a cluster matching to itself
  #membershipMatrix[diag(membershipMatrix)] <- 1
  commNums <- unique(finalNodeMatrix$community)
  
  for(c in 1:length(commNums)){
    
    indices <- na.omit(match(finalNodeMatrix$clust[which(finalNodeMatrix$community==commNums[c])],rownames(membershipMatrix)))
    if(length(indices)>0){
    
    membershipMatrix[indices,indices] <- 1
    
    }
    
  }
  
  #want to keep diagonals zero - don't count this as a cluster matching to itself
  diag(membershipMatrix) <- 0
  commNumsOrig <- commNums
  #for membership comparisons
  output_numMisMatches <- matrix(data=NA,ncol=(length(commNumsOrig)+1),nrow=numIter,
                          dimnames=list(c(1:numIter),c("total",commNumsOrig)))
  output_fractMisMatches  <- matrix(data=NA,ncol=(length(commNumsOrig)+1),nrow=numIter,
                                               dimnames=list(c(1:numIter),c("total",commNumsOrig)))
  
  for(i in 1:numIter){
    
    nodesRemove <- sample(x=uniqueNodes,size=numNodesRemove,replace=FALSE)
    #remove these nodes - aren't included in this analysis
    tmp_membershipMatrix <- membershipMatrix[-na.omit(match(nodesRemove,rownames(membershipMatrix))),
                                             -na.omit(match(nodesRemove,colnames(membershipMatrix))) ]
    tmp <-  removeNodes(nodesToRemove=nodesRemove, edgeMatrix=origEdgeMatrix, edgeWeightsMatrix=origEdgeWeightsMatrix)
    tmp_origEdgeMatrix <- tmp$edgeMatrix
    tmp_origEdgeWeightsMatrix <- tmp$edgeWeightsMatrix
    rm(list="tmp")

    commInfo <- findCommunities(edgeMatrix=tmp_origEdgeMatrix,edgeWeightMatrix= tmp_origEdgeWeightsMatrix ,clustIndexMatrix=clustIndexMatrix,fileTag="CoINcIDE_LeaveXOut_",
                                            saveDir=messageSaveDir,minNumUniqueStudiesPerCommunity=1,experimentName="exp",
                                            commMethod=commMethod,
                                            makePlots=FALSE,plotToScreen=FALSE,saveGraphData=TRUE,nodeFontSize=.7,nodePlotSize=10,
                                            findCommWithWeights=findCommWithWeights, plotSimilEdgeWeight = FALSE,fractEdgesInVsOutComm=0,
                                            fractEdgesInVsOutEdge=0)


 
    #create a membership matrix by community (meta-cluster)
    #keep same number of nodes as input nodes/tmp matrix.  want nodes that were removed by
    #community detection algorithms to still be identified here.
    new_membershipMatrix <- matrix(data=0,ncol=ncol(tmp_membershipMatrix), nrow=nrow(tmp_membershipMatrix),
                               dimnames = list(colnames(tmp_membershipMatrix),colnames(tmp_membershipMatrix)))
    #keep diagonal zero: don't want to count this in the Jaccard metric.
    #new_membershipMatrix[diag(new_membershipMatrix)] <- 1
    commInfo$attrDF$clust <- as.character(commInfo$attrDF$clust)
    commInfo$attrDF$community <- as.character(commInfo$attrDF$community)
    commNums <- unique(as.character(commInfo$attrDF$community))
    
    if(length(commNums)>0){
    for(c in 1:length(commNums)){
      #Note: if a node was removed/pruned by Girvan-Newman: then will just be zero for all values in here.
      indices <- na.omit(match(commInfo$attrDF$clust[which(commInfo$attrDF$community==commNums[c])],rownames(new_membershipMatrix)))
      if(length(indices)>0){
        
      new_membershipMatrix[indices,indices] <- 1
      
      }
      
    }
    }
    #    #keep diagonal zero: don't want to count this in the Jaccard metric.
    diag(new_membershipMatrix) <- 0
    
    #just to be on the safe side: make sure row, column indices of tmp, new membership matrix match up.
    tmp_membershipMatrix<-  tmp_membershipMatrix[na.omit(match(rownames(new_membershipMatrix),rownames(tmp_membershipMatrix))),
                                                 na.omit(match(rownames(new_membershipMatrix),rownames(tmp_membershipMatrix)))]
    if(nrow(new_membershipMatrix)!=nrow(tmp_membershipMatrix) ||
       !all(rownames(tmp_membershipMatrix)==rownames(new_membershipMatrix))
       ){
      
      stop("Dimensions of tmp membership and tmp original membership matrix not matching up.")
      
    }

    #only need the lower or upper part of symmetric matrix.
    #if difference is zero: either both started with zero, or both started with 1.
    #don't use diagonal - doesn't really count
    numMisMatches <- sum(abs(tmp_membershipMatrix[upper.tri(tmp_membershipMatrix,diag=FALSE)]-
                                                   new_membershipMatrix[upper.tri(new_membershipMatrix,diag=FALSE)]))
    
    fractMisMatches <- numMisMatches/length(upper.tri(tmp_membershipMatrix,diag=FALSE))
    

    output_numMisMatches[i,1] <- numMisMatches
      output_fractMisMatches[i,1] <- fractMisMatches 
    
  #use ORIGINAL communities. This is our ground truth level.
  for(c in 1:(length(commNumsOrig))){
    #split up by community.
    indices <- na.omit(match(finalNodeMatrix$clust[which(finalNodeMatrix$community==commNumsOrig[c])],rownames(tmp_membershipMatrix)))
    
    if(length(indices)>0){
      
    tmp1 <- new_membershipMatrix[indices,indices]
    tmp2 <- tmp_membershipMatrix[indices,indices]
    numMisMatches <- sum(abs(tmp2[upper.tri(tmp2,diag=FALSE)]-
                               tmp1[upper.tri(tmp1,diag=FALSE)]))
    fractMisMatches <- numMisMatches/length(upper.tri(tmp2,diag=FALSE))
    
    output_numMisMatches[i,(c+1)] <- numMisMatches
    output_fractMisMatches[i,(c+1)] <- fractMisMatches 
    
    }

    }
    #these will not necessarily be equal: only looking at mismatches within a single community
    #output_numMisMatches[,c(1:5)]!= output_numMisMatches[,1])
    #end of iteration i
  }
  
  warning("An NA means that randomly no nodes were selected from that original community in the leaveXOut analysis.")
  return(list(output_numMisMatches=output_numMisMatches,
              output_fractMisMatches=output_fractMisMatches,
              numNodesRemove=numNodesRemove,numStartNodes=uniqueNodes))
}
#plot
plotLeaveXOutAnalysis <- function(leaveXOutResults,
                                  experimentName,saveDir){
  
  for(i in 1:length(leaveXOutResults)){
    
    tmp <- data.frame(colnames(leaveXOutResults[[i]]$output_numMisMatches),names(leaveXOutResults)[i],colMeans(leaveXOutResults[[i]]$output_fractMisMatches,na.rm=TRUE),
                      colSds(leaveXOutResults[[i]]$output_fractMisMatches,na.rm=TRUE),
                      colMeans(leaveXOutResults[[i]]$output_numMisMatches,na.rm=TRUE),colSds(leaveXOutResults[[i]]$output_numMisMatches,na.rm=TRUE))
    
    colnames(tmp) <- c("comm","fractLO","fractMisMatch","sd_fractMM","numMisMatch","sd_numMM")
    
    if(i ==1){
      
      masterDF <- tmp
      
    }else{
      
      masterDF <- rbind(masterDF,tmp)
      
    }
    
    
  }
  
  for(i in 1:length(leaveXOutResults)){
    
    tmp <- data.frame(cbind(leaveXOutResults[[i]]$output_fractMisMatches,rep.int(names(leaveXOutResults)[i],times=nrow(leaveXOutResults[[i]]$output_fractMisMatches))))

    
    colnames(tmp) <- c(colnames(leaveXOutResults[[i]]$output_fractMisMatches),"fractLO")
    
    
    if(i ==1){
      
      fullMasterDF <- tmp
      
    }else{
      
      fullMasterDF <- rbind(fullMasterDF,tmp)
      
    }
 
  }
  

  #melt: want one column for community id
  #NEED reshape 2 library
  #number of levels are different bc different communities: may get a warning
  #http://stackoverflow.com/questions/25688897/reshape2-melt-warning-messages
  fullMasterDF <- melt(fullMasterDF,id.vars=c("fractLO"))
  colnames(fullMasterDF)[2] <- "comm"
  colnames(fullMasterDF)[3] <- "fractMisMatch"
  #make "Total" in caps
  levels(fullMasterDF[,2])[which( levels(fullMasterDF[,2])=="total")] <- "Total"
  #remove NAs - means all members of that meta-cluster were removed
  fullMasterDF <- fullMasterDF[-which(is.na(fullMasterDF[,3])), ]
  fullMasterDF[,3] <- round(as.numeric(fullMasterDF[,3]),digits=2)
  
LO_plot_average <- ggplot(data =masterDF,aes(x=fractLO,y=fractMisMatch,group=comm))  + geom_line(size=.4)+ geom_point(aes(shape=factor(comm),size=5))+
    labs(title = "",y="Mismatch Fraction",x="Fraction Left Out")+
    theme(panel.background = element_rect(fill='white', colour='black')) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.text.x = element_text(colour = "black",size=18),axis.title.x = element_text(colour = "black",size=20,vjust=0),
          axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1),
          plot.title=element_text(colour="black",size=22,vjust=1,hjust=.4))+
    theme(legend.title=element_text(colour="black",size=18),legend.text=element_text(colour="black",size=10))+
    labs(shape="Meta-cluster")
  
  
  options(bitmapType="cairo")
  png(filename=paste0(saveDir,"/LO_plot_",experimentName,"_",Sys.Date(),".png"),width=1500,height=1600,res=200);
  
  plot(LO_plot_average);
  
  dev.off();

  #quartile function taken from here: http://stackoverflow.com/questions/17319487/median-and-quartile-on-violin-plots-in-ggplot2
  median.quartile <- function(x){
    out <- quantile(x, probs = c(0.25,0.5,0.75))
    names(out) <- c("ymin","y","ymax")
    return(out) 
  }
  
  
  LO_plot_violin <- ggplot(data =  fullMasterDF ,aes(x=fractLO,y=fractMisMatch))  +
    geom_violin()+facet_grid(comm~.)+stat_summary(fun.y=median.quartile,geom='point')+
    labs(title = "",y="Mismatch Fraction",x="Fraction Left Out")+
   # theme(panel.background = element_rect(fill='white', colour='black')) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
   theme(axis.text.x = element_text(colour = "black",size=16),axis.title.x = element_text(colour = "black",size=20,vjust=0),
          axis.text.y = element_text(colour = "black",size=16),axis.title.y = element_text(colour = "black",size=20,vjust=1))

  
  options(bitmapType="cairo")
  png(filename=paste0(saveDir,"/LO_plot_violin_",experimentName,"_",Sys.Date(),".png"),width=1500,height=1600,res=200);
  
  plot(LO_plot_violin);
  
  dev.off();
  

#  dodge <- position_dodge(width=0.9)
#  limits <- aes(ymax = masterDF$fractMisMatch + masterDF$sd_fractMM, ymin=masterDF$fractMisMatch - masterDF$sd_fractMM)
  
  #eh...this is not working? no plotting along the points.
#   LO_plot_withErrorBars <-     LO_plot <- ggplot(data =masterDF,aes(x=fractLO,y=fractMisMatch,group=comm))  + geom_line(size=.4)+ geom_point(aes(shape=factor(comm),size=5))+
#     labs(title = "",y="Mismatch Fraction",x="Fraction Left Out")+
#     scale_color_manual(values=colorCodes)+
#     theme(panel.background = element_rect(fill='white', colour='black')) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#     theme(axis.text.x = element_text(colour = "black",size=18),axis.title.x = element_text(colour = "black",size=20,vjust=0),
#           axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1),
#           plot.title=element_text(colour="black",size=22,vjust=1,hjust=.4))+
#     theme(legend.title=element_text(colour="black",size=18),legend.text=element_text(colour="black",size=10))+
#     labs(shape="Fraction Left Out")+ 
#     geom_bar(position=dodge) + geom_errorbar(limits, position=dodge, width=0.25)
  
  return(list(masterDF=masterDF,fullMasterDF=fullMasterDF,LO_plot_average=LO_plot_average,LO_plot_violin=LO_plot_violin))
}

#edge output is edgeDF in findCommunities
networkVaryMinSimilAnalysis <- function(finalNodeMatrix, origEdgeMatrix,
                                    origEdgeWeightsMatrix,finalEdgeMarix,
                                    minSimilThreshVector=c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),
                                     commMethod="edgeBetween",saveDir="./",
                                    findCommWithWeights=TRUE,clustIndexMatrix,
                                    makePlots=TRUE,experimentName= "test"){
  
  summaryFile <- paste0(saveDir,"/similVaryTestsNetworkStats_",experimentName,".txt")
  #commMethod=c("edgeBetween","fastGreedy","walktrap","eigenvector","optimal","spinglass","multilevel")
  
  if(nrow(origEdgeWeightsMatrix) != nrow(origEdgeMatrix)){
    
    stop("Indices of your original weight and edge matrix aren't matching up.")
    
  }
  #careful! the append function does something weird with the factors, so just remove the levels here
  finalNodeMatrix$clust <- as.character(finalNodeMatrix$clust)
  finalNodeMatrix$igraph_id <- as.character(finalNodeMatrix$igraph_id)
  finalNodeMatrix$community <- as.character(finalNodeMatrix$community)
  origEdgeMatrix[,1] <- as.character(origEdgeMatrix[,1])
  origEdgeMatrix[,2] <- as.character(origEdgeMatrix[,2])
  finalEdgeMatrix[,1] <- as.character(finalEdgeMatrix[,1])
  finalEdgeMatrix[,2] <- as.character(finalEdgeMatrix[,2])
  
  
  #replace igraph clust/node ids with original clust ids
  
  if(nrow(finalEdgeMatrix) != length(na.omit(match(finalEdgeMatrix[,1],finalNodeMatrix$igraph_id))) 
     || nrow(finalEdgeMatrix) != length(na.omit(match(finalEdgeMatrix[,2],finalNodeMatrix$igraph_id))) ){
    
    stop("Indexing error going from igraph clust id to global clust id.")
    
  }
  
  finalEdgeMatrix[,1] <- finalNodeMatrix[na.omit(match(finalEdgeMatrix[,1],finalNodeMatrix$igraph_id)), "clust"]
  #do for both columns of clusters/nodes
  finalEdgeMatrix[,2] <- finalNodeMatrix[na.omit(match(finalEdgeMatrix[,2],finalNodeMatrix$igraph_id)), "clust"]
  
  
  
  uniqueNodes <- unique(append(finalEdgeMatrix[,1],finalEdgeMatrix[,2]))
  uniqueOrigNodes <-  unique(append(origEdgeMatrix[,1],origEdgeMatrix[,2]))
  #remove nodes (not edges) that are not in the final edge matrix
  #if Girvan-Newman removes edges connected to nodes that stayed in network:
  #want to keep these edges as input for the leave-X-out analysis.
  prunedNodes <- setdiff(uniqueOrigNodes,uniqueNodes)
  
  if(length(setdiff(uniqueNodes,uniqueOrigNodes)) !=0 ){
    
    stop("Error: somehow you have extra clusters (nodes) that were not in your original edge matrix.")
    
  }
  
  removeNodes <- function(nodesToRemove, edgeMatrix,edgeWeightsMatrix,
                          clustIndexMatrix){
    
    if(length(nodesToRemove)>0){
      
      #can have duplicated nodes, so must loop and not use match which only finds first match
      for(n in 1:length(nodesToRemove)){
        
        #look in first column
        removeRows <- which(edgeMatrix[,1]==nodesToRemove[n])
        if(length(removeRows)>0){
          
          edgeMatrix <- edgeMatrix[-removeRows, ,drop=FALSE]
          edgeWeightsMatrix <- edgeWeightsMatrix[-removeRows, ,drop=FALSE]
          
        }
        
      }
      
      
      if(nrow(edgeMatrix)>0){
        #now look in second column
        #can have duplicated nodes, so must loop and not use match which only finds first match
        for(n in 1:length(nodesToRemove)){
          
          removeRows <- which(edgeMatrix[,2]==nodesToRemove[n])
          
          if(length(removeRows)>0){
            
            edgeMatrix <- edgeMatrix[-removeRows, ,drop=FALSE]
            edgeWeightsMatrix <- edgeWeightsMatrix[-removeRows, ,drop=FALSE]
            
          }
          
          
        }
      }
      
    }
    return(list(edgeMatrix=edgeMatrix,edgeWeightsMatrix=edgeWeightsMatrix))
    
  }
  
  tmp <-  removeNodes(nodesToRemove=prunedNodes, edgeMatrix=origEdgeMatrix, 
                      edgeWeightsMatrix=origEdgeWeightsMatrix)
  origEdgeMatrix <- tmp$edgeMatrix
  origEdgeWeightsMatrix <- tmp$edgeWeightsMatrix
  rm(list="tmp")
  
  
  #create a membership matrix by community (meta-cluster)
  membershipMatrix <- matrix(data=0,ncol=length(unique(finalNodeMatrix$clust)), nrow=length(unique(finalNodeMatrix$clust)),
                             dimnames = list(unique(finalNodeMatrix$clust),unique(finalNodeMatrix$clust)))
  #want to keep diagonals zero - don't count this as a cluster matching to itself
  #membershipMatrix[diag(membershipMatrix)] <- 1
  commNums <- unique(finalNodeMatrix$community)
  
  for(c in 1:length(commNums)){
    
    indices <- na.omit(match(finalNodeMatrix$clust[which(finalNodeMatrix$community==commNums[c])],rownames(membershipMatrix)))
    if(length(indices)>0){
      
      membershipMatrix[indices,indices] <- 1
      
    }
    
  }
  
  #want to keep diagonals zero - don't count this as a cluster matching to itself
  diag(membershipMatrix) <- 0
  commNumsOrig <- commNums
  #for membership comparisons
  output_numMisMatches <- matrix(data=NA,ncol=(length(commNumsOrig)+1),nrow=length(minSimilThreshVector),
                                 dimnames=list(minSimilThreshVector,c("total",commNumsOrig)))
  output_fractMisMatches  <- matrix(data=NA,ncol=(length(commNumsOrig)+1),nrow=length(minSimilThreshVector),
                                    dimnames=list(minSimilThreshVector,c("total",commNumsOrig)))
  numNodesRemove <- list()
  commEdgeInfo <- list()
  
  for(i in 1:length(minSimilThreshVector)){
    
    #remove the nodes that are below the simil thresh.
    #weights matrix: same row indices as edge matrix.
    edgesRemove <- which(origEdgeWeightsMatrix[,"simil"] <= minSimilThreshVector[i])
    nodesRemove <- unique(append(origEdgeMatrix[edgesRemove,1], origEdgeMatrix[edgesRemove,1]))
    numNodesRemove[[i]] <- length(nodesRemove)
    #remove these nodes - aren't included in this analysis
    tmp_membershipMatrix <- membershipMatrix[-na.omit(match(nodesRemove,rownames(membershipMatrix))),
                                             -na.omit(match(nodesRemove,colnames(membershipMatrix))) ]
    tmp <-  removeNodes(nodesToRemove=nodesRemove, edgeMatrix=origEdgeMatrix, edgeWeightsMatrix=origEdgeWeightsMatrix)
    tmp_origEdgeMatrix <- tmp$edgeMatrix
    tmp_origEdgeWeightsMatrix <- tmp$edgeWeightsMatrix
    rm(list="tmp")
    
    commInfo <- findCommunities(edgeMatrix=tmp_origEdgeMatrix,edgeWeightMatrix= tmp_origEdgeWeightsMatrix ,clustIndexMatrix=clustIndexMatrix,fileTag="CoINcIDE_LeaveXOut_",
                                saveDir=saveDir,minNumUniqueStudiesPerCommunity=1,experimentName=paste0(experimentName,"_minSimil",minSimilThreshVector[i]),
                                commMethod=commMethod,
                                makePlots=makePlots,plotToScreen=FALSE,saveGraphData=TRUE,nodeFontSize=.7,nodePlotSize=10,
                                findCommWithWeights=findCommWithWeights, plotSimilEdgeWeight = FALSE,fractEdgesInVsOutComm=0,
                                fractEdgesInVsOutEdge=0)
    
    commEdgeInfo[[i]] <- commInfo$finalCommEdgeInfo
    
    if(i>1){
      
    textOut <- capture.output(commInfo$finalCommEdgeInfo)
    cat(paste0("Final community stats for simil ",minSimilThreshVector[i],":\n"),textOut,sep="\n",append=TRUE,file=summaryFile)
    cat(textOut,sep="\n",
        append=TRUE,file=summaryFile)
    
    }else{
      
      textOut <- capture.output(commInfo$finalCommEdgeInfo)
      cat(paste0("Final community stats for simil ",minSimilThreshVector[i],":\n"),textOut,sep="\n",append=FALSE,file=summaryFile)
      cat(textOut,sep="\n",
          append=FALSE,file=summaryFile)
      
      
    }
    
    if(commInfo$numCommunities>0){
      
    networkStats <- advancedNetworkPlots(communityMembership=commInfo,
                                         brewPal = "Set3",
                                         saveDir=saveDir,clustIndexMatrix=clustIndexMatrix,
                                         plotToScreen=FALSE,experimentName=paste0(experimentName,"_minSimil",minSimilThreshVector[i]),
                                         plotEdgeWeight=FALSE)$network_stats
    
      
    textOut <- capture.output(networkStats)
    cat(paste0("\nMore network stats for simil ",minSimilThreshVector[i],":\n"),textOut,sep="\n",
        file=summaryFile,append=TRUE)

    
    }
    #create a membership matrix by community (meta-cluster)
    #keep same number of nodes as input nodes/tmp matrix.  want nodes that were removed by
    #community detection algorithms to still be identified here.
    new_membershipMatrix <- matrix(data=0,ncol=ncol(tmp_membershipMatrix), nrow=nrow(tmp_membershipMatrix),
                                   dimnames = list(colnames(tmp_membershipMatrix),colnames(tmp_membershipMatrix)))
    #keep diagonal zero: don't want to count this in the Jaccard metric.
    #new_membershipMatrix[diag(new_membershipMatrix)] <- 1
    commInfo$attrDF$clust <- as.character(commInfo$attrDF$clust)
    commInfo$attrDF$community <- as.character(commInfo$attrDF$community)
    commNums <- unique(as.character(commInfo$attrDF$community))
    
    if(length(commNums)>0){
      for(c in 1:length(commNums)){
        #Note: if a node was removed/pruned by Girvan-Newman: then will just be zero for all values in here.
        indices <- na.omit(match(commInfo$attrDF$clust[which(commInfo$attrDF$community==commNums[c])],rownames(new_membershipMatrix)))
        if(length(indices)>0){
          
          new_membershipMatrix[indices,indices] <- 1
          
        }
        
      }
    }
    #    #keep diagonal zero: don't want to count this in the Jaccard metric.
    diag(new_membershipMatrix) <- 0
    
    #just to be on the safe side: make sure row, column indices of tmp, new membership matrix match up.
    tmp_membershipMatrix<-  tmp_membershipMatrix[na.omit(match(rownames(new_membershipMatrix),rownames(tmp_membershipMatrix))),
                                                 na.omit(match(rownames(new_membershipMatrix),rownames(tmp_membershipMatrix)))]
    if(nrow(new_membershipMatrix)!=nrow(tmp_membershipMatrix) ||
       !all(rownames(tmp_membershipMatrix)==rownames(new_membershipMatrix))
    ){
      
      stop("Dimensions of tmp membership and tmp original membership matrix not matching up.")
      
    }
    
    #only need the lower or upper part of symmetric matrix.
    #if difference is zero: either both started with zero, or both started with 1.
    #don't use diagonal - doesn't really count
    numMisMatches <- sum(abs(tmp_membershipMatrix[upper.tri(tmp_membershipMatrix,diag=FALSE)]-
                               new_membershipMatrix[upper.tri(new_membershipMatrix,diag=FALSE)]))
    
    fractMisMatches <- numMisMatches/length(upper.tri(tmp_membershipMatrix,diag=FALSE))
    
    
    output_numMisMatches[i,1] <- numMisMatches
    output_fractMisMatches[i,1] <- fractMisMatches 
    
    #use ORIGINAL communities. This is our ground truth level.
    for(c in 1:(length(commNumsOrig))){
      #split up by community.
      indices <- na.omit(match(finalNodeMatrix$clust[which(finalNodeMatrix$community==commNumsOrig[c])],rownames(tmp_membershipMatrix)))
      
      if(length(indices)>0){
        
        tmp1 <- new_membershipMatrix[indices,indices]
        tmp2 <- tmp_membershipMatrix[indices,indices]
        numMisMatches <- sum(abs(tmp2[upper.tri(tmp2,diag=FALSE)]-
                                   tmp1[upper.tri(tmp1,diag=FALSE)]))
        fractMisMatches <- numMisMatches/length(upper.tri(tmp2,diag=FALSE))
        
        output_numMisMatches[i,(c+1)] <- numMisMatches
        output_fractMisMatches[i,(c+1)] <- fractMisMatches 
        
      }
      
    }
    #these will not necessarily be equal: only looking at mismatches within a single community
    #output_numMisMatches[,c(1:5)]!= output_numMisMatches[,1])
    #end of iteration i
  }
  
  warning("An NA means no nodes passed the simil threshold for that original community.")
  
  output <- list(output_numMisMatches=output_numMisMatches,
                 output_fractMisMatches=output_fractMisMatches,
                 numNodesRemove=numNodesRemove,numStartNodes=uniqueNodes)

  
  output$output_fractMisMatches[which( output$output_numMisMatches[,"total"]==0), ] <- 0
  output$output_numMisMatches[which( output$output_numMisMatches[,"total"]==0), ] <- 0
  
  for(i in 1:ncol(output$output_fractMisMatches)){
    

    tmp <- data.frame(minSimilThreshVector,output$output_fractMisMatches[,i],output$output_numMisMatches[,i],
                    rep(colnames(output$output_fractMisMatches)[i]),length(minSimilThreshVector))
    
    colnames(tmp) <- c("minSimil","fractMisMatch","numMisMatch","comm")
    
    if(i ==1){
      
      masterDF <- tmp
      
    }else{
      
      masterDF <- rbind(masterDF,tmp)
      
    }
    
    
  }
  
  simil_plot <- ggplot(data =masterDF,aes(x=minSimil,y=fractMisMatch,point=comm))  + 
    geom_line(size=.4)+
    geom_point(aes(shape=factor(comm),size=5))+
    labs(title = "",y="Mismatch Fraction",x="Mininum Mean Similarity")+
    theme(panel.background = element_rect(fill='white', colour='black')) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.text.x = element_text(colour = "black",size=18),axis.title.x = element_text(colour = "black",size=20,vjust=0),
          axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1),
          plot.title=element_text(colour="black",size=22,vjust=1,hjust=.4))+
    theme(legend.title=element_text(colour="black",size=18),legend.text=element_text(colour="black",size=10))+
    labs(shape="Meta-cluster")
  
  
  options(bitmapType="cairo")
  png(filename=paste0(saveDir,"/varySimilMismatchFract_plot_",experimentName,"_",Sys.Date(),".png"),width=1500,height=1600,res=200);
  
  plot(simil_plot);
  
  dev.off();
  message("Summary file is :",summaryFile)
   return(list(output_numMisMatches=output_numMisMatches,commEdgeInfo=commEdgeInfo,
              output_fractMisMatches=output_fractMisMatches,
              numNodesRemove=numNodesRemove,numStartNodes=uniqueNodes,simil_plot=simil_plot))
}
# #plot
# plotLeaveXOutAnalysis <- function(leaveXOutResults,
#                                   experimentName,saveDir){
#   
#   for(i in 1:length(leaveXOutResults)){
#     
#     tmp <- data.frame(colnames(leaveXOutResults[[i]]$output_numMisMatches),names(leaveXOutResults)[i],colMeans(leaveXOutResults[[i]]$output_fractMisMatches,na.rm=TRUE),
#                       colSds(leaveXOutResults[[i]]$output_fractMisMatches,na.rm=TRUE),
#                       colMeans(leaveXOutResults[[i]]$output_numMisMatches,na.rm=TRUE),colSds(leaveXOutResults[[i]]$output_numMisMatches,na.rm=TRUE))
#     
#     colnames(tmp) <- c("comm","fractLO","fractMisMatch","sd_fractMM","numMisMatch","sd_numMM")
#     
#     if(i ==1){
#       
#       masterDF <- tmp
#       
#     }else{
#       
#       masterDF <- rbind(masterDF,tmp)
#       
#     }
#     
#     
#   }
#   LO_plot <- ggplot(data =masterDF,aes(x=fractLO,y=fractMisMatch,group=comm))  + geom_line(size=.4)+ geom_point(aes(shape=factor(comm),size=5))+
#     labs(title = "",y="Mismatch Fraction",x="Fraction Left Out")+
#     scale_color_manual(values=colorCodes)+
#     theme(panel.background = element_rect(fill='white', colour='black')) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#     theme(axis.text.x = element_text(colour = "black",size=18),axis.title.x = element_text(colour = "black",size=20,vjust=0),
#           axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1),
#           plot.title=element_text(colour="black",size=22,vjust=1,hjust=.4))+
#     theme(legend.title=element_text(colour="black",size=18),legend.text=element_text(colour="black",size=10))+
#     labs(shape="Meta-cluster")
#   
#   
#   options(bitmapType="cairo")
#   png(filename=paste0(saveDir,"/LO_plot_",experimentName,"_",Sys.Date(),".png"),width=1500,height=1600,res=200);
#   
#   plot(LO_plot);
#   
#   dev.off();
#   #  dodge <- position_dodge(width=0.9)
#   #  limits <- aes(ymax = masterDF$fractMisMatch + masterDF$sd_fractMM, ymin=masterDF$fractMisMatch - masterDF$sd_fractMM)
#   
#   #eh...this is not working? no plotting along the points.
#   #   LO_plot_withErrorBars <-     LO_plot <- ggplot(data =masterDF,aes(x=fractLO,y=fractMisMatch,group=comm))  + geom_line(size=.4)+ geom_point(aes(shape=factor(comm),size=5))+
#   #     labs(title = "",y="Mismatch Fraction",x="Fraction Left Out")+
#   #     scale_color_manual(values=colorCodes)+
#   #     theme(panel.background = element_rect(fill='white', colour='black')) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   #     theme(axis.text.x = element_text(colour = "black",size=18),axis.title.x = element_text(colour = "black",size=20,vjust=0),
#   #           axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1),
#   #           plot.title=element_text(colour="black",size=22,vjust=1,hjust=.4))+
#   #     theme(legend.title=element_text(colour="black",size=18),legend.text=element_text(colour="black",size=10))+
#   #     labs(shape="Fraction Left Out")+ 
#   #     geom_bar(position=dodge) + geom_errorbar(limits, position=dodge, width=0.25)
#   
#   return(list(masterDF=masterDF,LO_plot=LO_plot))
# }
# 
# 
# 
