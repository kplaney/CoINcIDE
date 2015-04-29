library("GSEABase")
library("plyr")
library("Biobase")
library("fdrtool")
library("rmeta")
library("RCurl")
library("RDAVIDWebService")
library("limma")
library("biomaRt")
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

#NOTE: groupings need names
#ovarian: use recurrence_status as binary variable
#for ovarian: not a big difference if neo or adjuvant
survivalAnalysis <- function(saveNameTag = "survivalAnalysis",outcomesVarBinary="os_event",outcomesVarCont = "os_days_to_event",
                             sampleClustCommPhenoData,CutoffPointYears=5,uniquePatientID="unique_patient_ID",
                             groupingTerm="community",addRX=FALSE){
  

  #only take samples with the groupingTerm you're looking at.
  sampleClustCommPhenoData <- sampleClustCommPhenoData[which(!is.na(sampleClustCommPhenoData[, groupingTerm])), ]
  #remove samples with NA values.
  groupings <- sampleClustCommPhenoData[, groupingTerm]
  
  message("Groupings total counts:\n")
  cat(table(groupings),"\n")
  
  groupingsBinary <- groupings[which(!is.na(sampleClustCommPhenoData[,outcomesVarBinary]))];
  
  
  message("Groupings total counts for binary outcomes variable:\n")
  cat(table(groupingsBinary),"\n")
  
  outcomesData <-  sampleClustCommPhenoData
  #outcomesData <-  outcomesData[which(!is.na( outcomesData[,outcomesVarBinary])),];
  
  #keep samples with NA days to event for now?
  groupingsCont <- as.numeric(as.factor(groupings[which(!is.na(outcomesData[,outcomesVarCont]))]))
  
  message("Groupings total counts for binary outcomes variable:\n")
  cat(table(groupingsCont),"\n")
  
  #outcomesData <- outcomesData[which(!is.na(outcomesData[,outcomesVarCont])),];
  
  #if binary is character string categories: make it a factor first, then numeric,
  #otherwise coxph function will throw errors.
  outcomesDataShort <- data.frame(as.numeric(as.factor(outcomesData[,outcomesVarBinary])),outcomesData[,outcomesVarCont],
  );
  
  #sometimes the names are duplicated across studies - remove this line
  #rownames(outcomesDataShort ) <- outcomesData[,uniquePatientID];
  colnames(outcomesDataShort) <- c("Censoring","TimeToLastContactOrEvent","groupings")
  
  nonCensoredTerm=1
  censoredTerm=0
  Survival <- outcomesDataShort
  #creating the survival objects with the time and censoring variables
  OverallSurvival <- Surv(Survival$TimeToLastContactOrEvent,Survival$Censoring==nonCensoredTerm);
  #KP:
  #creating a survival object cutoff at a certain point
  CutoffPoint <- CutoffPointYears*365;
  CutoffSamples=Survival$TimeToLastContactOrEvent>CutoffPoint & !is.na(Survival$TimeToLastContactOrEvent)
  SurvivalCutoff=Survival
  SurvivalCutoff$TimeToLastContactOrEvent[CutoffSamples]=CutoffPoint
  SurvivalCutoff$Censoring[CutoffSamples]=censoredTerm
  #"Surv" creates a survival object. really for binary outcomes data.
  OverallSurvivalCutoff=Surv(SurvivalCutoff$TimeToLastContactOrEvent,SurvivalCutoff$Censoring==nonCensoredTerm)
  
  results<-array(dim=c(1,13))
  #dataMatrix <- data.matrix(data);
  ########################################################
  # cox survival
  ########################################################
  
  ########################################################
  
  # STEP 1: Cox model on complete survival data + kaplan meier
  #Fits a Cox proportional hazards regression model.
  #***add in Rx info??
  #plot.survfit: can plot this data after using survfit.coxph
  #ex for plotting survival data: > leukemia.surv <- survfit(Surv(time, status) ~ x, data = aml) 
  #> plot(leukemia.surv, lty = 2:3) 
  #here: run test <- cox.zph(coxfit), plot.cox.zph(test)
  #COME BACK: include therapy too coxph( Surv(time, status) ~ therapy + ReceptorA + ReceptorB , data= sample.data) 
  #hmmm...not working?? is this because most patients died?
  #if groups are ill-balanced, may get a warnings: http://stats.stackexchange.com/questions/66591/coxph-ran-out-of-iterations-and-did-not-converge
  coxfit=coxph(OverallSurvival~groupings, data=Survival)
  plot.cox.zph(cox.zph(coxfit))
  
  tmp=summary(coxfit)
  #results[i,1]=tmp$logtest[3]  
  
  results[1,1]=tmp$waldtest[3]
  results[1,2]=exp(coxfit$coefficients)
  results[1,3]=tmp$conf.int[3]
  results[1,4]=tmp$conf.int[4]
  results[1,5]=tmp$coefficients[4] # the zscore
  
  #kaplan meier...is this correct?
  kmfit=survdiff(OverallSurvival ~ groupings)
  results[1,6]= 1 - pchisq(kmfit$chisq, length(kmfit$n) - 1)
  
  #calculate the sign of the survival relationship
  #do we save this somethere?? in 7?
  mfit=survfit(OverallSurvival ~ groupings)
  plot(mfit)
  mfit.summary=summary(mfit)
  
  ########################################################
  # STEP 2: Cox model on survival data with CUTOFF + kaplan meier
  
  coxfit=coxph(OverallSurvivalCutoff~groupings, data= SurvivalCutoff)
  
  tmp=summary(coxfit)
  results[,8]=tmp$waldtest[3]
  #what is this??
  results[,9]=exp(coxfit$coefficients)
  results[,10]=tmp$conf.int[3]
  results[,11]=tmp$conf.int[4]
  results[,12]=tmp$coefficients[4] # the zscore
  
  #The survfit function from the survival package computes the Kaplan-Meier estimator for truncated and/or censored data
  #survfit: This function creates survival curves from either a formula (e.g. the Kaplan-Meier), a previously fitted Cox model, or a previously fitted accelerated failure time model.
  kmfit=survdiff(OverallSurvivalCutoff ~ groupings)
  results[,13]= 1 - pchisq(kmfit$chisq, length(kmfit$n) - 1)  
  
  
  ########################################################
  # writing to file
  ########################################################
  #write.table(results,Args[7],sep="\t")
  
  #ADD FDR TESTING!
  ColumnLabels= c('Overall-WaldTest', 'Overall-HR','Overall-HRlower' ,'Overall-HRupper' ,'Overall-Zscore', 'Overall-KMtest' ,
                  'Cutoff-WaldTest', 'Cutoff-HR' ,'Cutoff-HRlower' ,'Cutoff-HRupper' ,'Cutoff-Zscore', 'Cutoff-KMtest' );
  
  colnames(results) <- ColumnLabels;
  rownames(results) <- colnames(dataMatrix);
  
  
  
  ########################################################
  #STEP 3:  covariate modeling....
  #COME BACK.
  
  if(addRX){
    
    
    #remove samples that don't have Rx info.
    sampleClustCommPhenoData <- sampleClustCommPhenoData[which(!is.na(sampleClustCommPhenoData[, c("chemo","estrogen_ther","her2_ther")])), ]
    #remove samples with NA values.
    groupings <- sampleClustCommPhenoData[, groupingTerm]
    groupings <- groupings[which(!is.na(sampleClustCommPhenoData[,outcomesVarBinary]))];
    outcomesData <- sampleClustCommPhenoData[which(!is.na(sampleClustCommPhenoData[,outcomesVarBinary])),];
    
    #keep samples with NA days to event for now?
    groupings <- as.numeric(as.factor(groupings[which(!is.na(outcomesData[,outcomesVarCont]))]))
    outcomesData <- outcomesData[which(!is.na(outcomesData[,outcomesVarCont])),];
    
    #if binary is character string categories: make it a factor first, then numeric,
    #otherwise coxph function will throw errors.
    outcomesDataShort <- data.frame(as.numeric(as.factor(outcomesData[,outcomesVarBinary])),outcomesData[,outcomesVarCont],
    );
    
    #sometimes the names are duplicated across studies - remove this line
    #rownames(outcomesDataShort ) <- outcomesData[,uniquePatientID];
    colnames(outcomesDataShort) <- c("Censoring","TimeToLastContactOrEvent","groupings")
    
    nonCensoredTerm=1
    censoredTerm=0
    Survival <- outcomesDataShort
    #creating the survival objects with the time and censoring variables
    OverallSurvival <- Surv(Survival$TimeToLastContactOrEvent,Survival$Censoring==nonCensoredTerm);
    #KP:
    #creating a survival object cutoff at a certain point
    CutoffPoint <- CutoffPointYears*365;
    CutoffSamples=Survival$TimeToLastContactOrEvent>CutoffPoint & !is.na(Survival$TimeToLastContactOrEvent)
    SurvivalCutoff=Survival
    SurvivalCutoff$TimeToLastContactOrEvent[CutoffSamples]=CutoffPoint
    SurvivalCutoff$Censoring[CutoffSamples]=censoredTerm
    #"Surv" creates a survival object. really for binary outcomes data.
    OverallSurvivalCutoff=Surv(SurvivalCutoff$TimeToLastContactOrEvent,SurvivalCutoff$Censoring==nonCensoredTerm)
    
    results_rx<-array(dim=c(1,13))
    #dataMatrix <- data.matrix(data);
    ########################################################
    # cox survival
    ########################################################
    
    ########################################################
    
    # STEP 1: Cox model on complete survival data + kaplan meier
    #Fits a Cox proportional hazards regression model.
    #***add in Rx info??
    #plot.survfit: can plot this data after using survfit.coxph
    #ex for plotting survival data: > leukemia.surv <- survfit(Surv(time, status) ~ x, data = aml) 
    #> plot(leukemia.surv, lty = 2:3) 
    #here: run test <- cox.zph(coxfit), plot.cox.zph(test)
    #COME BACK: include therapy too coxph( Surv(time, status) ~ therapy + ReceptorA + ReceptorB , data= sample.data) 
    #hmmm...not working?? is this because most patients died?
    #if groups are ill-balanced, may get a warnings: http://stats.stackexchange.com/questions/66591/coxph-ran-out-of-iterations-and-did-not-converge
    coxfit=coxph(OverallSurvival~groupings, data=Survival)
    plot.cox.zph(cox.zph(coxfit))
    
    tmp=summary(coxfit)
    #results[i,1]=tmp$logtest[3]  
    
    results_rx[1,1]=tmp$waldtest[3]
    results_rx[1,2]=exp(coxfit$coefficients)
    results_rx[1,3]=tmp$conf.int[3]
    results_rx[1,4]=tmp$conf.int[4]
    results_rx[1,5]=tmp$coefficients[4] # the zscore
    
    #kaplan meier...is this correct?
    kmfit=survdiff(OverallSurvival ~ groupings)
    results_rx[1,6]= 1 - pchisq(kmfit$chisq, length(kmfit$n) - 1)
    
    #calculate the sign of the survival relationship
    #do we save this somethere?? in 7?
    mfit=survfit(OverallSurvival ~ groupings)
    plot(mfit)
    mfit.summary=summary(mfit)
    
    ########################################################
    # STEP 2: Cox model on survival data with CUTOFF + kaplan meier
    
    coxfit=coxph(OverallSurvivalCutoff~groupings, data= SurvivalCutoff)
    
    tmp=summary(coxfit)
    results_rx=tmp$waldtest[3]
    #what is this??
    results_rx[,9]=exp(coxfit$coefficients)
    results_rx[,10]=tmp$conf.int[3]
    results_rx[,11]=tmp$conf.int[4]
    results_rx[,12]=tmp$coefficients[4] # the zscore
    
    #The survfit function from the survival package computes the Kaplan-Meier estimator for truncated and/or censored data
    #survfit: This function creates survival curves from either a formula (e.g. the Kaplan-Meier), a previously fitted Cox model, or a previously fitted accelerated failure time model.
    kmfit=survdiff(OverallSurvivalCutoff ~ groupings)
    results_rx[,13]= 1 - pchisq(kmfit$chisq, length(kmfit$n) - 1)  
    
    
    ########################################################
    # writing to file
    ########################################################
    #write.table(results,Args[7],sep="\t")
    
    #ADD FDR TESTING!
    ColumnLabels= c('Overall-WaldTest', 'Overall-HR','Overall-HRlower' ,'Overall-HRupper' ,'Overall-Zscore', 'Overall-KMtest' ,
                    'Cutoff-WaldTest', 'Cutoff-HR' ,'Cutoff-HRlower' ,'Cutoff-HRupper' ,'Cutoff-Zscore', 'Cutoff-KMtest' );
    
    colnames(results) <- ColumnLabels;
    rownames(results) <- colnames(dataMatrix);
    
    
    
  }else{
    
    results_rx <- NA
  }
  #COME BACK: save separate Survival objects for rx and no rx.
  return(list(results=results,results_rx=results_rx,Survival=Survival,SurvivalCutoff=SurvivalCutoff))
  
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
                                        featureNames,minOtherClass=5){
  
  
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
          
          #get p-value; compute/update to q-value after loop through all genes
          wilcoxon_qvalue[g,s] <- wilcox.test(t(matrix1[featureNames[g], ]), t(matrix2[featureNames[g],]),alternative = c("two.sided"),
                                                    paired=FALSE)$p.value
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
      #for this study: fdr correct all of the genes.
      wilcoxon_qvalue[which(!is.na(wilcoxon_qvalue[,s])), s] <- 
        
        p.adjust(wilcoxon_qvalue[which(!is.na(wilcoxon_qvalue[ ,s])), s],method="fdr")
      
      
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
    #weighted mean; inverse weigthing by standard deviation
    ES <- meta.summaries(d=d, se=se, method=c("random"),
                         logscale=TRUE,
                         conf.level=0.95);
    
    summHedgeG_ES[g,c] <- ES$summary;
    summHedgeG_ES_se[g,c] <- ES$se.summary;
    #summarize p-value across all studies (can use Fisher's because studies are independent.)
    summWilcoxon_qvalue[g,c] <-  fishersMethod(wilcoxon_qvalue[g, which(!is.na(wilcoxon_qvalue[g,]))])
    
    #end of looping over genes
    }
    #save all community-specific data matrices.
    wilcoxon_qvalueList[[commNames[c]]] <- wilcoxon_qvalue
    hedgeGList[[commNames[c]]] <- hedgeG
    hedgeG_seList[[commNames[c]]] <- hedgeG_se
    
    #end of looping over communities
  }
  
  #correct across the number of meta-clusters testing.
  #could do apply I guess here...but just sticking for for loop format throughout the function:
  for(g in 1:length(featureNames)){
    
    summWilcoxon_qvalue[g,which(!is.na(summWilcoxon_qvalue[g, ]))] <- p.adjust(summWilcoxon_qvalue[g,which(!is.na(summWilcoxon_qvalue[g, ]))],method="fdr");
  
  
  }
  
  output <- list(summWilcoxon_qvalue=summWilcoxon_qvalue,summHedgeG_ES=summHedgeG_ES,
                 summHedgeG_ES_se=summHedgeG_ES_se,wilcoxon_qvalueList=wilcoxon_qvalueList,
                 hedgeGList=hedgeGList,hedgeG_seList=hedgeG_seList)
  
  return(output)
  
  #EOF
}

#thresh is greater than.
selectMetaclustSigGenes <- function(computeMetaclustEffectSizesOutput,qvalueThresh=.1,
                                    ESthresh=0){
  
  commNames <- colnames(computeMetaclustEffectSizesOutput$summHedgeG_ES)
  
  sigMetaclustGenesWilcox <- list()
  sigMetaclustGenesES_pos <- list()
  sigMetaclustGenes_pos <- list()
  sigMetaclustGenesES_neg <- list()
  sigMetaclustGenes_neg <- list()
  
  for(c in 1:length(commNames)){
    
    sigMetaclustGenesWilcox[[commNames[c]]] <- rownames(computeMetaclustEffectSizesOutput$summWilcoxon_qvalue[
      which(computeMetaclustEffectSizesOutput$summWilcoxon_qvalue[,commNames[c]]<=qvalueThresh),commNames[c],drop=FALSE])
    
    sigMetaclustGenesES_pos[[commNames[c]]] <- rownames(computeMetaclustEffectSizesOutput$summHedgeG_ES[
      which(computeMetaclustEffectSizesOutput$summHedgeG_ES[,commNames[c]]>=ESthresh),commNames[c],drop=FALSE])
    
    sigMetaclustGenes_pos[[commNames[c]]] <- intersect(sigMetaclustGenesWilcox[[commNames[c]]],
                                                       sigMetaclustGenesES_pos[[commNames[c]]])
    
    
    sigMetaclustGenesES_neg[[commNames[c]]] <- rownames(computeMetaclustEffectSizesOutput$summHedgeG_ES[
      which(computeMetaclustEffectSizesOutput$summHedgeG_ES[,commNames[c]]<=ESthresh),commNames[c],drop=FALSE])
    
    sigMetaclustGenes_neg[[commNames[c]]] <- intersect(sigMetaclustGenesWilcox[[commNames[c]]],
                                                       sigMetaclustGenesES_neg[[commNames[c]]])
    
      
  }
  
  output <- list( sigMetaclustGenesWilcox= sigMetaclustGenesWilcox, sigMetaclustGenesES_pos= sigMetaclustGenesES_pos,
                  sigMetaclustGenes_pos=sigMetaclustGenes_pos,sigMetaclustGenesES_neg=sigMetaclustGenesES_neg,
                  sigMetaclustGenes_neg= sigMetaclustGenes_neg)
  
  return(output)
  
}

summarizePosMetaclustGenes <- function(selectMetaclustSigGenesOut,computeMetaclustEffectSizesOutput){
  
  featureNames <- c()
  commNames <- colnames(computeMetaclustEffectSizesOutput$summWilcoxon_qvalue)
  
  for(c in 1:length(commNames)){
    
    featureNames <- append(featureNames,selectMetaclustSigGenesOut$sigMetaclustGenes_pos[[commNames[c]]])
    
  }
  
  featureNames <- unique(featureNames)
  ESMatrix <- matrix(data=NA,ncol=length(commNames),nrow=length(featureNames),dimnames=list(featureNames,commNames))
  
  #all featureNames will be in this matrix, even if NA.
  for(c in 1:length(commNames)){
    
    ESMatrix[selectMetaclustSigGenesOut$sigMetaclustGenes_pos[[commNames[c]]],commNames[c]] <-
      computeMetaclustEffectSizesOutput$summHedgeG_ES[selectMetaclustSigGenesOut$sigMetaclustGenes_pos[[commNames[c]]], commNames[c]]
    
  }
  
  if(any(is.na(rowMeans(ESMatrix,na.rm=TRUE)))){
    
    stop("Not computing ESMatrix correctly; getting all NAs in some rows.")
  }
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
  
  
  

  sigFunct <- functAnnotChart[which(functAnnotChart[,"Benjamini"]<=qvalueThreshFunct) ,c("Category","Term")]
  
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
# mapGeneIDs <- function(geneticFeatures,bioMartSource="ensembl",bioMartDataset="hsapiens_gene_ensembl",
#                                           geneticFeatureSource = 'hgnc_symbol',returnOnlyOneSymbolPerFeature=TRUE){
#   
#   library("biomaRt")
#   mart <- useMart(bioMartSource,dataset=bioMartDataset)
#   
#   #can't have ' in a genetic features name for biomaRt to recognize it.
#   if(length(dim(geneticFeatures))>0){
#     
#     #take first column. otherwise gsub does weird stuff.
#     geneticFeatures <- geneticFeatures[,1];
#     
#   }
#   
#   if(length(geneFeatures)>450){
#     
#     stop("Your gene list is probably too long for biomaRt in R; try looping over chunks of the list.")
#     
#   }
#   #can't have ' in a genetic features name for biomaRt to recognize it.
#   geneticFeatures<- gsub("\'","",geneticFeatures);
#   #oftentimes this period is added to ensembl IDs - but it's not in biomart.
#   geneticFeatures <- strsplit2(geneticFeatures,"\\.")[,1];
#   
#   #some code alters geneticFeatures - we'll need the original list for some purposes later on down in the code
#   #(returning 1 row per original genetic feature source item.)
#   geneticFeaturesFull <- geneticFeatures;
#   warning("duplicate entries for a certain BiomaRt source may be returned.")
#   warning("if you're feeding in HGNC gene symbols,make sure these are the latest updated version, otherwise Biomart may not recognize them.")
#   #listAttributes(mart)
#   
#   #use BM, as BMList is slower and loops over each value you put in. BM looks at all values more efficiently.
#    
#       linkedData <- getBM(c(geneticFeatureSource,'hgnc_symbol','ensembl_gene_id','start_position','end_position','band','chromosome_name','gene_biotype','entrezgene'),
#                           filters=c(geneticFeatureSource),values=geneticFeatures,
#                           mart=mart,uniqueRows=TRUE); 
#   
# 
#   if(dim(linkedData)[1] < length(geneticFeatures)){
#     
#     warning("couldn't find a biomaRt record for all features.");
#     
#   }
#   
#   if(returnOnlyOneSymbolPerFeature){
#     
#     if(length(geneticFeaturesFull) != length(unique(geneticFeaturesFull))){
#       
#       warning("Some of your genetic feature inputs were duplicated. So you'll get duplicated output rows by choosing to have only 1 symbol returned per 
#               (unique) feature input.")
#       
#     }
#     
#     if(length(dim(linkedData))==0){
#       #so can grab column ID if it's a vector (shouldn't occur bc I put the feature source along with the gene symbol - so min 2 columns.)
#       linkedData <- as.matrix(linkedData);
#       
#     }
#     
#     featureInfo <- linkedData;
#     linkedData <- matrix(data=NA,nrow=length(geneticFeaturesFull),ncol=dim(linkedData)[2]);
#     colnames(linkedData) <- colnames(featureInfo);
#     
#     #looking for original column we fed in - want this many rows and to match on this.
#     #will only return 1 row per unique geneticFeatureSource
#     colIndex <- which(colnames(linkedData)==geneticFeatureSource)[1];
#     
#     for(g in 1:length(geneticFeaturesFull)){
#       
#       #may have duplicated entries per one feature.
#       #technically, if just did match function, would only return first match.
#       #but I liked doing this stuff deliberately so I can remember it better later...
#       matchIDs <- which(featureInfo[,colIndex ]==geneticFeaturesFull[g]);
#       
#       #may have returned NULL.
#       if(!is.na(all(matchIDs))){
#         
#         #just take the first match.
#         linkedData[g,] <- unlist(featureInfo[matchIDs[1],]);
#         
#       }else{
#         
#         cat("nothing for ", geneticFeaturesFull[g]);
#         #data is already NA. no need for this.
#         #linkedData[g,] <- NA;
#         
#       }
#       
#       #COME BACK: is this working OK?
#       if(length(matchIDs)>1){
#         warn <- TRUE;
#         cat("\n index",g," gene ",geneticFeaturesFull[g]," had multiple records returned.\n");
#         
#       }else{
#         warn <- FALSE;
#       }
#       
#       
#     }
#     
#     if(warn){
#       
#       warning("Found more than 1 gene symbol for some input features - printed them out above. The first symbol was chosen/kept each time - usually the best bet.")
#     }
#     
#   }
#   
#   
#   return(linkedData);
#   
#   }
# 
# #COME BACK: need to debug how compute this....
#  GSEA <- function(testGeneVector,refGeneLists=NULL,method=c("hypergeometric","fisher"),genomeSize,
#                   refGeneListDir="/home/data/MSigDB/GSEA_base_MSigDB_lists_merged.RData.gzip"){
#    
#    warning("\nThis GSEA test may overly conservative -Katie needs to make the genomeSize just all the genes across all of the refLists")
#    warning("\nThis code assumes that all of your genes in your test gene list and ref gene list are in the genome.");
#    
#    if(is.null(refGeneLists)){
#      
#    load(refGeneListDir)
#    #cat("\nUsing default MSigDB lists: MSigDB_onco_symbols, MSigDB_CanPath_symbols,MSigDB_TFT_symbols,MSigDB_immun_symbols,
#     #   and MSigDB_cancerNeigh_symbols.\n");
#    
#    }
#    
#   refGeneLists <- GSEA_base_MSigDB_lists_merged;
#   
#   pvalues <- list();
#   
#   for(r in 1:length(refGeneLists)){
#     
#     numIntersect <- length(intersect(testGeneVector,refGeneLists[[r]]));
#     
#     if(numIntersect>0){
#       
#       if(method=="hypergeometric"){
#         #we are taking the area under the curve to calculate a <= scenario.
#         #from 1 to the number of genes we've intersected - this is our vector of quantiles
#         #"representing the number of white balls drawn without replacement from an urn which contains both black and white balls"
#         #length(refGeneLists[[r]]) is the total number of balls we drew from the urn. so the urn has to be the entire genome here?
#         #tricky when filter down gene list first.  a larger genomeSize will make for more conservative p-values.
#         
#         #phyper(x, m, n, k) gives the probability of getting x or fewer
#         #phyper(numIntersect-1, length(geneList), genomeSize, length(refGeneLists[[r]]), lower.tail=FALSE)
#         #lower tail =FALSE same as 1-:
#         #what's the probability of have a larger tail than this?
#         pvalues[[r]] <- 1- phyper(numIntersect-1, length(testGeneVector), genomeSize, length(refGeneLists[[r]]));
#         
#         #this will also give the same answer:
#         #we want to know how probably it is by chance that we'd see MORE (or equal?) intersecting genes than we have.
#         #pvalue <- min(1-cumsum(dhyper(0:(numIntersect-1), length(testGeneVector), genomeSize-length(testGeneVector), length(refGeneLists[[r]])) ));
#       }else if(method=="fisher"){
#         
#         counts <-  matrix(data = c(numIntersect, length(refGeneLists[[r]])-numIntersect, length(testGeneVector),genomeSize-length(testGeneVector)),nrow = 2);
#         pvalues[[r]] <- fisher.test(counts,alternative="greater")$p.value;
#         
#       }else{
#         
#         stop("\nPlease input method=hypergeometric or method=fisher as the statistical method.");
#         
#       }
#       
#       
#       
#     }else{
#       
#       #set equal to one? no sure how to do this here....
#       pvalues[[r]] <- NA;
#       
#     }
#     
#   }
#   
#   nonNA_indices <- which(!is.na(unlist(unlist(pvalues))));
#   nonNA_pvalues <- unlist(pvalues)[which(!is.na(unlist(unlist(pvalues))))];
#   qvalues <- array(data=NA,dim=length(unlist(pvalues)));
#   #want to preserve NA values
#   #now do multiple hyp testing
#   qvalues[nonNA_indices] <- p.adjust(p=nonNA_pvalues,method="fdr");
#   
#   names(pvalues) <- names(refGeneLists)
#   names(qvalues) <- names(refGeneLists)
#   output <- list(qvalues=qvalues,pvalues=unlist(pvalues),refGeneListNames=names(refGeneLists));
#   
#   return(output);
#   
# }
# 
# # #for example: can create TCGA membership, overlay with ours.
# # visualizeTwoMetaclustMemberships <- function(memberMatrix1,memberMatrix2,returnSampleMemberMatrix()$fullMemberMatrix){
# #   
# #  totalColNames <- union(colnames(memberMatrix1),colnames(memberMatrix2))
# #  overlayMatrix <- matrix(data=NA,ncol=length(totalColNames),nrow=length(totalColNames),
# #                          dimnames=list(totalColNames,totalColNames))
# #  
# #  for(t in 1:length(totalColNames)){
# #    
# #    for(s in 1:length(totalColNames)){
# #    
# #      #have we computed this index yet?
# #      if(is.na(overlayMatrix[t,s])){
# #        
# #    if(!is.na(match(colnames(memberMatrix1),totalColNames[t]))){
# #      
# #      
# #      value1 <- memberMatrix1[totalColNames[t],totalColNames[s]] 
# #      
# #    }else{
# #      #want to be able to ID when a sample was simply not included in an analysis.
# #      value1  <- -2
# #      
# #    }
# #      
# #      if(!is.na(match(colnames(memberMatrix2),totalColNames[t]))){
# #        
# #        value2 <- memberMatrix2[totalColNames[t],totalColNames[s]]   
# #        
# #      }else{
# #        
# #        value2 <- -2
# #        
# #      }
# #      
# #          
# #          overlayMatrix[t,s] <- value1 + value 2
# #          
# #    #end of if NA      
# #    }
# #        
# #  
# #    }
# #    
# #  }
# #  
# #  if(any(is.na(overlayMatrix))){
# #    
# #    stop("Error: you are getting NAs in your overlayMatrix.")
# #    
# #  }
# #   return(overlayMatrix)
# #   #EOF
# # }