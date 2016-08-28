#library("sva")
#library("matrixStats")
batchNormalization <- function(countsMatrixNoNANoDup,outcomesAndCovariates,MinInBatch=4,combatModelFactorName=NULL,pvalueThresh=.05,batchColName="batch",outputFile="combatoutput.txt"){
  
  outputFile <- paste0(outputFile,"_",Sys.Date(),".txt");

    batchData <- outcomesAndCovariates[,batchColName,drop=FALSE]
    rownames(batchData) <- rownames(outcomesAndCovariates)
  
  batchNormCounts <- batchCorrection_MolecularData(GEN_Data=countsMatrixNoNANoDup,BatchData=batchData,MinInBatch=MinInBatch,
                                                                  combatModelFactorName=combatModelFactorName,pvalueThresh=pvalueThresh,batchColName=batchColName,
                                                                 outputFile=outputFile);
  
  return(batchNormCounts);
  
  
}

#ASSUMPTION: rownames of batchdata = colnames of matrix data
#rownames of batch data are column/sample IDs on GEN_ata
CheckBatchEffect <-function(GEN_Data,BatchData,batchColName) {
  # GEN_Data ;;;
  # Barch
  
  # PCA analysis
  # alternatively use fast.prcomp from package gmodels, but tests do not show this is faster
  PCAanalysis=prcomp(t(GEN_Data))
  PCdata=PCAanalysis$x
  
  batchEffectPlot <- (PCdata[,1]~as.numeric(BatchData[,batchColName]));
  
  if (length(unique(BatchData[,batchColName]))>1) {
    
    #AOV: anaysis of variance
    tmp=aov(PCdata[,1]~as.numeric(BatchData[,batchColName]));   
    
    return(list(Pvalue=summary(tmp)[[1]][["Pr(>F)"]][[1]],PCA=PCdata,BatchData=BatchData,batchEffectPlot=batchEffectPlot,model=tmp));
    
  } else {
    
    return(-1);
    
  }
}


#TO DO: describe mod.
BatchCorrection <- function(GEN_Data,BatchData,mod,batchColName,outputFile="./combatOutput.txt") {
  
  #KP: NEED THIS? JUST DID THIS IN WRAPPER FUNCTION? guess didn't do the second half...
  WithBatchSamples=is.element(colnames(GEN_Data),rownames(BatchData));
  
  if (length(which(WithBatchSamples==FALSE))>0) {
    
    GEN_Data=GEN_Data[,-which(WithBatchSamples==FALSE)];
    
  }
  
  # select only the batch data that is present in the current data set, remove others (remember, the batch data is for all of TCGA)
  PresentSamples <- is.element(rownames(BatchData),colnames(GEN_Data))
  
  if (sum(PresentSamples) != length(colnames(GEN_Data))) {
    
    BatchDataSelected=BatchData[-which(PresentSamples==FALSE), ,drop=FALSE]
    
  }else{
    
    
    BatchDataSelected <- BatchData;
  }
  
  BatchDataSelected[,batchColName] <- factor(BatchDataSelected[,batchColName])
 
  
  # reordening samples (not really necessary as Combat does this too)
  order <- match(colnames(GEN_Data),rownames(BatchDataSelected));
  BatchDataSelected <- BatchDataSelected[order, ,drop=FALSE];
  
  BatchDataSelected[,batchColName]<- factor(BatchDataSelected[,batchColName])
  
  # running combat
  #CombatResults=ComBat_NooutputFiles(GEN_Data,BatchDataSelected)
  #COME BACK: try filtering out low varying genes.
  #s long as don't put in numCovs argument, all variables will be treated as factors
  
  #final processing step: remove genes whose variance is near zero within any batch (not just globally)
  rowIndicesRemove <- c()
  for(b in 1:length(unique(BatchDataSelected[,batchColName]))){
    
    rowIndicesRemove <- append(rowIndicesRemove,which(rowVars(GEN_Data[,
                                                               which(BatchDataSelected[,batchColName]==unique(BatchDataSelected[,batchColName])[b])])<=.001))
    
  }
  
  if(length(rowIndicesRemove)>0){
    
    GEN_Data <- GEN_Data[-rowIndicesRemove, ]
    
  }
  
  #par.prior=TRUE indicates parametric adjustments will be used
  #need a numeric batch variable.
  CombatResults <- ComBat(dat=GEN_Data, batch=BatchDataSelected[,batchColName], mod=mod, par.prior=TRUE, prior.plots=FALSE);
  
  #GEN_Data_Corrected=CombatResults[,-1]
  #class(GEN_Data_Corrected) <- "numeric"
  return(CombatResults)
}

#end/final function you call.
batchCorrection_MolecularData <- function(GEN_Data,BatchData,MinInBatch,combatModelFactorName=NULL,pvalueThresh=.05,batchColName="bcr_batch",outputFile="./combatoutput.txt") {
  
  cat("running batch checks. Remember, they don't need to be in the same order",
      "\n but the row names of your batch data and column names of your expression match must contain the same entities.\n")
  # Remove samples with batch number 0
  #shouldn't this be column 1??? KP
  #IF there are none...this returns a zero matrix!
  #GEN_Data = GEN_Data[-which(BatchData[,batchColName]==0),]
  
  if(length(which(BatchData[,batchColName]==0))>0){
    
    message("Removing samples with batch ID=0")
    GEN_Data <- GEN_Data[,-which(BatchData[,batchColName]==0), drop=FALSE]
    BatchDataSelected <- BatchData[-which(BatchData[,batchColName]==0), , drop=FALSE]
    
  }else{
    
    BatchDataSelected <- BatchData;
  }
  
  # Remove samples with batch number NA
  if(any(is.na((BatchData[,batchColName])))){
    
    message("Removing samples with NA batch ID")
    GEN_Data <- GEN_Data[,-which(is.na(BatchData[,batchColName])), drop=FALSE]
    BatchDataSelected <- BatchData[-which(is.na(BatchData[,batchColName])), , drop=FALSE]
    
  }else{
    
    BatchDataSelected <- BatchData;
  }
  # remove batches that are too small     
  #MinInBatch=5
  #my colnames are NOT the GEN_data!
  
  PresentSamples <- is.element(rownames(BatchData),colnames(GEN_Data));
  
  if( sum(PresentSamples) != length(colnames(GEN_Data)) ) {
    
    cat("some of your batch data sample names weren't in the counts or expression matrix. Removinging these from the batch data.\n",append=TRUE,file=outputFile)
    BatchDataSelected <- BatchDataSelected[-which(PresentSamples==FALSE),]
    
    if(dim(as.matrix(BatchDataSelected))[1]==0){
      
      stop("none of the row names of your batch records were found in the column names of your counts or expression matrix! ")
    }
  }
  
  #COME BACK: also remove covariates that really don't vary much.
  BatchDataSelected[,batchColName] <- factor(BatchDataSelected[,batchColName]);
  
  NrPerBatch <- table(BatchDataSelected[,batchColName]);
  
  SmallBatches <- NrPerBatch<MinInBatch;
  
  BatchesToBeRemoved <- names(SmallBatches)[which(SmallBatches==TRUE)];
  
  SamplesToBeRemoved <- rownames(BatchDataSelected[which(BatchDataSelected[,batchColName] %in% BatchesToBeRemoved),]);
  
  if (length(colnames(GEN_Data))-length(which(colnames(GEN_Data) %in% SamplesToBeRemoved)) >5) { # just checking if we have enough samples after removing the too small batches
    
    if (length(which(colnames(GEN_Data) %in% SamplesToBeRemoved))>0) {
      
      cat("\n Removing",length(which(colnames(GEN_Data) %in% SamplesToBeRemoved)),"samples because their batches are too small.\n",append=TRUE,file=outputFile)
      
      GEN_Data <- GEN_Data[,-which(colnames(GEN_Data) %in% SamplesToBeRemoved)]
      #CHECK KP: also update batch data!
      BatchDataSelected <- BatchDataSelected[-which(BatchDataSelected[,batchColName] %in% BatchesToBeRemoved),];
      
    }          
    # batch correction with Combat, incorporate check for only 1 batch
    
    #now make sure batch effect row IDs and GEN_Data IDs are in order
    Order <- match(colnames(GEN_Data),rownames(BatchDataSelected));
    BatchDataSelected <- BatchDataSelected[Order, , drop=FALSE]
    colnames(BatchDataSelected) <- batchColName
    
    cat("\n original num samples was ",dim(BatchData)[1], " after removing batches with number of samples fewer than ",MinInBatch, " it's ",dim(BatchDataSelected)[1],"\n",append=TRUE,file=outputFile)
    #now check for batch effect.
    BatchCheck <- CheckBatchEffect(GEN_Data,BatchDataSelected,batchColName=batchColName);
    
    if (is.list(BatchCheck)) {
      
      #IF NEED BATCH CHECK: p-value thresh.
      beforeBatchEffectPlot <- BatchCheck$batchEffectPlot;
      beforePvalue <- BatchCheck$Pvalue;
      
      if(beforePvalue < pvalueThresh){
        
        #set up your model
        
        #Make ALL covariates factors!!! COME BACK.
        
        
        if(is.null(combatModelFactorName)){
          
          #
          mod <- NULL;
          #the result below is the same: with all the same outcomes, will give zero covariates.
          #mod=model.matrix(~ 1,data=as.data.frame(BatchDataSelected,stringsAsFactors=TRUE););
          
        }else{
          
          #did our pruning remove this variable?
          if(length(which(colnames(BatchDataSelected)==combatModelFactorName)) == 0){
            
            warning("your combat factor was pruned out, so using the NULL model.")
            
          }else{
            #data must be  d.f. for ~
            #we need BatchDataSelected to be a matrix later: so only make it a dataframe inside the model.
            mod<- model.matrix(~ combatModelFactorName,data=as.data.frame(BatchDataSelected,stringsAsFactors=TRUE));
            
          }
        }
        
        
        #ComBat returns a data frame
        GEN_Data_Corrected <- as.matrix(BatchCorrection(GEN_Data=GEN_Data,BatchData=BatchDataSelected,mod=mod,batchColName=batchColName));
        
        if(dim(GEN_Data_Corrected)[2] != dim(BatchDataSelected)[1]) {
          
          warning("Combat filtered out some samples. may affect your analyses");
          Order <- match(colnames(GEN_Data_Corrected),rownames(BatchDataSelected));
          BatchDataSelected <- BatchDataSelected[Order,]
          
        }
        #need matrices for this function
        BatchCheck <- CheckBatchEffect(GEN_Data_Corrected,BatchDataSelected,batchColName=batchColName)
        afterBatchEffectPlot <- BatchCheck$batchEffectPlot;
        afterPvalue <- BatchCheck$Pvalue;
        
        cat("\n batch effect detected. p-value on PCA test before was ",beforePvalue, " after combat correction it's ",afterPvalue," \n",append=TRUE,file=outputFile);
        
        return(list(GEN_Data_Corrected=GEN_Data_Corrected,BatchData=BatchDataSelected,afterBatchEffectPlot=afterBatchEffectPlot,beforeBatchEffectPlot=beforeBatchEffectPlot,beforePvalue=beforePvalue,afterPvalue=afterPvalue));
        
      }else{
        
        cat("\n no batch effect detected. Before p-value was ",beforePvalue," \n")
        cat("\n no batch effect detected. Before p-value was ",beforePvalue," \n",append=TRUE,file=outputFile);
        
        return(list(GEN_Data=GEN_Data,BatchData=BatchDataSelected,beforeBatchEffectPlot=beforeBatchEffectPlot,beforePvalue=beforePvalue));
      }
      
      
    } else {
      cat("\n only one batch after cleanup. Can't do batch correction.")
      cat("\n Only one batch, no batch correction possible.\n",append=TRUE,file=outputFile)
      return(NA)
    }
    
  } else {
    cat("\n The number of samples becomes to small, no batch correction possible.\n",append=TRUE,file=outputFile)
    return(NA)
  }
}
