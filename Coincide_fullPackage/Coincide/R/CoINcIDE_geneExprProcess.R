#library("impute")
#library("Biobase")


createS4exprSet <- function(expr,phenoData,featureData,featureDataFieldName="gene"){
  
  #search for class.
  if((!missing(phenoData))&&(!missing(featureData))){
    
    
    exprSet <- new("ExpressionSet", assayData = assayDataNew(exprs=new("matrix")), phenoData = new("AnnotatedDataFrame"), featureData = new("AnnotatedDataFrame"), experimentData = new("MIAME"), annotation = character(0))
    
    exprSet@assayData <- assayDataNew(exprs = expr)
    
    rownames(featureData) <- rownames(expr)
    featureData <-  new("AnnotatedDataFrame", data=as.data.frame(featureData))
    exprSet@featureData <- featureData
    
    rownames(phenoData) <- colnames(expr)
    
    phenoData <-  new("AnnotatedDataFrame", data=as.data.frame(phenoData))
    exprSet@phenoData <- phenoData
    #you can everything in your "new" implementation of exprSet, or just add the expression data (must add these this in new() or else can't add it later!)
    
    #sample names must be the same
    #how you would access the expression data: exprSet@assayData$expr
    #add the other data
    #weird...the way below didn't work! arrayQualityMetric couldn't read these slots..so just added them directly above.
    #phenoData <- as.data.frame(phenoData, row.names = colnames(expr),stringsAsFactors=FALSE)
    #exprSet@featureData <- AnnotatedDataFrame(as.data.frame(featureData,stringsAsFactors=FALSE))
    #want to actually name the column in phenodata - or else can't label covariates in arrayQualityMetrics
    #exprSet@phenoData$phenoData <- AnnotatedDataFrame(phenoData)
    
    
    
  }else if((missing(phenoData))&&(!missing(featureData))){
    
    rownames(featureData) <- rownames(expr)
    featureData <-  new("AnnotatedDataFrame", data=as.data.frame(featureData))
    #featureData <- as.data.frame(featureData,stringsAsFactors=FALSE)
    #exprSet <- new("ExpressionSet", assayData = assayDataNew(exprs=as.matrix(expr)), experimentData = new("MIAME"), annotation = character(0),featureData=annotatedDataFrameFrom(assayData, byrow=TRUE))
    #how you would access the expression data: exprSet@assayData$expr
    #exprSet <- ExpressionSet(assayData = as.matrix(expr),featureData=annotatedDataFrameFrom(assayData, byrow=TRUE), experimentData = new("MIAME"), annotation = character(0))
    
    exprSet <- new("ExpressionSet", assayData = assayDataNew(exprs=new("matrix")), phenoData = new("AnnotatedDataFrame"), featureData = new("AnnotatedDataFrame"), experimentData = new("MIAME"), annotation = character(0))
    
    exprSet@assayData <- assayDataNew(exprs = expr)
    
    exprSet@featureData <- AnnotatedDataFrame(featData)
    
  }else if ((!missing(phenoData))&&(missing(featureData))){
    
    rownames(phenoData) <- colnames(expr)
    
    phenoData <-  new("AnnotatedDataFrame", data=as.data.frame(phenoData))
  
    exprSet <- new("ExpressionSet", assayData = assayDataNew(exprs=new("matrix")), phenoData = new("AnnotatedDataFrame"), featureData = new("AnnotatedDataFrame"), experimentData = new("MIAME"), annotation = character(0))
    
    exprSet@assayData <- assayDataNew(exprs = expr)
    
    exprSet@phenoData <- phenoData
    
    featData <- data.frame(rownames(expr))
    colnames(featData) <-  featureDataFieldName
    exprSet@featureData <- AnnotatedDataFrame(featData)
    
  }else{
  exprSet <- new("ExpressionSet", assayData = assayDataNew(exprs=new("matrix")), phenoData = new("AnnotatedDataFrame"), featureData = new("AnnotatedDataFrame"), experimentData = new("MIAME"), annotation = character(0))
    
    exprSet@assayData <- assayDataNew(exprs = expr)      
    featData <- data.frame(rownames(expr))
    colnames(featData) <-  featureDataFieldName
    exprSet@featureData <- AnnotatedDataFrame(featData)
    
  }
  
  cat("\nExpressionSet successfull made. Dimensions of expr in AssayData is ",dim(exprSet@assayData$exprs),"\n")
  return(exprSet)
}

#####
#some of these variables are a bit noisy - remove
#phenoData <- patientDataFull[,c(1:112)]
createExpressionSetList <- function(exprMatrixList,masterPhenoData,patientKey="GEO_GSMID",
                                    featureDataFieldName){
  
  ExpressionSetList <- list()
  
  for(d in 1:length(exprMatrixList)){
    
    
    if(!is.matrix(exprMatrixList[[d]]$expr)){
      
      warning("This code assumes your expression data in the expr slot is a matrix, 
              as this allows for duplicated gene row names.")
      
    }
    #make an expression object for each exprMatrixList
    
    #add study Name,unique ID for phenoData
    phenoData <- rep(names(exprMatrixList)[d],ncol(exprMatrixList[[d]]$expr))
    
    if(!all(is.null(exprMatrixList[[d]]$phenoData))){
      
      pheno <- exprMatrixList[[d]]$phenoData
      
      #row # is number of patients here.
      phenoData <- cbind(phenoData,pheno)
      
      #user provided a master data list instead - pull out the relevant patients
      #for this loop.
    }else if(!missing(masterPhenoData)){
      
      #usually numeric vs. character not important here (if string is actually all numbers.)
      pData <- masterPhenoData[na.omit(match(colnames(exprMatrixList[[d]]$expr),masterPhenoData[,patientKey])), ]
      #keeps converting to a list!!! have do add in any potential NA patients in a roundabout way..
      if(any(is.na(match(colnames(exprMatrixList[[d]]$expr),masterPhenoData[,patientKey])))){
        
        stop(paste0("must have pheno data (even if NA) for each patient.\n",
                    "Error occured for dataset ",d,"\n"))
      }
      
      phenoData <- data.frame(phenoData,pData)
      
    }
    
    #Row names of phenoData must match column names of the matrix / matricies in expr.
    rownames(phenoData) <- colnames(exprMatrixList[[d]]$expr)
    #add a variable for study name.
    colnames(phenoData)[1] <- c("datasetName")
    
    #store row names in feature data - see below how will need to update row names of assay data to store featureData.
    featureData <- data.frame(rownames(exprMatrixList[[d]]$expr))

    if(!all(is.null(exprMatrixList[[d]]$featureData))){
      
      featureData <- data.frame(featureData,exprMatrixList[[d]]$featureData)
      
    }
        #want first column to be designated featureDataFieldName.
    colnames(featureData)[1] <- featureDataFieldName
    
    #Row names of featureData data frane must match row names of the matrix / matricies in expr.
    #unforutnately, can't have duplicated rownames in a dataframe...so just make unique IDs for assay data matrix.
    #even probes may have NAs, which will be interpreted as a duplicated row name.
    rownames(exprMatrixList[[d]]$expr) <- c(1:nrow(exprMatrixList[[d]]$expr))
    rownames(featureData) <- rownames(exprMatrixList[[d]]$expr)
    
    #   if(!all(is.null(exprMatrixList[[d]]$annotation))){
    #     
    #     annotation <- exprMatrixList[[d]]$annotation
    #     
    #   }
    #   
    #   
    #   if(!all(is.null(exprMatrixList[[d]]$procotolData))){
    #     
    #     protocolData <- exprMatrixList[[d]]$protocolData
    #     
    #   }
    
    expressionSetList[[d]] <- createS4exprSet(expr=exprMatrixList[[d]]$expr,phenoData=phenoData,featureData=featureData)
    
    
  }
  
  return(expressionSetList)
  
}
  
procExprSetList <- function(exprSetList,outputFileDirectory="./",
                                     numTopVarGenes,minVarPercentile,maxVarPercentile=1,minVar,
                                     featureDataFieldName="gene_symbol",uniquePDataID="unique_patient_ID"){
  
  outputFile = paste0(outputFileDirectory,
                      "/CoINcIDE_processExpressionSetMessages.txt")

  for(e in 1:length(exprSetList)){
    message("\nAnalyzing dataset ",e, " or dataset named ",
            names(exprSetList)[e]);
    
    if(!missing(numTopVarGenes)){
      
     exprSetList[[e]] <- procExprSet( exprSetList[[e]],
                        outputFileDirectory=outputFileDirectory,
                        numTopVarGenes=numTopVarGenes,featureDataFieldName=featureDataFieldName,uniquePDataID=uniquePDataID)
    
    }else if(!missing(minVarPercentile) && !missing(maxVarPercentile) 
             && missing(minVar)){
      
      exprSetList[[e]] <- procExprSet( exprSetList[[e]],
                          outputFileDirectory=outputFileDirectory,
                          maxVarPercentile=maxVarPercentile,
                          minVarPercentile=minVarPercentile,featureDataFieldName=featureDataFieldName,uniquePDataID=uniquePDataID)
      
    }else if(!missing(minVar)){
      
      exprSetList[[e]] <- procExprSet( exprSetList[[e]],
                          outputFileDirectory=outputFileDirectory,minVar=minVar,featureDataFieldName=featureDataFieldName,uniquePDataID=uniquePDataID)
      
      
   }else{
      
      exprSetList[[e]] <- procExprSet(exprSetList[[e]],
                          outputFileDirectory=outputFileDirectory,featureDataFieldName=featureDataFieldName,uniquePDataID=uniquePDataID)
    
    }
   
  }
  
  return(exprSetList)
  
}

#assumption:  your feature data has a column that labeled gene_symbol.
procExprSet <- function(exprSet,outputFileDirectory="./",
                                 numTopVarGenes,minVarPercentile,
                                 maxVarPercentile=1,minVar,featureDataFieldName="gene_symbol",uniquePDataID="unique_patient_ID"){
  
  if(length(na.omit(match(featureDataFieldName,colnames(fData(exprSet)))))!=1){
    
    stop(paste0("Your featureDataFieldName ",featureDataFieldName, " is not found as a unique column name in the fData feature data slot."))
 
  }
  study <- list(expr=exprSet@assayData$exprs,
                keys=fData(exprSet)[ ,featureDataFieldName],phenoData=pData(exprSet))
  tmp <- filterAndImputeSamples(study,studyName = "study",
         outputFile = paste0(outputFileDirectory,
         "/curatedBreastData_processExpressionSetMessages.txt"),impute=TRUE, 
          knnFractionSize=.01,fractionSampleNAcutoff=.005,
         fractionGeneNAcutoff = .01,exprIndex="expr",classIndex="phenoData",
         sampleCol=TRUE,returnErrorRate=FALSE)
  #NOTE: data.matrix() not put here in intial version - updated this bug.
  exprs(exprSet) <- data.matrix(tmp$exprFilterImpute)  
  featureNames(exprSet@assayData) <- c(1:length(tmp$keysFilterImpute))
   #need a data.frame and not cbind: otherwise coerces 
  #gene symbols into numeric data type.
  fData <- data.frame(c(1:length(tmp$keysFilterImpute)), tmp$keysFilterImpute)
  colnames(fData) <- c("number",featureDataFieldName)
    fData(exprSet) <- fData
  #only one set of classes data here - the phenoData, 
  #so will be in first list index.
    rownames(tmp$classesFilter[[1]]) <- colnames(exprs(exprSet))
    pData(exprSet) <- data.frame(tmp$classesFilter[[1]])
  #must also set global featureNames to updated keys, not just featureData.
  #otheriwse won't return a valid ExpressionObject.
  featureNames(exprSet@featureData) <- c(1:length(tmp$keysFilterImpute))
  #if samples were dropped: update protocolData
  #The number of rows and row names of protocolData must agree with the dimension and column names of assayData.
  #usually just an empty data frame - so just update the # of rows (samples.)
  protocolData(exprSet)@data <- data.frame(row.names=colnames(exprs(exprSet)))
  
  if(!validObject(exprSet)){
    
    warning("\nYour expression set is now not valid.
This is usually due to sample or feature names not matching up across
assay, pheno and feature slots.\n Proceed through data analysis with caution!")
    
  }
  
  #expand if multiple symbols on a single line
 # expand_multGeneSymbols(dataMatrix=exprSet@assayData$exprs, keys =fData(exprSet)[[2]], keys.sep = " /// ",checkForMissedMultSymbols=checkForMissedMultSymbols)
  #collapse duplicated probes (really gene symbols here)
   tmp <-  collapseDupProbes(expr=exprSet@assayData$exprs,
            sampleColNames=colnames(exprSet@assayData$exprs),
            keys=fData(exprSet)[[2]], method="highestVariance",
            debug=TRUE,removeNA_keys=TRUE,varMetric = "everything")

  exprs(exprSet) <- data.matrix(tmp$expr)
  fData <- data.frame(tmp$keys,stringsAsFactors=FALSE)
  colnames(fData) <- featureDataFieldName
  fData(exprSet) <- fData
  #must also set global featureNames to updated keys, not just featureData.
  #otheriwse won't return a valid ExpressionObject.
  featureNames(exprSet@featureData) <- tmp$keys
  featureNames(exprSet@assayData) <- tmp$keys
  
  exprMatrix <- exprSet@assayData$exprs
  #eset does not allow duplicated patient names - so have to
    #assign the unique IDs to this matrix.
  if(length(na.omit(match(uniquePDataID,colnames(pData(exprSet))))) != 1){
    
    stop(paste0("\nColumn name ", uniquePDataID, " not found in pData data frame\n"))
  
  }
  
  if(!all(is.na(pData(exprSet)[ ,uniquePDataID]))){
    #for certain datasets: this is actually all NA values. so don't use then!
    colnames(exprMatrix) <- pData(exprSet)[ ,uniquePDataID]
  
  }
  tmp <-  removeDupPatients(exprMatrix=exprMatrix, 
          outputFile=paste0(outputFileDirectory,
          "/curatedBreastData_processExpressionSetMessages.txt"), 
          varMetric = "everything")

  if(ncol(exprs(exprSet))!=ncol(tmp)){
    
    exprs(exprSet) <- tmp
    #lost some samples. need to subset the matrix
    pData(exprSet) <- pData(exprSet)[ ,na.omit(match(colnames(exprs(exprSet)),
                                    rownames(pData(exprSet))))]
    
  }

  study <- list(expr=exprSet@assayData$exprs,
                keys=fData(exprSet)[ ,featureDataFieldName])
  
  if(!missing(numTopVarGenes)){
    
    tmp <- filterGenesByVariance(study, plotSaveDir="~/",
          exprIndex = "expr", keysIndex = "keys", 
          outputFile=paste0(outputFileDirectory,
          "/curatedBreastData_processExpressionSetMessages.txt"),
          varMetric = c("everything"),sampleCol=TRUE,
          numTopVarGenes=numTopVarGenes,plotVarianceHist=FALSE);
    
    exprs(exprSet) <- data.matrix(tmp$filteredStudy$expr)
    fData <- data.frame(tmp$filteredStudy$keys,stringsAsFactors=FALSE)
    colnames(fData) <- featureDataFieldName
    fData(exprSet) <- fData
    #must also set global featureNames to updated keys, not just featureData.
    #otheriwse won't return a valid ExpressionObject.
    featureNames(exprSet@featureData) <- tmp$filteredStudy$keys
    featureNames(exprSet@assayData) <- tmp$filteredStudy$keys
    
    #do quantiles (percentage cut-offs.)
  }else if(!missing(minVarPercentile) && !missing(maxVarPercentile)
           && missing(minVar)){
    
    tmp <- filterGenesByVariance(study, plotSaveDir="~/",
          minVarPercentile=minVarPercentile,maxVarPercentile=maxVarPercentile,
          exprIndex = "expr", keysIndex = "keys", 
          outputFile=paste0(outputFileDirectory,
          "/curatedBreastData_processExpressionSetMessages.txt"),
          varMetric = c("everything"),sampleCol=TRUE,plotVarianceHist=FALSE);
    
    
    
    exprs(exprSet) <- data.matrix(tmp$filteredStudy$expr)
    fData <- data.frame(tmp$filteredStudy$keys,stringsAsFactors=FALSE)
    colnames(fData) <- featureDataFieldName
    fData(exprSet) <- fData
    #must also set global featureNames to updated keys, not just featureData.
    #otheriwse won't return a valid ExpressionObject.
    featureNames(exprSet@featureData) <- tmp$filteredStudy$keys
    featureNames(exprSet@assayData) <- tmp$filteredStudy$keys
    
  }else if(!missing(minVar)){
    
    tmp <- filterGenesByVariance(study, plotSaveDir="~/",minVar=minVar,
           exprIndex = "expr", keysIndex = "keys", 
           outputFile=paste0(outputFileDirectory,
            "/curatedBreastData_processExpressionSetMessages.txt"),
            varMetric = c("everything"),sampleCol=TRUE,plotVarianceHist=FALSE);
    
    
    exprs(exprSet) <- data.matrix(tmp$filteredStudy$expr)
    fData <- data.frame(tmp$filteredStudy$keys,stringsAsFactors=FALSE)
    colnames(fData) <- featureDataFieldName
    fData(exprSet) <- fData
    #must also set global featureNames to updated keys, not just featureData.
    #otheriwse won't return a valid ExpressionObject.
    featureNames(exprSet@featureData) <- tmp$filteredStudy$keys
    featureNames(exprSet@assayData) <- tmp$filteredStudy$keys
    
  }
  #test: still a valid expression set??
    if(!validObject(exprSet)){
      
      warning("\nYour expression set is now not valid.
This is usually due to sample or feature names not matching up across assay, 
pheno and feature slots.\n Proceed through data analysis with caution!")
      
    }

    return(exprSet)

}
#fractionGeneNAcutoff:  max fraction of NAs allowed for a certain gene across all samples.
#fractionSampleNAcutoff: max fraction of NAs allowed for a certain sample across all genes.
filterAndImputeSamples <- function(study,studyName = "study",
                                  outputFile = "createTestTrainSetsOutput.txt",
                                  impute=TRUE, knnFractionSize=.01,
                                  fractionSampleNAcutoff=.005,
                                  fractionGeneNAcutoff = .01,
                                  exprIndex="expr",classIndex,sampleCol=TRUE,
                                  returnErrorRate=TRUE){

  #   #knn impute format: An expression matrix with genes 
  #in the rows, samples in the columns
  if(sampleCol){
    #case were samples are the columns, not the rows.
    
    expr <- study[[exprIndex]]
    exprOrig <- study[[exprIndex]]
    
  }else{
    
    expr <- t(study[[exprIndex]])
    exprOrig <- t(study[[exprIndex]])
    warning("dimensions of expression study will be returned transposed: 
samples are now the columns for an pxn matrix.")
    
  }  

  #NOTE: must have a "keys" list in here!
  keysOrig <-study$keys
  numPatients <- dim(expr)[2]
  
  if(is.null(expr) || is.null(keysOrig) || all(is.na(expr)) || all(is.na(keysOrig))){
    
    stop("you didn't provide the right expr or keys keyword 
         for the filter genes by variance function.")
    
  }
  
  totalGen <- dim(expr)[1]
  
  warning("\nJust a warning: this function assumes your missing values
  are proper NAs, not \"null\",etc.\n")
  
  gene_fractionNAsamples <- apply(expr,MARGIN=1, 
                                  FUN = function(studyRow,numPatients){
    
    numNA <- length(which(is.na(studyRow)))
    fractionNAsamples <- numNA/numPatients
    return(fractionNAsamples)
    
    #pass in extra arguments to the function here
  },numPatients)
  
  #chances are, more cases where a ton of patients missing a 
  #certain gene than a ton of really bad patients with tons of NAs 
  #across lots of genes. so filter out these certain NA genes first.    
  goodGeneIndices <- which(gene_fractionNAsamples < fractionGeneNAcutoff)
  cat("there were ",length(goodGeneIndices), "genes that passed your max", 
      fractionGeneNAcutoff, "NA filter out of a total of ",totalGen, 
      "genes","\n",append=TRUE,file=outputFile)
  
  expr <- expr[goodGeneIndices,]
  #also need to update keys
  keysFilterImpute <- study$keys[goodGeneIndices]
  
  sample_fractionNAgenes <- apply(expr,MARGIN=2, 
                                  FUN = function(studyCol,totalGenes){
    
    numNA <- length(which(is.na(studyCol)))
    fractionNAgenes <- numNA/totalGenes
    return(fractionNAgenes)
    
    #pass in extra arguments to the function here
  },totalGen)
  
  goodIndices <- which(sample_fractionNAgenes < fractionSampleNAcutoff)   
  #now look at "bad" patient indices for training set
  
  cat("there were ",length(goodIndices), "patients that passed your max", 
      fractionSampleNAcutoff, "NA filter out of a total of ",numPatients, 
      "patients","\n",append=TRUE,file=outputFile)
  #cut out these patients from your already filtered gene list.
  expr <- expr[,goodIndices]
  
  #NOW run KNN impute
  #come back...what if it's a small sample size?
  knnFractionSize <- .1
  knnImpute <- function(knnFractionSize, expr,fractionSampleNAcutoff, 
                        fractionGeneNAcutoff,maxp=1500,rng.seed= 362436069){
    
    #k is the # samples use for neighborhood - try about 5% of patients?
    knnSize <- round(knnFractionSize*(dim(expr)[2]))
    #function below returns a few things - just keep them all, 
    #put out $data after function is run.
    exprFilterImpute <- impute.knn(expr,k=knnSize,colmax=fractionSampleNAcutoff,
                        rowmax=fractionGeneNAcutoff,maxp=maxp,rng.seed=rng.seed)
    return(exprFilterImpute)
    
  }
  
  if(impute){
    
    #what genes have at least 1 NA value?
    
    genesNA <- apply(expr,MARGIN=1, FUN = function(studyRow){
      
      return(any(is.na(studyRow)))
      
      #pass in extra arguments to the function here
    })
    
    #odd... if just try to index with -genesNA get nowhere...
    numNAgenes <- length(which(genesNA==TRUE))
    
    if(numNAgenes==0){
      
      cat("there are zero NA genes in this dataset 
          so returning original study \n",file=outputFile,append=TRUE)
      
      exprFilterImpute <- expr
      
      #there are some missing genes, continue on to impute!       
    }else{
      
      cat("there are ",numNAgenes, "genes with at least one NA value \n",
          file=outputFile,append=TRUE)
      
      #expr must NOT be a list!
      expr <- data.matrix(expr)
      exprFilterImpute <- knnImpute(k=knnFractionSize, expr=expr,
                          fractionSampleNAcutoff, fractionGeneNAcutoff,
                          maxp=1500,rng.seed= 362436069)
      exprFilterImpute <- exprFilterImpute$data
      
    }
    
    
    #don't bother with this if there aren't any NA genes.
    if(numNAgenes > 0 && returnErrorRate){
      
      noNAexpr <- expr[which(genesNA == FALSE),]
      
      #now try making some of these NA values, for each sample (randomly.)
      #mimic mean of actual NA genes missing in each sample.
      #average NA fraction
      fractionNA <- mean(sample_fractionNAgenes[
        which(sample_fractionNAgenes >0)])
      #how many NA genes per sample (column)?
      NAperCol <- round(dim(noNAexpr)[1]*fractionNA)
      
      
      fakeNAexpr <- apply(noNAexpr,MARGIN=2, FUN = function(studyCol, NAperCol){
        
        #make random genes (rows) for this patient NA
        NArows <- sample(1:length(studyCol),NAperCol,replace=FALSE)
        
        studyCol[NArows] <- NA
        
        return(studyCol)
        
        #pass in extra arguments to the function here.
        #match the limit for sample of fractions, or the average?
      },NAperCol)
      
      #now re-run imputation.
      #possible that randomly, for each sample, the same gene
      #was removed...so let that cutoff be 1.
      fakeNAindices <- which(is.na(fakeNAexpr))
      #add a .01 "buffer" to the fractionNA so the function doesn't 
      #throw an error, as our % column NA matches exactly fractionNA.
      fakeexprFilterImpute <- knnImpute(knnFractionSize, expr=fakeNAexpr,
                              fractionSampleNAcutoff=(fractionNA+.01), 
                              fractionGeneNAcutoff=1,maxp=1500,
                              rng.seed= 362436069)
      
      imputedNAs <- fakeexprFilterImpute$data[fakeNAindices]
      trueValues <- noNAexpr[fakeNAindices]
      errorRate <- sum(abs(trueValues-imputedNAs)/trueValues)/
        length(fakeNAindices)
      meanAbsDiff <- sum(abs(trueValues-imputedNAs))/length(fakeNAindices)
      #get NAs by comparing fakeNAexpr and noNAexpr.  compare output 
      #of imputation and noNAexpr.
      #now calculate error. take only non-NA genes. 
      cat("the average absolute error for imputation for study ", studyName,
" is ", errorRate, " with an average absolute difference of log-fold
          change expression of ", meanAbsDiff, "\n","\n", 
file=outputFile,append=TRUE)
      
      
    }else{
      
      cat("error rate not specified to be calculated or 
      there are no missing genes in study ", studyName, 
      " so imputation error rate not calculated. \n",
      file=outputFile,append=TRUE)
      errorRate <- NA
      meanAbsDiff <- NA
      
    }
    
    #end of gene filtering statements.
  }
  
  if(length(na.omit(match(colnames(exprFilterImpute),colnames(expr))))==0){
    
    stop("\nEither all of your samples were removed, or 
    somehow the sample names got changed after imputing data.\n")
    
  }
  #was a class index supplied? (i.e. is 
  #there a class list in this study object?)
  #then we need to also filter out the removed patients here.
  if(missing(classIndex)){
    
    cat("finished imputing study ",studyName, "\n",file=outputFile,append=TRUE)
    
    warning("no list index name for a class/outcomes given,
    so this will not be returned")
    
    study <- list(expr=exprOrig,exprFilterImpute = exprFilterImpute,
    keysFilterImpute=keysFilterImpute,keys=keysOrig,meanAbsDiff=meanAbsDiff,
    errorRate=errorRate)
    return(study)
    
  }else{
    
    classesFilter <- list()
    
    
    
    for(c in 1:length(classIndex)){
      
      if( !(classIndex[c] %in% names(study)) || 
            (length(which(names(study) ==classIndex[c])) >1) ){
        print(classIndex[c])
        stop("you either supplied an incorrect class/outcomes index name or
             there are duplicated names in your study list.")
        
      }
      #get the correctly subsetted classes/outcomes variables 
      #for all designated variables.
      classes <- study[[classIndex[c]]]
      cat("getting the re-ordered and trimmed down class labels.",
          append=TRUE,file=outputFile)
      
      classesFilter[[c]] <- study[[classIndex[c]]][
        na.omit(match(colnames(exprFilterImpute),colnames(study[[exprIndex]]))),
        ,drop=FALSE]
      
      
    }
    
    
    study <- list(expr=exprOrig,exprFilterImpute = exprFilterImpute, 
                  class=classes,classesFilter=classesFilter,
                  keysFilterImpute=keysFilterImpute,keys=keysOrig,
                  meanAbsDiff=meanAbsDiff,errorRate=errorRate)
    
    cat("finished imputing study ",studyName, "\n",file=outputFile,append=TRUE)
    return(study)
    
  }
  
  
  #end of function
}


#remove NAs is remove NA keys.
collapseDupProbes <- function(expr,sampleColNames=colnames(expr),keys, 
                              method=c("average","highestVariance"),debug=TRUE,
                              removeNA_keys=TRUE,
                              varMetric = c("everything", "all.obs", 
                                            "complete.obs", "na.or.complete", 
                                            "pairwise.complete.obs")){
  
  warning("It's best to impute NA values before running this function
otherwise it may set averages to NA if there is 1 NA present.
This function just removes any genes whose key is NA.")
  
  #make sure it's an expression matrix to start with
  #(otherwise with 1 sample won't work).
  #and make sure it's numeric! the data.matrix solves that...
  #but have to make it a data.frame first. weird!
  #if don't as stringsAsFactors=FALSE...get totally wrong values.
  expr <- data.matrix(as.data.frame(expr,stringsAsFactors=FALSE))
  
  if(length(colnames(expr))==0 && missing(sampleColNames)){
    
    stop("error: must have a sample (column) ID specified for this sample.")
    
  }
  
  if(dim(expr)[1] != length(keys)) 
    stop("error: length of keys doesn't match 
         number of rows in expression matrix.")
  
  if(ncol(expr)==1){
    
    warning("\nOnly 1 sample. May encounter edge cases 
            when collapsing duplicated probes.")

  }
  
  if(removeNA_keys){
    
    #catch both NA and "NA"
    #need to alternate between filtering expr and keys 
    #so that they keep the same indices
    #as.matrix() prevents it from turning into a vector after the subsetting. 
    #if you only have a 1-column matrix 
    expr <- as.matrix(expr[which(!is.na(keys)),]) 
    keys <- keys[which(!is.na(keys))]
    expr <- as.matrix(expr[which(keys != "NA"),])
    keys <- keys[which(keys != "NA")]
    #also remove "" keys just in case those 
    #somehow got coded up as "" instead of NA.
    expr <- as.matrix(expr[which(keys != ""),])
    keys <- keys[which(keys != "")]
    
  }
  
  singles.keys <- names(which(table(keys) == 1))
  singles.ind <- which(keys %in% singles.keys)
  
  warning("\nYou may get a warning here because key names are duplicated 
  so it can't use them as row names. That's OK.\n")
  gems <- data.frame(expr = expr, keys = keys,stringsAsFactors = FALSE)
  #hmm this isn't working...but aren't all the lengths the same now???
  #add data.frame() again to force it into dataframe unless it's one column.
  out <- data.frame(gems[ singles.ind, ])
  singleKeys <- out$keys
  rownames(out) <- singleKeys
  #now remove the (last) keys column so only have expression values
  #to convince yourself of this, just type in names(out) 
  #and you'll see what I mean.
  out <- out[,1:(dim(out)[2]-1)]
  
  #re-set to patient names (data.frame action put "expr." as prefix)
  #only 1 sample? that means it'll pop up as a numeric class, 
  #not a matrix. re-set indices.
  
  if(class(out)=="numeric"){
    
    #uggh so much workaround to keep row names for single columns!
    out <- as.matrix(out)
    rownames(out) <- singleKeys
    
  }
  
  #make sure can add sample ID for later reference.
  if(!missing(sampleColNames)){
    
    colnames(expr) <- sampleColNames
    
  }
  
  if(ncol(expr)==0){
    
    stop("\nBug in your code that removes all samples.")
  }
  
  colnames(out) <-colnames(expr)
  
  ## Next, deal with multiple keys within a study ##
  #don't include probes for right now - makes it tricky.
  multis <- gems[ -singles.ind, ]
  
  if(nrow(multis) > 0){  # skip if no multiple ID found
    
    #split will give a list of the length of all the duplicated keys 
    #(1 index per key name - THIS part isn't duplicated now!
    #duplicated values will be in the element's list for that specific key name.)
    multis$keys <- as.factor(multis$keys)
    #levels may be holding on to keys you removed earlier, 
    #so drop these that have no data attached to them now.
    tmp <- split(multis, multis$keys,drop=TRUE)
    
    if (method == "average") {
      
      out2 <- sapply(tmp, function(m) {
        
        #odd bc if you run this code outside of sapply: 
        #it will have flipped indices for m.
        #inside sapply: keys are in the ROWS. 
        #so colMean takes mean INTRA patient across the duplicated keys 
        #- what we want!
        
        #REMOVE LAST COLUMN: this is the key string character, 
        #not an expression value!! 
        m <- m[, (1:dim(m)[2]-1) ]
        #just average across all duplicated probes (for each patient/column: 
        #want a new collapsed value for all duplicated keys, for each patient.)
        
        #if only one sample: set as a matrix. if there's only 1 col but > 1 row:
        #means 1 GENE, not 1 sample.
        #we want to remedy the 1 sample, duplicated genes scenario. 
        #there should be no just 1 gene (row) scenario anyways:
        #we're only looking at multis here...
        
        #well have to re-set every time as a matrix anyways to make all numeric 
        #again so numbers don't end up weird!!
        m <- data.matrix(as.data.frame(m,stringsAsFactors=FALSE))
        #this colMeans is just the duplicates across each SEPARATE patient 
        #(not lumping patient expression values together)
        out2 <- colMeans(m,na.rm=TRUE)
        
      })
      
      #sapply flips it - puts the keys as columns.
      
      #but if it's one PATIENT (sample): that means we need to make it a matrix,
      #and don't transpose it after
      #...but then lose the row names...so make sure to add them back!
      
      
      if(class(out2)=="numeric"){
        
        out2 <- as.matrix(out2)
        
        
      }else{
        #sapply flips it- puts the keys as columns.
        out2 <- t(out2)
        
      }
      
      colnames(out2) <-colnames(expr)
      rownames(out2) <- names(tmp)
      
      
      if(debug){
        #randomly find a key that was duplicated, check its means.
        if(length(unique(multis$keys)) != dim(out2)[1]) stop ("function is not 
        collapsing duplicated keys to the expected number of unique keys")
        
        #try just the first duplicated key - we don't know if 
        #there will be more duplicated ones.
        keyInd <- grep(rownames(out2)[1],multis$keys)
        #is mean of these rows (excluding last column - 
        #the key characters) as expected?
        exprRows <- as.matrix(multis[keyInd,1:(dim(multis)[2]-1)])
        
        
        #if(dim(as.matrix(exprRows))[2]==1 && dim(as.matrix(out2))[1]>1){
        
        #if only 1 sample: have to make it a 2D matrix first.
        #means <- colMeans(as.matrix(exprRows)) 
        
        # }else{
        
        means <- as.numeric(colMeans(exprRows,na.rm=TRUE))
        
        #  }
        
        #ignore samples who actually just have an NA value for this key.
        #once I got a bug in this becuase I hadn't filtered out a few weird
        #keys with "" as names and not NA.
        if(any ((means==out2[1,]) == FALSE,na.rm=TRUE)) stop("function is not 
returning the expected mean of rows of expression values for duplicated keys ")
        
      }
      
      #end of average calculations.
    }
    
    
    ######
    
    if (method == "highestVariance") {
      
      #we need lapply here bc tmp is a list - can't just use apply.
      #sapply was too buggy with this variance function option - 
      #never quite figured out why.
      out2 <- lapply(tmp, function(m) {
        
        #just average across all duplicated probes (for each patient.)
        #var() computes the variance across columns...what we want - 
        #1 gene's variance across all patients
        #REMOVE LAST COLUMN: this is the key string character, 
        #not an expression value!!
        m <- m[, (1:dim(m)[2]-1) ]
        #make all numeric again so numbers don't end up weird!!
        m <- data.matrix(as.data.frame(m,stringsAsFactors=FALSE))
        #here, we want to look at the variance (diag of cov) of each 
        #row (i.e. across all the columns (patients))
        #remember covariance is t(X)X - want numrows*numPatients = X
        #m will have multiple rows...so use apply.
        geneVar <- apply(m, MARGIN=1, FUN = function(geneRow,varMetric){
          
          
          varGene <- var(geneRow,use=varMetric)
          
          
          return(varGene)  
        }, varMetric)
        
        #if there's a tie: just randomly pick the first one - 
        #otherwise will return double the patients!       
        rowIndex <- which(geneVar == max(geneVar))[1]
        
        #need to make this a vector - otherwise hard to format into a data.frame
        #when writing up to SQL table later...
        #and sapply tries to make m a list - not good here!
        #cat("\n",dim(m))
        out2 <- m[rowIndex, ]
        
      })
      
      
      #lapply flips it - puts the keys as columns.
      #also need to convert lists to 1 matrix
      out2 <- do.call(rbind,out2)
      colnames(out2) <-colnames(expr)
      rownames(out2) <- names(tmp)
      
      if(debug){
        #randomly find a key that was duplicated, check its means.
        if(length(unique(multis$keys)) != dim(out2)[1]) stop ("function is not 
        collapsing duplicated keys to the expected number of unique keys")
        
        #try just the first duplicated key - we don't know if there will
        #be more duplicated ones.
        keyInd <- grep(rownames(out2)[1],multis$keys)
        #is the row chosen out of this subset as expected?
        exprRows <- multis[keyInd,1:(dim(multis)[2]-1)]
        
        var <- diag(cov(t(exprRows),use=varMetric))
        
        rowIndex <- which(var == max(var))[1]          
        bestRow <- exprRows[rowIndex, ]
        
        #so do the values for each patient match up?
        #ignore samples who actually just have an NA value for this key.
        if( any((bestRow==out2[1,]) == FALSE ,na.rm=TRUE)) stop("function 
      is not returning the expected highest row variance for duplicated keys ")
        
      }
      
      #end of variance option. 
    }
    
    #save some original values to test this function does what we want.
    exprOut <- rbind(out, out2)
    #backstory: the m[row, ] part in method = "highVariance" makes like nested 
    #lists...and RmySQL goes crazy when you do this
    #tne entire expr table can be a list, but it shouldn't be a list of lists..
    #..found this trick on stack overflow:
    #http://stackoverflow.com/questions/8305735/error-writing-to-csv
    #print(dim(exprOut))
    #exprOutFirst <- exprOut
    #careful: below does weird stuff...also remove row names!
    #exprOut <- as.data.frame(lapply(exprOut,unlist))
    #add row names afterwards?
    
    #save(exprOutFirst, exprOut,file="exprOut.RData")
    #keys are in row names of out,out2.
    keysOut <- append(rownames(out),rownames(out2))
    #hmm getting matching probes may be harder...probesOut <- 
    #end of multis processing 
    
    #no duplicates found.
  }else{
    
    exprOut <- out
    keysOut <- rownames(out)
    
  }

  if(debug==TRUE){

    
    if(length(keysOut) != dim(exprOut)[1]) stop("function is not returning 
    the same number of expression rows as keys")
    
    if(dim(expr)[2] != dim(exprOut)[2]) stop("function is not returing the 
    orginal number of patients/samples (columns in expr matrix)")
    
    if(length(unique(keysOut)) != length(keysOut)) stop("function is not 
    returning the correct number of unique keys in final step.")
    
    
  }
  
  
  
  gems <- list(expr = exprOut, keys = keysOut)
  
  return(gems)
  
}

#picks patient's sample that has the highest variance across all features.
removeDupPatients <- function(exprMatrix, 
                                     outputFile="duplicatedPatientsOutput.txt", 
                                    varMetric = c("everything", "all.obs",
                                                  "complete.obs", 
                                                  "na.or.complete", 
                                                  "pairwise.complete.obs")){
  
  message("\nStarting with  ", ncol(exprMatrix) ," patients.")
  cat("\nStarting with ", ncol(exprMatrix) ,"patients. \n\n",
      append=TRUE,file=outputFile)
  
  
  if(length(varMetric)>1){
    
    warning("defaulting the everything variance metric.")
    varMetric = c("everything")
    
  }
  nSamples <- ncol(exprMatrix)
  
  if(all(is.na(exprMatrix))){
    
    stop("your counts matrix  is empty.")
  }
  #all.obs this means that missing observations will produce an error. 
  #we want that! there should be NO NAs.
  
  #take only the UNIQUE set of these (if there's 3 samples from 1 patient, 
  #will show up twice in duplicated list)
  duplicated_patients <- unique(colnames(exprMatrix)[
    which(duplicated(colnames(exprMatrix)))])
  
  
  
  if(length(duplicated_patients)>0){
    
    message("Found duplicated samples for ",length(duplicated_patients),
            "patients. There may be more than 2 samples per patient.")
    cat("Found duplicated samples for patients ",duplicated_patients,
        "patients. There may be more than 2 samples per patient.\n",
        append=TRUE,file=outputFile)
    
    #remember: duplicated only returns the index of the n-1 duplicated items, 
    #NOT the first time the item appeared!
    #why need the "which" in the for loop to grap all samples 
    #related to 1 patient.
    
    for(d in 1:length(duplicated_patients)){
      
      
      #how many NAs in each?
      dupIndexes <- which(colnames(exprMatrix)== duplicated_patients[d])
      
      sampleVars <- apply(exprMatrix[,dupIndexes], MARGIN=2, 
                          FUN = function(sampleCol,varMetric){
        
        
        varSample <- var(sampleCol,use=varMetric)
        
        
        return(varSample)  
      }, varMetric)
      
      
      dup_sample_highestVar <- which(sampleVars==max(sampleVars))
      
      if(length(dup_sample_highestVar)>1){
        
        #a tie. ha probably means a bug...that would be rare! 
        #or you literally just duplicated a data column but tacked
        #on a different aliquot name...
        #just take arbitrary first sample then.
        dup_sample_highestVar <- dup_sample_highestVar[1]
        
      }
      
      
      #remove the other duplicated samples besides the ith one.
      dupSampleRemove <- dupIndexes[-dup_sample_highestVar]
      
      if(length(dupSampleRemove)==0){
        
        stop("bug in your remove duplicated patients code...no indices 
        left to remove the duplicated patient sample.")
        
      }
      
      exprMatrix <- exprMatrix[,-dupSampleRemove]
      
      
      
      #end of loop
    }
    
    
    #how many samples did we lose?
    numSamplesRemove <- nSamples - ncol(exprMatrix)
    message("removed ",numSamplesRemove, " samples from patients with 
            multiple samples. \n There are ", ncol(exprMatrix) ,"left.")
    cat("removed ",numSamplesRemove, " samples from patients with multiple 
        samples. \n There are ", ncol(exprMatrix) ,"left. \n\n",
        append=TRUE,file=outputFile)
    
    #end of if statement.
  }else{
    
    message("found no multiple samples from the same patient(s)")
    cat("found no multiple samples from the same patient(s) \n\n",
        append=TRUE,file=outputFile)
    
  }
  
  #return study (processed or not depending on if found duplicated)
  return(exprMatrix)
}


#note: this can handle duplicated keys: it never actually uses the key names, 
#so they don't need to be unique.
#this filters by variance, not coefficient of variance.
#if not using logged values, sometimes filtering by coeff of var is better.
#to-do: replace missing() with each variance option as NULL to make running code 
#in another function less cumbersome
filterGenesByVariance <- function(study, plotSaveDir="~/",minVarPercentile,
                                  maxVarPercentile=1,maxVar,minVar,
                                  exprIndex = "expr", keysIndex = "keys", 
                                  outputFile="varCal.txt",
                                  plotVarianceHist=FALSE,
                                  varMetric = c("everything", "all.obs", 
                                                "complete.obs", 
                                                "na.or.complete", 
                                                "pairwise.complete.obs"),
                                  sampleCol=TRUE,
                                  numTopVarGenes){
  
  
  if(sampleCol){
    #case were samples are the columns, not the rows.
    
    expr <- study[[exprIndex]];
    
  }else{
    
    expr <- t(study[[exprIndex]]);
    
    warning("dimensions of expression data will be returned transposed: 
    samples are now the columns for a pxn matrix.");
    
  }
  
  
  keys <- study[[keysIndex]];
  
  if(all(is.na(expr)) || all(is.na(keys))){
    
    stop("you didn't provide the right expr and/or keys keyword for the 
         filter genes by variance function.")
    
  }
  
  cat("calculating variance across all genes for plotting","\n",
      file=outputFile,append=TRUE);
  
  
  #we want a variance calculation for each gene (row here)
  geneVar <- apply(expr, MARGIN=1, FUN = function(geneRow,varMetric){
    
    
    varGene <- var(geneRow,use=varMetric);
    
    
    return(varGene);  
  }, varMetric);
  
  names(geneVar) <- keys;
  cat("finished calculating variance","\n",file=outputFile,append=TRUE);
  
  #remove any genes with NA variance - cannot estimate these...
  
  geneVarNA <- which(is.na(geneVar));
  
  if(length(geneVarNA)>0){
    
    expr <- expr[-geneVarNA,];
    
    keys <- keys[-geneVarNA];
    
  }
  
  cat("found ",length(geneVarNA), "genes with NA variance. Removing them. \n",
      file=outputFile,append=TRUE);
  
  if(plotVarianceHist){
    
    cat("plotting a histogram of the variance for each gene 
    across all patients.","\n",file=outputFile,append=TRUE)
    p <- ggplot(as.data.frame(geneVar),aes(x=geneVar)) + 
      geom_histogram(binwidth=round((max(geneVar,na.rm=TRUE)-
                                       min(geneVar,na.rm=TRUE))/10)) +
      labs(title = "Histogram of gene expression variances",
           y ="number of genes", x="variance of expression");
    
    fileSave <- paste(plotSaveDir,"geneVarianceHist.png",sep="")
    ggsave(filename=fileSave);
    
  }else{
    
    p <- "no plot requested from function inputs";
    
  }
  
  
  if(missing(numTopVarGenes)){
    
    if(!missing(minVar)){
      #use an arbitrary variance cutoff
      
      medianGeneVar <- median(geneVar);
      
      
      geneVarMedianThresholded <- geneVar[which(geneVar > medianGeneVar)];
      cat("there were ",length(geneVar)," genes with no NA values 
     and the median variance of these genes was ", medianGeneVar, 
      "\n","there were", length(geneVarMedianThresholded) ,
      "genes above this median","\n" ,append=TRUE,file=outputFile);
      
      
      geneVarThreshholded <- geneVar[which(geneVar > minVar)];
      cat("there were ",length(geneVar),
          " genes with no NA values and there were ", 
          length(geneVarThreshholded) ,"genes above the mininum variance of ",
          minVar,"\n" ,append=TRUE,file=outputFile);
      
      expr <- expr[which(geneVar > minVar),];
      
      keys <- keys[which(geneVar > minVar)];
      
      filteredStudy <- list(expr=expr,keys=keys);
      
    }
    
    if(!missing(maxVar)){
      
      #only take genes UNDER a certain variance
      medianGeneVar <- median(geneVar);
      
      
      geneVarMedianThresholded <- geneVar[which(geneVar < medianGeneVar)];
      cat("there were ",length(geneVar)," genes with no NA values 
      and the median variance of these genes was ", medianGeneVar, "\n",
          "there were", length(geneVarMedianThresholded) ,
          "genes below this median","\n" ,append=TRUE,file=outputFile);
      
      
      geneVarThreshholded <- geneVar[which(geneVar > maxVar)];
      cat("there were ",length(geneVar),
          " genes with no NA values and there were ", 
          length(geneVarThreshholded) ,"genes below the maximum variance of ",
          minVar,"\n" ,append=TRUE,file=outputFile);
      
      expr <- expr[which(geneVar < maxVar),];
      
      keys <- keys[which(geneVar < maxVar)];
      
      filteredStudy <- list(expr=expr,keys=keys);
      
      
    }
    
    #percentiles will get messed up if don't do this together
    #default to min/max Var if that was also given.
    if(!missing(minVarPercentile) && !missing(maxVarPercentile)
       && missing(minVar) && missing(maxVar)){
      
      #do quantiles instead
      #if want .001: then times *1000+1 indices.
      varQuantiles <- quantile(geneVar,probs=seq(0,1,.01),na.rm=TRUE);
      #want to match up with varQuantiles index (.95 is 96th index)
      minVarPercentile <- minVarPercentile*100+1;
      minVar <- varQuantiles[minVarPercentile];
      
      
      filteredStudy <- list(expr=expr,keys=keys);
      #want to match up with varQuantiles index (.95 is 96th index)
      maxVarPercentile <- maxVarPercentile*100+1;
      maxVar <- varQuantiles[maxVarPercentile];
      
      #now filter out these genes. must do intersect: often, 
      #a lot of the genes will overlap!
      
      expr <- expr[intersect(which(geneVar > minVar),which(geneVar < maxVar)),]; 
      keys <- keys[intersect(which(geneVar > minVar),which(geneVar < maxVar))];
      filteredStudy <- list(expr=expr,keys=keys);
      
    }else if(missing(minVar) && missing(maxVar)){
      
      stop("you didn't specify the number of genes to threshold but didn't 
           #give a variance threshold metric either.")
    }
    
  }else{
    
    #just take top X varying genes.
    if(any(is.na(geneVar))){
      
      warning("You have NA values in your gene variances. 
      variances with NA values will not be considered.")
      
    }
    #na.last=TRUE not allowed if returning an index like we are.
    geneVarSorted <- sort.int(geneVar,na.last=NA,decreasing=TRUE,
                              index.return=TRUE);
    #this way, can actually handle duplicated gene names. 
    #never sees the gene names, just the indices.
    if(length(geneVarSorted$ix)>= numTopVarGenes){
      
      topGeneIndices <- geneVarSorted$ix[1:numTopVarGenes];
      
    }else{
      #just take the rest of the genes
      cat("there weren't ",numTopVarGenes, " to choose from, 
      so just took all of the remaining genes.\n",append=TRUE,file=outputFile)
      topGeneIndices <- geneVarSorted$ix;
      
    }
    
    expr <- expr[topGeneIndices, ];
    keys <- keys[topGeneIndices];
    
    #come back and TEST:
    #take first index - there may be a tie
    topGeneIndex <- which(geneVar ==max(geneVar,na.rm=TRUE))[1];
    
    if(! (topGeneIndex %in% topGeneIndices) ){
      
      stop("looks like there's a bug in how you're 
    sorting the top gene indices.")
    }
    
    filteredStudy <- list(expr=expr,keys=keys);
    
    output <- list(study=study,filteredStudy=filteredStudy,p=p);
    
    return(output);
    
  }
  
  output <- list(study=study,filteredStudy=filteredStudy,p=p);
  
  return(output);
  
}


exprSetListToMatrixList <- function(esetList,featureDataFieldName="gene_symbol"){
  
  dataMatrixList <- list()
  for(e in 1:length(esetList)){
    
  dataMatrixList[[e]] <- exprs(esetList[[e]])
  rownames(dataMatrixList[[e]]) <- fData(esetList[[e]])[,featureDataFieldName]
  
  }
  
  return(dataMatrixList)
  
}

mergeDatasetList <- function(datasetList,minNumGenes = 10000, minNumPatients = 40,batchNormalize = c('BMC','none','combat'),NA_genesRemove=TRUE,
                              outputFile="./mergeDatasets",combatPvalueThresh=.05){
  
  dataMatrixList <- list();
  numDatasets <- 0;
  removeDatasetIndices <- c()
  for(d in 1:length(datasetList)){
    

    if(nrow(datasetList[[d]]) >= minNumGenes){
      
      if(ncol(datasetList[[d]])>= minNumPatients){
        
        numDatasets <- numDatasets + 1;
        
        if(batchNormalize=='BMC'){
          
          dataMatrixList[[numDatasets]] <- data.matrix(datasetList[[d]]) - rowMeans(data.matrix(datasetList[[d]]),na.rm=TRUE,dims=1)
          
        }else if(batchNormalize=='none' || batchNormalize=='combat'){
          
          dataMatrixList[[numDatasets]] <-  data.matrix(datasetList[[d]]);
          
        }
        
        names(dataMatrixList)[numDatasets] <- names(datasetList)[d];

        rownames(dataMatrixList[[numDatasets]]) <- rownames(datasetList[[d]])
        
      }else{
        
        removeDatasetIndices <- append (removeDatasetIndices,d)
      }
      
    }else{
      
      removeDatasetIndices <- append (removeDatasetIndices,d)
    }
    
    
  }
  
  rm(list="datasetList");
  study <- c();
  sampleName <- c();
  
  cat("\nFound ",numDatasets," passing the input thresholds.\n")
  
  for(d in 1:length(dataMatrixList)){
    study <- append(study,rep(names(dataMatrixList)[d],ncol(dataMatrixList[[d]])));
    #study <- c(1:length(datasetList[[d]]$class));
    sampleName <- append(sampleName,colnames(dataMatrixList[[d]]));
    
    
    if(d>=2){
      
      #I tried the regular old merge() function and it crashed the server!
      #merge() also orders your patients/samples/column in a weird way?
      #besides, it takes only 1 extra line of code to do it myself...
      tmp <- data.matrix(dataMatrixList[[d]]);
      intersectGenes <- intersect(rownames(tmp),rownames(mergedExprMatrix));
      
      if(NA_genesRemove){
        
        if(length(which(is.na(intersectGenes)))>0){
          #don't want to match up on NA genes!
          intersectGenes <- intersectGenes[-which(is.na(intersectGenes))];
          
        }
        
      }
      
      mergedExprMatrix <- cbind(mergedExprMatrix[na.omit(match(intersectGenes,rownames(mergedExprMatrix))), ],
                                tmp[na.omit(match(intersectGenes,rownames(tmp))), ]);
      
      
    }else if(d==1){
      
      mergedExprMatrix <- data.matrix(dataMatrixList[[d]]);
      
    }
    
  }
  #factor easier for combat to interpret.
  study <- factor(study)
  sampleData <-  data.frame(study,sampleName)
  #make study names the numeric factors for combat.
  colnames(sampleData) <- c('batch','sampleName')
  sampleDataOrig <- sampleData
  sampleData$batch <- as.numeric(study)
  #may get a warning here if sample names aren't unique
  if(any(duplicated(sampleName))){
    
    stop("Must have unique sample names across the entire dataset (if this is regularly a pain, please provide feedback to CoINcIDE author.)")
    
  }
  row.names(sampleData) <- sampleName
  if(batchNormalize=='combat'){

    combatOut <- batchNormalization(countsMatrixNoNANoDup=mergedExprMatrix,outcomesAndCovariates=sampleData,
                       MinInBatch=4,combatModelFactorName=NULL,pvalueThresh=combatPvalueThresh,batchColName="batch",outputFile=outputFile)
    
    if(!is.null(combatOut$GEN_Data_Corrected)){
      
      mergedExprMatrix <- combatOut$GEN_Data_Corrected
      sampleDataPostCombat <- sampleData[na.omit(match(combatOut$BatchData[,"batch"],sampleData[,"batch"])), ]
      
    }
  }else{
    
    combatOut <- NULL
    sampleDataPostCombat <- NULL
  }
  
  output <- list(mergedExprMatrix=mergedExprMatrix,sampleData=sampleDataOrig, removeDatasetIndices= removeDatasetIndices,
                 combatOut = combatOut,sampleDataPostCombat=sampleDataPostCombat);
  return(output);
  
}

expand_multGeneSymbols <-function (dataMatrix=NULL, keys, keys.sep = " /// ",checkForMissedMultSymbols=TRUE){

  keys <- as.character(keys);
  
  if(is.null(dataMatrix)){
    
    #just make a dummy expr list then.
    dataMatrix <- matrix(rep.int(1,times=length(keys)));
    
  }
  
  #print(df$probes[1:10])
  print("old number of gene symbols is:");
  print(length(keys));
  #gsub("///", ",", keys);
  skey <- strsplit(keys, split = keys.sep);
  times <- sapply(skey,length);
  #need a 1 where times = 0...want to replicate this gene symbol still once!
  zeroIndices <- which(times==0, arr.ind=TRUE);
  times[zeroIndices] <- 1;
  #remove indices that are only "" i.e. blank.
  #look for empty values! otherwise "unlist" will remove these
  #empty does not equal NULL or NA in R.
  #can't assign NULL values, use NA.
  skey[zeroIndices]<-NA;
  #print(skey[1:30]);
  #string out the symbols
  skeyList <-unlist(skey);
  
  #pass over Gene Symbol check...already did this!
  keys <-as.character(as.matrix(skeyList));
  print("new number of gene symbols is:");
  print(length(keys));
  #just print some out to gut-check
  #print(keys[1:30])
  #print(keys[1000:1020])
  
  #we need more rows to string out the duplicate probes.
  numRow <- length(keys);
  #still have same number of samples
  numCol <- ncol(dataMatrix);
  expr <-matrix(data=NA, nrow = numRow, ncol = numCol);
  #error check: any weird symbols in your probe table?? stop then!
  #margin must =2 (2 is for columns)
  keysWithPunct<- grep("[[:punct:]]",keys);
  #now take this subset. will still have "-" and this is accepted.
  #why unique: sometimes have symbols like "RP11-35N6.1" with punctuation marks
  #want unique indices - could still mean some actual symbols are duplicated
  keysWithPunct <- unique(keysWithPunct);
  keysWithPunct <- keys[keysWithPunct];
  keysWithPunctPeriod <- grep("\\.",keysWithPunct);
  
  if(length( keysWithPunctPeriod )>0){
    
    warning("your gene symbols have . version tags on them.\nThis could mean distinct gene types (common for lncRNAs) or actually just versions.\nLook into this.")
    
  }
  
  #don't run this if you KNOW that your code has a few different mult symbols - like / and ///. Then just run this function twice.
  if(checkForMissedMultSymbols){
  #hopefully all the probes with punctation were all the probes with "-" or "@"...
  #(yes, a double gene symbol that is say HEXA-AS1///HLA-F won't be caught...but that would be an extreme outlier)
  #why a union:  sometimes have symbols like "RP11-35N6.1" with two punctuation marks.
  #don't want to count this twice
    #look for standard punctuation found in gene symbols that doesn't link to mult symbls on one line.
    keysWithPunctHyphen <- grep("-",keysWithPunct);
    keysWithPunctAt <- grep("@",keysWithPunct);
    keysWithPunctUnder <- grep("_",keysWithPunct);
    keysWithPunctQuote <- grep("\\'",keysWithPunct);
    
  if( length(keysWithPunct)!=length(union(keysWithPunctHyphen,union(keysWithPunctAt,union( keysWithPunctUnder,union(keysWithPunctQuote,keysWithPunctPeriod))))) ){
    
    stop("Error: your gene symbols contain extra punctuation beyond the standard -, @, and your keys.sep you entered. check this!")
    
  }
  
  }
  dataMatrix <- as.matrix(dataMatrix)
  for(i in 1:numCol){    
    #rep(x, time): replicates values in x
    #times: number of times to repeat each element of x
    #perfect for replicating gene expression rows based on
    #replicated Gene Symbols.
    #want 1:nrow(expr) becuase we want to repeat all columns at once
    #sapply(skey, lenght): produces a vector, and each unit
    #is the number of duplicate probes for a row (because strplit parses
    #the probes into separate columns)  
    expr[,i]   <- rep(dataMatrix[,i], times)
    
  }
  row.names(expr)<-keys;
  colnames(expr) <-colnames(dataMatrix)
  

  return(expr)
}


