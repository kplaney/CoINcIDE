
library("curatedOvarianData")

datasetNames <- data(package="curatedOvarianData")
datasetNames <- datasetNames[3]$results[,"Item"]



#odd way to get a dataset list...I just downloaded the source file
#and pull out their code.
expandProbesets <- function (eset, sep = "///"){
    x <- lapply(featureNames(eset), function(x) strsplit(x, sep)[[1]])
    eset <- eset[order(sapply(x, length)), ]
    x <- lapply(featureNames(eset), function(x) strsplit(x, sep)[[1]])
    idx <- unlist(sapply(1:length(x), function(i) rep(i, length(x[[i]]))))
    xx <- !duplicated(unlist(x))
    idx <- idx[xx]
    x <- unlist(x)[xx]
    eset <- eset[idx, ]
    featureNames(eset) <- x
    eset
}


add.surv.y <- function(X) Surv(X$days_to_death, X$vital_status=="deceased")
## -----------------------------------------------------------------------------
##load the esets
## -----------------------------------------------------------------------------

delim <- ":" 

data(list=data(package=package.name)[[3]][,3])

strEsets <- ls(pattern="^.*_eset$")
esets <- list()
for (strEset in strEsets){
    eset <- get(strEset)
    ##Split out rows without unique gene name
    eset <- expandProbesets(eset)
    esets[[strEset]] <- eset

}
#pData(esets) already contains survival data.

#remove the extra TCGA datasets
#TCGA.RNASeqV2_eset                           Integrated genomic analyses of ovarian carcinoma.
#TCGA.mirna.8x15kv2_eset                      Integrated genomic analyses of ovarian carcinoma.
#TCGA_eset                                    Integrated genomic analyses of ovarian carcinoma.
indicesRemove <- na.omit(match(c("TCGA.RNASeqV2_eset","TCGA.mirna.8x15kv2_eset"),names(esets)))
esets <- esets[-indicesRemove]


#look for batches - break these apart?
#but it only looks like TCGA has actual site batches -
#the rest look like the date the samples were collected?
#in the R documentation it says 
#"Hybridization date when Affymetrix CEL files are available"
#perhaps TOO granular for our purposes.

for(e in 1:length(esets)){
  
  if(!all(is.na(pData(esets[[e]])[,"batch"]))){
    
    cat("\n",names(esets[[e]]))
    
  }
  
}

#only split up TCGA into batches
TCGAbatches <- unique(pData(esets[["TCGA_eset"]])[,"batch"])
#any samples with NA batch? these will be thrown out.
which(is.na(pData(esets[["TCGA_eset"]])[,"batch"]))

for(t in 1:length(TCGAbatches)){
  
  expr <- exprs(esets[["TCGA_eset"]])[ ,
                                which(pData(esets[["TCGA_eset"]])[,"batch"]==TCGAbatches[t])]
  
  phenoData <-  pData(esets[["TCGA_eset"]])[
                                which(pData(esets[["TCGA_eset"]])[,"batch"]==TCGAbatches[t]), ]
  
  esets[[paste0("TCGA_eset_batch_",t)]] <-  createS4exprSet(expr=expr,phenoData=phenoData)
  
}

#remove original entire TCGA matrix
esets <- esets[-which(names(esets)=="TCGA_eset")]

#remove datasets that are too small? allow them for now.
#now process, then grab meta-features
#sounds like there SHOULD be some duplicated samples  in here.
esets <- processExpressionSetList(exprSetList=esets,outputFileDirectory="~/Desktop/",
                                     minVar=0,featureDataFieldName="gene")

#now format just as a list of data matrices.
dataMatrixList <- exprSetListToMatrixList(esets,featureDataFieldName="gene")

names(dataMatrixList) <- names(esets)
#run meta-feature analysis.
metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),numFeatSelectByGlobalRank=1000,
                                       numTopFeatFromEachDataset=10,fractNATopFeatAllowedPerDataset=.2,selectMethod=c("median"),
                                       outputFile="~/Desktop//selectFeaturesMetaVarianceOut.txt")
#remove datasets with too many missing top gene features
dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]

#STOPPED HERE.
clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- metaFeatures$finalFeatures
  
}
clusterOut <- clusterMatrixListHclustGap(dataMatrixList=dataMatrixList,clustFeaturesList=clustFeaturesList,maxNumClusters=20,algorithm="ward.D",
                              distMethod=c("euclidean"),outputFile="~/Desktop//cluster_hclustGap_output.txt",
                              corUse=c("everything"),numSims=100)


  
  
  
  
  
  
  
  
  
####################
#extra functions from curatedBreatData code.
  esetToSurv <- function(eset, output, months=TRUE){
  if (sum(eset$vital_status =="deceased", na.rm=TRUE) == 0) return(NULL)
  eset <- eset[ ,!is.na(eset$vital_status) & !is.na(eset$days_to_death)]
  time <- eset$days_to_death
  if (months)
      time <- time / 30.4
  cens <- rep(NA,length(time))
  if(output == "event"){
      cens[eset$vital_status == "living"] <- 0
      cens[eset$vital_status == "deceased"] <- 1
  }else if(output == "reverseevent"){
      cens[eset$vital_status == "living"] <- 1
      cens[eset$vital_status == "deceased"] <- 0
  }
  Surv(time,cens)
}

medianSurvival <- function(eset){
    ##make surv objects
    y <- esetToSurv(eset, "event")
    if (is.null(y)) {
        binary.table <- paste(table(eset$os_binary), collapse="/")
        if(length(binary.table) == 0 || binary.table == "")
            binary.table <- NA
        output <- c(rep(NA,3), binary.table )
    } else {
        y.rev <- esetToSurv(eset, "reverseevent")
        ##kaplan-meier fits
        km.fit <- survfit(y ~ 1)
        reverse.km.fit <- survfit(y.rev ~ 1)
        ##get median survival and follow-up times
        surv.median.position <- min( which(km.fit$surv < 0.5) )
        surv.median.time <- km.fit$time[ surv.median.position ]
        followup.median.position <- min( which(reverse.km.fit$surv < 0.5) )
        followup.median.time <- reverse.km.fit$time[ followup.median.position ]
        ## % censoring
        percent.censoring <- round(sum(y[,2] == 0) / nrow(y) * 100)
        output <- round( c(surv.median.time, followup.median.time, percent.censoring, NA) )
    }
    names(output) <- c("median.survival", "median follow-up", "percent censoring", "binarized OS (long/short)")
    return(output)
}
getEsetData <- function(eset, 
    possible.stages = c("early", "late"),
    possible.histological_types = c("ser", "clearcell", "endo", "mucinous", "other", "unknown")
){
    makeSummary <- function(field, possible.values){
        my.table <- sapply(possible.values, function(x) sum(eset[[field]] == x, na.rm=TRUE))
        my.table <- c(my.table, ncol(eset) - sum(my.table))
        my.table <- paste(my.table, collapse="/")
        return(my.table)
    }
    survdata <- medianSurvival(eset)
    output <- c(experimentData(eset)@lab,  #name/year
                experimentData(eset)@pubMedIds,
                ncol(eset),  # number of samples
                makeSummary("summarystage", possible.stages),
                makeSummary("histological_type", possible.histological_types),
                sum(eset$histological_type == possible.histological_types[1] &
                eset$summarystage == possible.stages[2], na.rm=TRUE), ##number of serous samples
                ifelse(length(experimentData(eset)@other$platform_shorttitle), experimentData(eset)@other$platform_shorttitle,""),  #platform
                survdata)
    names(output) <- c("Study",
                       "PMID",
                       "N samples",
                       "stage",
                       "histology",
                       "known late-stage serous",
                       "Platform",
                       names(survdata))
    return(output)
}
