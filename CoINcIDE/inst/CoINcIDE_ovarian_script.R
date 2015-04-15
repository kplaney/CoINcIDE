library("Biobase")
library("curatedOvarianData")
datasetNames <- data(package="curatedOvarianData")
datasetNames <- datasetNames[3]$results[,"Item"]



#odd way to get a dataset list...I just downloaded the source file
#and pull out their code.
#hmmm...is this the case for all genes?? WTF?
expandProbesets <- function (eset, sep = "///"){
  x <- lapply(featureNames(eset), function(x) strsplit(x, sep)[[1]])
  eset <- eset[order(sapply(x, length)), ]
  x <- lapply(featureNames(eset), function(x) strsplit(x, sep)[[1]])
  idx <- unlist(sapply(1:length(x), function(i) rep(i, length(x[[i]]))))
  xx <- !duplicated(unlist(x))
  idx <- idx[xx]
  x <- unlist(x)[xx]
  eset <- eset[idx, ]
  #Katie's notes: featureNames is from Biobase, for eset objects.
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
#ON SERVER: just download the tarbell to get package, code
#cd /home/kplaney/R/x86_64-redhat-linux-gnu-library
#sudo wget http://www.bioconductor.org/packages/release/data/experiment/src/contrib/curatedOvarianData_1.3.5.tar.gz
#tar xvzf curatedOvarianData_1.3.5.tar.gz
# strEsets <- list.files("/home/kplaney/R/x86_64-redhat-linux-gnu-library/curatedOvarianData/data/")
# #datalist is a text file...shouldn't be any other text files in this directory, only .rda files.
# strEsets <- strEsets[-which(strEsets=="datalist")]
# library("limma")
# names(strEsets) <- strsplit2(strEsets,split=".rda")[,1]
# #now load all of this data. equivalent of data() call above if can install package.
# #is having trouble reading from a connection?? what's up...
# #WTF..may have to manually step through this one..
# #weird...now suddently dataset 16 isn't working? why only THIS one?
# #makes me wonder if there is a space issue on the server.
# #I had to re-download the datasets...
# for (strEset in strEsets){
#   
#   load(paste0("/home/kplaney/R/x86_64-redhat-linux-gnu-library/curatedOvarianData/data/",strEset))
#   
# }
# #now call strEset again like above:
# strEsets <- ls(pattern="^.*_eset$")

##now back to common code for server or local computer:
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

#save for easier access next time:
#save(esets, file="/home/kplaney/curatedOvarianData_esetList.RData.gzip",compress="gzip")
#CHECK: remove any healthy samples. do this for ALL studies in case some have health tissue.
#for example: TCGA has some normal samples
table(pData(esets[["TCGA_eset"]])[,"batch"],pData(esets[["TCGA_eset"]])[,"sample_type"])
indicesRemove <- c()

for(e in 1:length(esets)){
  
  #are these all primary tumor? I wish they had kept the full TCGA labels.
  pDat <- pData(esets[[e]])[which(pData(esets[[e]])["sample_type"]=="tumor"),]
  expr <- exprs(esets[[e]])[ ,which(pData(esets[[e]])["sample_type"]=="tumor")]
  
  #if no samples left: remove
  if(ncol(expr)==0){
    
    indicesRemove <- append(indicesRemove,e)
    
  }else{
    
    exprs(esets[[e]]) <- expr
    pData(esets[[e]]) <- pDat
    protocolData(esets[[e]])@data <- data.frame(row.names=colnames(exprs(esets[[e]])))
    #check: make sure is still a valid object
    validObject(esets[[e]])
    
  }
  
}
message(paste("dataset indices to remove : ",indicesRemove))
#gut check: in TCGA, should be only tumor samples now!
if(any(unique(pData(esets[["TCGA_eset"]])[,"sample_type"]) != "tumor")){
  
  stop("\nNormal filtering appeared to not work.")
  
}

#check? any duplicated samples?
#it looks like unique_patient_ID ids duplicated samples, NOT the actual expression colnames?
#a little annoying - re-label the colnames to be this unique ID, if duplicated samples are found
#but esets can't have duplicated column names.
#in the end, looks like no samples were duplicated? 
for(e in 1:length(esets)){
  
  if(!all(is.na(pData(esets[[e]])[ ,"unique_patient_ID"]))){
    #for certain datasets: this is actually all NA values. so don't use then!
  if(any(duplicated(pData(esets[[e]])[,"unique_patient_ID"]))){
    
    #if(any(duplicated(colnames(exprs(esets[[e]]))))){
  
    
    message(paste0("Duplicated samples in study ",names(esets)[e]))
  
  }
  
}

}
#}
#remove datasets that are too small? allow them for now.
#now process, then grab meta-features
#sounds like there SHOULD be some duplicated samples  in here.
#want very lowly varying genes to be removed for Combat.
#COME BACK: do variance filtering by batch for Combat.
esets <- processExpressionSetList(exprSetList=esets,outputFileDirectory="/home/kplaney/ovarian_analysis/",
                                  minVar=.001,featureDataFieldName="gene",uniquePDataID="unique_patient_ID")

#just strictly limit TCGA dataset now - bc of combat
#esets[["TCGA_eset"]] <- processExpressionSet(exprSet,esets[["TCGA_eset"]],
#                                 outputFileDirectory="/home/kplaney/ovarian_analysis/",
#                                  minVar=.1,featureDataFieldName="gene",uniquePDataID="unique_patient_ID")
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
#na.omit: if any NA batch variables, then unique will return one "NA" 
TCGAbatches <- na.omit(unique(pData(esets[["TCGA_eset"]])[,"batch"]))
#any samples with NA batch? these will be thrown out.
which(is.na(pData(esets[["TCGA_eset"]])[,"batch"]))
#let's see if batches seem to have uneven distributions of known subtypes
#and/or grade...that will make me hesitant to use batch correction.
#hmm...so most are serious?
table(pData(esets[["TCGA_eset"]])[,"batch"],pData(esets[["TCGA_eset"]])[,"histological_type"])
unique(pData(esets[["TCGA_eset"]])[,"histological_type"])
#yep...only 10 a NA, rest are serous
length(which(is.na(pData(esets[["TCGA_eset"]])[,"histological_type"])))
#this matches the Nature paper saying the samples are only serous: http://www.nature.com/nature/journal/v474/n7353/full/nature10166.html

#and predominantly late-stage for all batches
table(pData(esets[["TCGA_eset"]])[,"batch"],pData(esets[["TCGA_eset"]])[,"summarystage"])
#same trend for substage
table(pData(esets[["TCGA_eset"]])[,"batch"],pData(esets[["TCGA_eset"]])[,"substage"])


#let's see: is there even a discernable batch effect?

splitOutTCGABatches <- FALSE


#if(splitOutTCGABatches){
  
  for(t in 1:length(TCGAbatches)){
    
    expr <- exprs(esets[["TCGA_eset"]])[ ,
                                        which(pData(esets[["TCGA_eset"]])[,"batch"]==TCGAbatches[t])]
    
    phenoData <-  pData(esets[["TCGA_eset"]])[
      which(pData(esets[["TCGA_eset"]])[,"batch"]==TCGAbatches[t]), ]
    
    esets[[paste0("TCGA_eset_batch_",t)]] <-  createS4exprSet(expr=expr,phenoData=phenoData,
                                                              featureDataFieldName="gene")
    
  }
  
  #remove original entire TCGA matrix? perhaps keep for now? then can do a comparison between the split-up groups and the non split-up groups.
  #esets <- esets[-which(names(esets)=="TCGA_eset")]
  
#}else{
  #just do batch correction
  #do this after removely very lowly varying genes.
#DID NOT RUN: test later....
  outcomesAndCovariates <- pData(esets[["TCGA_eset"]])[,"batch",drop=FALSE]
  rownames(outcomesAndCovariates) <- colnames(exprs(esets[["TCGA_eset"]]))
  exprMatrix <- exprs(esets[["TCGA_eset"]])
  rownames(exprMatrix) <- featureNames(esets[["TCGA_eset"]])
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_batchCorrection.R")
  batchCorrect <- batchNormalization(countsMatrixNoNANoDup=exprMatrix,outcomesAndCovariates=outcomesAndCovariates,MinInBatch=4,combatModelFactorName=NULL,pvalueThresh=.05,batchColName="batch",outputFile="combatoutput.txt")
  
esets[["TCGA_eset"]] <- esets[["TCGA_eset"]][rownames(batchCorrect$GEN_Data_Corrected), colnames(batchCorrect$GEN_Data_Corrected)]
featureNames(esets[["TCGA_eset"]] ) <- rownames(batchCorrect$GEN_Data_Corrected)
  validObject(esets[["TCGA_eset"]])
    
#}


#now format just as a list of data matrices.
dataMatrixList <- exprSetListToMatrixList(esets,featureDataFieldName="gene")

names(dataMatrixList) <- names(esets)
#run meta-feature analysis.
metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),numFeatSelectByGlobalRank=1000,
                                                   numTopFeatFromEachDataset=20,fractNATopFeatAllowedPerDataset=.2,selectMethod=c("median"),
                                                   outputFile="/home/kplaney/ovarian_analysis//selectFeaturesMetaVarianceOut.txt")

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))
#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]

}


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- metaFeatures$finalFeatures
  
}

output <- list(dataMatrixList=dataMatrixList,clustFeaturesList=clustFeaturesList)
save(output,file="/home/kplaney/ovarian_analysis/curatedBreastData_proc.RData.gzip",compress="gzip")


load("/home/kplaney/ovarian_analysis/curatedBreastData_proc.RData.gzip")
#remove end TCGA clusters
dataMatrixList <- output$dataMatrixList[-c(24:37)]
clustFeaturesList <- output$clustFeaturesList[c(24:37)]
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")
clusterOut_hclust <- clusterMatrixListHclustGap(dataMatrixList=dataMatrixList,clustFeaturesList=clustFeaturesList,maxNumClusters=15,algorithm="ward.D",
                                         distMethod=c("euclidean"),outputFile="/home/kplaney/ovarian_analysis//cluster_hclustGap_output.txt",
                                         corUse=c("everything"),numSims=100)

save(clusterOut_hclust,file="/home/kplaney/ovarian_analysis/clusterOut_hclust_noTCGAcombat.RData.gzip",compress="gzip")


clusterOut_kmeans <- clustMatrixListKmeansGap(dataMatrixList,clustFeaturesList,maxNumClusters=15,iter.max=30,nstart=25,numSims=100,algorithm="Hartigan-Wong",
                                              outputFile="/home/kplaney/ovarian_analysis//cluster_kmeansGap_output.txt"
                               )

clusterOut_kmeansConsensus <- consensusClusterMatrixList(dataMatrixList,clustFeaturesList,maxNumClusters = 15, numSims=100, pItem=0.8, pFeature=1, clusterAlg=c("km"),
                                                                     hclustAlgorithm="ward.D", consensusHclustAlgorithm="ward.D",
                                                                     corUse=c("everything"),
                                                                     distMethod=c("euclidean"),
                                                                     minClustConsensus=.7,bestKmethod=c("highestMinConsensus"))
 

save(clusterOut_kmeansConsensus,file="/home/kplaney/ovarian_analysis/clusterOut_kmeansConsensus_noTCGAcombat.RData.gzip",compress="gzip")


clusterOut_hclustConsensus <- consensusClusterMatrixList(dataMatrixList,clustFeaturesList,maxNumClusters = 15, numSims=100, pItem=0.8, pFeature=1, clusterAlg=c("hc"),
                                                         hclustAlgorithm="ward.D", consensusHclustAlgorithm="ward.D",
                                                                     corUse=c("everything"),
                                                                     distMethod=c("euclidean"),
                                                                     minClustConsensus=.7,bestKmethod=c("highestMinConsensus"))


save(clusterOut_hclustConsensus,file="/home/kplaney/ovarian_analysis/clusterOut_hclustConsensus_noTCGAcombat.RData.gzip",compress="gzip")


#compute edges
fractFeatIntersectThresh=.7
numFeatIntersectThresh=500
clustSizeThresh=3
clustSizeFractThresh=.05
numParallelCores=8
minTrueSimilThresh=.25
maxNullFractSize=.2
maxTrueSimilThresh=Inf
includeRefClustInNull=TRUE
numSims=100
sigMethod=c("centroid")
edgeMethod=c("pearson")


indEdgePvalueThresh=.3
meanEdgePairPvalueThresh=.2
minTrueSimilThresh = .25
maxTrueSimilThresh=Inf


###hclust consensus
load("/home/kplaney/ovarian_analysis/curatedBreastData_proc_TCGAnoBatch.RData.gzip")
dataMatrixList <- output$dataMatrixList
load("/home/kplaney/ovarian_analysis/clusterOut_hclustConsensus_noTCGAcombat.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")
#remove full TCGA dataset
indexRemove <- which(names(dataMatrixList)=="TCGA_eset")
dataMatrixList <- dataMatrixList[-indexRemove]
clustSampleIndexList <- clusterOut_hclustConsensus$clustSampleIndexList[-indexRemove]
clustFeatureIndexList <- clusterOut_hclustConsensus$clustFeatureIndexList[-indexRemove]


test <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                    edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                    sigMethod=sigMethod,maxNullFractSize=maxNullFractSize,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                    
                                    outputFile="/home/kplaney/ovarian_analysis//CoINcIDE_hclustConsensus_messages.txt",fractFeatIntersectThresh=fractFeatIntersectThresh,
                                numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, clustSizeFractThresh=clustSizeFractThresh)

hclustConsensus_getAdjMatricesOut_centroid <- test
save(hclustConsensus_getAdjMatricesOut_centroid,file="/home/kplaney/ovarian_analysis/hclustConsensus_getAdjMatricesOut_centroid.RData.gzip",compress="gzip")

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_communityDetection.R")
load("/home/kplaney/ovarian_analysis/hclustConsensus_getAdjMatricesOut_centroid.RData.gzip")
testEdges <- assignFinalEdges(computeTrueSimilOutput=hclustConsensus_getAdjMatricesOut_centroid$computeTrueSimilOutput,
                              pvalueMatrix=hclustConsensus_getAdjMatricesOut_centroid$pvalueMatrix,indEdgePvalueThresh=indEdgePvalueThresh,
                              meanEdgePairPvalueThresh=meanEdgePairPvalueThresh,
                              minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                              minFractFeatureIntersect=minFractFeatureIntersect,fractFeatIntersectThresh=fractFeatIntersectThresh,numFeatIntersectThresh=numFeatIntersectThresh ,
                              clustSizeThresh=clustSizeThresh, clustSizeFractThresh= clustSizeFractThresh,saveDir="/home/kplaney/ovarian_analysis/",fileTag="CoINcIDE_edges_hclustConsensus"
)

communityInfo <- findCommunities(edgeMatrix=testEdges$filterEdgeOutput$edgeMatrix,edgeWeightMatrix=testEdges$filterEdgeOutput$edgeWeightMatrix,hclustConsensus_getAdjMatricesOut_centroid$clustIndexMatrix,fileTag="hclustConsensus_centroid",
                                 saveDir="/home/kplaney/ovarian_analysis//",minNumUniqueStudiesPerCommunity=2,clustMethodName="hclustConsensus_centroid",
                                 commMethod=c("edgeBetween"),
                                 makePlots=TRUE,saveGraphData=FALSE)

##kmeans consensus
load("/home/kplaney/ovarian_analysis/clusterOut_kmeansConsensus_noTCGAcombat.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")
load("/home/kplaney/ovarian_analysis/curatedBreastData_proc_TCGAnoBatch.RData.gzip")
dataMatrixList <- output$dataMatrixList
indexRemove <- which(names(dataMatrixList)=="TCGA_eset")
dataMatrixList <- dataMatrixList[-indexRemove]
clustSampleIndexList <- clusterOut_kmeansConsensus$clustSampleIndexList[-indexRemove]
clustFeatureIndexList <- clusterOut_kmeansConsensus$clustFeatureIndexList[-indexRemove]

test <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                sigMethod=sigMethod,maxNullFractSize=maxNullFractSize,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                
                                outputFile="/home/kplaney/ovarian_analysis//CoINcIDE_kmeansConsensus_messages.txt",fractFeatIntersectThresh=fractFeatIntersectThresh,
                                numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, clustSizeFractThresh=clustSizeFractThresh)
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_communityDetection.R")

testEdges <- assignFinalEdges(computeTrueSimilOutput=test$computeTrueSimilOutput,pvalueMatrix=test$pvalueMatrix,indEdgePvalueThresh=indEdgePvalueThresh,
                              meanEdgePairPvalueThresh=meanEdgePairPvalueThresh,
                              minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                              minFractFeatureIntersect=minFractFeatureIntersect,fractFeatIntersectThresh=fractFeatIntersectThresh,numFeatIntersectThresh=numFeatIntersectThresh ,
                              clustSizeThresh=clustSizeThresh, clustSizeFractThresh= clustSizeFractThresh,saveDir="/home/kplaney/ovarian_analysis/",fileTag="CoINcIDE_edges_kmeansConsensus"
)


##regular old hclust
#####have not run yet:
load("/home/kplaney/ovarian_analysis/clusterOut_hclust_noTCGAcombat.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")
load("/home/kplaney/ovarian_analysis/curatedBreastData_proc_TCGAnoBatch.RData.gzip")
dataMatrixList <- output$dataMatrixList
indexRemove <- which(names(dataMatrixList)=="TCGA_eset")
dataMatrixList <- dataMatrixList[-indexRemove]
clustSampleIndexList <- clusterOut_hclust$clustSampleIndexList[-indexRemove]
clustFeatureIndexList <- clusterOut_hclust$clustFeatureIndexList[-indexRemove]

test <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                sigMethod=sigMethod,maxNullFractSize=maxNullFractSize,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                
                                outputFile="/home/kplaney/ovarian_analysis//CoINcIDE_hclust_messages.txt",fractFeatIntersectThresh=fractFeatIntersectThresh,
                                numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, clustSizeFractThresh=clustSizeFractThresh)
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_communityDetection.R")

testEdges <- assignFinalEdges(computeTrueSimilOutput=test$computeTrueSimilOutput,pvalueMatrix=test$pvalueMatrix,indEdgePvalueThresh=indEdgePvalueThresh,
                              meanEdgePairPvalueThresh=meanEdgePairPvalueThresh,
                              minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                              minFractFeatureIntersect=minFractFeatureIntersect,fractFeatIntersectThresh=fractFeatIntersectThresh,numFeatIntersectThresh=numFeatIntersectThresh ,
                              clustSizeThresh=clustSizeThresh, clustSizeFractThresh= clustSizeFractThresh,saveDir="/home/kplaney/ovarian_analysis/",fileTag="CoINcIDE_edges_hclust"
)



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
