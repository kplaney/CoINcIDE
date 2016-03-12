#!/usr/bin/Rscript --default-packages=utils

library("Coincide")
#options(device=NULL) 

#grab arguments
args <- commandArgs(TRUE)
edgeMethod <- args[1]
centroidMethod <- args[2]
#these are full paths to rds objects (saved using saveRDS)
clusterData <- args[3]
dataMatrixList <- args[4]

experimentName <- args[5]
saveDir <- args[6]

message("edge method is ",edgeMethod)
message("centroid method is ", centroidMethod)
message("cluster data object is ", clusterData)
message("dataMatrixList object is ", dataMatrixList)
message("save directory is ",saveDir)
message("experiment name is ",experimentName)
message("running CoINcIDE for null datasets")

clusterData <- readRDS(clusterData)
dataMatrixList <- readRDS(dataMatrixList)
 
 message("Loaded data")

clustSampleIndexList <- clusterData$clustSampleIndexList_PACR
clustFeatureIndexList <- clusterData$clustFeatureIndexList_PACR

#some datasets may have not been used if not enough genes. so remove them from master data list
#so indices match up:
indicesKeep <- na.omit(match(names(clustFeatureIndexList),names(dataMatrixList)))
dataMatrixList <- dataMatrixList[indicesKeep]

if((length(dataMatrixList) != length(clustFeatureIndexList)) || (length(dataMatrixList) != length(clustSampleIndexList))){
  
  stop("Not indexing correctly.")
  
}
#just make these default for now.
outputFile <- "~/CoINcIDE_messages.txt"
numSims <- 500
numNullIter <- 10

CoINcIDE_nullOutput <- computeAdjMatricesNullMatrixList(dataMatrixList,numIter=numNullIter,
                                          clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                          edgeMethod=edgeMethod,
                                          numSims=numSims,
                                          outputFile=outputFile,
                                          centroidMethod=centroidMethod)


saveRDS(CoINcIDE_nullOutput,file=paste0(saveDir,"/CoINcIDE_NullOutput_",experimentName,"_",edgeMethod,"edgeMethod_",centroidMethod,"_centroidMethod",Sys.Date(),".rds"),compress=TRUE)
  
 message("saved these files:")
# message(paste0(saveDir,"/CoINcIDE_results_",experimentName,"_",edgeMethod,"_edgeMethod_",centroidMethod,"_centroidMethod",Sys.Date(),".rds"))
 message(paste0(saveDir,"/CoINcIDE_Nullresults_",experimentName,"_",edgeMethod,"_edgeMethod_",centroidMethod,"_centroidMethod",Sys.Date(),".rds"))
