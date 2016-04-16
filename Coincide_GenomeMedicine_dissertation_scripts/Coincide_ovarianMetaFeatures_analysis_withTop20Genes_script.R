#!/usr/bin/Rscript --default-packages=utils

library("CoINcIDE")
source("/home/ywrfc09/CoINcIDE/coincide/GenomeBiology/GenomeBiology_metaFeatures_analysis_wrapper.R")   
outputFile <- "~/CoINcIDE_messages.txt"
saveDir <- "/home/ywrfc09/ovarian_analysis/metaRankWithTop20Genes"
globalSaveDir <- "/home/ywrfc09/ovarian_analysis/"
#200 features, pearson:
load(paste0(saveDir,"/metaFeatures_200.RData.gzip"))
clusterCoINcIDE_output =  readRDS("/home/ywrfc09/ovarian_analysis/metaRankWithTop20Genes/curatedOvarianData_kmeansConsensus_nstart1_200Features_2015-04-28.rds")


clustSampleIndexList <-  clusterCoINcIDE_output$clustSampleIndexList_PACR
clustFeatureIndexList <- clusterCoINcIDE_output$clustFeatureIndexList_PACR

load(paste0(globalSaveDir,"/esets_proc_TCGAcombat.RData.gzip"))
esets = esets

CoINcIDE_output = readRDS("/home/ywrfc09/ovarian_analysis/metaRankWithTop20Genes/CoINcIDE_results_ovarian240F_pearson_edgeMethod_mean_centroidMethod2015-07-08.rds")
experimentName <- "240F_pear_meanCent_MM5"
eset_featureDataFieldName="gene"
networkColors = "Set2"
outcomesVarBinary="vital_status"
outcomesVarCont = "days_to_death"
ovarian <- TRUE
fisherTestVariables <- c("histological_type","tumorstage","recurrence_status","grade","age_at_initial_pathologic_diagnosis")
fisherTestVariableLegendNames <- c("hist\ntype","tumor\nstage","recurrence","tumor\ngrade","age")
fisherTestVariableTitleNames <- c("histological type","tumor stage", "recurrence status","tumor grade","age at diagnosis")

#breast: must also input dataMatrixList
#brewPal = c("Set3","Paired","Spectral","BrBG","PiYG","RdYlGn","RdYlBu","RdBu","PiYG","Set2"),


#looks like simil thresh of .5 appropriate
densityPlot <- meanMetricDensityPlot(CoINcIDE_output$meanMetricMatrix,saveDir=saveDir,
                                     experimentName=experimentName,savePlot=TRUE,yLimit=2.5)




ovarian_240genes_pearson_meanCentroid_analysis <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                                              meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .01, minTrueSimilThresh = .5, maxTrueSimilThresh = Inf,minFractNN=.8,
                                                                              clustSizeThresh = 0,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                                              commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=3, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                                              survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                                              CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                                              ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames,
                                                                              GSEAanalysis=TRUE,clinVarPlots=TRUE, fractFeatIntersectThresh=.6,numFeatIntersectThresh =0,clustSizeFractThresh =0,
                                                                              findCommWithWeights=TRUE, plotSimilEdgeWeight = TRUE,plotToScreen=FALSE,fractEdgesInVsOutEdge=0, fractEdgesInVsOutComm=0,minNumEdgesForCluster=1)





CoINcIDE_nullOutput <- readRDS("/home/ywrfc09/ovarian_analysis/metaRankWithTop20Genes/CoINcIDE_NullOutput_ovarian240F_pearsonedgeMethod_mean_centroidMethod2015-07-09.rds")
globalFDR_results <- globalFDR(CoINcIDE_outputList=CoINcIDE_nullOutput$CoINcIDE_NullOutputList,
                               edgeMethod="pearson",minTrueSimilThresh=0.5,maxTrueSimilThresh=Inf,
                               outputFile=outputFile,fractFeatIntersectThresh=0,numFeatIntersectThresh=0 ,clustSizeThresh=0, clustSizeFractThresh=0,
                               meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .01, 
                               saveDir = saveDir,experimentName = "nullTest",
                               commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=3,minFractNN =.8,
                               findCommWithWeights=TRUE,minNumEdgesForCluster=1,fractEdgesInVsOutComm=0,fractEdgesInVsOutEdge=0)

saveRDS(globalFDR_results,file=paste0(saveDir,"/CoINcIDE_globalFDRresults_",experimentName,"_",edgeMethod,"edgeMethod_",centroidMethod,"_centroidMethod_minTrueSimil",minTrueSimilThresh,"_",Sys.Date(),".rds"),compress=TRUE)

options(bitmapType="cairo")
png(filename=paste0(saveDir,"/",experimentName,"_",Sys.Date(),"/dendogram_",Sys.Date(),".png"),width=1800,height=1000,res=160)

plot(ovarian_240genes_pearson_meanCentroid_analysis$commInfo$unPrunedDendogram)

dev.off()

finalNodeMatrix <-  ovarian_240genes_pearson_meanCentroid_analysis$commInfo$attrDF
origEdgeMatrix <- ovarian_240genes_pearson_meanCentroid_analysis$finalEdgeInfo$filterEdgeOutput$edgeMatrix
origEdgeWeightsMatrix <- ovarian_240genes_pearson_meanCentroid_analysis$finalEdgeInfo$filterEdgeOutput$edgeWeightMatrix
#right now these are igraph ids: will need to change that.
finalEdgeMatrix <- ovarian_240genes_pearson_meanCentroid_analysis$commInfo$edgeDF[,c(1:2)]
clustIndexMatrix <- ovarian_240genes_pearson_meanCentroid_analysis$CoINcIDE_computeEdgesObject$clustIndexMatrix


fractLeaveOutVector <- seq(0.1,0.9,by=.1)  
leaveXOutResults <- list()
numIter <- 100

for(f in 1:length( fractLeaveOutVector)){
  message("noiseLevel ",f)
  leaveXOutResults[[f]] <- networkLeaveOutAnalysis(finalNodeMatrix=finalNodeMatrix, 
                                                   origEdgeMatrix=origEdgeMatrix,
                                                   origEdgeWeightsMatrix=origEdgeWeightsMatrix,
                                                   finalEdgeMarix=finalEdgeMarix,fractLeaveOut=fractLeaveOutVector[f],
                                                   numIter=numIter,commMethod="edgeBetween",
                                                   findCommWithWeights=TRUE,clustIndexMatrix=clustIndexMatrix)
  
}

names(leaveXOutResults) <- as.character(  fractLeaveOutVector)
saveDir <- paste0(saveDir,"/",experimentName,"_",Sys.Date())
#date was 2015-08-16
saveRDS(leaveXOutResults,file=paste0(saveDir,"/CoINcIDE_LeaveXOutAnalysis_",experimentName,"_",Sys.Date(),".rds"),compress=TRUE)
source('~/CoINcIDE/coincide/CoINcIDE_packageVersion/CoINcIDE/R/CoINcIDE_metaClusterAnalysis.R')
LO_analysis <- plotLeaveXOutAnalysis(leaveXOutResults,
                                     experimentName=experimentName,saveDir=saveDir)


#also look at when vary simil ranges
varySimil <- networkVaryMinSimilAnalysis(finalNodeMatrix, origEdgeMatrix,
                                         origEdgeWeightsMatrix,finalEdgeMarix,
                                         minSimilThreshVector=c(0.6,0.7,0.8,0.9,1.0),
                                         commMethod="edgeBetween",saveDir=saveDir,
                                         findCommWithWeights=TRUE,clustIndexMatrix,
                                         makePlots=TRUE,experimentName= experimentName)

saveRDS(varySimil,file=paste0(saveDir,"/varySimilResults.rds"),compress=TRUE)


#2014
library("CoINcIDE")
source("/home/ywrfc09/CoINcIDE/coincide/GenomeBiology/GenomeBiology_metaFeatures_analysis_wrapper.R")   
outputFile <- "~/CoINcIDE_messages.txt"
saveDir <- "/home/ywrfc09/ovarian_analysis/metaRankWithTop20Genes"
globalSaveDir <- "/home/ywrfc09/ovarian_analysis/"
#200 features, pearson:
load(paste0(saveDir,"/metaFeatures_2000.RData.gzip"))
clusterCoINcIDE_output =  readRDS("/home/ywrfc09/ovarian_analysis/metaRankWithTop20Genes/curatedOvarianData_kmeansConsensus_nstart1_2000Features_2015-04-30.rds")


clustSampleIndexList <-  clusterCoINcIDE_output$clustSampleIndexList_PACR
clustFeatureIndexList <- clusterCoINcIDE_output$clustFeatureIndexList_PACR

load(paste0(globalSaveDir,"/esets_proc_TCGAcombat.RData.gzip"))
esets = esets

CoINcIDE_output = readRDS("/home/ywrfc09/ovarian_analysis/metaRankWithTop20Genes/CoINcIDE_results_ovarian2014F_pearson_edgeMethod_mean_centroidMethod2015-07-08.rds")
experimentName <- "2014F_pear_meanCent_MM5"
eset_featureDataFieldName="gene"
networkColors = "Set2"
outcomesVarBinary="vital_status"
outcomesVarCont = "days_to_death"
ovarian <- TRUE
fisherTestVariables <- c("histological_type","tumorstage","recurrence_status","grade","age_at_initial_pathologic_diagnosis")
fisherTestVariableLegendNames <- c("hist\ntype","tumor\nstage","recurrence","tumor\ngrade","age")
fisherTestVariableTitleNames <- c("histological type","tumor stage", "recurrence status","tumor grade","age at diagnosis")

#breast: must also input dataMatrixList
#brewPal = c("Set3","Paired","Spectral","BrBG","PiYG","RdYlGn","RdYlBu","RdBu","PiYG","Set2"),

#looks like .5 or .7 appropriate
densityPlot <- meanMetricDensityPlot(CoINcIDE_output$meanMetricMatrix,saveDir=saveDir,
                                     experimentName=experimentName,savePlot=TRUE,yLimit=2.5)



ovarian_2014genes_pearson_meanCentroid_analysis <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                                              meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .01, minTrueSimilThresh = .5, maxTrueSimilThresh = Inf,minFractNN=.8,
                                                                              clustSizeThresh = 0,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                                              commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=3, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                                              survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                                              CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                                              ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames,
                                                                              GSEAanalysis=TRUE,clinVarPlots=TRUE, fractFeatIntersectThresh=.6,numFeatIntersectThresh =0,clustSizeFractThresh =0,
                                                                              findCommWithWeights=TRUE, plotSimilEdgeWeight = TRUE,plotToScreen=FALSE,fractEdgesInVsOutEdge=0, fractEdgesInVsOutComm=0,minNumEdgesForCluster=1)





CoINcIDE_nullOutput <- readRDS("/home/ywrfc09/ovarian_analysis/metaRankWithTop20Genes/CoINcIDE_NullOutput_ovarian2014F_pearsonedgeMethod_mean_centroidMethod2015-07-10.rds")

globalFDR_results <- globalFDR(CoINcIDE_outputList=CoINcIDE_nullOutput$CoINcIDE_NullOutputList,
                               edgeMethod="pearson",minTrueSimilThresh=0.5,maxTrueSimilThresh=Inf,
                               outputFile=outputFile,fractFeatIntersectThresh=0,numFeatIntersectThresh=0 ,clustSizeThresh=0, clustSizeFractThresh=0,
                               meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .01, 
                               saveDir = saveDir,experimentName = "nullTest",
                               commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=3,minFractNN =.8,
                               findCommWithWeights=TRUE,minNumEdgesForCluster=1,fractEdgesInVsOutComm=0,fractEdgesInVsOutEdge=0)

saveRDS(globalFDR_results,file=paste0(saveDir,"/CoINcIDE_globalFDRresults_",experimentName,"_",edgeMethod,"edgeMethod_",centroidMethod,"_centroidMethod_minTrueSimil",minTrueSimilThresh,"_",Sys.Date(),".rds"),compress=TRUE)

####Julia's feedback: remove TCGA, only keep serous samples.
table(ovarian_2014genes_pearson_meanCentroid_analysis$sampleClustCommPhenoData$histological_type,ovarian_2014genes_pearson_meanCentroid_analysis$sampleClustCommPhenoData$studyNum)
#we see that the studies below are 100% serous:
#NO #23
#serousDataset <- c(1,2,3,5,6,12,14,15,20,21,22,23)
#so what cluster numbers does this translate to?
# serousClusters <- c()
# 
# for(s in 1:length(serousDataset)){
#   
# tmp <- ovarian_2014genes_pearson_meanCentroid_analysis$sampleClustCommPhenoData$globalClustNum[
#   which(ovarian_2014genes_pearson_meanCentroid_analysis$sampleClustCommPhenoData$studyNum==serousDataset[s])]
# serousClusters  <- append(tmp,serousClusters )
# }
#serousClusters <- unique(serousClusters)

TCGA_clust <- ovarian_2014genes_pearson_meanCentroid_analysis$sampleClustCommPhenoData$globalClustNum[
     which(ovarian_2014genes_pearson_meanCentroid_analysis$sampleClustCommPhenoData$studyNum==23)]

TCGA_clust <- unique(TCGA_clust)
#subset:
tmp_CoINcIDE_output=CoINcIDE_output
#set to NA
tmp_CoINcIDE_output$pvalueMatrix[TCGA_clust,] <- NA
tmp_CoINcIDE_output$pvalueMatrix[,TCGA_clust] <- NA
noTCGA_ovarian_2014genes_pearson_meanCentroid_analysis <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=tmp_CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                                               meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .01, minTrueSimilThresh = .5, maxTrueSimilThresh = Inf,minFractNN=.8,
                                                                               clustSizeThresh = 0,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                                               commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=3, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                                               survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                                               CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                                               ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames,
                                                                               GSEAanalysis=FALSE,clinVarPlots=TRUE, fractFeatIntersectThresh=.6,numFeatIntersectThresh =0,clustSizeFractThresh =0,
                                                                               findCommWithWeights=TRUE, plotSimilEdgeWeight = TRUE,plotToScreen=FALSE,fractEdgesInVsOutEdge=0, fractEdgesInVsOutComm=0,minNumEdgesForCluster=1)



#how did meta-cluster status of serous clusters change?
noTCGA_membership <- noTCGA_ovarian_2014genes_pearson_meanCentroid_analysis$sampleClustCommPhenoData[,c("sampleName","globalClustNum","community","studyNum")]
noTCGA_membership <- noTCGA_membership[-which(noTCGA_membership$community==4),]
noTCGA_membership <- noTCGA_membership[-which(noTCGA_membership$community==5),]
noTCGA_membership <- noTCGA_membership[-which(noTCGA_membership$community==6),]
noTCGA_membership <- noTCGA_membership[-which(is.na(noTCGA_membership$community)),]

yesTCGA_membership <- ovarian_2014genes_pearson_meanCentroid_analysis$sampleClustCommPhenoData[,c("sampleName","globalClustNum","community","studyNum")]
yesTCGA_membership <- yesTCGA_membership[-which(yesTCGA_membership$community==4),]
yesTCGA_membership <- yesTCGA_membership[-which(yesTCGA_membership$community==5),]
yesTCGA_membership <- yesTCGA_membership[-which(yesTCGA_membership$community==6),]
yesTCGA_membership <- yesTCGA_membership[-which(is.na(yesTCGA_membership$community)),]
yesTCGA_membership <- yesTCGA_membership[-which(yesTCGA_membership$studyNum==23),]

if(any(is.na(match(noTCGA_membership$sampleName,yesTCGA_membership$sampleName)))){
  
  stop("Error - some clusters are not in both datasets")
}

noTCGA_membership <- noTCGA_membership[match(noTCGA_membership$sampleName,yesTCGA_membership$sampleName), ]
all(noTCGA_membership$sampleName==yesTCGA_membership$sampleName)
#yay - all the same!
which(noTCGA_membership$community!=yesTCGA_membership$community)

#also try full set with 0.7 (other global max)
experimentName <- "2014F_pear_meanCent_MM7"
ovarian_2014genes_pearson_meanCentroid_analysis7 <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                                               meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .01, minTrueSimilThresh = .7, maxTrueSimilThresh = Inf,minFractNN=.8,
                                                                               clustSizeThresh = 0,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                                               commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=3, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                                               survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                                               CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                                               ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames,
                                                                               GSEAanalysis=TRUE,clinVarPlots=TRUE, fractFeatIntersectThresh=.6,numFeatIntersectThresh =0,clustSizeFractThresh =0,
                                                                               findCommWithWeights=TRUE, plotSimilEdgeWeight = TRUE,plotToScreen=FALSE,fractEdgesInVsOutEdge=0, fractEdgesInVsOutComm=0,minNumEdgesForCluster=1)







options(bitmapType="cairo")
png(filename=paste0(saveDir,"/",experimentName,"_",Sys.Date(),"/dendogram_",Sys.Date(),".png"),width=1800,height=1000,res=160)

plot(ovarian_2014genes_pearson_meanCentroid_analysis$commInfo$unPrunedDendogram)

dev.off()

finalNodeMatrix <-  ovarian_2014genes_pearson_meanCentroid_analysis$commInfo$attrDF
origEdgeMatrix <- ovarian_2014genes_pearson_meanCentroid_analysis$finalEdgeInfo$filterEdgeOutput$edgeMatrix
origEdgeWeightsMatrix <- ovarian_2014genes_pearson_meanCentroid_analysis$finalEdgeInfo$filterEdgeOutput$edgeWeightMatrix
#right now these are igraph ids: will need to change that.
finalEdgeMatrix <- ovarian_2014genes_pearson_meanCentroid_analysis$commInfo$edgeDF[,c(1:2)]
clustIndexMatrix <- ovarian_2014genes_pearson_meanCentroid_analysis$CoINcIDE_computeEdgesObject$clustIndexMatrix


fractLeaveOutVector <- seq(0.1,0.9,by=.1)  
leaveXOutResults <- list()
numIter <- 100

for(f in 1:length( fractLeaveOutVector)){
  message("noiseLevel ",f)
  leaveXOutResults[[f]] <- networkLeaveOutAnalysis(finalNodeMatrix=finalNodeMatrix, 
                                                   origEdgeMatrix=origEdgeMatrix,
                                                   origEdgeWeightsMatrix=origEdgeWeightsMatrix,
                                                   finalEdgeMarix=finalEdgeMarix,fractLeaveOut=fractLeaveOutVector[f],
                                                   numIter=numIter,commMethod="edgeBetween",
                                                   findCommWithWeights=TRUE,clustIndexMatrix=clustIndexMatrix)
  
}

names(leaveXOutResults) <- as.character(  fractLeaveOutVector)

saveDir <- paste0(saveDir,"/",experimentName,"_",Sys.Date())
saveRDS(leaveXOutResults,file=paste0(saveDir,"/CoINcIDE_LeaveXOutAnalysis_",experimentName,"_",Sys.Date(),".rds"),compress=TRUE)
source('~/CoINcIDE/coincide/CoINcIDE_packageVersion/CoINcIDE/R/CoINcIDE_metaClusterAnalysis.R')
LO_analysis <- plotLeaveXOutAnalysis(leaveXOutResults,
                                     experimentName=experimentName,saveDir=saveDir)


#also look at when vary simil ranges
varySimil <- networkVaryMinSimilAnalysis(finalNodeMatrix, origEdgeMatrix,
                                         origEdgeWeightsMatrix,finalEdgeMarix,
                                         minSimilThreshVector=c(0.6,0.7,0.8,0.9,1.0),
                                         commMethod="edgeBetween",saveDir=saveDir,
                                         findCommWithWeights=TRUE,clustIndexMatrix,
                                         makePlots=TRUE,experimentName= experimentName)

saveRDS(varySimil,file=paste0(saveDir,"/varySimilResults.rds"),compress=TRUE)


############
##500plus features

source("/home/ywrfc09/CoINcIDE/coincide/CoINcIDE_packageVersion//CoINcIDE/R/CoINcIDE_geneExprProcess.R")
source("/home/ywrfc09/CoINcIDE/coincide/CoINcIDE_packageVersion//CoINcIDE/R/CoINcIDE_computeEdges.R")
source("/home/ywrfc09/CoINcIDE/coincide/CoINcIDE_packageVersion/CoINcIDE/R/CoINcIDE_metaClusterAnalysis.R")
source("/home/ywrfc09/CoINcIDE/coincide/GenomeBiology/GenomeBiology_metaFeatures_analysis_wrapper.R")   
saveDir <- "/home/ywrfc09/ovarian_analysis/metaRankWithTop20Genes"
globalSaveDir <- "/home/ywrfc09/ovarian_analysis/"
#200 features, pearson:
edgeMethod <- "pearson"
load(paste0(saveDir,"/curatedOvarianData_kmeansConsensus_nstart1_500Features_2015-04-29.RData.gzip"))
load(paste0(saveDir,"/metaFeatures_500.RData.gzip"))

load(paste0(globalSaveDir,"/esets_proc_TCGAcombat.RData.gzip"))

#now format just as a list of data matrices.
dataMatrixList <- exprSetListToMatrixList(esets,featureDataFieldName="gene")

names(dataMatrixList) <- names(esets)

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

clustSampleIndexList <-  kmeansConsensus$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensus$clustFeatureIndexList_PACR

clusterCoINcIDE_output =  kmeansConsensus

edgeMethod <- "pearson"
centroidMethod <- "mean"
outputFile <- "~/CoINcIDE_messages.txt"
numSims <- 500
source("/home/ywrfc09/CoINcIDE/coincide/CoINcIDE_packageVersion/CoINcIDE/R/CoINcIDE_computeEdges.R")
ovarian_500genes_pearson_meanCentroid <-CoINcIDE_getAdjMatrices(dataMatrixList,clustSampleIndexList,clustFeatureIndexList,
                                                                 edgeMethod=edgeMethod,centroidMethod=centroidMethod,
                                                                 numSims=500,
                                                                 outputFile=outputFile)


save(ovarian_500genes_pearson_meanCentroid,file=paste0(saveDir,"ovarian_500genes_pearson_meanCentroid.RData.gzip"),compress="gzip")

###null tests
CoINcIDE_nullOutputList <- create_nullMatrixList_results(dataMatrixList=dataMatrixList,numIter=5,
                                                         clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                         edgeMethod=edgeMethod,
                                                         numSims=numSims,
                                                         outputFile=outputFile,
                                                         centroidMethod=centroidMethod)





FDRstats_500Features <- global_FDR(CoINcIDE_outputList=CoINcIDE_nullOutputList,
                                    minTrueSimilThresh=.4,maxTrueSimilThresh=Inf,
                                    outputFile=outputFile,fractFeatIntersectThresh=.8,numFeatIntersectThresh=0 ,clustSizeThresh=0, clustSizeFractThresh=0,
                                    meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, 
                                    saveDir = saveDir,experimentName = "nullTest240Features",
                                    commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=3,minFractNN =.8,findCommWithWeights=TRUE)


##now analyze results:
load(paste0(globalSaveDir,"/esets_proc_TCGAcombat.RData.gzip"))
esets = esets
#dataMatrixList
load(paste0(saveDir,"ovarian_500genes_pearson_meanCentroid.RData.gzip"))
CoINcIDE_output = ovarian_500genes_pearson_meanCentroid
experimentName <- "ovarian_500_featuresWithGSEA"
eset_featureDataFieldName="gene"
networkColors = "Set2"
outcomesVarBinary="vital_status"
outcomesVarCont = "days_to_death"
ovarian <- TRUE
fisherTestVariables <- c("histological_type","tumorstage","recurrence_status","grade","age_at_initial_pathologic_diagnosis")
fisherTestVariableLegendNames <- c("hist\ntype","tumor\nstage","recurrence","tumor\ngrade","age")
fisherTestVariableTitleNames <- c("histological type","tumor stage", "recurrence status","tumor grade","age at diagnosis")

#breast: must also input dataMatrixList
#brewPal = c("Set3","Paired","Spectral","BrBG","PiYG","RdYlGn","RdYlBu","RdBu","PiYG","Set2"),



ovarian_500genes_pearson_meanCentroid_analysis <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                                               meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .01, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,minFractNN=.8,
                                                                               clustSizeThresh = 0,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                                               commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=3, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                                               survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                                               CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                                               ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames,
                                                                               GSEAanalysis=FALSE,clinVarPlots=FALSE, fractFeatIntersectThresh=.8,numFeatIntersectThresh =0,clustSizeFractThresh =0,
                                                                               findCommWithWeights=TRUE, plotSimilEdgeWeight = TRUE,plotToScreen=TRUE,minMedianNumEdgesPerNodeInCommunity=1)








###1000 plus features

source("/home/ywrfc09/CoINcIDE/coincide/CoINcIDE_packageVersion//CoINcIDE/R/CoINcIDE_geneExprProcess.R")
source("/home/ywrfc09/CoINcIDE/coincide/CoINcIDE_packageVersion//CoINcIDE/R/CoINcIDE_computeEdges.R")
source("/home/ywrfc09/CoINcIDE/coincide/CoINcIDE_packageVersion/CoINcIDE/R/CoINcIDE_metaClusterAnalysis.R")
source("/home/ywrfc09/CoINcIDE/coincide/GenomeBiology/GenomeBiology_metaFeatures_analysis_wrapper.R")   
saveDir <- "/home/ywrfc09/ovarian_analysis/metaRankWithTop20Genes"
globalSaveDir <- "/home/ywrfc09/ovarian_analysis/"
#200 features, pearson:
edgeMethod <- "pearson"
load(paste0(saveDir,"/curatedOvarianData_kmeansConsensus_nstart1_1000Features_2015-04-29.RData.gzip"))
load(paste0(saveDir,"/metaFeatures_1000.RData.gzip"))

load(paste0(globalSaveDir,"/esets_proc_TCGAcombat.RData.gzip"))

#now format just as a list of data matrices.
dataMatrixList <- exprSetListToMatrixList(esets,featureDataFieldName="gene")

names(dataMatrixList) <- names(esets)

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

clustSampleIndexList <-  kmeansConsensus$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensus$clustFeatureIndexList_PACR

clusterCoINcIDE_output =  kmeansConsensus

edgeMethod <- "pearson"
centroidMethod <- "mean"
outputFile <- "~/CoINcIDE_messages.txt"
numSims <- 500
source("/home/ywrfc09/CoINcIDE/coincide/CoINcIDE_packageVersion/CoINcIDE/R/CoINcIDE_computeEdges.R")
ovarian_1000genes_pearson_meanCentroid <-CoINcIDE_getAdjMatrices(dataMatrixList,clustSampleIndexList,clustFeatureIndexList,
                                                                edgeMethod=edgeMethod,centroidMethod=centroidMethod,
                                                                numSims=500,
                                                                outputFile=outputFile)


save(ovarian_1000genes_pearson_meanCentroid,file=paste0(saveDir,"ovarian_1000genes_pearson_meanCentroid.RData.gzip"),compress="gzip")

###null tests
CoINcIDE_nullOutputList <- create_nullMatrixList_results(dataMatrixList=dataMatrixList,numIter=5,
                                                         clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                         edgeMethod=edgeMethod,
                                                         numSims=numSims,
                                                         outputFile=outputFile,
                                                         centroidMethod=centroidMethod)





FDRstats_1000Features <- global_FDR(CoINcIDE_outputList=CoINcIDE_nullOutputList,
                                   minTrueSimilThresh=.4,maxTrueSimilThresh=Inf,
                                   outputFile=outputFile,fractFeatIntersectThresh=.8,numFeatIntersectThresh=0 ,clustSizeThresh=0, clustSizeFractThresh=0,
                                   meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, 
                                   saveDir = saveDir,experimentName = "nullTest240Features",
                                   commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=3,minFractNN =.8,findCommWithWeights=TRUE)


##now analyze results:
load(paste0(globalSaveDir,"/esets_proc_TCGAcombat.RData.gzip"))
esets = esets
#dataMatrixList
load(paste0(saveDir,"ovarian_1000genes_pearson_meanCentroid.RData.gzip"))
CoINcIDE_output = ovarian_1000genes_pearson_meanCentroid
experimentName <- "ovarian_1000_features"
eset_featureDataFieldName="gene"
networkColors = "Set2"
outcomesVarBinary="vital_status"
outcomesVarCont = "days_to_death"
ovarian <- TRUE
fisherTestVariables <- c("histological_type","tumorstage","recurrence_status","grade","age_at_initial_pathologic_diagnosis")
fisherTestVariableLegendNames <- c("hist\ntype","tumor\nstage","recurrence","tumor\ngrade","age")
fisherTestVariableTitleNames <- c("histological type","tumor stage", "recurrence status","tumor grade","age at diagnosis")

#breast: must also input dataMatrixList
#brewPal = c("Set3","Paired","Spectral","BrBG","PiYG","RdYlGn","RdYlBu","RdBu","PiYG","Set2"),



ovarian_1000genes_pearson_meanCentroid_analysis <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                                              meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .01, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,minFractNN=.8,
                                                                              clustSizeThresh = 0,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                                              commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=3, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                                              survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                                              CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                                              ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames,
                                                                              GSEAanalysis=FALSE,clinVarPlots=TRUE, fractFeatIntersectThresh=.8,numFeatIntersectThresh =0,clustSizeFractThresh =0,
                                                                              findCommWithWeights=TRUE, plotSimilEdgeWeight =TRUE,plotToScreen=TRUE)











###2014 features
source("/home/ywrfc09/CoINcIDE/coincide/CoINcIDE_packageVersion//CoINcIDE/R/CoINcIDE_geneExprProcess.R")
source("/home/ywrfc09/CoINcIDE/coincide/CoINcIDE_packageVersion//CoINcIDE/R/CoINcIDE_computeEdges.R")
source("/home/ywrfc09/CoINcIDE/coincide/CoINcIDE_packageVersion/CoINcIDE/R/CoINcIDE_metaClusterAnalysis.R")
source("/home/ywrfc09/CoINcIDE/coincide/GenomeBiology/GenomeBiology_metaFeatures_analysis_wrapper.R")   
saveDir <- "/home/ywrfc09/ovarian_analysis/metaRankWithTop20Genes"
globalSaveDir <- "/home/ywrfc09/ovarian_analysis/"
#200 features, pearson:
edgeMethod <- "pearson"
load(paste0(saveDir,"/curatedOvarianData_kmeansConsensus_nstart1_2000Features_2015-04-30.RData.gzip"))
load(paste0(saveDir,"/metaFeatures_2000.RData.gzip"))

load(paste0(globalSaveDir,"/esets_proc_TCGAcombat.RData.gzip"))

#now format just as a list of data matrices.
dataMatrixList <- exprSetListToMatrixList(esets,featureDataFieldName="gene")

names(dataMatrixList) <- names(esets)

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

clustSampleIndexList <-  kmeansConsensus$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensus$clustFeatureIndexList_PACR

clusterCoINcIDE_output =  kmeansConsensus

edgeMethod <- "pearson"
centroidMethod <- "mean"
outputFile <- "~/CoINcIDE_messages.txt"
numSims <- 500
source("/home/ywrfc09/CoINcIDE/coincide/CoINcIDE_packageVersion/CoINcIDE/R/CoINcIDE_computeEdges.R")
ovarian_2014genes_pearson_meanCentroid <-CoINcIDE_getAdjMatrices(dataMatrixList,clustSampleIndexList,clustFeatureIndexList,
                                                                edgeMethod=edgeMethod,centroidMethod=centroidMethod,
                                                                numSims=500,
                                                                outputFile=outputFile)


save(ovarian_2014genes_pearson_meanCentroid,file=paste0(saveDir,"ovarian_2014genes_pearson_meanCentroid.RData.gzip"),compress="gzip")

###null tests
CoINcIDE_nullOutputList <- create_nullMatrixList_results(dataMatrixList=dataMatrixList,numIter=5,
                                                         clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                         edgeMethod=edgeMethod,
                                                         numSims=numSims,
                                                         outputFile=outputFile,
                                                         centroidMethod=centroidMethod)





FDRstats_2014Features <- global_FDR(CoINcIDE_outputList=CoINcIDE_nullOutputList,
                                   minTrueSimilThresh=.4,maxTrueSimilThresh=Inf,
                                   outputFile=outputFile,fractFeatIntersectThresh=.8,numFeatIntersectThresh=0 ,clustSizeThresh=0, clustSizeFractThresh=0,
                                   meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, 
                                   saveDir = saveDir,experimentName = "nullTest240Features",
                                   commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=3,minFractNN =.8,findCommWithWeights=TRUE)


##now analyze results:
load(paste0(globalSaveDir,"/esets_proc_TCGAcombat.RData.gzip"))
esets = esets
#dataMatrixList
load(paste0(saveDir,"ovarian_2014genes_pearson_meanCentroid.RData.gzip"))
CoINcIDE_output = ovarian_2014genes_pearson_meanCentroid
experimentName <- "ovarian_2014_pear_meanCentroid"
eset_featureDataFieldName="gene"
networkColors = "Set2"
outcomesVarBinary="vital_status"
outcomesVarCont = "days_to_death"
ovarian <- TRUE
fisherTestVariables <- c("histological_type","tumorstage","recurrence_status","grade","age_at_initial_pathologic_diagnosis")
fisherTestVariableLegendNames <- c("hist\ntype","tumor\nstage","recurrence","tumor\ngrade","age")
fisherTestVariableTitleNames <- c("histological type","tumor stage", "recurrence status","tumor grade","age at diagnosis")

#breast: must also input dataMatrixList
#brewPal = c("Set3","Paired","Spectral","BrBG","PiYG","RdYlGn","RdYlBu","RdBu","PiYG","Set2"),


#.5 gives clearer clusters:
ovarian_2014genes_pearson_meanCentroid_analysis <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                                              meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .01, minTrueSimilThresh = .5, maxTrueSimilThresh = Inf,minFractNN=.8,
                                                                              clustSizeThresh = 0,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                                              commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=3, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                                              survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                                              CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                                              ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames,
                                                                              GSEAanalysis=FALSE,clinVarPlots=FALSE, fractFeatIntersectThresh=.8,numFeatIntersectThresh =0,clustSizeFractThresh =0,
                                                                              findCommWithWeights=TRUE, plotSimilEdgeWeight = TRUE,plotToScreen=TRUE)









#now do  spearman
edgeMethod <- "spearman"
load("/home/kplaney/ovarian_analysis/curatedOvarianData_kmeansConsensus_nstart1_200Features_2015-04-28.RData.gzip")
load("/home/kplaney/ovarian_analysis/metaFeatures_200.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load("/home/kplaney/ovarian_analysis/esets_proc_TCGAcombat.RData.gzip")


#now format just as a list of data matrices.
dataMatrixList <- exprSetListToMatrixList(esets,featureDataFieldName="gene")

names(dataMatrixList) <- names(esets)

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

clustSampleIndexList <-  kmeansConsensus$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensus$clustFeatureIndexList_PACR



source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


ov_200F_spearman_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                      edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                      sigMethod=sigMethod,maxNullFractSize=maxNullFractSize,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                      outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                      numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                      clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


save(ov_200F_spearman_meanMatrix ,file=
       paste0("/home/kplaney/ovarian_analysis/adjMatrices_200F_spearman_meanMatrix_",Sys.Date(),"RData.gzip"),compress="gzip")


########500:
fractFeatIntersectThresh=.8
#we know these all share at least 35 genes.
numFeatIntersectThresh=425
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=3
#this is "softer" than breast cancer.
maxNullFractSize=.2
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/ovarian_analysis//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"

# 500features, pearson:
edgeMethod <- "pearson"
load("/home/kplaney/ovarian_analysis//curatedOvarianData_kmeansConsensus_nstart1_500Features_2015-04-29.RData.gzip")
load("/home/kplaney/ovarian_analysis/metaFeatures_500.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load("/home/kplaney/ovarian_analysis/esets_proc_TCGAcombat.RData.gzip")


#now format just as a list of data matrices.
dataMatrixList <- exprSetListToMatrixList(esets,featureDataFieldName="gene")

names(dataMatrixList) <- names(esets)

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

clustSampleIndexList <-  kmeansConsensus$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensus$clustFeatureIndexList_PACR

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


ov_500F_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                      edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                      sigMethod=sigMethod,maxNullFractSize=maxNullFractSize,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                      outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                      numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                      clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


#had a bug in paste:
save(ov_500F_pearson_meanMatrix ,file=
       paste0("/home/kplaney/ovarian_analysis/adjMatrices_500F_pearson_meanMatrix_updated",Sys.Date(),"RData.gzip"),compress="gzip")


##############
#now do  spearman
edgeMethod <- "spearman"
load("/home/kplaney/ovarian_analysis//curatedOvarianData_kmeansConsensus_nstart1_500Features_2015-04-29.RData.gzip")
load("/home/kplaney/ovarian_analysis/metaFeatures_500.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load("/home/kplaney/ovarian_analysis/esets_proc_TCGAcombat.RData.gzip")


#now format just as a list of data matrices.
dataMatrixList <- exprSetListToMatrixList(esets,featureDataFieldName="gene")

names(dataMatrixList) <- names(esets)

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

clustSampleIndexList <-  kmeansConsensus$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensus$clustFeatureIndexList_PACR



source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


ov_500F_spearman_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                       edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                       sigMethod=sigMethod,maxNullFractSize=maxNullFractSize,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                       outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                       numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                       clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


save(ov_500F_spearman_meanMatrix ,file=
       paste0("/home/kplaney/ovarian_analysis/adjMatrices_500F_spearman_meanMatrix_",Sys.Date(),"RData.gzip"),compress="gzip")

########1000
fractFeatIntersectThresh=.8
#we know these all share at least 35 genes.
numFeatIntersectThresh=800
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=3
#this is "softer" than breast cancer.
maxNullFractSize=.2
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/ovarian_analysis//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"

# 500features, pearson:
edgeMethod <- "pearson"
load("/home/kplaney/ovarian_analysis//curatedOvarianData_kmeansConsensus_nstart1_1000Features_2015-04-29.RData.gzip")
load("/home/kplaney/ovarian_analysis/metaFeatures_1000.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load("/home/kplaney/ovarian_analysis/esets_proc_TCGAcombat.RData.gzip")


#now format just as a list of data matrices.
dataMatrixList <- exprSetListToMatrixList(esets,featureDataFieldName="gene")

names(dataMatrixList) <- names(esets)

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

clustSampleIndexList <-  kmeansConsensus$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensus$clustFeatureIndexList_PACR

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


ov_1000F_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                      edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                      sigMethod=sigMethod,maxNullFractSize=maxNullFractSize,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                      outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                      numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                      clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


#had a bug in paste:
save(ov_1000F_pearson_meanMatrix ,file=
       paste0("/home/kplaney/ovarian_analysis/adjMatrices_1000F_pearson_meanMatrix_updated",Sys.Date(),"RData.gzip"),compress="gzip")



###spearman
########1000
fractFeatIntersectThresh=.8
#we know these all share at least 35 genes.
numFeatIntersectThresh=800
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=3
#this is "softer" than breast cancer.
maxNullFractSize=.2
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/ovarian_analysis//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"

edgeMethod <- "spearman"
load("/home/kplaney/ovarian_analysis_withTop20Genes///curatedOvarianData_kmeansConsensus_nstart1_1000Features_2015-04-29.RData.gzip")
load("/home/kplaney/ovarian_analysis_withTop20Genes//metaFeatures_1000.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load("/home/kplaney/ovarian_analysis/esets_proc_TCGAcombat.RData.gzip")


#now format just as a list of data matrices.
dataMatrixList <- exprSetListToMatrixList(esets,featureDataFieldName="gene")

names(dataMatrixList) <- names(esets)

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

clustSampleIndexList <-  kmeansConsensus$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensus$clustFeatureIndexList_PACR

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


ov_1000F_spearman_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                       edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                       sigMethod=sigMethod,maxNullFractSize=maxNullFractSize,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                       outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                       numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                       clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


#had a bug in paste:
save(ov_1000F_spearman_meanMatrix ,file=
       paste0("/home/kplaney/ovarian_analysis_withTop20Genes//adjMatrices_1000F_pearson_meanMatrix_updated",Sys.Date(),".RData.gzip"),compress="gzip")


####2000:
fractFeatIntersectThresh=.85
#we know these all share at least 35 genes.
numFeatIntersectThresh=1700
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=3
#this is "softer" than breast cancer.
maxNullFractSize=.2
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/ovarian_analysis//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"

# 1000features, pearson:
edgeMethod <- "pearson"
load("/home/kplaney/ovarian_analysis//curatedOvarianData_kmeansConsensus_nstart1_2000Features_2015-04-30.RData.gzip")
load("/home/kplaney/ovarian_analysis/metaFeatures_2000.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load("/home/kplaney/ovarian_analysis/esets_proc_TCGAcombat.RData.gzip")


#now format just as a list of data matrices.
dataMatrixList <- exprSetListToMatrixList(esets,featureDataFieldName="gene")

names(dataMatrixList) <- names(esets)

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

clustSampleIndexList <-  kmeansConsensus$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensus$clustFeatureIndexList_PACR

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


ov_2000F_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                       edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                       sigMethod=sigMethod,maxNullFractSize=maxNullFractSize,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                       outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                       numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                       clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


#had a bug in paste:
save(ov_2000F_pearson_meanMatrix ,file=
       paste0("/home/kplaney/ovarian_analysis/adjMatrices_2000F_pearson_meanMatrix_updated",Sys.Date(),"RData.gzip"),compress="gzip")

##spearman
fractFeatIntersectThresh=.85
#we know these all share at least 35 genes.
numFeatIntersectThresh=1700
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=3
#this is "softer" than breast cancer.
maxNullFractSize=.2
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/ovarian_analysis//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"

# 1000features, pearson:
edgeMethod <- "spearman"
load("/home/kplaney/ovarian_analysis_withTop20Genes///curatedOvarianData_kmeansConsensus_nstart1_2000Features_2015-04-30.RData.gzip")
load("/home/kplaney/ovarian_analysis_withTop20Genes//metaFeatures_2000.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load("/home/kplaney/ovarian_analysis_withTop20Genes//esets_proc_TCGAcombat.RData.gzip")


#now format just as a list of data matrices.
dataMatrixList <- exprSetListToMatrixList(esets,featureDataFieldName="gene")

names(dataMatrixList) <- names(esets)

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

clustSampleIndexList <-  kmeansConsensus$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensus$clustFeatureIndexList_PACR

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


ov_2000F_spearman_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                       edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                       sigMethod=sigMethod,maxNullFractSize=maxNullFractSize,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                       outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                       numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                       clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


#had a bug in paste:
save(ov_2000F_spearman_meanMatrix ,file=
       paste0("/home/kplaney/ovarian_analysis_withTop20Genes//adjMatrices_2000F_spearman_meanMatrix",Sys.Date(),".RData.gzip"),compress="gzip")
