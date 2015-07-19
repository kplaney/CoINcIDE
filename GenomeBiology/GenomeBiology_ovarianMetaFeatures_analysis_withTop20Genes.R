#!/usr/bin/Rscript --default-packages=utils

library("CoINcIDE")
source("/home/ywrfc09/CoINcIDE/coincide/oldCode/CoINcIDE_metaFeatures_analysis_wrapper.R")   
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
                                     experimentName=experimentName,savePlot=TRUE)




ovarian_240genes_pearson_meanCentroid_analysis <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                                              meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .01, minTrueSimilThresh = .5, maxTrueSimilThresh = Inf,minFractNN=.8,
                                                                              clustSizeThresh = 0,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                                              commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=3, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                                              survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                                              CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                                              ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames,
                                                                              GSEAanalysis=TRUE,clinVarPlots=TRUE, fractFeatIntersectThresh=.6,numFeatIntersectThresh =0,clustSizeFractThresh =0,
                                                                              findCommWithWeights=TRUE, plotSimilEdgeWeight = TRUE,plotToScreen=FALSE,fractEdgesInVsOutEdge=0, fractEdgesInVsOutComm=0,minNumEdgesForCluster=1)







#2000
library("CoINcIDE")
source("/home/ywrfc09/CoINcIDE/coincide/oldCode/CoINcIDE_metaFeatures_analysis_wrapper.R")   
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
                                     experimentName=experimentName,savePlot=TRUE)




ovarian_2014genes_pearson_meanCentroid_analysis <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                                              meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .01, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,minFractNN=.8,
                                                                              clustSizeThresh = 0,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                                              commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=3, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                                              survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                                              CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                                              ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames,
                                                                              GSEAanalysis=TRUE,clinVarPlots=TRUE, fractFeatIntersectThresh=.6,numFeatIntersectThresh =0,clustSizeFractThresh =0,
                                                                              findCommWithWeights=TRUE, plotSimilEdgeWeight = TRUE,plotToScreen=FALSE,fractEdgesInVsOutEdge=0, fractEdgesInVsOutComm=0,minNumEdgesForCluster=1)












############
##500plus features

source("/home/ywrfc09/CoINcIDE/coincide/CoINcIDE_packageVersion//CoINcIDE/R/CoINcIDE_geneExprProcess.R")
source("/home/ywrfc09/CoINcIDE/coincide/CoINcIDE_packageVersion//CoINcIDE/R/CoINcIDE_computeEdges.R")
source("/home/ywrfc09/CoINcIDE/coincide/CoINcIDE_packageVersion/CoINcIDE/R/CoINcIDE_metaClusterAnalysis.R")
source("/home/ywrfc09/CoINcIDE/coincide/oldCode/CoINcIDE_metaFeatures_analysis_wrapper.R")   
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
source("/home/ywrfc09/CoINcIDE/coincide/oldCode/CoINcIDE_metaFeatures_analysis_wrapper.R")   
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
source("/home/ywrfc09/CoINcIDE/coincide/oldCode/CoINcIDE_metaFeatures_analysis_wrapper.R")   
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
