#!/usr/bin/Rscript --default-packages=utils

#200 features, pearson:
library("CoINcIDE")
saveDir <- "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/"
globalSaveDir <-  "/home/ywrfc09/breast_analysis"
load(paste0(saveDir,"/metaFeatures_200.RData.gzip"))

clusterCoINcIDE_output <- readRDS("/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/curatedbreastData_kmeansConsensus_nstart1pItem9200Features_2015-05-04.rds")
clustSampleIndexList <-  clusterCoINcIDE_output$clustSampleIndexList_PACR
clustFeatureIndexList <- clusterCoINcIDE_output$clustFeatureIndexList_PACR


outputFile <- "~/CoINcIDE_messages.txt"

 source("/home/ywrfc09/CoINcIDE/coincide/GenomeBiology/GenomeBiology_metaFeatures_analysis_wrapper.R")
CoINcIDE_output <- readRDS("/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/CoINcIDE_results_breast278F_pearson_edgeMethod_mean_centroidMethod2015-07-08.rds")

saveDirNew <- "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/"


load(paste0(globalSaveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies

experimentName <- "breast_278genes_pear_meanCent_MM3"
eset_featureDataFieldName="gene_symbol"
networkColors = "Set3"
outcomesVarBinary=NA
#not enough continuous variables for breast data
outcomesVarCont = NA
ovarian <- FALSE
eset_uniquePatientID="dbUniquePatientID"
fisherTestVariables <- c("DFS","RFS","metastasis","pCR",
                         "tumor_stage_preTrt" ,
                         "preTrt_lymph_node_status" , 
                         "neoadjuvant_or_adjuvant","ER_preTrt","HER2_preTrt","hist_grade","chemotherapyClass",
                         "anti_HER2","anti_estrogen")

fisherTestVariableLegendNames <- c("DFS","RFS","metastasis","pCR","pretreat\ntumor stage",
                                   "pretreat\nlymph node","neoadjuvant\nvs adjuvant","pretreat\nER status",
                                   "pretreat\nHER2 status","hist grade","chemotherapy","anti-HER2\ntreat",
                                   "anti-ER\ntreat")


fisherTestVariableTitleNames <- c("DFS","RFS","metastasis","pCR","pretreatment tumor stage",
                                  "pretreatment lymph node","neoadjuvant vs adjuvant","pretreatment ER status",
                                  "pretreatment HER2 status","histological grade","chemotherapy","anti-HER2 treatment",
                                  "anti-ER treatment")

 source("/home/ywrfc09/gitRepos/CoINcIDE/coincide/GenomeBiology/GenomeBiology_metaFeatures_analysis_wrapper.R")

#looks like simil thresh of .3 or .5 appropriate
densityPlot <- meanMetricDensityPlot(CoINcIDE_output$meanMetricMatrix,saveDir=saveDir,
                                     experimentName=experimentName,savePlot=TRUE)


#match up Y axis limit for both 278 and 2052
breast_278genes_pearson_meanCentroid_analysis <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                                    meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .01, minTrueSimilThresh = .3, maxTrueSimilThresh = Inf,minFractNN=.8,
                                                                    clustSizeThresh = 0,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                                    commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=3, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                                    survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                                    CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                                    ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames,
                                                                    GSEAanalysis=TRUE,clinVarPlots=TRUE, fractFeatIntersectThresh=.6,numFeatIntersectThresh =0,clustSizeFractThresh =0,
                                                                    findCommWithWeights=TRUE, plotSimilEdgeWeight = TRUE,plotToScreen=FALSE,fractEdgesInVsOutEdge=0, fractEdgesInVsOutComm=0,minNumEdgesForCluster=1,plotStackedYLimit=800)




CoINcIDE_nullOutput <- readRDS("/home/ywrfc09/breast_analysis/metaRankWithTop20Genes///CoINcIDE_NullOutput_breast278F_pearsonedgeMethod_mean_centroidMethod2015-07-09.rds")
globalFDR_results <- globalFDR(CoINcIDE_outputList=CoINcIDE_nullOutput$CoINcIDE_NullOutputList,
                               edgeMethod="pearson",minTrueSimilThresh=0.3,maxTrueSimilThresh=Inf,
                               outputFile=outputFile,fractFeatIntersectThresh=0,numFeatIntersectThresh=0 ,clustSizeThresh=0, clustSizeFractThresh=0,
                               meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .01, 
                               saveDir = saveDir,experimentName = "nullTest",
                               commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=3,minFractNN =.8,
                               findCommWithWeights=TRUE,minNumEdgesForCluster=1,fractEdgesInVsOutComm=0,fractEdgesInVsOutEdge=0)

saveRDS(globalFDR_results,file=paste0(saveDir,"/CoINcIDE_globalFDRresults_",experimentName,"_",edgeMethod,"edgeMethod_",centroidMethod,"_centroidMethod_minTrueSimil",minTrueSimilThresh,"_",Sys.Date(),".rds"),compress=TRUE)


options(bitmapType="cairo")
png(filename=paste0(saveDir,"/",experimentName,"_","dendogram_",Sys.Date(),".png"),width=1800,height=1000,res=160)

plot(breast_278genes_pearson_meanCentroid_analysis$commInfo$unPrunedDendogram)

dev.off()
finalNodeMatrix <-  breast_278genes_pearson_meanCentroid_analysis$commInfo$attrDF
origEdgeMatrix <- breast_278genes_pearson_meanCentroid_analysis$finalEdgeInfo$filterEdgeOutput$edgeMatrix
origEdgeWeightsMatrix <- breast_278genes_pearson_meanCentroid_analysis$finalEdgeInfo$filterEdgeOutput$edgeWeightMatrix
#right now these are igraph ids: will need to change that.
finalEdgeMatrix <- breast_278genes_pearson_meanCentroid_analysis$commInfo$edgeDF[,c(1:2)]
clustIndexMatrix <- breast_278genes_pearson_meanCentroid_analysis$CoINcIDE_computeEdgesObject$clustIndexMatrix


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
saveRDS(leaveXOutResults,file=paste0(saveDir,"/CoINcIDE_LeaveXOutAnalysis_",experimentName,"_",Sys.Date(),".rds"),compress=TRUE)

LO_analysis <- plotLeaveXOutAnalysis(leaveXOutResults,
                                     experimentName=experimentName,saveDir=saveDir)


##2000 genes:
library("CoINcIDE")
saveDir <- "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/"
globalSaveDir <-  "/home/ywrfc09/breast_analysis"
load(paste0(saveDir,"/metaFeatures_2000.RData.gzip"))

clusterCoINcIDE_output <- readRDS("/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/curatedbreastData_kmeansConsensus_nstart1pItem92000Features_2015-05-05.rds")
clustSampleIndexList <-  clusterCoINcIDE_output$clustSampleIndexList_PACR
clustFeatureIndexList <- clusterCoINcIDE_output$clustFeatureIndexList_PACR


outputFile <- "~/CoINcIDE_messages.txt"

 source("/home/ywrfc09/CoINcIDE/coincide/GenomeBiology/GenomeBiology_metaFeatures_analysis_wrapper.R")
CoINcIDE_output <- readRDS("/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/CoINcIDE_results_breast2052F_pearson_edgeMethod_mean_centroidMethod2015-07-08.rds")

saveDirNew <- "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/"


load(paste0(globalSaveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies

experimentName <- "breast_2052genes_pear_meanCent_MM4"
eset_featureDataFieldName="gene_symbol"
networkColors = "Set3"
outcomesVarBinary=NA
#not enough continuous variables for breast data
outcomesVarCont = NA
ovarian <- FALSE
eset_uniquePatientID="dbUniquePatientID"
fisherTestVariables <- c("DFS","RFS","metastasis","pCR",
                         "tumor_stage_preTrt" ,
                         "preTrt_lymph_node_status" , 
                         "neoadjuvant_or_adjuvant","ER_preTrt","HER2_preTrt","hist_grade","chemotherapyClass",
                         "anti_HER2","anti_estrogen")

fisherTestVariableLegendNames <- c("DFS","RFS","metastasis","pCR","pretreat\ntumor stage",
                                   "pretreat\nlymph node","neoadjuvant\nvs adjuvant","pretreat\nER status",
                                   "pretreat\nHER2 status","hist grade","chemotherapy","anti-HER2\ntreat",
                                   "anti-ER\ntreat")


fisherTestVariableTitleNames <- c("DFS","RFS","metastasis","pCR","pretreatment tumor stage",
                                  "pretreatment lymph node","neoadjuvant vs adjuvant","pretreatment ER status",
                                  "pretreatment HER2 status","histological grade","chemotherapy","anti-HER2 treatment",
                                  "anti-ER treatment")

 source("/home/ywrfc09/CoINcIDE/coincide/GenomeBiology/GenomeBiology_metaFeatures_analysis_wrapper.R")

#looks like simil thresh of .4 appropriate
densityPlot <- meanMetricDensityPlot(CoINcIDE_output$meanMetricMatrix,saveDir=saveDir,
                                     experimentName=experimentName,savePlot=TRUE)




breast_2052genes_pearson_meanCentroid_analysis <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                                                                                               meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .01, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,minFractNN=.8,
                                                                                                                               clustSizeThresh = 0,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                                                                                               commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=3, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                                                                                               survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                                                                                               CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                                                                                               ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames,
                                                                                                                               GSEAanalysis=TRUE,clinVarPlots=TRUE, fractFeatIntersectThresh=.6,numFeatIntersectThresh =0,clustSizeFractThresh =0,
                                                                                                                               findCommWithWeights=TRUE, plotSimilEdgeWeight = TRUE,plotToScreen=FALSE,fractEdgesInVsOutEdge=0, fractEdgesInVsOutComm=0,minNumEdgesForCluster=1,plotStackedYLimit=800)








CoINcIDE_nullOutput <- readRDS("/home/ywrfc09/breast_analysis/metaRankWithTop20Genes///CoINcIDE_NullOutput_breast2052F_pearsonedgeMethod_mean_centroidMethod2015-07-09.rds")
globalFDR_results <- globalFDR(CoINcIDE_outputList=CoINcIDE_nullOutput$CoINcIDE_NullOutputList,
                               edgeMethod="pearson",minTrueSimilThresh=0.4,maxTrueSimilThresh=Inf,
                               outputFile=outputFile,fractFeatIntersectThresh=0,numFeatIntersectThresh=0 ,clustSizeThresh=0, clustSizeFractThresh=0,
                               meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .01, 
                               saveDir = saveDir,experimentName = "nullTest",
                               commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=3,minFractNN =.8,
                               findCommWithWeights=TRUE,minNumEdgesForCluster=1,fractEdgesInVsOutComm=0,fractEdgesInVsOutEdge=0)

saveRDS(globalFDR_results,file=paste0(saveDir,"/CoINcIDE_globalFDRresults_",experimentName,"_",edgeMethod,"edgeMethod_",centroidMethod,"_centroidMethod_minTrueSimil",minTrueSimilThresh,"_",Sys.Date(),".rds"),compress=TRUE)


options(bitmapType="cairo")
png(filename=paste0(saveDir,"/",experimentName,"_","dendogram_",Sys.Date(),".png"),width=1800,height=1000,res=160)

plot(breast_2052genes_pearson_meanCentroid_analysis$commInfo$unPrunedDendogram)

dev.off()
finalNodeMatrix <-  breast_2052genes_pearson_meanCentroid_analysis$commInfo$attrDF
origEdgeMatrix <- breast_2052genes_pearson_meanCentroid_analysis$finalEdgeInfo$filterEdgeOutput$edgeMatrix
origEdgeWeightsMatrix <- breast_2052genes_pearson_meanCentroid_analysis$finalEdgeInfo$filterEdgeOutput$edgeWeightMatrix
#right now these are igraph ids: will need to change that.
finalEdgeMatrix <- breast_2052genes_pearson_meanCentroid_analysis$commInfo$edgeDF[,c(1:2)]
clustIndexMatrix <- breast_2052genes_pearson_meanCentroid_analysis$CoINcIDE_computeEdgesObject$clustIndexMatrix


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
saveRDS(leaveXOutResults,file=paste0(saveDir,"/CoINcIDE_LeaveXOutAnalysis_",experimentName,"_",Sys.Date(),".rds"),compress=TRUE)

LO_analysis <- plotLeaveXOutAnalysis(leaveXOutResults,
                                     experimentName=experimentName,saveDir=saveDir)





#####50 but NO PAM50 (no top 20 genes either)
numFeatures <- 50
saveDirPAM50 <- "/home/ywrfc09/breast_analysis/PAM50_analyses"
load(paste0(saveDirPAM50,"/pam50_centroids_updatedSymbols.RData"))
PAM50genes <- centroidMatrix[,"geneSymbol"]
##load(paste0(saveDirPAM50,"/pam50_centroids_updatedSymbols.RData"))
library("CoINcIDE")
saveDir <- "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/"
globalSaveDir <-  "/home/ywrfc09/breast_analysis"
load(paste0(globalSaveDir,"/metaFeatures_50_NOPAM50_NOTOP20GENES.RData.gzip")


noPAM50_50MetaF <- setdiff(metaFeatures$finalFeatures,PAM50genes)

metaFeatures$finalFeatures <- noPAM50_50MetaF

#careful - sometimes saved as RDATA, sometimes rds.
clusterCoINcIDE_output  <-readRDS(paste0("/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/curatedbreastData_kmeansConsensus_nstart1pItem950Features_NOPAM50_50genesNOTOP20_2015-12-17.rds"))
clustSampleIndexList <-  clusterCoINcIDE_output$clustSampleIndexList_PACR
clustFeatureIndexList <- clusterCoINcIDE_output$clustFeatureIndexList_PACR


outputFile <- "~/CoINcIDE_messages.txt"

 source("/home/ywrfc09/CoINcIDE/coincide/GenomeBiology/GenomeBiology_metaFeatures_analysis_wrapper.R")
CoINcIDE_output <- readRDS("/home/ywrfc09/breast_analysis/metaRankWithTop20Genes//CoINcIDE_results_50F_NOPAM50_pearson_edgeMethod_mean_centroidMethod2015-12-17.rds")

saveDirNew <- "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/"


load(paste0(globalSaveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies

eset_featureDataFieldName="gene_symbol"
networkColors = "Set3"
outcomesVarBinary=NA
#not enough continuous variables for breast data
outcomesVarCont = NA
ovarian <- FALSE
eset_uniquePatientID="dbUniquePatientID"
fisherTestVariables <- c("DFS","RFS","metastasis","pCR",
                         "tumor_stage_preTrt" ,
                         "preTrt_lymph_node_status" , 
                         "neoadjuvant_or_adjuvant","ER_preTrt","HER2_preTrt","hist_grade","chemotherapyClass",
                         "anti_HER2","anti_estrogen")

fisherTestVariableLegendNames <- c("DFS","RFS","metastasis","pCR","pretreat\ntumor stage",
                                   "pretreat\nlymph node","neoadjuvant\nvs adjuvant","pretreat\nER status",
                                   "pretreat\nHER2 status","hist grade","chemotherapy","anti-HER2\ntreat",
                                   "anti-ER\ntreat")


fisherTestVariableTitleNames <- c("DFS","RFS","metastasis","pCR","pretreatment tumor stage",
                                  "pretreatment lymph node","neoadjuvant vs adjuvant","pretreatment ER status",
                                  "pretreatment HER2 status","histological grade","chemotherapy","anti-HER2 treatment",
                                  "anti-ER treatment")

experimentName <- paste0("breast_50genesNoPAM50NoTop20_pear_meanCent_",Sys.Date())
densityPlot <- meanMetricDensityPlot(CoINcIDE_output$meanMetricMatrix,saveDir=saveDir,
                                     experimentName=experimentName,savePlot=TRUE)

#looks like simil thresh of just under .4 appropriate so round up
experimentName <- paste0("breast_50genesNoPAM50NoTop20_pear_meanCent_MM4",Sys.Date())
minTrueSimil <- 0.4
 source("/home/ywrfc09/CoINcIDE/coincide/GenomeBiology/GenomeBiology_metaFeatures_analysis_wrapper.R")


breast_50genes_NoPAM50_pearson_meanCentroid_analysis <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                                                                                               meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .01, minTrueSimilThresh = minTrueSimil, maxTrueSimilThresh = Inf,minFractNN=.8,
                                                                                                                               clustSizeThresh = 0,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                                                                                               commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=3, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                                                                                               survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                                                                                               CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                                                                                               ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames,
                                                                                                                               GSEAanalysis=TRUE,clinVarPlots=TRUE, fractFeatIntersectThresh=.6,numFeatIntersectThresh =0,clustSizeFractThresh =0,
                                                                                                                               findCommWithWeights=TRUE, plotSimilEdgeWeight = TRUE,plotToScreen=FALSE,fractEdgesInVsOutEdge=0, fractEdgesInVsOutComm=0,minNumEdgesForCluster=1,plotStackedYLimit=800)


###200 but NO PAM50
numFeatures <- 200
saveDirPAM50 <- "/home/ywrfc09/breast_analysis/PAM50_analyses"
load(paste0(saveDirPAM50,"/pam50_centroids_updatedSymbols.RData"))
##load(paste0(saveDirPAM50,"/pam50_centroids_updatedSymbols.RData"))
library("CoINcIDE")
saveDir <- "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/"
globalSaveDir <-  "/home/ywrfc09/breast_analysis"
load(paste0(saveDir,"/metaFeatures_",numFeatures,".RData.gzip"))

PAM50genes <- centroidMatrix[,"geneSymbol"]

noPAM50_200MetaF <- setdiff(metaFeatures$finalFeatures,PAM50genes)

metaFeatures$finalFeatures <- noPAM50_200MetaF

#careful - sometimes saved as RDATA, sometimes rds.
clusterCoINcIDE_output  <-readRDS(paste0("/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/curatedbreastData_kmeansConsensus_nstart1pItem9",numFeatures,"Features_NOPAM50_2015-12-17.rds"))
clustSampleIndexList <-  clusterCoINcIDE_output$clustSampleIndexList_PACR
clustFeatureIndexList <- clusterCoINcIDE_output$clustFeatureIndexList_PACR

outputFile <- "~/CoINcIDE_messages.txt"

CoINcIDE_output <- readRDS(paste0("/home/ywrfc09/breast_analysis/metaRankWithTop20Genes//CoINcIDE_results_",numFeatures,"F_NOPAM50_pearson_edgeMethod_mean_centroidMethod2015-12-17.rds"))

saveDirNew <- "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/"


load(paste0(globalSaveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies

eset_featureDataFieldName="gene_symbol"
networkColors = "Set3"
outcomesVarBinary=NA
#not enough continuous variables for breast data
outcomesVarCont = NA
ovarian <- FALSE
eset_uniquePatientID="dbUniquePatientID"
fisherTestVariables <- c("DFS","RFS","metastasis","pCR",
                         "tumor_stage_preTrt" ,
                         "preTrt_lymph_node_status" , 
                         "neoadjuvant_or_adjuvant","ER_preTrt","HER2_preTrt","hist_grade","chemotherapyClass",
                         "anti_HER2","anti_estrogen")

fisherTestVariableLegendNames <- c("DFS","RFS","metastasis","pCR","pretreat\ntumor stage",
                                   "pretreat\nlymph node","neoadjuvant\nvs adjuvant","pretreat\nER status",
                                   "pretreat\nHER2 status","hist grade","chemotherapy","anti-HER2\ntreat",
                                   "anti-ER\ntreat")


fisherTestVariableTitleNames <- c("DFS","RFS","metastasis","pCR","pretreatment tumor stage",
                                  "pretreatment lymph node","neoadjuvant vs adjuvant","pretreatment ER status",
                                  "pretreatment HER2 status","histological grade","chemotherapy","anti-HER2 treatment",
                                  "anti-ER treatment")

 source("/home/ywrfc09/CoINcIDE/coincide/GenomeBiology/GenomeBiology_metaFeatures_analysis_wrapper.R")

experimentName <- paste0("breast_",numFeatures,"FNoPAM50_pear_meanCent_",Sys.Date())
densityPlot <- meanMetricDensityPlot(CoINcIDE_output$meanMetricMatrix,saveDir=saveDir,
                                     experimentName=experimentName,savePlot=TRUE)

#looks like simil thresh of around 0.2, or 0.5 - tried 0.2
experimentName <- paste0("breast_",numFeatures",FNoPAM50_pear_meanCent_MM2")
minTrueSimil <- 0.2
breast_200genes_NoPAM50_pearson_meanCentroid_analysis <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                                                                                               meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .01, minTrueSimilThresh = minTrueSimil, maxTrueSimilThresh = Inf,minFractNN=.8,
                                                                                                                               clustSizeThresh = 0,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                                                                                               commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=3, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                                                                                               survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                                                                                               CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                                                                                               ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames,
                                                                                                                               GSEAanalysis=TRUE,clinVarPlots=TRUE, fractFeatIntersectThresh=.6,numFeatIntersectThresh =0,clustSizeFractThresh =0,
                                                                                                                               findCommWithWeights=TRUE, plotSimilEdgeWeight = TRUE,plotToScreen=FALSE,fractEdgesInVsOutEdge=0, fractEdgesInVsOutComm=0,minNumEdgesForCluster=1,plotStackedYLimit=800)

##also try 0.5
#looks like simil thresh of around 0.2, or 0.5 - tried 0.2
minTrueSimil <- 0.5
experimentName <- paste0("breast_",numFeatures,"FNoPAM50_pear_meanCent_MM",minTrueSimil)
breast_200genes_NoPAM50_pearson_meanCentroid_analysis <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                                                                                               meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .01, minTrueSimilThresh = minTrueSimil, maxTrueSimilThresh = Inf,minFractNN=.8,
                                                                                                                               clustSizeThresh = 0,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                                                                                               commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=3, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                                                                                               survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                                                                                               CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                                                                                               ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames,
                                                                                                                               GSEAanalysis=TRUE,clinVarPlots=TRUE, fractFeatIntersectThresh=.6,numFeatIntersectThresh =0,clustSizeFractThresh =0,
                                                                                                                               findCommWithWeights=TRUE, plotSimilEdgeWeight = TRUE,plotToScreen=FALSE,fractEdgesInVsOutEdge=0, fractEdgesInVsOutComm=0,minNumEdgesForCluster=1,plotStackedYLimit=800)


###500 features, No PAM50
numFeatures <- 500
saveDirPAM50 <- "/home/ywrfc09/breast_analysis/PAM50_analyses"
load(paste0(saveDirPAM50,"/pam50_centroids_updatedSymbols.RData"))
##load(paste0(saveDirPAM50,"/pam50_centroids_updatedSymbols.RData"))
library("CoINcIDE")
saveDir <- "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/"
globalSaveDir <-  "/home/ywrfc09/breast_analysis"
load(paste0(saveDir,"/metaFeatures_",numFeatures,".RData.gzip"))

PAM50genes <- centroidMatrix[,"geneSymbol"]

noPAM50_500MetaF <- setdiff(metaFeatures$finalFeatures,PAM50genes)

metaFeatures$finalFeatures <- noPAM50_500MetaF

#careful - sometimes saved as RDATA, sometimes rds.
clusterCoINcIDE_output  <-readRDS(paste0("/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/curatedbreastData_kmeansConsensus_nstart1pItem9",numFeatures,"Features_NOPAM50_2015-12-17.rds"))
clustSampleIndexList <-  clusterCoINcIDE_output$clustSampleIndexList_PACR
clustFeatureIndexList <- clusterCoINcIDE_output$clustFeatureIndexList_PACR

outputFile <- "~/CoINcIDE_messages.txt"

CoINcIDE_output <- readRDS(paste0("/home/ywrfc09/breast_analysis/metaRankWithTop20Genes//CoINcIDE_results_",numFeatures,"F_NOPAM50_pearson_edgeMethod_mean_centroidMethod2015-12-17.rds"))

saveDirNew <- "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/"


load(paste0(globalSaveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies

eset_featureDataFieldName="gene_symbol"
networkColors = "Set3"
outcomesVarBinary=NA
#not enough continuous variables for breast data
outcomesVarCont = NA
ovarian <- FALSE
eset_uniquePatientID="dbUniquePatientID"
fisherTestVariables <- c("DFS","RFS","metastasis","pCR",
                         "tumor_stage_preTrt" ,
                         "preTrt_lymph_node_status" , 
                         "neoadjuvant_or_adjuvant","ER_preTrt","HER2_preTrt","hist_grade","chemotherapyClass",
                         "anti_HER2","anti_estrogen")

fisherTestVariableLegendNames <- c("DFS","RFS","metastasis","pCR","pretreat\ntumor stage",
                                   "pretreat\nlymph node","neoadjuvant\nvs adjuvant","pretreat\nER status",
                                   "pretreat\nHER2 status","hist grade","chemotherapy","anti-HER2\ntreat",
                                   "anti-ER\ntreat")


fisherTestVariableTitleNames <- c("DFS","RFS","metastasis","pCR","pretreatment tumor stage",
                                  "pretreatment lymph node","neoadjuvant vs adjuvant","pretreatment ER status",
                                  "pretreatment HER2 status","histological grade","chemotherapy","anti-HER2 treatment",
                                  "anti-ER treatment")

 source("/home/ywrfc09/CoINcIDE/coincide/GenomeBiology/GenomeBiology_metaFeatures_analysis_wrapper.R")

experimentName <- paste0("breast_",numFeatures,"FNoPAM50_pear_meanCent_",Sys.Date())
densityPlot <- meanMetricDensityPlot(CoINcIDE_output$meanMetricMatrix,saveDir=saveDir,
                                     experimentName=experimentName,savePlot=TRUE)

#looks like simil thresh of around 0.5 best.
minTrueSimil <- 0.5
experimentName <- paste0("breast_",numFeatures,"FNoPAM50_pear_meanCent_MM",minTrueSimil)
breast_500genes_NoPAM50_pearson_meanCentroid_analysis <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                                                                                               meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .01, minTrueSimilThresh = minTrueSimil, maxTrueSimilThresh = Inf,minFractNN=.8,
                                                                                                                               clustSizeThresh = 0,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                                                                                               commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=3, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                                                                                               survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                                                                                               CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                                                                                               ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames,
                                                                                                                               GSEAanalysis=TRUE,clinVarPlots=TRUE, fractFeatIntersectThresh=.6,numFeatIntersectThresh =0,clustSizeFractThresh =0,
                                                                                                                               findCommWithWeights=TRUE, plotSimilEdgeWeight = TRUE,plotToScreen=FALSE,fractEdgesInVsOutEdge=0, fractEdgesInVsOutComm=0,minNumEdgesForCluster=1,plotStackedYLimit=800)


##1000 features, no PAM50
numFeatures <- 1000
saveDirPAM50 <- "/home/ywrfc09/breast_analysis/PAM50_analyses"
load(paste0(saveDirPAM50,"/pam50_centroids_updatedSymbols.RData"))
##load(paste0(saveDirPAM50,"/pam50_centroids_updatedSymbols.RData"))
library("CoINcIDE")
saveDir <- "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/"
globalSaveDir <-  "/home/ywrfc09/breast_analysis"
load(paste0(saveDir,"/metaFeatures_",numFeatures,".RData.gzip"))

PAM50genes <- centroidMatrix[,"geneSymbol"]

noPAM50_1000MetaF <- setdiff(metaFeatures$finalFeatures,PAM50genes)

metaFeatures$finalFeatures <- noPAM50_1000MetaF

#careful - sometimes saved as RDATA, sometimes rds.
clusterCoINcIDE_output  <-readRDS(paste0("/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/curatedbreastData_kmeansConsensus_nstart1pItem9",numFeatures,"Features_NOPAM50_2015-12-17.rds"))
clustSampleIndexList <-  clusterCoINcIDE_output$clustSampleIndexList_PACR
clustFeatureIndexList <- clusterCoINcIDE_output$clustFeatureIndexList_PACR

outputFile <- "~/CoINcIDE_messages.txt"

CoINcIDE_output <- readRDS(paste0("/home/ywrfc09/breast_analysis/metaRankWithTop20Genes//CoINcIDE_results_",numFeatures,"F_NOPAM50_pearson_edgeMethod_mean_centroidMethod2015-12-17.rds"))

saveDirNew <- "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/"


load(paste0(globalSaveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies

eset_featureDataFieldName="gene_symbol"
networkColors = "Set3"
outcomesVarBinary=NA
#not enough continuous variables for breast data
outcomesVarCont = NA
ovarian <- FALSE
eset_uniquePatientID="dbUniquePatientID"
fisherTestVariables <- c("DFS","RFS","metastasis","pCR",
                         "tumor_stage_preTrt" ,
                         "preTrt_lymph_node_status" , 
                         "neoadjuvant_or_adjuvant","ER_preTrt","HER2_preTrt","hist_grade","chemotherapyClass",
                         "anti_HER2","anti_estrogen")

fisherTestVariableLegendNames <- c("DFS","RFS","metastasis","pCR","pretreat\ntumor stage",
                                   "pretreat\nlymph node","neoadjuvant\nvs adjuvant","pretreat\nER status",
                                   "pretreat\nHER2 status","hist grade","chemotherapy","anti-HER2\ntreat",
                                   "anti-ER\ntreat")


fisherTestVariableTitleNames <- c("DFS","RFS","metastasis","pCR","pretreatment tumor stage",
                                  "pretreatment lymph node","neoadjuvant vs adjuvant","pretreatment ER status",
                                  "pretreatment HER2 status","histological grade","chemotherapy","anti-HER2 treatment",
                                  "anti-ER treatment")

 source("/home/ywrfc09/CoINcIDE/coincide/GenomeBiology/GenomeBiology_metaFeatures_analysis_wrapper.R")

experimentName <- paste0("breast_",numFeatures,"FNoPAM50_pear_meanCent_",Sys.Date())
densityPlot <- meanMetricDensityPlot(CoINcIDE_output$meanMetricMatrix,saveDir=saveDir,
                                     experimentName=experimentName,savePlot=TRUE)

#looks like simil thresh of around 0.5 best.
minTrueSimil <- 0.5
experimentName <- paste0("breast_",numFeatures,"FNoPAM50_pear_meanCent_MM",minTrueSimil)
breast_1000genes_NoPAM50_pearson_meanCentroid_analysis <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                                                                                               meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .01, minTrueSimilThresh = minTrueSimil, maxTrueSimilThresh = Inf,minFractNN=.8,
                                                                                                                               clustSizeThresh = 0,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                                                                                               commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=3, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                                                                                               survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                                                                                               CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                                                                                               ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames,
                                                                                                                               GSEAanalysis=TRUE,clinVarPlots=TRUE, fractFeatIntersectThresh=.6,numFeatIntersectThresh =0,clustSizeFractThresh =0,
                                                                                                                               findCommWithWeights=TRUE, plotSimilEdgeWeight = TRUE,plotToScreen=FALSE,fractEdgesInVsOutEdge=0, fractEdgesInVsOutComm=0,minNumEdgesForCluster=1,plotStackedYLimit=800)


##2000 features, no PAM50
numFeatures <- 2000
saveDirPAM50 <- "/home/ywrfc09/breast_analysis/PAM50_analyses"
load(paste0(saveDirPAM50,"/pam50_centroids_updatedSymbols.RData"))
##load(paste0(saveDirPAM50,"/pam50_centroids_updatedSymbols.RData"))
library("CoINcIDE")
saveDir <- "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/"
globalSaveDir <-  "/home/ywrfc09/breast_analysis"
load(paste0(saveDir,"/metaFeatures_",numFeatures,".RData.gzip"))

PAM50genes <- centroidMatrix[,"geneSymbol"]

noPAM50_2000MetaF <- setdiff(metaFeatures$finalFeatures,PAM50genes)

metaFeatures$finalFeatures <- noPAM50_2000MetaF

#careful - sometimes saved as RDATA, sometimes rds.
clusterCoINcIDE_output  <-readRDS(paste0("/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/curatedbreastData_kmeansConsensus_nstart1pItem9",numFeatures,"Features_NOPAM50_2015-12-17.rds"))
clustSampleIndexList <-  clusterCoINcIDE_output$clustSampleIndexList_PACR
clustFeatureIndexList <- clusterCoINcIDE_output$clustFeatureIndexList_PACR

outputFile <- "~/CoINcIDE_messages.txt"

CoINcIDE_output <- readRDS(paste0("/home/ywrfc09/breast_analysis/metaRankWithTop20Genes//CoINcIDE_results_",numFeatures,"F_NOPAM50_pearson_edgeMethod_mean_centroidMethod2015-12-17.rds"))

saveDirNew <- "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/"


load(paste0(globalSaveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies

eset_featureDataFieldName="gene_symbol"
networkColors = "Set3"
outcomesVarBinary=NA
#not enough continuous variables for breast data
outcomesVarCont = NA
ovarian <- FALSE
eset_uniquePatientID="dbUniquePatientID"
fisherTestVariables <- c("DFS","RFS","metastasis","pCR",
                         "tumor_stage_preTrt" ,
                         "preTrt_lymph_node_status" , 
                         "neoadjuvant_or_adjuvant","ER_preTrt","HER2_preTrt","hist_grade","chemotherapyClass",
                         "anti_HER2","anti_estrogen")

fisherTestVariableLegendNames <- c("DFS","RFS","metastasis","pCR","pretreat\ntumor stage",
                                   "pretreat\nlymph node","neoadjuvant\nvs adjuvant","pretreat\nER status",
                                   "pretreat\nHER2 status","hist grade","chemotherapy","anti-HER2\ntreat",
                                   "anti-ER\ntreat")


fisherTestVariableTitleNames <- c("DFS","RFS","metastasis","pCR","pretreatment tumor stage",
                                  "pretreatment lymph node","neoadjuvant vs adjuvant","pretreatment ER status",
                                  "pretreatment HER2 status","histological grade","chemotherapy","anti-HER2 treatment",
                                  "anti-ER treatment")

 source("/home/ywrfc09/CoINcIDE/coincide/GenomeBiology/GenomeBiology_metaFeatures_analysis_wrapper.R")

experimentName <- paste0("breast_",numFeatures,"FNoPAM50_pear_meanCent_",Sys.Date())
densityPlot <- meanMetricDensityPlot(CoINcIDE_output$meanMetricMatrix,saveDir=saveDir,
                                     experimentName=experimentName,savePlot=TRUE)

#looks like simil thresh of around 0.5 best.
minTrueSimil <- 0.5
experimentName <- paste0("breast_",numFeatures,"FNoPAM50_pear_meanCent_MM",minTrueSimil)
breast_2000genes_NoPAM50_pearson_meanCentroid_analysis <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                                                                                               meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .01, minTrueSimilThresh = minTrueSimil, maxTrueSimilThresh = Inf,minFractNN=.8,
                                                                                                                               clustSizeThresh = 0,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                                                                                               commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=3, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                                                                                               survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                                                                                               CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                                                                                               ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames,
                                                                                                                               GSEAanalysis=TRUE,clinVarPlots=TRUE, fractFeatIntersectThresh=.6,numFeatIntersectThresh =0,clustSizeFractThresh =0,
                                                                                                                               findCommWithWeights=TRUE, plotSimilEdgeWeight = TRUE,plotToScreen=FALSE,fractEdgesInVsOutEdge=0, fractEdgesInVsOutComm=0,minNumEdgesForCluster=1,plotStackedYLimit=800)













#######did not run for most recent analyses; was playing around with effects of
#parameters here.












#intersect: p-value, fract, meanMetric


#####now do  spearman
edgeMethod <- "spearman"
load("/home/kplaney/breast_analysis//curatedbreastData_kmeansConsensus_nstart1pItem9200Features_2015-05-04.RData.gzip")
load("/home/kplaney/breast_analysis/metaFeatures_200.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")


load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

clustSampleIndexList <-  kmeansConsensus$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensus$clustFeatureIndexList_PACR



source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast200F_spearman_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                       edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                       sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                       outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                       numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                       clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


save(breast200F_spearman_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_200F_spearman_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")

###now: pearson, but with centroid method
fractFeatIntersectThresh=.8
#we know these all share at least 35 genes.
numFeatIntersectThresh=150
clustSizeThresh=5
clustSizeFractThresh=.05
#more cores: centroid methods run slower!
numParallelCores=4
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "centroid"

#200 features, pearson:
edgeMethod <- "pearson"
load("/home/kplaney/breast_analysis//curatedbreastData_kmeansConsensus_nstart1pItem9200Features_2015-05-04.RData.gzip")
load("/home/kplaney/breast_analysis/metaFeatures_200.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
#load data matrix list
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

clustSampleIndexList <-  kmeansConsensus$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensus$clustFeatureIndexList_PACR

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast200F_pearson_centroid <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                         edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                         sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                         outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                         numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                         clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)



save(breast200F_pearson_centroid ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_200F_pearson_centroid_",Sys.Date(),".RData.gzip"),compress="gzip")

########500:
fractFeatIntersectThresh=.8
#we know these all share at least 35 genes.
numFeatIntersectThresh=425
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=3
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"

# 500features, pearson:
edgeMethod <- "pearson"
load("/home/kplaney/breast_analysis//curatedbreastData_kmeansConsensus_nstart1pItem9500Features_2015-05-04.RData.gzip")
load("/home/kplaney/breast_analysis/metaFeatures_500.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

clustSampleIndexList <-  kmeansConsensus$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensus$clustFeatureIndexList_PACR

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast500F_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                      edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                      sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                      outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                      numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                      clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)



save(breast500F_pearson_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_500F_pearson_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")


##############
#now do  spearman
edgeMethod <- "spearman"
load("/home/kplaney/breast_analysis//curatedbreastData_kmeansConsensus_nstart1pItem9500Features_2015-05-04.RData.gzip")
load("/home/kplaney/breast_analysis/metaFeatures_500.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")

load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}
clustSampleIndexList <-  kmeansConsensus$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensus$clustFeatureIndexList_PACR



source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast500F_spearman_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                       edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                       sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                       outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                       numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                       clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


save(breast500F_spearman_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_500F_spearman_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")

###now pearson, but with centroid:
fractFeatIntersectThresh=.8
#we know these all share at least 35 genes.
numFeatIntersectThresh=425
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=4
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "centroid"

# 500features, pearson:
edgeMethod <- "pearson"
load("/home/kplaney/breast_analysis//curatedbreastData_kmeansConsensus_nstart1pItem9500Features_2015-05-04.RData.gzip")
load("/home/kplaney/breast_analysis/metaFeatures_500.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

clustSampleIndexList <-  kmeansConsensus$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensus$clustFeatureIndexList_PACR

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast500F_pearson_centroid <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                         edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                         sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                         outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                         numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                         clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


#had a bug in paste:
save(breast500F_pearson_centroid,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_500F_pearson_centroid_",Sys.Date(),".RData.gzip"),compress="gzip")


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
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"

# 500features, pearson:
edgeMethod <- "pearson"
load("/home/kplaney/breast_analysis//curatedbreastData_kmeansConsensus_nstart1pItem91000Features_2015-05-05.RData.gzip")
load("/home/kplaney/breast_analysis/metaFeatures_1000.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

clustSampleIndexList <-  kmeansConsensus$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensus$clustFeatureIndexList_PACR

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast1000F_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                       edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                       sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                       outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                       numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                       clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)



save(breast1000F_pearson_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_1000F_pearson_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")



###1000 spearman:

edgeMethod <- "spearman"
load("/home/kplaney/breast_analysis//curatedbreastData_kmeansConsensus_nstart1pItem91000Features_2015-05-05.RData.gzip")
load("/home/kplaney/breast_analysis/metaFeatures_1000.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

clustSampleIndexList <-  kmeansConsensus$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensus$clustFeatureIndexList_PACR

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast1000F_spearman_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                          edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                          sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                          outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                          numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                          clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)



save(breast1000F_spearman_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_1000F_spearman_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")



###now pearson, but with centroid:
fractFeatIntersectThresh=.8
#we know these all share at least 35 genes.
numFeatIntersectThresh=800
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=4
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "centroid"

# 500features, pearson:
edgeMethod <- "pearson"
load("/home/kplaney/breast_analysis//curatedbreastData_kmeansConsensus_nstart1pItem91000Features_2015-05-05.RData.gzip")
load("/home/kplaney/breast_analysis/metaFeatures_1000.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

clustSampleIndexList <-  kmeansConsensus$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensus$clustFeatureIndexList_PACR

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast1000F_pearson_centroid <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                          edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                          sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                          outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                          numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                          clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)



save(breast1000F_pearson_centroid ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_1000F_pearson_centroid_",Sys.Date(),".RData.gzip"),compress="gzip")


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
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"

# 2000features, pearson:
edgeMethod <- "pearson"
load("/home/kplaney/breast_analysis//curatedbreastData_kmeansConsensus_nstart1pItem92000Features_2015-05-05.RData.gzip")
load("/home/kplaney/breast_analysis/metaFeatures_2000.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
#load("/home/kplaney/breast_analysis/esets_proc_TCGAcombat.RData.gzip")

load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}
clustSampleIndexList <-  kmeansConsensus$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensus$clustFeatureIndexList_PACR

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast2000F_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                       edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                       sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                       outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                       numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                       clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


#had a bug in paste:
save(breast2000F_pearson_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_2000F_pearson_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")


###spearman 2000

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
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"

# 2000features, pearson:
edgeMethod <- "spearman"
load("/home/kplaney/breast_analysis//curatedbreastData_kmeansConsensus_nstart1pItem92000Features_2015-05-05.RData.gzip")
load("/home/kplaney/breast_analysis/metaFeatures_2000.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
#load("/home/kplaney/breast_analysis/esets_proc_TCGAcombat.RData.gzip")

load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}
clustSampleIndexList <-  kmeansConsensus$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensus$clustFeatureIndexList_PACR

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast2000F_spearman_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                          edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                          sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                          outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                          numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                          clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)



save(breast2000F_spearman_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_2000F_spearman_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")

###now with pearson, but centroid
fractFeatIntersectThresh=.85
numFeatIntersectThresh=1700
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=4
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "centroid"

# 2000features, pearson:
edgeMethod <- "pearson"
load("/home/kplaney/breast_analysis//curatedbreastData_kmeansConsensus_nstart1pItem92000Features_2015-05-05.RData.gzip")
load("/home/kplaney/breast_analysis/metaFeatures_2000.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
#load("/home/kplaney/breast_analysis/esets_proc_TCGAcombat.RData.gzip")

load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}
clustSampleIndexList <-  kmeansConsensus$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensus$clustFeatureIndexList_PACR

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast2000F_pearson_centroid <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                          edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                          sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                          outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                          numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                          clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)



save(breast2000F_pearson_centroid,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_2000F_pearson_meanMatrix_centroid_",Sys.Date(),".RData.gzip"),compress="gzip")

#########
####gap test results as opposed to k-means consensus.
##gap test with k-means

###200 features
fractFeatIntersectThresh=.8
#we know these all share at least 35 genes.
numFeatIntersectThresh=150
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=3
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"

#200 features, pearson:
edgeMethod <- "pearson"
load("/home/kplaney/breast_analysis/curatedbreastData_kmeansGap_nstart25_200_features_2015-05-18.RData.gzip")
load("/home/kplaney/breast_analysis/metaFeatures_200.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
#load data matrix list
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

clustSampleIndexList <-  kmeansGap$clustSampleIndexList
clustFeatureIndexList <- kmeansGap$clustFeatureIndexList

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast200F_kmeansGap_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                                   edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                                   sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                                   outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                                   numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                                   clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)



save(breast200F_kmeansGap_pearson_meanMatrix,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_breast200F_kmeansGap_pearson_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")


###500 features
fractFeatIntersectThresh=.8
#we know these all share at least 35 genes.
numFeatIntersectThresh=425
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=3
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"

# 500features, pearson:
edgeMethod <- "pearson"
load("/home/kplaney/breast_analysis/curatedbreastData_kmeansGap_nstart25_500_features_2015-05-18.RData.gzip")
load("/home/kplaney/breast_analysis/metaFeatures_500.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

clustSampleIndexList <-  kmeansGap$clustSampleIndexList
clustFeatureIndexList <- kmeansGap$clustFeatureIndexList

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast500F_kmeansGap_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                                   edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                                   sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                                   outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                                   numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                                   clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)



save(breast500F_kmeansGap_pearson_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_breast500F_kmeansGap_pearson_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")


#TO DO: 1000,2000 have not finished running yet.
#1000

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
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"


edgeMethod <- "pearson"
#TO DO: this hasn't finished running yet:
load("/home/kplaney/breast_analysis//curatedbreastData_kmeansGap_nstart25_1000_features_2015-05-18.RData.gzip")
load("/home/kplaney/breast_analysis/metaFeatures_1000.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

clustSampleIndexList <-  kmeansGap$clustSampleIndexList
clustFeatureIndexList <- kmeansGap$clustFeatureIndexList

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast1000F_kmeansGap_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                                    edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                                    sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                                    outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                                    numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                                    clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)



save(breast1000F_kmeansGap_pearson_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_breast1000F_kmeansGap_pearson_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")


###2000
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
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"

# 2000features, pearson:
edgeMethod <- "pearson"
load("/home/kplaney/breast_analysis/curatedbreastData_kmeansGap_nstart25_2000_features_")
load("/home/kplaney/breast_analysis/metaFeatures_2000.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
#load("/home/kplaney/breast_analysis/esets_proc_TCGAcombat.RData.gzip")

load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}
clustSampleIndexList <-  hclustOut$clustSampleIndexList
clustFeatureIndexList <- hclustOut$clustFeatureIndexList

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast2000F_kmeansGap_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                                    edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                                    sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                                    outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                                    numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                                    clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)



save(breast2000F_kmeansGap_pearson_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_2000F_kmeansGap_pearson_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")


##########
##gap test using hierarchical clustering:
###200 features
fractFeatIntersectThresh=.8
#we know these all share at least 35 genes.
numFeatIntersectThresh=150
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=3
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"

#200 features, pearson:
edgeMethod <- "pearson"
load("/home/kplaney/breast_analysis/curatedbreastData_hclust_200Features_2015-05-15.RData.gzip")
load("/home/kplaney/breast_analysis/metaFeatures_200.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
#load data matrix list
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){

dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]

}

clustSampleIndexList <-  hclustOut$clustSampleIndexList
clustFeatureIndexList <- hclustOut$clustFeatureIndexList

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast200F_hclustGap_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)



save(breast200F_hclustGap_pearson_meanMatrix,file=
paste0("/home/kplaney/breast_analysis/adjMatrices_breast200F_hclustGap_pearson_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")


###500 features
fractFeatIntersectThresh=.8
#we know these all share at least 35 genes.
numFeatIntersectThresh=425
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=3
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"

# 500features, pearson:
edgeMethod <- "pearson"
load("/home/kplaney/breast_analysis/curatedbreastData_hclust_500Features_2015-05-15.RData.gzip")
load("/home/kplaney/breast_analysis/metaFeatures_500.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

clustSampleIndexList <-  hclustOut$clustSampleIndexList
clustFeatureIndexList <- hclustOut$clustFeatureIndexList

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast500F_hclustGap_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                         edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                         sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                         outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                         numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                         clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)



save(breast500F_hclustGap_pearson_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_breast500F_hclustGap_pearson_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")


#
#10000

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
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"


edgeMethod <- "pearson"
load("/home/kplaney/breast_analysis/curatedbreastData_hclust_1000Features_2015-05-15.RData.gzip")
load("/home/kplaney/breast_analysis/metaFeatures_1000.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

clustSampleIndexList <-  hclustOut$clustSampleIndexList
clustFeatureIndexList <- hclustOut$clustFeatureIndexList

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast1000F_hclustGap_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                          edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                          sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                          outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                          numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                          clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)



save(breast1000F_hclustGap_pearson_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_breast1000F_hclustGap_pearson_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")


###2000
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
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"

# 2000features, pearson:
edgeMethod <- "pearson"
load("/home/kplaney/breast_analysis/curatedbreastData_hclust_2000Features_2015-05-16.RData.gzip")
load("/home/kplaney/breast_analysis/metaFeatures_2000.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
#load("/home/kplaney/breast_analysis/esets_proc_TCGAcombat.RData.gzip")

load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}
clustSampleIndexList <-  hclustOut$clustSampleIndexList
clustFeatureIndexList <- hclustOut$clustFeatureIndexList

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast2000F_hclustGap_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                          edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                          sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                          outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                          numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                          clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)



save(breast2000F_hclustGap_pearson_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_2000F_hclustGap_pearson_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")



