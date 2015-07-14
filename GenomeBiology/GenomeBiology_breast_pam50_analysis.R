library("CoINcIDE")
saveDirPAM50 <- "/home/ywrfc09/breast_analysis/PAM50_analyses"
saveDir <- "/home/ywrfc09/breast_analysis/"

outputFile <- "~/CoINcIDE_messages.txt"


####analyze
  source("/home/ywrfc09/CoINcIDE/coincide/oldCode/CoINcIDE_metaFeatures_analysis_wrapper.R")
  #grab data matrix list, clust features list
  load(paste0(saveDirPAM50,"/pam50_centroids_updatedSymbols.RData"))
  pam50Genes <- centroidMatrix[,1]
  metaFeatures <- list(finalFeatures=pam50Genes,datasetListIndicesToRemove=NULL)
  load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
  esets=esets_minVar001_17_studies
  load(paste0(saveDirPAM50,"/pam50Full_centroidCluster.RData.gzip"))
  clusterCoINcIDE_output =  readRDS("/home/ywrfc09/breast_analysis/PAM50_analyses/pam50Full_centroidCluster.rds")
  #need to add fake PACR variable
  clusterCoINcIDE_output$bestK_PACR <- unlist(sapply(clusterCoINcIDE_output$clustSampleIndexList,FUN=function(unit){
    
    length(unit)
    
  })) 
  
  names(clusterCoINcIDE_output$bestK_PACR) <- names(esets)
  clusterCoINcIDE_output$clustSampleIndexList_PACR <- clusterCoINcIDE_output$clustSampleIndexList
  clusterCoINcIDE_output$clustFeatureIndexList_PACR <- clusterCoINcIDE_output$clustFeatureIndexList
  
  CoINcIDE_output = readRDS("/home/ywrfc09/breast_analysis/PAM50_analyses/CoINcIDE_results_PAM50centroidCluster_pearson_edgeMethod_mean_centroidMethod2015-07-08.rds")
  experimentName <- "pam50Centroid_pear_meanCent"
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
  
  #options(device=NULL)
  breast_pam50FullCentroidsOut <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                              meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .01, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,minFractNN=.8,
                                                              clustSizeThresh = 0,saveDir =saveDirPAM50,experimentName = experimentName,networkColors = networkColors,
                                                              commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=3, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                              survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                              CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                              ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames,
                                                              GSEAanalysis=TRUE, fractFeatIntersectThresh=.6,numFeatIntersectThresh =0,clustSizeFractThresh =0,
                                                              findCommWithWeights=TRUE, plotSimilEdgeWeight = TRUE,plotToScreen=FALSE,fractEdgesInVsOutEdge=.5, fractEdgesInVsOutComm=0)
  
  
  
 
  
  



  
####PAM50 short
  library("CoINcIDE")
saveDirPAM50 <- "/home/ywrfc09/breast_analysis/PAM50_analyses"
saveDir <- "/home/ywrfc09/breast_analysis/"
outputFile <- "~/CoINcIDE_messages.txt"

load("~/pam50Short_kmeansConsClusters_pearson_meanCentroid.RData.gzip")
CoINcIDE_output = readRDS("/home/ywrfc09/breast_analysis/PAM50_analyses//CoINcIDE_results_PAM50kmeansShort_pearson_edgeMethod_mean_centroidMethod2015-07-08.rds")

#load esets
load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
#35/50 is 70%
load(paste0(saveDirPAM50,"/pam50Short_genes.RData"))
metaFeatures <- list(finalFeatures=pam50Short,datasetListIndicesToRemove=NULL)
clusterCoINcIDE_output <- readRDS("/home/ywrfc09/breast_analysis/PAM50_analyses//kmeansConsensuspam50_short_Nstart1pItem9_pam50ShortFeatures_2015-04-28.rds")
     
     clustSampleIndexList <- clusterCoINcIDE_output$clustSampleIndexList_PACR
     clustFeatureIndexList <- clusterCoINcIDE_output$clustFeatureIndexList_PACR
     
     load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
     esets=esets_minVar001_17_studies


 
     experimentName <- "pam50ShortKmeans_pear_meanCent"
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
     
     source("/home/ywrfc09/CoINcIDE/coincide/oldCode/CoINcIDE_metaFeatures_analysis_wrapper.R")   
breast_pam50Short_pearson_centroidMean <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                            meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .01, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,minFractNN=.8,
                                                            clustSizeThresh = 0,saveDir =saveDirPAM50,experimentName = experimentName,networkColors = networkColors,
                                                            commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                            survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                            CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                            ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames,
                                                            GSEAanalysis=TRUE,clinVarPlots=TRUE, fractFeatIntersectThresh=.6,numFeatIntersectThresh =0,clustSizeFractThresh =0,
                                                            findCommWithWeights=TRUE, plotSimilEdgeWeight = TRUE,plotToScreen=FALSE,fractEdgesInVsOutEdge=.5, fractEdgesInVsOutComm=0)


###PAM50 full
library("CoINcIDE")
saveDirPAM50 <- "/home/ywrfc09/breast_analysis/PAM50_analyses"
saveDir <- "/home/ywrfc09/breast_analysis/"
outputFile <- "~/CoINcIDE_messages.txt"
     #grab data matrix list, clust features list
     load(paste0(saveDirPAM50,"/pam50_centroids_updatedSymbols.RData"))
     pam50Genes <- centroidMatrix[,1]
     metaFeatures <- list(finalFeatures=pam50Genes,datasetListIndicesToRemove=NULL)
     load("~/pam50Short_kmeansConsClusters_pearson_meanCentroid.RData.gzip")
     CoINcIDE_output = readRDS("/home/ywrfc09/breast_analysis/PAM50_analyses//CoINcIDE_results_PAM50kmeansFull_pearson_edgeMethod_mean_centroidMethod2015-07-08.rds")
     
     #load esets
     load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
     esets=esets_minVar001_17_studies
     #35/50 is 70%
     load(paste0(saveDirPAM50,"/pam50Short_genes.RData"))
     metaFeatures <- list(finalFeatures=pam50Short,datasetListIndicesToRemove=NULL)
     clusterCoINcIDE_output <- readRDS("/home/ywrfc09/breast_analysis/PAM50_analyses//kmeansConsensuspam50_full_Nstart1pItem9_pam50FullFeatures_2015-04-28.rds")
     
     clustSampleIndexList <- clusterCoINcIDE_output$clustSampleIndexList_PACR
     clustFeatureIndexList <- clusterCoINcIDE_output$clustFeatureIndexList_PACR
     
     load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
     esets=esets_minVar001_17_studies
     
     
     
     
     experimentName <- "pam50FullKmeans_pear_meanCent"
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
     
     source("/home/ywrfc09/CoINcIDE/coincide/oldCode/CoINcIDE_metaFeatures_analysis_wrapper.R")   
     
  breast_pam50Full_pearsonMeanCentroid <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                                         meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .01, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,minFractNN=.8,
                                                                         clustSizeThresh = 0,saveDir =saveDirPAM50,experimentName = experimentName,networkColors = networkColors,
                                                                         commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=3, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                                         survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                                         CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                                         ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames,
                                                                         GSEAanalysis=TRUE,clinVarPlots=TRUE, fractFeatIntersectThresh=.6,numFeatIntersectThresh =0,clustSizeFractThresh =0,
                                                                         findCommWithWeights=TRUE, plotSimilEdgeWeight = TRUE,plotToScreen=FALSE,fractEdgesInVsOutEdge=.5, fractEdgesInVsOutComm=0)
     
     
   
 ######################## 
  ####have not run....
  
    #####pearson, median centroid
     
     saveDirPAM50 <- "/home/ywrfc09/breast_analysis/PAM50_analyses"
     saveDir <- "/home/ywrfc09/breast_analysis/"
     ####pam50 centroids clustering
     load(paste0(saveDirPAM50,"/pam50Full_centroidCluster.RData.gzip"))
     #load data matrix list
     load(paste0(saveDir,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
     
     clustSampleIndexList <- pam50Full_centroidCluster$clustSampleIndexList
     clustFeatureIndexList <- pam50Full_centroidCluster$clustFeatureIndexList
     
     edgeMethod <- "pearson"
     centroidMethod <- "mean"
     outputFile <- "~/CoINcIDE_messages.txt"
     numSims <- 500
     source("/home/ywrfc09/CoINcIDE/coincide/CoINcIDE_packageVersion/CoINcIDE/R/CoINcIDE_computeEdges.R")
     pam50Full_centroidCluster_pearson_meanCentroid <-CoINcIDE_getAdjMatrices(dataMatrixList,clustSampleIndexList,clustFeatureIndexList,
                                                                              edgeMethod=edgeMethod,centroidMethod=centroidMethod,
                                                                              numSims=500,
                                                                              outputFile=outputFile)
     
     
     
     
     create_nullMatrixList_results <- function(dataMatrixList,numIter=5,numParallelCores=1,
                                               clustSampleIndexList,clustFeatureIndexList,
                                               edgeMethod=c("distCor","spearman","pearson","kendall","Euclidean","cosine",
                                                            "Manhattan","Minkowski","Mahalanobis"),
                                               numSims=500,
                                               outputFile="./CoINcIDE_messages.txt",
                                               centroidMethod=c("mean","median"))
       
       global_FDR <- function(CoINcIDE_outputList,
                              edgeMethod=c("distCor","spearman","pearson","kendall","Euclidean","cosine",
                                           "Manhattan","Minkowski","Mahalanobis"),numParallelCores=1,minTrueSimilThresh=-Inf,maxTrueSimilThresh=Inf,
                              sigMethod=c("meanMatrix","centroid"),numSims=500,
                              outputFile="./CoINcIDE_messages.txt",fractFeatIntersectThresh=0,numFeatIntersectThresh=0 ,clustSizeThresh=0, clustSizeFractThresh=0,
                              checkNA=FALSE,centroidMethod=c("mean","median"), 
                              meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, 
                              saveDir = "/home/kplaney/ovarian_analysis/",experimentName = "nullTest",networkColors = "Set3",
                              commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4,clustIndexMatrix,minFractNN =.7,findCommWithWeights=FALSE
                              
                              
       )
     
     
     ####analyze
     source("/home/ywrfc09/CoINcIDE/coincide/oldCode/CoINcIDE_metaFeatures_analysis_wrapper.R")
     #grab data matrix list, clust features list
     load(paste0(saveDirPAM50,"/pam50_centroids_updatedSymbols.RData"))
     pam50Genes <- centroidMatrix[,1]
     metaFeatures <- list(finalFeatures=pam50Genes,datasetListIndicesToRemove=NULL)
     
     load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
     esets=esets_minVar001_17_studies
     load(paste0(saveDirPAM50,"/pam50Full_centroidCluster.RData.gzip"))
     clusterCoINcIDE_output =  pam50Full_centroidCluster
     #need to add fake PACR variable
     clusterCoINcIDE_output$bestK_PACR <- unlist(sapply(clusterCoINcIDE_output$clustSampleIndexList,FUN=function(unit){
       
       length(unit)
       
     })) 
     
     names(clusterCoINcIDE_output$bestK_PACR) <- names(esets)
     clusterCoINcIDE_output$clustSampleIndexList_PACR <- clusterCoINcIDE_output$clustSampleIndexList
     clusterCoINcIDE_output$clustFeatureIndexList_PACR <- clusterCoINcIDE_output$clustFeatureIndexList
     
     #test run for now
     load("~/test_PAM50_pearson_centroidMean_extraParam.RData.gzip")
     CoINcIDE_output = PAM50_pearson_centroidMean_extraParam
     experimentName <- "breast_pam50FullCentroidMean_features"
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
     
     breast_pam50FullCentroidsOut <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                                 meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .01, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,minFractNN=.8,
                                                                 clustSizeThresh = 0,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                                 commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=3, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                                 survivalAnalysis=FALSE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                                 CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                                 ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames,
                                                                 GSEAanalysis=FALSE,clinVarPlots=TRUE, fractFeatIntersectThresh=.6,numFeatIntersectThresh =0,clustSizeFractThresh =0,
                                                                 findCommWithWeights=TRUE, plotSimilEdgeWeight = TRUE,plotToScreen=TRUE)
     
     
     
     
     
     
     
     ####PAM50 short
     saveDirPAM50 <- "/home/ywrfc09/breast_analysis/PAM50_analyses"
     saveDir <- "/home/ywrfc09/breast_analysis/"
     
     load(paste0(saveDir,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
     
     load(paste0(saveDirPAM50,"/kmeansConsensuspam50_short_Nstart1pItem9_pam50ShortFeatures_2015-04-28.RData.gzip"))
     
     clustSampleIndexList <- kmeansConsensuspam50_short_Nstart1pItem9$clustSampleIndexList_PACR
     clustFeatureIndexList <- kmeansConsensuspam50_short_Nstart1pItem9$clustFeatureIndexList_PACR
     
     edgeMethod <- "pearson"
     centroidMethod <- "mean"
     outputFile <- "~/CoINcIDE_messages.txt"
     numSims <- 500
     source("/home/ywrfc09/CoINcIDE/coincide/CoINcIDE_packageVersion/CoINcIDE/R/CoINcIDE_computeEdges.R")
     pam50Short_kmeansConsClusters_pearson_meanCentroid <-CoINcIDE_getAdjMatrices(dataMatrixList,clustSampleIndexList,clustFeatureIndexList,
                                                                                  edgeMethod=edgeMethod,centroidMethod=centroidMethod,
                                                                                  numSims=500,
                                                                                  outputFile=outputFile)
     
     
     
     
     load("~/pam50Short_kmeansConsClusters_pearson_meanCentroid.RData.gzip")
     CoINcIDE_output = pam50Short_kmeansConsClusters_pearson_meanCentroid
     
     #load data matrix list
     load(paste0(saveDir,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
     #35/50 is 70%
     load(paste0(saveDirPAM50,"/pam50Short_genes.RData"))
     metaFeatures <- list(finalFeatures=pam50Short,datasetListIndicesToRemove=NULL)
     load(paste0(saveDirPAM50,"//kmeansConsensuspam50_short_Nstart1pItem9_pam50ShortFeatures_2015-04-28.RData.gzip"))
     
     clustSampleIndexList <- kmeansConsensuspam50_short_Nstart1pItem9$clustSampleIndexList_PACR
     clustFeatureIndexList <- kmeansConsensuspam50_short_Nstart1pItem9$clustFeatureIndexList_PACR
     clusterCoINcIDE_output =  kmeansConsensuspam50_short_Nstart1pItem9
     
     load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
     esets=esets_minVar001_17_studies
     
     saveDirNew <- "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/"
     
     
     experimentName <- "breast_pam50ShortKmeans_centroidMean"
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
     
     source("/home/ywrfc09/CoINcIDE/coincide/oldCode/CoINcIDE_metaFeatures_analysis_wrapper.R")   
     
     breast_pam50Short_pearson_centroidMean <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                                           meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .5, maxTrueSimilThresh = Inf,minFractNN=.8,
                                                                           clustSizeThresh = 0,saveDir =saveDirNew,experimentName = experimentName,networkColors = networkColors,
                                                                           commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=3, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                                           survivalAnalysis=FALSE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                                           CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                                           ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames,
                                                                           GSEAanalysis=FALSE,clinVarPlots=TRUE, fractFeatIntersectThresh=.6,numFeatIntersectThresh =0,clustSizeFractThresh =0,
                                                                           findCommWithWeights=TRUE, plotSimilEdgeWeight = TRUE,plotToScreen=TRUE)
     
     
     ###PAM50 full
     saveDirPAM50 <- "/home/ywrfc09/breast_analysis/PAM50_analyses"
     saveDir <- "/home/ywrfc09/breast_analysis/"
     #grab data matrix list, clust features list
     load(paste0(saveDirPAM50,"/pam50_centroids_updatedSymbols.RData"))
     pam50Genes <- centroidMatrix[,1]
     metaFeatures <- list(finalFeatures=pam50Genes,datasetListIndicesToRemove=NULL)
     
     load(paste0(saveDir,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
     
     edgeMethod <- "pearson"
     load(paste0(saveDirPAM50,"//kmeansConsensuspam50_full_Nstart1pItem9_pam50FullFeatures_2015-04-28.RData.gzip"))
     
     clustSampleIndexList <- kmeansConsensuspam50_full_Nstart1pItem9$clustSampleIndexList_PACR
     clustFeatureIndexList <- kmeansConsensuspam50_full_Nstart1pItem9$clustFeatureIndexList_PACR
     clusterCoINcIDE_output =  kmeansConsensuspam50_full_Nstart1pItem9
     
     edgeMethod <- "pearson"
     centroidMethod <- "mean"
     outputFile <- "~/CoINcIDE_messages.txt"
     numSims <- 500
     source("/home/ywrfc09/CoINcIDE/coincide/CoINcIDE_packageVersion/CoINcIDE/R/CoINcIDE_computeEdges.R")
     pam50Full_kmeansConsClusters_pearson_meanCentroid <-CoINcIDE_getAdjMatrices(dataMatrixList,clustSampleIndexList,clustFeatureIndexList,
                                                                                 edgeMethod=edgeMethod,centroidMethod=centroidMethod,
                                                                                 numSims=500,
                                                                                 outputFile=outputFile)
     
     
     
     
     
     load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
     esets=esets_minVar001_17_studies
     
     saveDirNew <- "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/"
     
     
     experimentName <- "breast_pam50FullKmeans_centroidMean"
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
     
     source("/home/ywrfc09/CoINcIDE/coincide/oldCode/CoINcIDE_metaFeatures_analysis_wrapper.R")   
     
     
     load("~/pam50Full_kmeansConsClusters_pearson_meanCentroid.RData.gzip")
     CoINcIDE_output <- pam50Full_kmeansConsClusters_pearson_meanCentroid
     
     saveDirNew <- "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/"
     breast_pam50Full_pearsonMeanCentroid <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                                         meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .5, maxTrueSimilThresh = Inf,minFractNN=.8,
                                                                         clustSizeThresh = 0,saveDir =saveDirNew,experimentName = experimentName,networkColors = networkColors,
                                                                         commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=3, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                                         survivalAnalysis=FALSE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                                         CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                                         ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames,
                                                                         GSEAanalysis=FALSE,clinVarPlots=TRUE, fractFeatIntersectThresh=.6,numFeatIntersectThresh =0,clustSizeFractThresh =0,
                                                                         findCommWithWeights=FALSE, plotSimilEdgeWeight = TRUE,plotToScreen=TRUE)
     
     
     
     
     
     
     
     
     
     
     ##########   
####pam50 centroids clustering
load(paste0(saveDir,"/pam50Full_centroidCluster.RData.gzip"))
#load data matrix list
load(paste0(saveDir,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))

clustSampleIndexList <- pam50Full_centroidCluster$clustSampleIndexList
clustFeatureIndexList <- pam50Full_centroidCluster$clustFeatureIndexList
#35/50 is 70%
fractFeatIntersectThresh=.69
#we know these all share at least 35 genes.
numFeatIntersectThresh=34
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=3
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile=paste0(saveDir,"/CoINcIDE_messages.txt")
centroidMethod <- "mean"
edgeMethod <- "pearson"
#not used here
maxNullFractSize=.2

#test FDR first:

source("/home/ywrfc09CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")
pam50Full_centroidCluster_pearson_meanMatrix <-CoINcIDE_getAdjMatrices(dataMatrixList,clustSampleIndexList,clustFeatureIndexList,
                                                                                   edgeMethod=edgeMethod,centroidMethod=centroidMethod,
                                                                                   numSims=500,
                                                                                   outputFile="~/CoINcIDE_messages.txt")

                                          
##test median centroids now:                                                                       
load(paste0(saveDirPAM50,"//pam50Full_centroidCluster.RData.gzip")
#load data matrix list
load(paste0(saveDir,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))

clustSampleIndexList <- pam50Full_centroidCluster$clustSampleIndexList
clustFeatureIndexList <- pam50Full_centroidCluster$clustFeatureIndexList
#35/50 is 70%
fractFeatIntersectThresh=.69
#we know these all share at least 35 genes.
numFeatIntersectThresh=34
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=4
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile=paste0(savePAM50,"breast_analysis//CoINcIDE_messages.txt")
centroidMethod <- "median"
edgeMethod <- "pearson"


source("/home/ywrfc09/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")
pam50Full_centroidCluster_pearson_centroid <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                                        edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                                        sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                                        outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                                        numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                                        clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


save(pam50Full_centroidCluster_pearson_centroid,file=paste0("/home/kplaney/breast_analysis/adjMatrices_pam50Full_centroidCluster_pearson_centroid_",Sys.Date(),".RData.gzip"),compress="gzip")



###and run with spearman
load(paste0(saveDirPAM50,"//pam50Full_centroidCluster.RData.gzip"))
#load data matrix list
load(paste0(saveDir,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))

clustSampleIndexList <- pam50Full_centroidCluster$clustSampleIndexList
clustFeatureIndexList <- pam50Full_centroidCluster$clustFeatureIndexList
#35/50 is 70%
fractFeatIntersectThresh=.69
#we know these all share at least 35 genes.
numFeatIntersectThresh=34
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=4
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500

outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
centroidMethod <- "mean"
edgeMethod <- "spearman"


#FDR first?
source("/home/ywrfc09CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")
pam50Full_centroidCluster_spearman_meanMatrix <- source("/home/ywrfc09CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")
pam50Full_centroidCluster_pearson_meanMatrix <-CoINcIDE_getAdjMatrices(dataMatrixList,clustSampleIndexList,clustFeatureIndexList,
                                                                       edgeMethod=edgeMethod,centroidMethod=centroidMethod,
                                                                       numSims=500,
                                                                       outputFile="~/CoINcIDE_messages.txt")



save(pam50Full_centroidCluster_spearman_meanMatrix,file=paste0("/home/kplaney/breast_analysis/adjMatrices_pam50Full_centroidCluster_spearman_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")

###STOPPED HERE
####spearman, median
load(paste0(saveDirPAM50,"//pam50Full_centroidCluster.RData.gzip"))
#load data matrix list
load(paste0(saveDir,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))

clustSampleIndexList <- pam50Full_centroidCluster$clustSampleIndexList
clustFeatureIndexList <- pam50Full_centroidCluster$clustFeatureIndexList
#35/50 is 70%
fractFeatIntersectThresh=.69
#we know these all share at least 35 genes.
numFeatIntersectThresh=34
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=4
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"

centroidMethod <- "median"
edgeMethod <- "spearman"


#FDR first?
create_nullMatrixList_results <- function(dataMatrixList,numIter=5,numParallelCores=1,
                                          clustSampleIndexList,clustFeatureIndexList,
                                          edgeMethod=c("distCor","spearman","pearson","kendall","Euclidean","cosine",
                                                       "Manhattan","Minkowski","Mahalanobis"),
                                          numSims=500,
                                          outputFile="./CoINcIDE_messages.txt",
                                          centroidMethod=c("mean","median")){
  
  global_FDR <- function(CoINcIDE_outputList,
                         edgeMethod=c("distCor","spearman","pearson","kendall","Euclidean","cosine",
                                      "Manhattan","Minkowski","Mahalanobis"),numParallelCores=1,minTrueSimilThresh=-Inf,maxTrueSimilThresh=Inf,
                         sigMethod=c("meanMatrix","centroid"),numSims=500,
                         outputFile="./CoINcIDE_messages.txt",fractFeatIntersectThresh=0,numFeatIntersectThresh=0 ,clustSizeThresh=0, clustSizeFractThresh=0,
                         checkNA=FALSE,centroidMethod=c("mean","median"), 
                         meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, 
                         saveDir = "/home/kplaney/ovarian_analysis/",experimentName = "nullTest",networkColors = "Set3",
                         commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4,clustIndexMatrix,minFractNN =.7,findCommWithWeights=FALSE
                         
                         
  )
source("/home/ywrfc09/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")
pam50Full_centroidCluster_spearman_medianCentroid <-CoINcIDE_getAdjMatrices(dataMatrixList,clustSampleIndexList,clustFeatureIndexList,
                                                                       edgeMethod=edgeMethod,centroidMethod=centroidMethod,
                                                                       numSims=500,
                                                                       outputFile="~/CoINcIDE_messages.txt")


#######
####now do k-means pam50 results:
#load data matrix list
load(paste0(saveDir,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
#35/50 is 70%
fractFeatIntersectThresh=.69
#we know these all share at least 35 genes.
numFeatIntersectThresh=34
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=3
#not used here
maxNullFractSize=.2
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"


#full features, pearson
edgeMethod <- "pearson"
load(paste0(saveDir,"//kmeansConsensuspam50_full_Nstart1pItem9_pam50FullFeatures_2015-04-28.RData.gzip")

clustSampleIndexList <- kmeansConsensuspam50_full_Nstart1pItem9$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensuspam50_full_Nstart1pItem9$clustFeatureIndexList_PACR

source("/home/ywrfc09CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")

pam50Full_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                        edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                        sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                        outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                        numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                        clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)

save(pam50Full_pearson_meanMatrix,file=paste0("/home/kplaney/breast_analysis/adjMatrices_pam50Full_pearson_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")


#now spearman:
edgeMethod <- "spearman"
load(paste0(saveDir,"//hclustConsensuspam50_full_pItem9_pam50FullFeatures_2015-04-28.RData.gzip")

clustSampleIndexList <- kmeansConsensuspam50_full_Nstart1pItem9$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensuspam50_full_Nstart1pItem9$clustFeatureIndexList_PACR

source("/home/ywrfc09CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")

pam50Full_spearman_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                        edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                        sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                        outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                        numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                        clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


save(pam50Full_spearman_meanMatrix,file=paste0("/home/kplaney/breast_analysis/adjMatrices_pam50Full_spearman_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")
############
#with pearson, but centroid
#load data matrix list
load(paste0(saveDir,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
#35/50 is 70%
fractFeatIntersectThresh=.69
#we know these all share at least 35 genes.
numFeatIntersectThresh=34
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


#full features, pearson
edgeMethod <- "pearson"
load(paste0(saveDir,"//kmeansConsensuspam50_full_Nstart1pItem9_pam50FullFeatures_2015-04-28.RData.gzip")

clustSampleIndexList <- kmeansConsensuspam50_full_Nstart1pItem9$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensuspam50_full_Nstart1pItem9$clustFeatureIndexList_PACR

source("/home/ywrfc09CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")

pam50Full_pearson_centroid <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                        edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                        sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                        outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                        numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                        clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)

save(pam50Full_pearson_centroid,file=paste0("/home/kplaney/breast_analysis/adjMatrices_pam50Full_pearson_centroid_",Sys.Date(),".RData.gzip"),compress="gzip")



####now pam50 short:

#load data matrix list
load(paste0(saveDir,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
#35/50 is 70%
fractFeatIntersectThresh=.69
#we know these all share at least 35 genes.
numFeatIntersectThresh=34
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


#full features, pearson
edgeMethod <- "pearson"
load(paste0(saveDir,"//kmeansConsensuspam50_short_Nstart1pItem9_pam50ShortFeatures_2015-04-28.RData.gzip")

clustSampleIndexList <- kmeansConsensuspam50_short_Nstart1pItem9$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensuspam50_short_Nstart1pItem9$clustFeatureIndexList_PACR

source("/home/ywrfc09CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")

pam50Short_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                         edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                         sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                         outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                         numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                         clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


save(pam50Short_pearson_meanMatrix="/home/kplaney/breast_analysis/adjMatrices_pam50Short_pearson_meanMatrix_",Sys.Date(),".RData.gzip",compress="gzip")

######
load(paste0(saveDir,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))


#35/50 is 70%
fractFeatIntersectThresh=.69
#we know these all share at least 35 genes.
numFeatIntersectThresh=34
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


#full features, pearson
edgeMethod <- "spearman"
load(paste0(saveDir,"//kmeansConsensuspam50_short_Nstart1pItem9_pam50ShortFeatures_2015-04-28.RData.gzip")

clustSampleIndexList <- kmeansConsensuspam50_short_Nstart1pItem9$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensuspam50_short_Nstart1pItem9$clustFeatureIndexList_PACR

source("/home/ywrfc09CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")

pam50Short_spearman_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                         edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                         sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                         outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                         numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                         clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


save(pam50Short_spearman_meanMatrix,file=paste0("/home/kplaney/breast_analysis/adjMatrices_pam50Short_spearman_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")

#####pearson, but with centroid
#load data matrix list
load(paste0(saveDir,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
#35/50 is 70%
fractFeatIntersectThresh=.69
#we know these all share at least 35 genes.
numFeatIntersectThresh=34
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=3
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "centroid"


#full features, pearson
edgeMethod <- "pearson"
load(paste0(saveDir,"//kmeansConsensuspam50_short_Nstart1pItem9_pam50ShortFeatures_2015-04-28.RData.gzip")

clustSampleIndexList <- kmeansConsensuspam50_short_Nstart1pItem9$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensuspam50_short_Nstart1pItem9$clustFeatureIndexList_PACR

source("/home/ywrfc09CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")

pam50Short_pearson_centroid <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                         edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                         sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                         outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                         numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                         clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


save(pam50Short_pearson_centroid="/home/kplaney/breast_analysis/adjMatrices_pam50Short_pearson_centroid_",Sys.Date(),".RData.gzip",compress="gzip")


##########gap test versions
###k-means
##pam50Short
load(paste0(saveDir,"//gapTestKmeans_pam50Short_nstart252015-05-04.RData.gzip")
load(paste0(saveDir,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
#35/50 is 70%
fractFeatIntersectThresh=.69
#we know these all share at least 35 genes.
numFeatIntersectThresh=34
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=3
#not used here
maxNullFractSize=.2
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"


#full features, pearson
edgeMethod <- "pearson"

clustSampleIndexList <- kmeansGapTest_pam50_short_Nstart25$clustSampleIndexList
clustFeatureIndexList <- kmeansGapTest_pam50_short_Nstart25$clustFeatureIndexList

source("/home/ywrfc09CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")

pam50ShortGapTest_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                        edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                        sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                        outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                        numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                        clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


save(pam50ShortGapTest_pearson_meanMatrix="/home/kplaney/breast_analysis/adjMatrices_pam50ShortGapTest_pearson_meanMatrix_",Sys.Date(),".RData.gzip",compress="gzip")

###pam50Full
load(paste0(saveDir,"//gapTestKmeans_pam50Full_nstart252015-05-04.RData.gzip")
load(paste0(saveDir,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
#35/50 is 70%
fractFeatIntersectThresh=.69
#we know these all share at least 35 genes.
numFeatIntersectThresh=34
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=3
#not used here
maxNullFractSize=.2
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"


#full features, pearson
edgeMethod <- "pearson"

clustSampleIndexList <- kmeansGapTest_pam50_full_Nstart25$clustSampleIndexList
clustFeatureIndexList <- kmeansGapTest_pam50_full_Nstart25$clustFeatureIndexList

source("/home/ywrfc09CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")

pam50FullGapTest_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                                edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                                sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                                outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                                numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                                clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


save(pam50FullGapTest_pearson_meanMatrix="/home/kplaney/breast_analysis/adjMatrices_pam50FullGapTest_pearson_meanMatrix_",Sys.Date(),".RData.gzip",compress="gzip")



####with hierarchical clustering, gap test.
load(paste0(saveDir,"//hclustGapTest_pam50_short_2015-05-15.RData.gzip")
load(paste0(saveDir,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
#35/50 is 70%
fractFeatIntersectThresh=.69
#we know these all share at least 35 genes.
numFeatIntersectThresh=34
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=3
#not used here
maxNullFractSize=.2
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"


#full features, pearson
edgeMethod <- "pearson"

clustSampleIndexList <- hclustGapTest_pam50_short$clustSampleIndexList
clustFeatureIndexList <- hclustGapTest_pam50_short$clustFeatureIndexList

source("/home/ywrfc09CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")

pam50ShortHclustGapTest_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                                edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                                sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                                outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                                numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                                clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


save(pam50ShortHclustGapTest_pearson_meanMatrix,file=paste0("/home/kplaney/breast_analysis/adjMatrices_pam50ShortHclustGapTest_pearson_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")

###pam50Full
load(paste0(saveDir,"//gapTesthclust_pam50Full_2015-05-15.RData.gzip")
load(paste0(saveDir,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
#35/50 is 70%
fractFeatIntersectThresh=.69
#we know these all share at least 35 genes.
numFeatIntersectThresh=34
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=3
#not used here
maxNullFractSize=.2
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"


#full features, pearson
edgeMethod <- "pearson"

clustSampleIndexList <- hclustGapTest_pam50_full_Nstart25$clustSampleIndexList
clustFeatureIndexList <- hclustGapTest_pam50_full_Nstart25$clustFeatureIndexList

source("/home/ywrfc09CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")

pam50FullHclustGapTest_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                               edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                               sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                               outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                               numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                               clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


save(pam50FullHclustGapTest_pearson_meanMatrix,file=paste0("/home/kplaney/breast_analysis/adjMatrices_pam50FullHclustGapTest_pearson_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")

