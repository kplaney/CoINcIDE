#/home/kplaney/breast_analysis_withTop20Genes/updatedSurvivalAnalyses



###pam50 full centroids
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis_withTop20Genes/updatedSurvivalAnalyses"
load("/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData")
pam50Genes <- centroidMatrix[,1]
metaFeatures <- list(finalFeatures=pam50Genes,datasetListIndicesToRemove=NULL)

load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/pam50Full_centroidCluster.RData.gzip"))
clusterCoINcIDE_output =  pam50Full_centroidCluster
#need to add fake PACR variable
clusterCoINcIDE_output$bestK_PACR <- unlist(sapply(clusterCoINcIDE_output$clustSampleIndexList,FUN=function(unit){
  
  length(unit)
  
})) 

names(clusterCoINcIDE_output$bestK_PACR) <- names(esets)
clusterCoINcIDE_output$clustSampleIndexList_PACR <- clusterCoINcIDE_output$clustSampleIndexList
clusterCoINcIDE_output$clustFeatureIndexList_PACR <- clusterCoINcIDE_output$clustFeatureIndexList
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/pam50Full_centroidCluster_pearson_meanMatrix_2015-05-03.RData.gzip"))
CoINcIDE_output = pam50Full_centroidCluster_pearson_meanMatrix
experimentName <- "breast_pam50FullCentroid_features"
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
                                                            meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                                            clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                            commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                            survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                            CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                            ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)


save(breast_pam50FullCentroidsOut,file="/home/kplaney/breast_analysis_withTop20Genes/breast_pam50FullCentroid_features_2015-05-07//breast_pam50FullCentroidsOut.RData.gzip",compress="gzip")

###all other metrics (mean matrix, spearman, etc. returned exactly the same meta-cluster.)

#####pam50 full kmeans
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis_withTop20Genes/updatedSurvivalAnalyses"
load("/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData")
pam50Genes <- centroidMatrix[,1]
metaFeatures <- list(finalFeatures=pam50Genes,datasetListIndicesToRemove=NULL)
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/kmeansConsensuspam50_full_Nstart1pItem9_pam50FullFeatures_2015-04-28.RData.gzip"))

clusterCoINcIDE_output =  kmeansConsensuspam50_full_Nstart1pItem9
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
#from the testing K runs:
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/adjMatrices_kmeansConsensuspam50_full_Nstart1pItem9.RData.gzip"))
CoINcIDE_output = pam50Full_pearson_meanMatrix
experimentName <- "pam50Full_meanMatrixPearson"
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

breast_pam50FullOut <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                   meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                                   clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                   commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                   survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                   CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                   ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)

save(breast_pam50FullOut,file="/home/kplaney/breast_analysis_withTop20Genes/breast_pam50Full_features_2015-05-08//breast_pam50FullOut.RData.gzip",compress="gzip")


####pam50 short kmeans
load("/home/kplaney/breast_analysis_withTop20Genes/pam50Short_genes.RData")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis_withTop20Genes/updatedSurvivalAnalyses"
metaFeatures <- list(finalFeatures=pam50Short,datasetListIndicesToRemove=NULL)
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/kmeansConsensuspam50_short_Nstart1pItem9_pam50ShortFeatures_2015-04-28.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensuspam50_short_Nstart1pItem9
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
#from the testing K runs:
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/adjMatrices_pam50Short_pearson_meanMatrix_2015-04-29.RData.gzip"))
CoINcIDE_output = pam50Short_pearson_meanMatrix
experimentName <- "breast_pam50Short_features"
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

breast_pam50ShortOut <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                    meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                                    clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                    commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                    survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                    CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                    ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)

save(breast_pam50ShortOut,file="/home/kplaney/breast_analysis_withTop20Genes/breast_pam50Short_features_",Sys.Date(),"/breast_pam50ShortOut.RData.gzip",compress="gzip")


###200
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis_withTop20Genes/updatedSurvivalAnalyses"
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/metaFeatures_200.RData.gzip"))
metaFeatures =metaFeatures
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/curatedbreastData_kmeansConsensus_nstart1pItem9200Features_2015-05-04.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensus
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/adjMatrices_200F_pearson_meanMatrix_2015-05-05.RData.gzip"))
CoINcIDE_output = breast200F_pearson_meanMatrix
experimentName <- "200_pearson_meanMatrix"
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

breast_200Out <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                             meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                             clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                             commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,
                                             eset_featureDataFieldName=eset_featureDataFieldName,
                                             survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                             CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                             ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)



save(breast_200Out,file="/home/kplaney/breast_analysis_withTop20Genes/breast_200_features_",Sys.Date(),"/breast_200_features_analysis_withTop20GenesOut.RData.gzip",compress="gzip")

####500
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis_withTop20Genes/updatedSurvivalAnalyses"
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/metaFeatures_500.RData.gzip"))
metaFeatures =metaFeatures
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/curatedbreastData_kmeansConsensus_nstart1pItem9500Features_2015-05-04.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensus
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/adjMatrices_500F_pearson_meanMatrix_2015-05-05.RData.gzip"))
CoINcIDE_output = breast500F_pearson_meanMatrix
experimentName <- "breast_500_features"
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

breast_500Out <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                             meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                             clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                             commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,
                                             eset_featureDataFieldName=eset_featureDataFieldName,
                                             survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                             CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                             ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)



save(breast_500Out,file="/home/kplaney/breast_analysis_withTop20Genes/breast_500_features_",Sys.Date(),"/breast_500_features_analysisOut.RData.gzip",compress="gzip")


####1000
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis_withTop20Genes/updatedSurvivalAnalyses"
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/metaFeatures_1000.RData.gzip"))
metaFeatures =metaFeatures
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/curatedbreastData_kmeansConsensus_nstart1pItem91000Features_2015-05-05.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensus
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/adjMatrices_1000F_pearson_meanMatrix_2015-05-05.RData.gzip"))
CoINcIDE_output = breast1000F_pearson_meanMatrix
experimentName <- "1000_pearson_meanMatrix"
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

breast_1000Out <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                              meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                              clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                              commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,
                                              eset_featureDataFieldName=eset_featureDataFieldName,
                                              survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                              CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                              ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)



save(breast_1000Out,file="/home/kplaney/breast_analysis_withTop20Genes/breast_1000_features_",Sys.Date(),"/breast_1000_features_analysisOut.RData.gzip",compress="gzip")

###2000
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list

saveDir <- "/home/kplaney/breast_analysis_withTop20Genes/updatedSurvivalAnalyses"
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/metaFeatures_2000.RData.gzip"))
metaFeatures =metaFeatures
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/curatedbreastData_kmeansConsensus_nstart1pItem92000Features_2015-05-05.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensus
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/adjMatrices_2000F_pearson_meanMatrix_2015-05-05.RData.gzip"))
CoINcIDE_output = breast2000F_pearson_meanMatrix
experimentName <- "2000_pearson_meanMatrix"

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

breast_2000Out <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                              meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                              clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                              commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,
                                              eset_featureDataFieldName=eset_featureDataFieldName,
                                              survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                              CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                              ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)



save(breast_2000Out,file="/home/kplaney/breast_analysis_withTop20Genes/breast_2000_features_",Sys.Date(),"/breast_2000_features_analysisOut.RData.gzip",compress="gzip")



########breast spearman analyses

#####pam50 full kmeans
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis_withTop20Genes/updatedSurvivalAnalyses"
load("/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData")
pam50Genes <- centroidMatrix[,1]
metaFeatures <- list(finalFeatures=pam50Genes,datasetListIndicesToRemove=NULL)
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/kmeansConsensuspam50_full_Nstart1pItem9_pam50FullFeatures_2015-04-28.RData.gzip"))

clusterCoINcIDE_output =  kmeansConsensuspam50_full_Nstart1pItem9
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
#from the testing K runs:
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/adjMatrices_pam50Full_spearman_meanMatrix_2015-04-29.RData.gzip"))
CoINcIDE_output = pam50Full_spearman_meanMatrix
experimentName <- "breast_pam50Full_spearman"
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

breast_spearman_pam50FullOut <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                            meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                                            clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                            commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                            survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                            CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                            ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)

save(breast_spearman_pam50FullOut,file="/home/kplaney/breast_analysis_withTop20Genes/breast_spearman_pam50Full_features_2015-05-08//breast_spearman_pam50FullOut.RData.gzip",compress="gzip")


####pam50 short kmeans
load("/home/kplaney/breast_analysis_withTop20Genes/pam50Short_genes.RData")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis_withTop20Genes/updatedSurvivalAnalyses"
metaFeatures <- list(finalFeatures=pam50Short,datasetListIndicesToRemove=NULL)
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/kmeansConsensuspam50_short_Nstart1pItem9_pam50ShortFeatures_2015-04-28.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensuspam50_short_Nstart1pItem9
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
#from the testing K runs:
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/adjMatrices_pam50Short_spearman_meanMatrix_2015-04-29.RData.gzip"))
CoINcIDE_output = pam50Short_spearman_meanMatrix
experimentName <- "pam50Short_spearman_meanMatrix"
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

breast_spearman_pam50ShortOut <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                             meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                                             clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                             commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                             survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                             CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                             ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)

save(breast_spearman_pam50ShortOut,file="/home/kplaney/breast_analysis_withTop20Genes/breast_pam50Short_spearman_",Sys.Date()),"/breast_spearman_pam50ShortOut.RData.gzip",compress="gzip")


###200
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis_withTop20Genes/updatedSurvivalAnalyses"
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/metaFeatures_200.RData.gzip"))
metaFeatures =metaFeatures
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/curatedbreastData_kmeansConsensus_nstart1pItem9200Features_2015-05-04.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensus
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/adjMatrices_200F_spearman_meanMatrix_2015-05-15.RData.gzip"))
CoINcIDE_output = breast200F_spearman_meanMatrix
experimentName <- "200_spearman_meanMatrix"
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

breast_spearman_200Out <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                      meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                                      clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                      commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,
                                                      eset_featureDataFieldName=eset_featureDataFieldName,
                                                      survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                      CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                      ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)



save(breast_spearman_200Out,file="/home/kplaney/breast_analysis_withTop20Genes/breast_200_spearman_",Sys.Date(),"/breast_200_features_analysisOut.RData.gzip",compress="gzip")

####500
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis_withTop20Genes/updatedSurvivalAnalyses"
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/metaFeatures_500.RData.gzip"))
metaFeatures =metaFeatures
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/curatedbreastData_kmeansConsensus_nstart1pItem9500Features_2015-05-04.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensus
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/adjMatrices_500F_spearman_meanMatrix_2015-05-15.RData.gzip"))
CoINcIDE_output = breast500F_spearman_meanMatrix
experimentName <- "500_spearman_meanMatrix"
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

breast_spearman_500Out <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                      meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                                      clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                      commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,
                                                      eset_featureDataFieldName=eset_featureDataFieldName,
                                                      survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                      CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                      ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)



save(breast_spearman_500Out,file="/home/kplaney/breast_analysis_withTop20Genes/breast_500_spearman_",Sys.Date(),"/breast_500_features_analysisOut.RData.gzip",compress="gzip")


####1000
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis_withTop20Genes/updatedSurvivalAnalyses"
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/metaFeatures_1000.RData.gzip"))
metaFeatures =metaFeatures
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/curatedbreastData_kmeansConsensus_nstart1pItem91000Features_2015-05-05.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensus
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/adjMatrices_1000F_spearman_meanMatrix_2015-05-15.RData.gzip"))
CoINcIDE_output = breast1000F_spearman_meanMatrix
experimentName <- "1000_spearman_meanMatrix"
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

breast_spearman_1000Out <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                       meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                                       clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                       commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,
                                                       eset_featureDataFieldName=eset_featureDataFieldName,
                                                       survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                       CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                       ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)



save(breast_spearman_1000Out,file="/home/kplaney/breast_analysis_withTop20Genes/breast_1000_spearman_",Sys.Date(),"/breast_1000_features_analysisOut.RData.gzip",compress="gzip")

###2000
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list

saveDir <- "/home/kplaney/breast_analysis_withTop20Genes/updatedSurvivalAnalyses"
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/metaFeatures_2000.RData.gzip"))
metaFeatures =metaFeatures
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/curatedbreastData_kmeansConsensus_nstart1pItem92000Features_2015-05-05.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensus
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/adjMatrices_2000F_spearman_meanMatrix_2015-05-15.RData.gzip"))
CoINcIDE_output = breast2000F_spearman_meanMatrix
experimentName <- "2000_spearman_meanMatrix"

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

breast_spearman_2000Out <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                       meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                                       clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                       commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,
                                                       eset_featureDataFieldName=eset_featureDataFieldName,
                                                       survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                       CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                       ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)



save(breast_spearman_2000Out,file="/home/kplaney/breast_analysis_withTop20Genes/breast_2000_spearman_",Sys.Date(),"/breast_2000_features_analysisOut.RData.gzip",compress="gzip")

####pam50 with centroid:

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis_withTop20Genes/updatedSurvivalAnalyses/"
load("/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData")
pam50Genes <- centroidMatrix[,1]
metaFeatures <- list(finalFeatures=pam50Genes,datasetListIndicesToRemove=NULL)
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/kmeansConsensuspam50_full_Nstart1pItem9_pam50FullFeatures_2015-04-28.RData.gzip"))

clusterCoINcIDE_output =  kmeansConsensuspam50_full_Nstart1pItem9
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
#from the testing K runs:
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/adjMatrices_pam50Full_pearson_centroid_2015-05-18.RData.gzip"))
CoINcIDE_output = pam50Full_pearson_centroid
#use the .3 instead of the .1 version:
CoINcIDE_output$pvalueMatrix <- CoINcIDE_output$pvalueMatrix3
experimentName <- "pam50Full_pearson_centroid"
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

breast_centroid_pam50FullOut <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                            meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                                            clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                            commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                            survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                            CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                            ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)

save(breast_centroid_pam50FullOut,file=paste0("/home/kplaney/breast_analysis_withTop20Genes/",experimentName,"_",Sys.Date(),
                                              "/",experimentName,"Out.RData.gzip"),compress="gzip")


####pam50 short kmeans
load("/home/kplaney/breast_analysis_withTop20Genes/pam50Short_genes.RData")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis_withTop20Genes/updatedSurvivalAnalyses/"
metaFeatures <- list(finalFeatures=pam50Short,datasetListIndicesToRemove=NULL)
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/kmeansConsensuspam50_short_Nstart1pItem9_pam50ShortFeatures_2015-04-28.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensuspam50_short_Nstart1pItem9
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
#from the testing K runs:
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/adjMatrices_pam50Short_pearson_centroid_2015-05-18.RData.gzip"))
CoINcIDE_output = pam50Short_pearson_centroid
#use the .3 instead of the .1 version:
CoINcIDE_output$pvalueMatrix <- CoINcIDE_output$pvalueMatrix3
experimentName <- "pam50Short_pearson_centroid"
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

breast_centroid_pam50ShortOut <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                             meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                                             clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                             commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                             survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                             CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                             ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)

save(breast_centroid_pam50ShortOut,file=paste0("/home/kplaney/breast_analysis_withTop20Genes/",experimentName,"_",Sys.Date(),
                                               "/",experimentName,"Out.RData.gzip"),compress="gzip")

####more centroid
###200
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis_withTop20Genes/updatedSurvivalAnalyses/"
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/metaFeatures_200.RData.gzip"))
metaFeatures =metaFeatures
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/curatedbreastData_kmeansConsensus_nstart1pItem9200Features_2015-05-04.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensus
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/adjMatrices_200F_pearson_centroid_2015-05-18.RData.gzip"))
CoINcIDE_output = breast200F_pearson_centroid
#use the .3 instead of the .1 version:
CoINcIDE_output$pvalueMatrix <- CoINcIDE_output$pvalueMatrix3
experimentName <- "breast_200_centroid"
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

breast_centroid_200Out <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                      meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                                      clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                      commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,
                                                      eset_featureDataFieldName=eset_featureDataFieldName,
                                                      survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                      CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                      ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)



save(breast_centroid_200Out,file=paste0("/home/kplaney/breast_analysis_withTop20Genes/",experimentName,"_",Sys.Date(),
                                        "/",experimentName,"Out.RData.gzip"),compress="gzip")
####500
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis_withTop20Genes/updatedSurvivalAnalyses/"
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/metaFeatures_500.RData.gzip"))
metaFeatures =metaFeatures
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/curatedbreastData_kmeansConsensus_nstart1pItem9500Features_2015-05-04.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensus
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/adjMatrices_500F_pearson_centroid_2015-05-18.RData.gzip"))
CoINcIDE_output = breast500F_pearson_centroid
#use the .3 instead of the .1 version:
CoINcIDE_output$pvalueMatrix <- CoINcIDE_output$pvalueMatrix3
experimentName <- "breast_500_centroid"
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

breast_centroid_500Out <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                      meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                                      clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                      commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,
                                                      eset_featureDataFieldName=eset_featureDataFieldName,
                                                      survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                      CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                      ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)



save(breast_centroid_500Out,file=paste0("/home/kplaney/breast_analysis_withTop20Genes/",experimentName,"_",Sys.Date(),
                                        "/",experimentName,"Out.RData.gzip"),compress="gzip")

####1000
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis_withTop20Genes/updatedSurvivalAnalyses/"
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/metaFeatures_1000.RData.gzip"))
metaFeatures =metaFeatures
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/curatedbreastData_kmeansConsensus_nstart1pItem91000Features_2015-05-05.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensus
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/adjMatrices_1000F_pearson_centroid_2015-05-18.RData.gzip"))
CoINcIDE_output = breast1000F_pearson_centroid
#use the .3 instead of the .1 version:
CoINcIDE_output$pvalueMatrix <- CoINcIDE_output$pvalueMatrix3
experimentName <- "breast_1000_centroid"
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

breast_centroid_1000Out <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                       meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                                       clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                       commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,
                                                       eset_featureDataFieldName=eset_featureDataFieldName,
                                                       survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                       CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                       ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)



save(breast_centroid_1000Out,file=paste0("/home/kplaney/breast_analysis_withTop20Genes/",experimentName,"_",Sys.Date(),
                                         "/",experimentName,"Out.RData.gzip"),compress="gzip")
###2000
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list

saveDir <- "/home/kplaney/breast_analysis_withTop20Genes/updatedSurvivalAnalyses/"
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/metaFeatures_2000.RData.gzip"))
metaFeatures =metaFeatures
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/curatedbreastData_kmeansConsensus_nstart1pItem92000Features_2015-05-05.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensus
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0("/home/kplaney/breast_analysis_withTop20Genes/","/adjMatrices_2000F_pearson_meanMatrix_centroid_2015-05-18.RData.gzip"))
CoINcIDE_output = breast2000F_pearson_centroid
#use the .3 instead of the .1 version:
CoINcIDE_output$pvalueMatrix <- CoINcIDE_output$pvalueMatrix3
experimentName <- "breast_2000_centroid"

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

breast_centroid_2000Out <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                       meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                                       clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                       commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,
                                                       eset_featureDataFieldName=eset_featureDataFieldName,
                                                       survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                       CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                       ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)



save(breast_centroid_2000Out,file=paste0("/home/kplaney/breast_analysis_withTop20Genes/",experimentName,"_",Sys.Date(),
                                         "/",experimentName,"Out.RData.gzip"),compress="gzip")

