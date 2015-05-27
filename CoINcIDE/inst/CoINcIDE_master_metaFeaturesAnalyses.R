

######ovarian 200
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
load("/home/kplaney/ovarian_analysis/metaFeatures_200.RData.gzip")
metaFeatures =metaFeatures
load("/home/kplaney/ovarian_analysis//curatedOvarianData_kmeansConsensus_nstart1_200Features_2015-05-19.RData.gzip")
clusterCoINcIDE_output =  kmeansConsensus
load("/home/kplaney/ovarian_analysis/esets_proc_TCGAcombat.RData.gzip")
esets = esets
#dataMatrixList
load("/home/kplaney/ovarian_analysis/adjMatrices_200F_pearson_meanMatrix_updated_2015-05-19RData.gzip")
CoINcIDE_output = ov_200F_pearson_meanMatrix
experimentName <- "ovarian_200_features"
eset_featureDataFieldName="gene"
networkColors = "Set2"
saveDir <- "/home/kplaney/ovarian_analysis/"
outcomesVarBinary="vital_status"
outcomesVarCont = "days_to_death"
ovarian <- TRUE
fisherTestVariables <- c("histological_type","tumorstage","recurrence_status","grade","age_at_initial_pathologic_diagnosis")
fisherTestVariableLegendNames <- c("hist\ntype","tumor\nstage","recurrence","tumor\ngrade","age")
fisherTestVariableTitleNames <- c("histological type","tumor stage", "recurrence status","tumor grade","age at diagnosis")

#breast: must also input dataMatrixList
#brewPal = c("Set3","Paired","Spectral","BrBG","PiYG","RdYlGn","RdYlBu","RdBu","PiYG","Set2"),
ov_200Out <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                        meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                        clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                        commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                        survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                        CutoffPointYears=5, eset_uniquePatientID="unique_patient_ID", fisherTestVariables = fisherTestVariables,fisherTestVariableLegendNames=fisherTestVariableLegendNames,
                                        fisherTestVariableTitleNames=fisherTestVariableTitleNames)

  
save(ov_200Out,file=paste0(saveDir,"/",experimentName,"_",Sys.Date(),"/ov200Out.RData.gzip"),compress="gzip")
######ovarian 500
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
load("/home/kplaney/ovarian_analysis/metaFeatures_500.RData.gzip")
metaFeatures =metaFeatures
load("/home/kplaney/ovarian_analysis/curatedOvarianData_kmeansConsensus_nstart1_500Features_2015-05-19.RData.gzip")
clusterCoINcIDE_output =  kmeansConsensus
load("/home/kplaney/ovarian_analysis/esets_proc_TCGAcombat.RData.gzip")
esets = esets
#dataMatrixList
load("/home/kplaney/ovarian_analysis/adjMatrices_500F_pearson_meanMatrix_updated2015-05-20RData.gzip")
CoINcIDE_output = ov_500F_pearson_meanMatrix
experimentName <- "ovarian_500_features"
fisherTestVariables <- c("histological_type","tumorstage","recurrence_status","grade","age_at_initial_pathologic_diagnosis")
fisherTestVariableLegendNames <- c("hist\ntype","tumor\nstage","recurrence","tumor\ngrade","age")
fisherTestVariableTitleNames <- c("histological type","tumor stage", "recurrence status","tumor grade","age at diagnosis")
eset_featureDataFieldName="gene"
networkColors = "Set2"
saveDir <- "/home/kplaney/ovarian_analysis/"
outcomesVarBinary="vital_status"
outcomesVarCont = "days_to_death"
ovarian <- TRUE
eset_uniquePatientID="unique_patient_ID"
#breast: must also input dataMatrixList
#brewPal = c("Set3","Paired","Spectral","BrBG","PiYG","RdYlGn","RdYlBu","RdBu","PiYG","Set2"),
ov_500Out <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                         meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                         clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                         commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                         survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                         CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                         fisherTestVariableTitleNames=fisherTestVariableTitleNames,fisherTestVariableLegendNames =fisherTestVariableLegendNames )

save(ov_500Out,file=paste0(saveDir,"/",experimentName,"_",Sys.Date(),"/ov500Out.RData.gzip"),compress="gzip")

##########ovarian 1000

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
load("/home/kplaney/ovarian_analysis/metaFeatures_1000.RData.gzip")
metaFeatures =metaFeatures
load("/home/kplaney/ovarian_analysis/curatedOvarianData_kmeansConsensus_nstart1_1000Features_2015-04-29.RData.gzip")
clusterCoINcIDE_output =  kmeansConsensus
load("/home/kplaney/ovarian_analysis/esets_proc_TCGAcombat.RData.gzip")
esets = esets
#dataMatrixList
load("/home/kplaney/ovarian_analysis/adjMatrices_1000F_pearson_meanMatrix_2015-04-29RData.gzip")
CoINcIDE_output = ov_1000F_pearson_meanMatrix
experimentName <- "ovarian_1000_features"
fisherTestVariables <- c("histological_type","tumorstage","recurrence_status","grade","age_at_initial_pathologic_diagnosis")
fisherTestVariableLegendNames <- c("hist\ntype","tumor\nstage","recurrence","tumor\ngrade","age")
fisherTestVariableTitleNames <- c("histological type","tumor stage", "recurrence status","tumor grade","age at diagnosis")
eset_featureDataFieldName="gene"
networkColors = "Set2"
saveDir <- "/home/kplaney/ovarian_analysis/"
outcomesVarBinary="vital_status"
outcomesVarCont = "days_to_death"
ovarian <- TRUE
eset_uniquePatientID="unique_patient_ID"
#breast: must also input dataMatrixList
#brewPal = c("Set3","Paired","Spectral","BrBG","PiYG","RdYlGn","RdYlBu","RdBu","PiYG","Set2"),
ov_1000Out <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                         meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                         clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                         commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                         survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                         CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                         fisherTestVariableTitleNames=fisherTestVariableTitleNames,fisherTestVariableLegendNames=fisherTestVariableLegendNames)




save(ov_1000Out,file=paste0(saveDir,"/",experimentName,"_",Sys.Date(),"/ov1000Out.RData.gzip"),compress="gzip")
###########


source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
load("/home/kplaney/ovarian_analysis/metaFeatures_2000.RData.gzip")
metaFeatures =metaFeatures
load("/home/kplaney/ovarian_analysis/curatedOvarianData_kmeansConsensus_nstart1_2000Features_2015-04-30.RData.gzip")
clusterCoINcIDE_output =  kmeansConsensus
load("/home/kplaney/ovarian_analysis/esets_proc_TCGAcombat.RData.gzip")
esets = esets
#dataMatrixList
load("/home/kplaney/ovarian_analysis/adjMatrices_2000F_pearson_meanMatrix_updated2015-05-01RData.gzip")
CoINcIDE_output = ov_2000F_pearson_meanMatrix
experimentName <- "ovarian_2000_features"
fisherTestVariables <- c("histological_type","tumorstage","recurrence_status","grade","age_at_initial_pathologic_diagnosis")
fisherTestVariableLegendNames <- c("hist\ntype","tumor\nstage","recurrence","tumor\ngrade","age")
fisherTestVariableTitleNames <- c("histological type","tumor stage", "recurrence status","tumor grade","age at diagnosis")

eset_featureDataFieldName="gene"
networkColors = "Set2"
saveDir <- "/home/kplaney/ovarian_analysis/"
outcomesVarBinary="vital_status"
outcomesVarCont = "days_to_death"
ovarian <- TRUE
eset_uniquePatientID="unique_patient_ID"
#breast: must also input dataMatrixList
#brewPal = c("Set3","Paired","Spectral","BrBG","PiYG","RdYlGn","RdYlBu","RdBu","PiYG","Set2"),
ov_2000Out <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                          meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                          clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                          commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                          survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                          CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                          fisherTestVariableTitleNames =fisherTestVariableTitleNames ,fisherTestVariableLegendNames=fisherTestVariableLegendNames)


save(ov_2000Out,file=paste0(saveDir,"/",experimentName,"_",Sys.Date(),"/ov2000Out.RData.gzip"),compress="gzip")
##############
##############
####breast
###pam50 full centroids
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
load("/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData")
pam50Genes <- centroidMatrix[,1]
metaFeatures <- list(finalFeatures=pam50Genes,datasetListIndicesToRemove=NULL)

load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0(saveDir,"/pam50Full_centroidCluster.RData.gzip"))
clusterCoINcIDE_output =  pam50Full_centroidCluster
#need to add fake PACR variable
clusterCoINcIDE_output$bestK_PACR <- unlist(sapply(clusterCoINcIDE_output$clustSampleIndexList,FUN=function(unit){
  
  length(unit)
  
})) 

names(clusterCoINcIDE_output$bestK_PACR) <- names(esets)
clusterCoINcIDE_output$clustSampleIndexList_PACR <- clusterCoINcIDE_output$clustSampleIndexList
clusterCoINcIDE_output$clustFeatureIndexList_PACR <- clusterCoINcIDE_output$clustFeatureIndexList
load(paste0(saveDir,"/pam50Full_centroidCluster_pearson_meanMatrix_2015-05-03.RData.gzip"))
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


save(breast_pam50FullCentroidsOut,file="/home/kplaney/breast_analysis/breast_pam50FullCentroid_features_2015-05-07//breast_pam50FullCentroidsOut.RData.gzip",compress="gzip")



###pam50 centroid clustering pearson, centroid sim
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
load("/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData")
pam50Genes <- centroidMatrix[,1]
metaFeatures <- list(finalFeatures=pam50Genes,datasetListIndicesToRemove=NULL)

load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0(saveDir,"/adjMatrices_pam50Full_centroidCluster_pearson_centroid_2015-05-20.RData.gzip"))
clusterCoINcIDE_output =  pam50Full_centroidCluster
#need to add fake PACR variable
clusterCoINcIDE_output$bestK_PACR <- unlist(sapply(clusterCoINcIDE_output$clustSampleIndexList,FUN=function(unit){
  
  length(unit)
  
})) 

names(clusterCoINcIDE_output$bestK_PACR) <- names(esets)
clusterCoINcIDE_output$clustSampleIndexList_PACR <- clusterCoINcIDE_output$clustSampleIndexList
clusterCoINcIDE_output$clustFeatureIndexList_PACR <- clusterCoINcIDE_output$clustFeatureIndexList
load(paste0(saveDir,"/pam50Full_centroidCluster_pearson_meanMatrix_2015-05-03.RData.gzip"))
CoINcIDE_output = pam50Full_centroidCluster_pearson_meanMatrix
experimentName <- "pam50FullCentroid_pearson_centroid"
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

breast_pam50FullCentroidsOut_centroid <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                            meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                                            clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                            commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                            survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                            CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                            ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)


#save(breast_pam50FullCentroidsOut_centroid,file=paste0(saveDir,"/",experimentName,/


###pam50 centroid clustering spearman sim
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
load("/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData")
pam50Genes <- centroidMatrix[,1]
metaFeatures <- list(finalFeatures=pam50Genes,datasetListIndicesToRemove=NULL)

load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0(saveDir,"/adjMatrices_pam50Full_centroidCluster_spearman_meanMatrix_2015-05-20.RData.gzip"))
clusterCoINcIDE_output =  pam50Full_centroidCluster
#need to add fake PACR variable
clusterCoINcIDE_output$bestK_PACR <- unlist(sapply(clusterCoINcIDE_output$clustSampleIndexList,FUN=function(unit){
  
  length(unit)
  
})) 

names(clusterCoINcIDE_output$bestK_PACR) <- names(esets)
clusterCoINcIDE_output$clustSampleIndexList_PACR <- clusterCoINcIDE_output$clustSampleIndexList
clusterCoINcIDE_output$clustFeatureIndexList_PACR <- clusterCoINcIDE_output$clustFeatureIndexList
load(paste0(saveDir,"/pam50Full_centroidCluster_pearson_meanMatrix_2015-05-03.RData.gzip"))
CoINcIDE_output = pam50Full_centroidCluster_pearson_meanMatrix
experimentName <- "pam50FullCentroid_spearman_meanMatrix"
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

breast_pam50FullCentroidsOut_spearmanMM <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                                     meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                                                     clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                                     commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                                     survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                                     CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                                     ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)


###pam50 spearman, centroid
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
load("/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData")
pam50Genes <- centroidMatrix[,1]
metaFeatures <- list(finalFeatures=pam50Genes,datasetListIndicesToRemove=NULL)

load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0(saveDir,"/adjMatrices_pam50Full_centroidCluster_spearman_centroid_2015-05-20.RData.gzip"))
clusterCoINcIDE_output =  pam50Full_centroidCluster
#need to add fake PACR variable
clusterCoINcIDE_output$bestK_PACR <- unlist(sapply(clusterCoINcIDE_output$clustSampleIndexList,FUN=function(unit){
  
  length(unit)
  
})) 

names(clusterCoINcIDE_output$bestK_PACR) <- names(esets)
clusterCoINcIDE_output$clustSampleIndexList_PACR <- clusterCoINcIDE_output$clustSampleIndexList
clusterCoINcIDE_output$clustFeatureIndexList_PACR <- clusterCoINcIDE_output$clustFeatureIndexList
load(paste0(saveDir,"/pam50Full_centroidCluster_pearson_meanMatrix_2015-05-03.RData.gzip"))
CoINcIDE_output = pam50Full_centroidCluster_pearson_meanMatrix
experimentName <- "pam50FullCentroid_spearman_centroid"
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

breast_pam50FullCentroidsOut_spearman_centroid <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                                       meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                                                       clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                                       commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                                       survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                                       CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                                       ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)




#####pam50 full kmeans
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
load("/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData")
pam50Genes <- centroidMatrix[,1]
metaFeatures <- list(finalFeatures=pam50Genes,datasetListIndicesToRemove=NULL)
load(paste0(saveDir,"/kmeansConsensuspam50_full_Nstart1pItem9_pam50FullFeatures_2015-04-28.RData.gzip"))

clusterCoINcIDE_output =  kmeansConsensuspam50_full_Nstart1pItem9
load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
#from the testing K runs:
load(paste0(saveDir,"/adjMatrices_kmeansConsensuspam50_full_Nstart1pItem9.RData.gzip"))
CoINcIDE_output = pam50Full_pearson_meanMatrix
experimentName <- "breast_pam50Full_features"
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

save(breast_pam50FullOut,file="/home/kplaney/breast_analysis/breast_pam50Full_features_2015-05-08//breast_pam50FullOut.RData.gzip",compress="gzip")


####pam50 short kmeans
load("/home/kplaney/breast_analysis/pam50Short_genes.RData")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
metaFeatures <- list(finalFeatures=pam50Short,datasetListIndicesToRemove=NULL)
load(paste0(saveDir,"/kmeansConsensuspam50_short_Nstart1pItem9_pam50ShortFeatures_2015-04-28.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensuspam50_short_Nstart1pItem9
load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
#from the testing K runs:
load(paste0(saveDir,"/adjMatrices_pam50Short_pearson_meanMatrix_2015-04-29.RData.gzip"))
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
                                                   meanEdgePairPvalueThresh = .05,indEdgePvalueThresh = .1, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                                   clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                   commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                   survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                   CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                   ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)

save(breast_pam50ShortOut,file="/home/kplaney/breast_analysis/breast_pam50Short_features_",Sys.Date(),"/breast_pam50ShortOut.RData.gzip",compress="gzip")

##250
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
load(paste0(saveDir,"/metaFeatures_250.RData.gzip"))
metaFeatures =metaFeatures
load(paste0(saveDir,"/curatedbreastData_kmeansConsensus_nstart1pItem9250Features_2015-05-20.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensus
load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0(saveDir,"/adjMatrices_250F_pearson_meanMatrix_2015-05-20.RData.gzip"))
CoINcIDE_output = breast250F_pearson_meanMatrix
experimentName <- "breast_250_features"
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

breast_250Out <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                             meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .5, maxTrueSimilThresh = Inf,
                                             clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                             commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,
                                             eset_featureDataFieldName=eset_featureDataFieldName,
                                             survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                             CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                             ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)



save(breast_250Out,file=paste0("/home/kplaney/breast_analysis/",experimentName,"_",Sys.Date(),"/breast_250_features_analysisOut.RData.gzip"),compress="gzip")

####500

###200
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
load(paste0(saveDir,"/metaFeatures_200.RData.gzip"))
metaFeatures =metaFeatures
load(paste0(saveDir,"/curatedbreastData_kmeansConsensus_nstart1pItem9200Features_2015-05-19.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensus
load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0(saveDir,"/adjMatrices_200F_pearson_meanMatrix_2015-05-19.RData.gzip"))
CoINcIDE_output = breast200F_pearson_meanMatrix
experimentName <- "breast_200_features"
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
                                          meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .5, maxTrueSimilThresh = Inf,
                                          clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                          commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,
                                          eset_featureDataFieldName=eset_featureDataFieldName,
                                          survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                          CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                          ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)



save(breast_200Out,file=paste0("/home/kplaney/breast_analysis/",experimentName,"_",Sys.Date(),"/breast_200_features_analysisOut.RData.gzip"),compress="gzip")

####500
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
load(paste0(saveDir,"/metaFeatures_500.RData.gzip"))
metaFeatures =metaFeatures
load(paste0(saveDir,"/curatedbreastData_kmeansConsensus_nstart1pItem9500Features_2015-05-19.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensus
load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0(saveDir,"/adjMatrices_500F_pearson_meanMatrix_2015-05-20.RData.gzip"))
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
     
     
     
     save(breast_500Out,file=paste0(saveDir,"/",experimentName,"_",Sys.Date(),"/breast_500_features_analysisOut.RData.gzip"),compress="gzip")
     
     
####1000
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
load(paste0(saveDir,"/metaFeatures_1000.RData.gzip"))
metaFeatures =metaFeatures
load(paste0(saveDir,"/curatedbreastData_kmeansConsensus_nstart1pItem91000Features_2015-05-20.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensus
load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0(saveDir,"/adjMatrices_1000F_pearson_meanMatrix_2015-05-20.RData.gzip"))
CoINcIDE_output = breast1000F_pearson_meanMatrix
experimentName <- "breast_1000_features"
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
     
     
     
     save(breast_1000Out,file="/home/kplaney/breast_analysis/breast_1000_features_",Sys.Date(),"/breast_1000_features_analysisOut.RData.gzip",compress="gzip")
     
###2000
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list

saveDir <- "/home/kplaney/breast_analysis/"
load(paste0(saveDir,"/metaFeatures_2000.RData.gzip"))
metaFeatures =metaFeatures
load(paste0(saveDir,"/curatedbreastData_kmeansConsensus_nstart1pItem92000Features_2015-05-20.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensus
load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0(saveDir,"/adjMatrices_2000F_pearson_meanMatrix_2015-05-27.RData.gzip"))
CoINcIDE_output = breast2000F_pearson_meanMatrix
experimentName <- "breast_2000_features"

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
     
     
     
     save(breast_2000Out,file="/home/kplaney/breast_analysis/breast_2000_features_",Sys.Date(),"/breast_2000_features_analysisOut.RData.gzip",compress="gzip")



########breast spearman analyses

#####pam50 full kmeans
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
load("/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData")
pam50Genes <- centroidMatrix[,1]
metaFeatures <- list(finalFeatures=pam50Genes,datasetListIndicesToRemove=NULL)
load(paste0(saveDir,"/kmeansConsensuspam50_full_Nstart1pItem9_pam50FullFeatures_2015-04-28.RData.gzip"))

clusterCoINcIDE_output =  kmeansConsensuspam50_full_Nstart1pItem9
load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
#from the testing K runs:
load(paste0(saveDir,"/adjMatrices_pam50Full_spearman_meanMatrix_2015-04-29.RData.gzip"))
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

save(breast_spearman_pam50FullOut,file="/home/kplaney/breast_analysis/breast_spearman_pam50Full_features_2015-05-08//breast_spearman_pam50FullOut.RData.gzip",compress="gzip")


####pam50 short kmeans
load("/home/kplaney/breast_analysis/pam50Short_genes.RData")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
metaFeatures <- list(finalFeatures=pam50Short,datasetListIndicesToRemove=NULL)
load(paste0(saveDir,"/kmeansConsensuspam50_short_Nstart1pItem9_pam50ShortFeatures_2015-04-28.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensuspam50_short_Nstart1pItem9
load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
#from the testing K runs:
load(paste0(saveDir,"/adjMatrices_pam50Short_spearman_meanMatrix_2015-04-29.RData.gzip"))
CoINcIDE_output = pam50Short_spearman_meanMatrix
experimentName <- "breast_pam50Short_spearman"
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

save(breast_spearman_pam50ShortOut,file="/home/kplaney/breast_analysis/breast_pam50Short_spearman_",Sys.Date()),"/breast_spearman_pam50ShortOut.RData.gzip",compress="gzip")


###200
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
load(paste0(saveDir,"/metaFeatures_200.RData.gzip"))
metaFeatures =metaFeatures
load(paste0(saveDir,"/curatedbreastData_kmeansConsensus_nstart1pItem9200Features_2015-05-19.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensus
load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0(saveDir,"/adjMatrices_200F_spearman_meanMatrix_2015-05-15.RData.gzip"))
CoINcIDE_output = breast200F_spearman_meanMatrix
experimentName <- "breast_200_spearman"
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



save(breast_spearman_200Out,file="/home/kplaney/breast_analysis/breast_200_spearman_",Sys.Date(),"/breast_200_features_analysisOut.RData.gzip",compress="gzip")

####500
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
load(paste0(saveDir,"/metaFeatures_500.RData.gzip"))
metaFeatures =metaFeatures
load(paste0(saveDir,"/curatedbreastData_kmeansConsensus_nstart1pItem9500Features_2015-05-04.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensus
load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0(saveDir,"/adjMatrices_500F_spearman_meanMatrix_2015-05-15.RData.gzip"))
CoINcIDE_output = breast500F_spearman_meanMatrix
experimentName <- "breast_500_spearman"
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



save(breast_spearman_500Out,file="/home/kplaney/breast_analysis/breast_500_spearman_",Sys.Date(),"/breast_500_features_analysisOut.RData.gzip",compress="gzip")


####1000
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
load(paste0(saveDir,"/metaFeatures_1000.RData.gzip"))
metaFeatures =metaFeatures
load(paste0(saveDir,"/curatedbreastData_kmeansConsensus_nstart1pItem91000Features_2015-05-05.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensus
load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0(saveDir,"/adjMatrices_1000F_spearman_meanMatrix_2015-05-15.RData.gzip"))
CoINcIDE_output = breast1000F_spearman_meanMatrix
experimentName <- "breast_1000_spearman"
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



save(breast_spearman_1000Out,file="/home/kplaney/breast_analysis/breast_1000_spearman_",Sys.Date(),"/breast_1000_features_analysisOut.RData.gzip",compress="gzip")

###2000
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list

saveDir <- "/home/kplaney/breast_analysis/"
load(paste0(saveDir,"/metaFeatures_2000.RData.gzip"))
metaFeatures =metaFeatures
load(paste0(saveDir,"/curatedbreastData_kmeansConsensus_nstart1pItem92000Features_2015-05-05.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensus
load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0(saveDir,"/adjMatrices_2000F_spearman_meanMatrix_2015-05-15.RData.gzip"))
CoINcIDE_output = breast2000F_spearman_meanMatrix
experimentName <- "breast_2000_spearman"

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



save(breast_spearman_2000Out,file="/home/kplaney/breast_analysis/breast_2000_spearman_",Sys.Date(),"/breast_2000_features_analysisOut.RData.gzip",compress="gzip")

#############
#TO DO: run all of this.need to edit!
##gap test, k-means
####pam50Short, full gap test (k-means)

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
load("/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData")
pam50Genes <- centroidMatrix[,1]
metaFeatures <- list(finalFeatures=pam50Genes,datasetListIndicesToRemove=NULL)
load(paste0(saveDir,"/gapTestKmeans_pam50Full_nstart12015-05-04.RData.gzip"))

load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
#nstart is a misnomer - k-means only uses nstart.
clusterCoINcIDE_output =  hclustGapTest_pam50_full_Nstart25
#need to add fake PACR variable
clusterCoINcIDE_output$bestK_PACR <- unlist(sapply(clusterCoINcIDE_output$clustSampleIndexList,FUN=function(unit){
  
  length(unit)
  
})) 

datasetNames <- names(esets)
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  datasetNames <- datasetNames[-metaFeatures$datasetListIndicesToRemove]
}
names(clusterCoINcIDE_output$bestK_PACR) <- datasetNames
clusterCoINcIDE_output$clustSampleIndexList_PACR <- clusterCoINcIDE_output$clustSampleIndexList
clusterCoINcIDE_output$clustFeatureIndexList_PACR <- clusterCoINcIDE_output$clustFeatureIndexList
names(clusterCoINcIDE_output$clustSampleIndexList_PACR)<- datasetNames
names(clusterCoINcIDE_output$clustFeatureIndexList_PACR )<- datasetNames
#load CoINcIDE output
load(paste0(saveDir,"/adjMatrices_pam50FullGapTest_pearson_meanMatrix_2015-05-15.RData.gzip"))
CoINcIDE_output = pam50FullHclustGapTest_pearson_meanMatrix
experimentName <- "breast_pam50Full_hclustGap"
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

breast_hclustGap_pam50FullOut <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                             meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                                             clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                             commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                             survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                             CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                             ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)

save(breast_hclustGap_pam50FullOut,file=paste0("/home/kplaney/breast_analysis/",experimentName,"_",Sys.Date(),"/breast_hclustGap_pam50FullOut.RData.gzip"),compress="gzip")


###pam50 short
load("/home/kplaney/breast_analysis/pam50Short_genes.RData")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
metaFeatures <- list(finalFeatures=pam50Short,datasetListIndicesToRemove=NULL)
load(paste0(saveDir,"/gapTestKmeans_pam50Short_nstart12015-05-04.RData.gzip"))
clusterCoINcIDE_output =  kmeansGapTest_pam50_short_Nstart1
#need to add fake PACR variable
clusterCoINcIDE_output$bestK_PACR <- unlist(sapply(clusterCoINcIDE_output$clustSampleIndexList,FUN=function(unit){
  
  length(unit)
  
})) 

load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
datasetNames <- names(esets)
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  datasetNames <- datasetNames[-metaFeatures$datasetListIndicesToRemove]
  
}
names(clusterCoINcIDE_output$bestK_PACR) <- datasetNames
clusterCoINcIDE_output$clustSampleIndexList_PACR <- clusterCoINcIDE_output$clustSampleIndexList
clusterCoINcIDE_output$clustFeatureIndexList_PACR <- clusterCoINcIDE_output$clustFeatureIndexList
names(clusterCoINcIDE_output$clustSampleIndexList_PACR) <- datasetNames
names(clusterCoINcIDE_output$clustFeatureIndexList_PACR) <- datasetNames

#from the testing K runs:
#RE-DO this CLUSTERING: IT APPEARS TO HAVE USED THE WRONG CLUSTERING INPUTS.
load(paste0(saveDir,"/adjMatrices_pam50ShortGapTest_pearson_meanMatrix_2015-05-15.RData.gzip"))
CoINcIDE_output = pam50ShortGapTest_pearson_meanMatrix
experimentName <- "breast_pam50Short_kmeansGap"
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

breast_kmeanGap_pam50ShortOut <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                              meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                                              clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                              commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                              survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                              CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                              ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)

save(breast_kmeansGap_pam50ShortOut,file=paste0("/home/kplaney/breast_analysis/",experimentName,"_",Sys.Date(),"/breast_hclustGap_pam50FullOut.RData.gzip"),compress="gzip")

####200 features
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
load(paste0(saveDir,"/metaFeatures_200.RData.gzip"))
metaFeatures =metaFeatures
load(paste0(saveDir,"/curatedbreastData_kmeansGap_nstart25_200_features_2015-05-19.RData.gzip"))
clusterCoINcIDE_output =  kmeansGap
clusterCoINcIDE_output$bestK_PACR <- unlist(sapply(clusterCoINcIDE_output$clustSampleIndexList,FUN=function(unit){
  
  length(unit)
  
})) 
load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
datasetNames <- names(esets)
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  datasetNames <- datasetNames[-metaFeatures$datasetListIndicesToRemove]
}
names(clusterCoINcIDE_output$bestK_PACR) <- datasetNames
clusterCoINcIDE_output$clustSampleIndexList_PACR <- clusterCoINcIDE_output$clustSampleIndexList
clusterCoINcIDE_output$clustFeatureIndexList_PACR <- clusterCoINcIDE_output$clustFeatureIndexList
names(clusterCoINcIDE_output$clustSampleIndexList_PACR)<- datasetNames
names(clusterCoINcIDE_output$clustFeatureIndexList_PACR )<- datasetNames
load(paste0(saveDir,"/adjMatrices_breast200F_kmeansGap_pearson_meanMatrix_2015-05-19.RData.gzip"))
CoINcIDE_output =  breast200F_kmeansGap_pearson_meanMatrix
experimentName <- "breast_200_kmeansGap"
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

breast_kmeansGap_200Out <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                             meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .5, maxTrueSimilThresh = Inf,
                                             clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                             commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,
                                             eset_featureDataFieldName=eset_featureDataFieldName,
                                             survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                             CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                             ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)



save(breast_kmeansGap_200Out,file="/home/kplaney/breast_analysis/",experimentName,"_",Sys.Date(),"/breast_200_features_analysisOut.RData.gzip",compress="gzip")



#######
"/home/kplaney/breast_analysis//curatedbreastData_kmeansGap_nstart25_200_features_2015-05-18.RData.gzip"
"/home/kplaney/breast_analysis/adjMatrices_breast200F_kmeansGap_pearson_meanMatrix_2015-05-18.RData.gzip"
#500 features
"/home/kplaney/breast_analysis//curatedbreastData_kmeansGap_nstart25_500_features_2015-05-18.RData.gzip"
"/home/kplaney/breast_analysis/adjMatrices_breast500F_kmeansGap_pearson_meanMatrix_2015-05-18.RData.gzip"
#1000 features
"/home/kplaney/breast_analysis//curatedbreastData_kmeansGap_nstart25_1000_features_2015-05-18.RData.gzip"

######
###gap test, hclust
##pam50 full

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
load("/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData")
pam50Genes <- centroidMatrix[,1]
metaFeatures <- list(finalFeatures=pam50Genes,datasetListIndicesToRemove=NULL)
load(paste0(saveDir,"/gapTesthclust_pam50Full_2015-05-15.RData.gzip"))

load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
#nstart is a misnomer - k-means only uses nstart.
clusterCoINcIDE_output =  hclustGapTest_pam50_full_Nstart25
#need to add fake PACR variable
clusterCoINcIDE_output$bestK_PACR <- unlist(sapply(clusterCoINcIDE_output$clustSampleIndexList,FUN=function(unit){
  
  length(unit)
  
})) 

datasetNames <- names(esets)
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  datasetNames <- datasetNames[-metaFeatures$datasetListIndicesToRemove]
}
names(clusterCoINcIDE_output$bestK_PACR) <- datasetNames
clusterCoINcIDE_output$clustSampleIndexList_PACR <- clusterCoINcIDE_output$clustSampleIndexList
clusterCoINcIDE_output$clustFeatureIndexList_PACR <- clusterCoINcIDE_output$clustFeatureIndexList
names(clusterCoINcIDE_output$clustSampleIndexList_PACR)<- datasetNames
names(clusterCoINcIDE_output$clustFeatureIndexList_PACR )<- datasetNames
#load CoINcIDE output
load(paste0(saveDir,"/adjMatrices_pam50FullHclustGapTest_pearson_meanMatrix_2015-05-18.RData.gzip"))
CoINcIDE_output = pam50FullHclustGapTest_pearson_meanMatrix
experimentName <- "breast_pam50Full_hclustGap"
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

breast_hclustGap_pam50FullOut <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                   meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                                   clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                   commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                   survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                   CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                   ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)

save(breast_hclustGap_pam50FullOut,file=paste0("/home/kplaney/breast_analysis/",experimentName,"_",Sys.Date(),"/breast_hclustGap_pam50FullOut.RData.gzip"),compress="gzip")


###pam50 short
load("/home/kplaney/breast_analysis/pam50Short_genes.RData")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
metaFeatures <- list(finalFeatures=pam50Short,datasetListIndicesToRemove=NULL)
load(paste0(saveDir,"/hclustGapTest_pam50_short_2015-05-15.RData.gzip"))
clusterCoINcIDE_output =  hclustGapTest_pam50_short
#need to add fake PACR variable
clusterCoINcIDE_output$bestK_PACR <- unlist(sapply(clusterCoINcIDE_output$clustSampleIndexList,FUN=function(unit){
  
  length(unit)
  
})) 

load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
datasetNames <- names(esets)
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  datasetNames <- datasetNames[-metaFeatures$datasetListIndicesToRemove]
  
}
names(clusterCoINcIDE_output$bestK_PACR) <- datasetNames
clusterCoINcIDE_output$clustSampleIndexList_PACR <- clusterCoINcIDE_output$clustSampleIndexList
clusterCoINcIDE_output$clustFeatureIndexList_PACR <- clusterCoINcIDE_output$clustFeatureIndexList
names(clusterCoINcIDE_output$clustSampleIndexList_PACR) <- datasetNames
names(clusterCoINcIDE_output$clustFeatureIndexList_PACR) <- datasetNames

#from the testing K runs:
load(paste0(saveDir,"/adjMatrices_pam50ShortHclustGapTest_pearson_meanMatrix_2015-05-18.RData.gzip"))
CoINcIDE_output = pam50ShortHclustGapTest_pearson_meanMatrix
experimentName <- "breast_pam50Short_hclustGap"
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

breast_hclustGap_pam50ShortOut <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                    meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                                    clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                    commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                    survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                    CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                    ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)

save(breast_hclustGap_pam50ShortOut,file=paste0("/home/kplaney/breast_analysis/",experimentName,"_",Sys.Date(),"/breast_hclustGap_pam50FullOut.RData.gzip"),compress="gzip")


####200 genes
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
load(paste0(saveDir,"/metaFeatures_200.RData.gzip"))
metaFeatures =metaFeatures
load(paste0(saveDir,"/curatedbreastData_hclust_200Features_2015-05-15.RData.gzip"))
clusterCoINcIDE_output =  hclustOut

#need to add fake PACR variable
clusterCoINcIDE_output$bestK_PACR <- unlist(sapply(clusterCoINcIDE_output$clustSampleIndexList,FUN=function(unit){
  
  length(unit)
  
})) 

load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
datasetNames <- names(esets)
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  datasetNames <- datasetNames[-metaFeatures$datasetListIndicesToRemove]
  
}
names(clusterCoINcIDE_output$bestK_PACR) <- datasetNames
clusterCoINcIDE_output$clustSampleIndexList_PACR <- clusterCoINcIDE_output$clustSampleIndexList
clusterCoINcIDE_output$clustFeatureIndexList_PACR <- clusterCoINcIDE_output$clustFeatureIndexList
names(clusterCoINcIDE_output$clustSampleIndexList_PACR) <- datasetNames
names(clusterCoINcIDE_output$clustFeatureIndexList_PACR) <- datasetNames

load(paste0(saveDir,"/adjMatrices_breast200F_hclustGap_pearson_meanMatrix_2015-05-17.RData.gzip"))
CoINcIDE_output = breast200F_hclustGap_pearson_meanMatrix
experimentName <- "breast_200_hclustGap"
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

breast_hclustGap_200Out <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                             meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                             clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                             commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,
                                             eset_featureDataFieldName=eset_featureDataFieldName,
                                             survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                             CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                             ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)




save(breast_hclustGap_200Out,file=paste0("/home/kplaney/breast_analysis/",experimentName,"_",Sys.Date(),"/breast_hclustGap_200Out.RData.gzip"),compress="gzip")


###500 features
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
load(paste0(saveDir,"/metaFeatures_500.RData.gzip"))
metaFeatures =metaFeatures
load(paste0(saveDir,"/curatedbreastData_hclust_500Features_2015-05-15.RData.gzip"))
clusterCoINcIDE_output =  hclustOut

#need to add fake PACR variable
clusterCoINcIDE_output$bestK_PACR <- unlist(sapply(clusterCoINcIDE_output$clustSampleIndexList,FUN=function(unit){
  
  length(unit)
  
})) 

load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
datasetNames <- names(esets)
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  datasetNames <- datasetNames[-metaFeatures$datasetListIndicesToRemove]
  
}
names(clusterCoINcIDE_output$bestK_PACR) <- datasetNames
clusterCoINcIDE_output$clustSampleIndexList_PACR <- clusterCoINcIDE_output$clustSampleIndexList
clusterCoINcIDE_output$clustFeatureIndexList_PACR <- clusterCoINcIDE_output$clustFeatureIndexList
names(clusterCoINcIDE_output$clustSampleIndexList_PACR) <- datasetNames
names(clusterCoINcIDE_output$clustFeatureIndexList_PACR) <- datasetNames

load(paste0(saveDir,"/adjMatrices_breast500F_hclustGap_pearson_meanMatrix_2015-05-17.RData.gzip"))
CoINcIDE_output = breast500F_hclustGap_pearson_meanMatrix
experimentName <- "breast_500_hclustGap"
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

breast_hclustGap_500Out <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                       meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                                       clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                       commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,
                                                       eset_featureDataFieldName=eset_featureDataFieldName,
                                                       survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                       CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                       ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)




save(breast_hclustGap_500Out,file=paste0("/home/kplaney/breast_analysis/",experimentName,"_",Sys.Date(),"/breast_hclustGap_500Out.RData.gzip"),compress="gzip")



####1000 genes
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
load(paste0(saveDir,"/metaFeatures_1000.RData.gzip"))
metaFeatures =metaFeatures
load(paste0(saveDir,"/curatedbreastData_hclust_1000Features_2015-05-15.RData.gzip"))
clusterCoINcIDE_output =  hclustOut

#need to add fake PACR variable
clusterCoINcIDE_output$bestK_PACR <- unlist(sapply(clusterCoINcIDE_output$clustSampleIndexList,FUN=function(unit){
  
  length(unit)
  
})) 

load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
datasetNames <- names(esets)
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  datasetNames <- datasetNames[-metaFeatures$datasetListIndicesToRemove]
  
}
names(clusterCoINcIDE_output$bestK_PACR) <- datasetNames
clusterCoINcIDE_output$clustSampleIndexList_PACR <- clusterCoINcIDE_output$clustSampleIndexList
clusterCoINcIDE_output$clustFeatureIndexList_PACR <- clusterCoINcIDE_output$clustFeatureIndexList
names(clusterCoINcIDE_output$clustSampleIndexList_PACR) <- datasetNames
names(clusterCoINcIDE_output$clustFeatureIndexList_PACR) <- datasetNames

load(paste0(saveDir,"/adjMatrices_breast1000F_hclustGap_pearson_meanMatrix_2015-05-17.RData.gzip"))
CoINcIDE_output = breast1000F_hclustGap_pearson_meanMatrix
experimentName <- "breast_1000_hclustGap"
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

breast_hclustGap_1000Out <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                       meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                                       clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                       commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,
                                                       eset_featureDataFieldName=eset_featureDataFieldName,
                                                       survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                       CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                       ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)




save(breast_hclustGap_1000Out,file=paste0("/home/kplaney/breast_analysis/",experimentName,"_",Sys.Date(),"/breast_hclustGap_1000Out.RData.gzip"),compress="gzip")


###2000 genes
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
load(paste0(saveDir,"/metaFeatures_2000.RData.gzip"))
metaFeatures =metaFeatures
load(paste0(saveDir,"/curatedbreastData_hclust_2000Features_2015-05-16.RData.gzip"))
clusterCoINcIDE_output =  hclustOut

#need to add fake PACR variable
clusterCoINcIDE_output$bestK_PACR <- unlist(sapply(clusterCoINcIDE_output$clustSampleIndexList,FUN=function(unit){
  
  length(unit)
  
})) 

load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
datasetNames <- names(esets)
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  datasetNames <- datasetNames[-metaFeatures$datasetListIndicesToRemove]
  
}
names(clusterCoINcIDE_output$bestK_PACR) <- datasetNames
clusterCoINcIDE_output$clustSampleIndexList_PACR <- clusterCoINcIDE_output$clustSampleIndexList
clusterCoINcIDE_output$clustFeatureIndexList_PACR <- clusterCoINcIDE_output$clustFeatureIndexList
names(clusterCoINcIDE_output$clustSampleIndexList_PACR) <- datasetNames
names(clusterCoINcIDE_output$clustFeatureIndexList_PACR) <- datasetNames

load(paste0(saveDir,"/adjMatrices_2000F_hclustGap_pearson_meanMatrix_2015-05-17.RData.gzip"))
CoINcIDE_output = breast2000F_hclustGap_pearson_meanMatrix
experimentName <- "breast_2000_hclustGap"
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

breast_hclustGap_2000Out <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                        meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                                        clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                        commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,
                                                        eset_featureDataFieldName=eset_featureDataFieldName,
                                                        survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                        CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                        ovarian=ovarian,fisherTestVariableLegendNames=fisherTestVariableLegendNames,fisherTestVariableTitleNames=fisherTestVariableTitleNames)




save(breast_hclustGap_2000Out,file=paste0("/home/kplaney/breast_analysis/",experimentName,"_",Sys.Date(),"/breast_hclustGap_2000Out.RData.gzip"),compress="gzip")



############
#THEN: pearson + k-means consensus, using centroid method

####pam50Full
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
load("/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData")
pam50Genes <- centroidMatrix[,1]
metaFeatures <- list(finalFeatures=pam50Genes,datasetListIndicesToRemove=NULL)
load(paste0(saveDir,"/kmeansConsensuspam50_full_Nstart1pItem9_pam50FullFeatures_2015-04-28.RData.gzip"))

clusterCoINcIDE_output =  kmeansConsensuspam50_full_Nstart1pItem9
load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
#from the testing K runs:
load(paste0(saveDir,"/adjMatrices_pam50Full_pearson_centroid_2015-05-18.RData.gzip"))
CoINcIDE_output = pam50Full_pearson_centroid
#use the .3 instead of the .1 version:
CoINcIDE_output$pvalueMatrix <- CoINcIDE_output$pvalueMatrix3
experimentName <- "breast_pam50Full_centroid"
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

save(breast_centroid_pam50FullOut,file=paste0("/home/kplaney/breast_analysis/",experimentName,"_",Sys.Date(),
                                              "/",experimentName,"Out.RData.gzip"),compress="gzip")


####pam50 short kmeans
load("/home/kplaney/breast_analysis/pam50Short_genes.RData")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
metaFeatures <- list(finalFeatures=pam50Short,datasetListIndicesToRemove=NULL)
load(paste0(saveDir,"/kmeansConsensuspam50_short_Nstart1pItem9_pam50ShortFeatures_2015-04-28.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensuspam50_short_Nstart1pItem9
load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
#from the testing K runs:
load(paste0(saveDir,"/adjMatrices_pam50Short_pearson_centroid_2015-05-18.RData.gzip"))
CoINcIDE_output = pam50Short_pearson_centroid
#use the .3 instead of the .1 version:
CoINcIDE_output$pvalueMatrix <- CoINcIDE_output$pvalueMatrix3
experimentName <- "breast_pam50Short_centroid"
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

save(breast_centroid_pam50ShortOut,file=paste0("/home/kplaney/breast_analysis/",experimentName,"_",Sys.Date(),
                                      "/",experimentName,"Out.RData.gzip"),compress="gzip")

###200
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
load(paste0(saveDir,"/metaFeatures_200.RData.gzip"))
metaFeatures =metaFeatures
load(paste0(saveDir,"/curatedbreastData_kmeansConsensus_nstart1pItem9200Features_2015-05-19.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensus
load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0(saveDir,"/adjMatrices_200F_pearson_centroid_2015-05-20.RData.gzip"))
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



save(breast_centroid_200Out,file=paste0("/home/kplaney/breast_analysis/",experimentName,"_",Sys.Date(),
                               "/",experimentName,"Out.RData.gzip"),compress="gzip")
####500
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
load(paste0(saveDir,"/metaFeatures_500.RData.gzip"))
metaFeatures =metaFeatures
load(paste0(saveDir,"/curatedbreastData_kmeansConsensus_nstart1pItem9500Features_2015-05-04.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensus
load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0(saveDir,"/adjMatrices_500F_pearson_centroid_2015-05-18.RData.gzip"))
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



save(breast_centroid_500Out,file=paste0("/home/kplaney/breast_analysis/",experimentName,"_",Sys.Date(),
                               "/",experimentName,"Out.RData.gzip"),compress="gzip")

####1000
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
load(paste0(saveDir,"/metaFeatures_1000.RData.gzip"))
metaFeatures =metaFeatures
load(paste0(saveDir,"/curatedbreastData_kmeansConsensus_nstart1pItem91000Features_2015-05-05.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensus
load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0(saveDir,"/adjMatrices_1000F_pearson_centroid_2015-05-18.RData.gzip"))
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



save(breast_centroid_1000Out,file=paste0("/home/kplaney/breast_analysis/",experimentName,"_",Sys.Date(),
                                "/",experimentName,"Out.RData.gzip"),compress="gzip")
###2000
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list

saveDir <- "/home/kplaney/breast_analysis/"
load(paste0(saveDir,"/metaFeatures_2000.RData.gzip"))
metaFeatures =metaFeatures
load(paste0(saveDir,"/curatedbreastData_kmeansConsensus_nstart1pItem92000Features_2015-05-05.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensus
load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0(saveDir,"/adjMatrices_2000F_pearson_meanMatrix_centroid_2015-05-18.RData.gzip"))
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



save(breast_centroid_2000Out,file=paste0("/home/kplaney/breast_analysis/",experimentName,"_",Sys.Date(),
                                "/",experimentName,"Out.RData.gzip"),compress="gzip")

