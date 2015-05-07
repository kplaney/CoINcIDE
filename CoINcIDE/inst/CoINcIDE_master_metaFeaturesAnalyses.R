
######ovarian 200
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
load("/home/kplaney/ovarian_analysis/metaFeatures_200.RData.gzip")
metaFeatures =metaFeatures
load("/home/kplaney/ovarian_analysis/curatedOvarianData_kmeansConsensus_nstart1_200Features_2015-04-28.RData.gzip")
clusterCoINcIDE_output =  kmeansConsensus
load("/home/kplaney/ovarian_analysis/esets_proc_TCGAcombat.RData.gzip")
esets = esets
#dataMatrixList
load("/home/kplaney/ovarian_analysis/adjMatrices_200F_pearson_meanMatrix_2015-04-29RData.gzip")
CoINcIDE_output = ov_200F_pearson_meanMatrix
experimentName <- "ovarian_200_features"
fisherTestVariables <- c("histological_type","tumorstage","recurrence_status","grade","age_at_initial_pathologic_diagnosis")
eset_featureDataFieldName="gene"
networkColors = "Set2"
saveDir <- "/home/kplaney/ovarian_analysis/"
outcomesVarBinary="vital_status"
outcomesVarCont = "days_to_death"
ovarian <- TRUE
#breast: must also input dataMatrixList
#brewPal = c("Set3","Paired","Spectral","BrBG","PiYG","RdYlGn","RdYlBu","RdBu","PiYG","Set2"),
ov_200Out <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                        meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                        clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                        commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                        survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                        CutoffPointYears=5, eset_uniquePatientID="unique_patient_ID", fisherTestVariables = fisherTestVariables)

  

######ovarian 500

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
load("/home/kplaney/ovarian_analysis/metaFeatures_500.RData.gzip")
metaFeatures =metaFeatures
load("/home/kplaney/ovarian_analysis/curatedOvarianData_kmeansConsensus_nstart1_500Features_2015-04-29.RData.gzip")
clusterCoINcIDE_output =  kmeansConsensus
load("/home/kplaney/ovarian_analysis/esets_proc_TCGAcombat.RData.gzip")
esets = esets
#dataMatrixList
load("/home/kplaney/ovarian_analysis/adjMatrices_500F_pearson_meanMatrix_2015-04-29RData.gzip")
CoINcIDE_output = ov_500F_pearson_meanMatrix
experimentName <- "ovarian_500_features"
fisherTestVariables <- c("histological_type","tumorstage","recurrence_status","grade","age_at_initial_pathologic_diagnosis")
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
                                         CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables)


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
                                         CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables)




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
                                          CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables)



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
fisherTestVariables <- c("OS","DFS","RFS","RCB","metastasis","pCR","near_pCR","clinical_AJCC_stage" ,"tumor_stage_preTrt",
                         "preTrt_lymph_node_status" , "preTrt_numPosLymphNodes","family_history","race","nationality",
                         "neoadjuvant_or_adjuvant","PR_preTrt","HER2_preTrt","hist_grade","path_diagnosis","path")

breast_pam50FullCentroidsOut <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                             meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                             clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                             commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                             survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                             CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                             ovarian=ovarian)


save(breast_pam50FullCentroidsOut,file="/home/kplaney/breast_analysis/breast_pam50FullCentroid_features_2015-05-06//breast_pam50FullCentroidsOut.RData.gzip",compress="gzip")

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
fisherTestVariables <- c("OS","DFS","RFS","RCB","metastasis","pCR","near_pCR","clinical_AJCC_stage" ,"tumor_stage_preTrt",
                         "preTrt_lymph_node_status" , "preTrt_numPosLymphNodes","family_history","race","nationality",
                         "neoadjuvant_or_adjuvant","PR_preTrt","HER2_preTrt","hist_grade","path_diagnosis","path")

breast_pam50FullOut <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                            meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                                            clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                            commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                            survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                            CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                            ovarian=ovarian)

save(breast_pam50FullOut,file="/home/kplaney/breast_analysis/breast_pam50Full_features_2015-05-06//breast_pam50FullOut.RData.gzip",compress="gzip")


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
fisherTestVariables <- c("OS","DFS","RFS","RCB","metastasis","pCR","near_pCR","clinical_AJCC_stage" ,"tumor_stage_preTrt",
                         "preTrt_lymph_node_status" , "preTrt_numPosLymphNodes","family_history","race","nationality",
                         "neoadjuvant_or_adjuvant","PR_preTrt","HER2_preTrt","hist_grade","path_diagnosis","path")

breast_pam50ShortOut <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                   meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                                   clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                   commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName=eset_featureDataFieldName,
                                                   survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                   CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                   ovarian=ovarian)

save(breast_pam50ShortOut,file="/home/kplaney/breast_analysis/breast_pam50Short_features_2015-05-06//breast_pam50ShortOut.RData.gzip",compress="gzip")


###200
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
load(paste0(saveDir,"/metaFeatures_200.RData.gzip"))
metaFeatures =metaFeatures
load(paste0(saveDir,"/curatedbreastData_kmeansConsensus_nstart1pItem9200Features_2015-05-04.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensus
load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0(saveDir,"/adjMatrices_pam50Short_pearson_meanMatrix_2015-04-29.RData.gzip")
CoINcIDE_output = breast200F_pearson_meanMatrix
experimentName <- "breast_200_features"
eset_featureDataFieldName="gene_symbol"
networkColors = "Set3"
outcomesVarBinary=NA
#not enough continuous variables for breast data
outcomesVarCont = NA
ovarian <- FALSE
eset_uniquePatientID="dbUniquePatientID"
fisherTestVariables <- c("OS","DFS","RFS","RCB","metastasis","pCR","near_pCR","clinical_AJCC_stage" ,"tumor_stage_preTrt",
  "preTrt_lymph_node_status" , "preTrt_numPosLymphNodes","family_history","race","nationality",
  "neoadjuvant_or_adjuvant","PR_preTrt","HER2_preTrt","hist_grade","path_diagnosis","path")

breast_200Out <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                          meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                          clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                          commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,
                                          eset_featureDataFieldName=eset_featureDataFieldName,
                                          survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                          CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                          ovarian=ovarian)



save(breast_200Out,file="/home/kplaney/breast_analysis/breast_200_features_2015-05-06/breast_200_features_analysisOut.RData.gzip",compress="gzip")

####500
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
load(paste0(saveDir,"/metaFeatures_200.RData.gzip"))
metaFeatures =metaFeatures
load(paste0(saveDir,"/curatedbreastData_kmeansConsensus_nstart1pItem9200Features_2015-05-04.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensus
load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0(saveDir,"/adjMatrices_pam50Short_pearson_meanMatrix_2015-04-29.RData.gzip")
     CoINcIDE_output = breast200F_pearson_meanMatrix
     experimentName <- "breast_200_features"
     eset_featureDataFieldName="gene_symbol"
     networkColors = "Set3"
     outcomesVarBinary=NA
     #not enough continuous variables for breast data
     outcomesVarCont = NA
     ovarian <- FALSE
     eset_uniquePatientID="dbUniquePatientID"
     fisherTestVariables <- c("OS","DFS","RFS","RCB","metastasis","pCR","near_pCR","clinical_AJCC_stage" ,"tumor_stage_preTrt",
                              "preTrt_lymph_node_status" , "preTrt_numPosLymphNodes","family_history","race","nationality",
                              "neoadjuvant_or_adjuvant","PR_preTrt","HER2_preTrt","hist_grade","path_diagnosis","path")
     
     breast_200Out <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                  meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                                  clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                  commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,
                                                  eset_featureDataFieldName=eset_featureDataFieldName,
                                                  survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                  CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                  ovarian=ovarian)
     
     
     
     save(breast_200Out,file="/home/kplaney/breast_analysis/breast_200_features_2015-05-06/breast_200_features_analysisOut.RData.gzip",compress="gzip")
     
     
####1000
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
load(paste0(saveDir,"/metaFeatures_200.RData.gzip"))
metaFeatures =metaFeatures
load(paste0(saveDir,"/curatedbreastData_kmeansConsensus_nstart1pItem9200Features_2015-05-04.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensus
load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0(saveDir,"/adjMatrices_pam50Short_pearson_meanMatrix_2015-04-29.RData.gzip")
     CoINcIDE_output = breast200F_pearson_meanMatrix
     experimentName <- "breast_200_features"
     eset_featureDataFieldName="gene_symbol"
     networkColors = "Set3"
     outcomesVarBinary=NA
     #not enough continuous variables for breast data
     outcomesVarCont = NA
     ovarian <- FALSE
     eset_uniquePatientID="dbUniquePatientID"
     fisherTestVariables <- c("OS","DFS","RFS","RCB","metastasis","pCR","near_pCR","clinical_AJCC_stage" ,"tumor_stage_preTrt",
                              "preTrt_lymph_node_status" , "preTrt_numPosLymphNodes","family_history","race","nationality",
                              "neoadjuvant_or_adjuvant","PR_preTrt","HER2_preTrt","hist_grade","path_diagnosis","path")
     
     breast_200Out <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                  meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                                  clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                  commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,
                                                  eset_featureDataFieldName=eset_featureDataFieldName,
                                                  survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                  CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                  ovarian=ovarian)
     
     
     
     save(breast_200Out,file="/home/kplaney/breast_analysis/breast_200_features_2015-05-06/breast_200_features_analysisOut.RData.gzip",compress="gzip")
     
###2000
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/inst/CoINcIDE_metaFeatures_analysis_wrapper.R")
#grab data matrix list, clust features list
saveDir <- "/home/kplaney/breast_analysis/"
load(paste0(saveDir,"/metaFeatures_200.RData.gzip"))
metaFeatures =metaFeatures
load(paste0(saveDir,"/curatedbreastData_kmeansConsensus_nstart1pItem9200Features_2015-05-04.RData.gzip"))
clusterCoINcIDE_output =  kmeansConsensus
load(paste0(saveDir,"/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
esets=esets_minVar001_17_studies
load(paste0(saveDir,"/adjMatrices_pam50Short_pearson_meanMatrix_2015-04-29.RData.gzip")
     CoINcIDE_output = breast200F_pearson_meanMatrix
     experimentName <- "breast_200_features"
     eset_featureDataFieldName="gene_symbol"
     networkColors = "Set3"
     outcomesVarBinary=NA
     #not enough continuous variables for breast data
     outcomesVarCont = NA
     ovarian <- FALSE
     eset_uniquePatientID="dbUniquePatientID"
     fisherTestVariables <- c("OS","DFS","RFS","RCB","metastasis","pCR","near_pCR","clinical_AJCC_stage" ,"tumor_stage_preTrt",
                              "preTrt_lymph_node_status" , "preTrt_numPosLymphNodes","family_history","race","nationality",
                              "neoadjuvant_or_adjuvant","PR_preTrt","HER2_preTrt","hist_grade","path_diagnosis","path")
     
     breast_200Out <- metaFeaturesAnalysisWrapper(metaFeatures=metaFeatures,esets=esets,CoINcIDE_output=CoINcIDE_output , clusterCoINcIDE_output=clusterCoINcIDE_output,
                                                  meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
                                                  clustSizeThresh = 5,saveDir =saveDir,experimentName = experimentName,networkColors = networkColors,
                                                  commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,
                                                  eset_featureDataFieldName=eset_featureDataFieldName,
                                                  survivalAnalysis=TRUE,outcomesVarBinary=outcomesVarBinary,outcomesVarCont = outcomesVarCont,
                                                  CutoffPointYears=5, eset_uniquePatientID=eset_uniquePatientID, fisherTestVariables = fisherTestVariables,
                                                  ovarian=ovarian)
     
     
     
     save(breast_200Out,file="/home/kplaney/breast_analysis/breast_200_features_2015-05-06/breast_200_features_analysisOut.RData.gzip",compress="gzip")
     