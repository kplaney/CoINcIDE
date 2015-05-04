
load("/home/kplaney/breast_analysis/curatedBreastData_esets_proc.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_featureSelection.R")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
#combat code:
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_batchCorrection.R")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_visualization.R")
#now format just as a list of data matrices.
dataMatrixList <- exprSetListToMatrixList(esets,featureDataFieldName="gene_symbol")

names(dataMatrixList) <- names(esets)

##ALSO: merge matrices first to help decide which studies to keep overall

#also merge this one
output <- merge_datasetList(datasetList=dataMatrixList,minNumGenes = 10000, minNumPatients = 40,batchNormalize = c('none'));
save(output,file="/home/kplaney/breast_analysis//mergedExprMatrix_minVar001_17_studies_no_norm.RData.gzip",compress="gzip")

load("/home/kplaney/breast_analysis/pam50Short_genes.RData")

dataMatrixList <- list(mergedNoNorm=output$mergedExprMatrix)

clustFeaturesList <- list(pam50Short=pam50Short)

#we know these are strong clusters. have  minMeanClustConsensus=.8
kmeansConsensus <-  clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                             pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                             numSims=500,maxNumClusters=20,
                                             outputFile="/home/kplaney/breast_analysis/test.txt",distMethod=c("euclidean"),
                                             hclustAlgorithm=c("average"),
                                             consensusHclustAlgorithm=c("average"),
                                             minClustConsensus=.7, minMeanClustConsensus=.85,
                                             corUse="everything",pItem=.9,maxPAC=.15)


save(kmeansConsensus,file=paste0("/home/kplaney/breast_analysis/curatedBreastData_kmeansConsensus_mergedNoNorm_pam50Short_nstart1",Sys.Date(),".RData.gzip"),compress="gzip")

#how do the centroids look?
load("/home/kplaney/breast_analysis/curatedBreastData_kmeansConsensus_mergedNoNorm_pam50Short_nstart12015-05-01.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_visualization.R")
# combat must be run on intersecting genes only; so limited pam50 gene set is a confounding factor.
#but still, if subtypes changes between no norm and combat with limited pam50 gene set, this is still an issue:
##also: look at effects of BMC normalization:
load("/home/kplaney/breast_analysis//mergedExprMatrix_minVar001_17_studies_no_norm.RData.gzip")
mergedNoNorm <- output$mergedExprMatrix
Pam50_subtypes_noNorm <- assignCentroidSubtype(t(mergedNoNorm),minNumGenes=30,centroidRData="/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData");


#ooof..these are really mixed up!
table(Pam50_subtypes_noNorm$subtype[kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]] ,"subtypeLabels"])
table(Pam50_subtypes_noNorm$subtype[kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]] ,"subtypeLabels"])
#and this cluster is super tiny!
table(Pam50_subtypes_noNorm$subtype[kmeansConsensus$clustSampleIndexList_PACR[[1]][[3]] ,"subtypeLabels"])

attrDF <- data.frame(c(rep.int(1,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]])),
                       rep.int(2,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]])),
                       rep.int(3,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[3]]))
                       ),
                     c(Pam50_subtypes_noNorm$subtype[kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]] ,"subtypeLabels"],
                       Pam50_subtypes_noNorm$subtype[kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]] ,"subtypeLabels"],
                       Pam50_subtypes_noNorm$subtype[kmeansConsensus$clustSampleIndexList_PACR[[1]][[3]] ,"subtypeLabels"]))
library("ggplot2")
colnames(attrDF) <- c("clustNum","subtype")
#....not good
ggplot(data=attrDF,aes(subtype))+geom_bar() + facet_grid(.~clustNum)

#try with the pull pam50 subtypes, too?
load("/home/kplaney/breast_analysis/pam50Full_subtypeDF.RData.gzip")
subtypes <- subtypeDF$subtype[na.omit(match(colnames(mergedNoNorm),subtypeDF$sampleName))]

#ooof..these are really mixed up!
table(subtypes[kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]]])
table(subtypes[kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]]])
#and this cluster is super tiny!
table(subtypes[kmeansConsensus$clustSampleIndexList_PACR[[1]][[3]]])

attrDF_pam50Full <- data.frame(c(rep.int(1,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]])),
                       rep.int(2,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]])),
                       rep.int(3,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[3]]))
),
c(subtypes[kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]]],
  subtypes[kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]]],
  subtypes[kmeansConsensus$clustSampleIndexList_PACR[[1]][[3]]]))
library("ggplot2")
colnames(attrDF_pam50Full) <- c("clustNum","subtype")
#....not good
ggplot(data=attrDF_pam50Full,aes(subtype))+geom_bar() + facet_grid(.~clustNum)

#########
#BMC
output <- merge_datasetList(datasetList=dataMatrixList,minNumGenes = 10000, minNumPatients = 40,batchNormalize = c('BMC'));
save(output,file="/home/kplaney/breast_analysis//mergedExprMatrix_minVar001_17_studies_BMC_norm.RData.gzip",compress="gzip")


load("/home/kplaney/breast_analysis/pam50Short_genes.RData")

dataMatrixList <- list(mergedBMC=output$mergedExprMatrix)

clustFeaturesList <- list(pam50Short=pam50Short)

#we know these are strong clusters. have  minMeanClustConsensus=.8
kmeansConsensus <-  clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                           pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                           numSims=500,maxNumClusters=20,
                                           outputFile="/home/kplaney/breast_analysis/test.txt",distMethod=c("euclidean"),
                                           hclustAlgorithm=c("average"),
                                           consensusHclustAlgorithm=c("average"),
                                           minClustConsensus=.7, minMeanClustConsensus=.85,
                                           corUse="everything",pItem=.9,maxPAC=.15)


save(kmeansConsensus,file=paste0("/home/kplaney/breast_analysis/curatedBreastData_kmeansConsensus_mergedBMCNorm_pam50Short_nstart1",Sys.Date(),".RData.gzip"),compress="gzip")

#########
#combat
###did not detect a batch effect!! p-value was .055....so try .1:
output <- merge_datasetList(datasetList=dataMatrixList,minNumGenes = 10000, minNumPatients = 40,batchNormalize = c('combat'),combatPvalueThresh=.1);
save(output,file="/home/kplaney/breast_analysis//mergedExprMatrix_minVar001_17_studies_combat_norm.RData.gzip",compress="gzip")
message("Done with Combat")

load("/home/kplaney/breast_analysis/pam50Short_genes.RData")

dataMatrixList <- list(mergedCombat=output$mergedExprMatrix)

clustFeaturesList <- list(pam50Short=pam50Short)

#we know these are strong clusters. have  minMeanClustConsensus=.8
kmeansConsensus <-  clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                           pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                           numSims=500,maxNumClusters=20,
                                           outputFile="/home/kplaney/breast_analysis/test.txt",distMethod=c("euclidean"),
                                           hclustAlgorithm=c("average"),
                                           consensusHclustAlgorithm=c("average"),
                                           minClustConsensus=.7, minMeanClustConsensus=.85,
                                           corUse="everything",pItem=.9,maxPAC=.15)


save(kmeansConsensus,file=paste0("/home/kplaney/breast_analysis/curatedBreastData_kmeansConsensus_mergedCombatNorm_pam50Short_nstart1",Sys.Date(),".RData.gzip"),compress="gzip")


###look at Pam50 breakdown
load("/home/kplaney/breast_analysis/curatedBreastData_kmeansConsensus_mergedCombatNorm_pam50Short_nstart12015-05-02.RData.gzip")
load("/home/kplaney/breast_analysis//mergedExprMatrix_minVar001_17_studies_combat_norm.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_visualization.R")
mergedCombat <- output$mergedExprMatrix
Pam50_subtypes_noNorm <- assignCentroidSubtype(t(mergedCombat),minNumGenes=30,centroidRData="/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData");


#combat does make it better:
table(Pam50_subtypes_noNorm$subtype[kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]] ,"subtypeLabels"])
table(Pam50_subtypes_noNorm$subtype[kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]] ,"subtypeLabels"])
#and this cluster is super tiny!
table(Pam50_subtypes_noNorm$subtype[kmeansConsensus$clustSampleIndexList_PACR[[1]][[3]] ,"subtypeLabels"])

attrDF <- data.frame(c(rep.int(1,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]])),
                       rep.int(2,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]])),
                       rep.int(3,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[3]]))
),
c(Pam50_subtypes_noNorm$subtype[kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]] ,"subtypeLabels"],
  Pam50_subtypes_noNorm$subtype[kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]] ,"subtypeLabels"],
  Pam50_subtypes_noNorm$subtype[kmeansConsensus$clustSampleIndexList_PACR[[1]][[3]] ,"subtypeLabels"]))
library("ggplot2")
colnames(attrDF) <- c("clustNum","subtype")
#....not good
ggplot(data=attrDF,aes(subtype))+geom_bar() + facet_grid(.~clustNum)

#try with the pull pam50 subtypes, too?
load("/home/kplaney/breast_analysis/pam50Full_subtypeDF.RData.gzip")
subtypes <- subtypeDF$subtype[na.omit(match(colnames(mergedCombat),subtypeDF$sampleName))]

#ooof..these are really mixed up!
table(subtypes[kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]]])
table(subtypes[kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]]])
#and this cluster is super tiny!
table(subtypes[kmeansConsensus$clustSampleIndexList_PACR[[1]][[3]]])

attrDF_pam50Full <- data.frame(c(rep.int(1,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]])),
                                 rep.int(2,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]])),
                                 rep.int(3,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[3]]))
),
c(subtypes[kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]]],
  subtypes[kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]]],
  subtypes[kmeansConsensus$clustSampleIndexList_PACR[[1]][[3]]]))
library("ggplot2")
colnames(attrDF_pam50Full) <- c("clustNum","subtype")
#....not good
ggplot(data=attrDF_pam50Full,aes(subtype))+geom_bar() + facet_grid(.~clustNum)

########inspecting how normalization affects pam50 groups:
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_visualization.R")
# combat must be run on intersecting genes only; so limited pam50 gene set is a confounding factor.
#but still, if subtypes changes between no norm and combat with limited pam50 gene set, this is still an issue:
##also: look at effects of BMC normalization:
load("/home/kplaney/breast_analysis//mergedExprMatrix_minVar001_17_studies_BMC_norm.RData.gzip")
mergedBMC <- output$mergedExprMatrix
load("/home/kplaney/breast_analysis//mergedExprMatrix_minVar001_17_studies_no_norm.RData.gzip")
mergedNoNorm <- output$mergedExprMatrix
Pam50_subtypes_noNorm <- assignCentroidSubtype(t(mergedNoNorm),minNumGenes=30,centroidRData="/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData");
Pam50_subtypes_BMCNorm <- assignCentroidSubtype(t(mergedBMC),minNumGenes=30,centroidRData="/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData");
#wow...a lot changed!!! of course, if all pam50 genes were in all datasets, these numbers may slightly changed, but unlikely by much.
length(which(Pam50_subtypes_noNorm$subtype[,"subtypeLabels"] != Pam50_subtypes_BMCNorm$subtypes[,"subtypeLabels"]))

load("/home/kplaney/breast_analysis//mergedExprMatrix_minVar001_17_studies_combat_norm.RData.gzip")
mergedCombat <- output$mergedExprMatrix
Pam50_subtypes_combatNorm <- assignCentroidSubtype(t(mergedCombat),minNumGenes=30,centroidRData="/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData");


