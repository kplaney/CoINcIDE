
#this script follows from: breastProcessAndGeneFeature_script.R
#It analyses the breast data using concatenated clustering on a merged 
#dataset.

#it assumes you have these files in your path (they are in the same folder 
#as this script in the github repo):
#pam50_centroids_updatedSymbols.RData
#pam50Short_genes.RData

#and note the libraries called later down: plyr, ggplot2

#CHANGE these paths to match ones have used in scripts that 
#feed into this one
saveDirGlobal <- "/home/ywrfc09/breast_analysis"
saveDirMerged <- "/home/ywrfc09/breast_analysis/mergedMatrix/"
outputFile <- "/home/kplaney/breast_analysis/clust_test_outMessages.txt"
library("Coincide")
#need to source a function tailored to analyzing specifically the curatedBreastData
#binary survival outcomes:
source("Coincide_metaFeaturesAnalysisWrapper.R")


####concatenated/merged clustering on no normalization, BMC and ComBat normalization:

dataMatrixList <- readRDS(paste0(saveDirGlobal,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.rds"))
##this code was already run in GenomeBiology_breastCluster_script.R to help with selection of final datasets.
output <- merge_datasetList(datasetList=dataMatrixList,minNumGenes = 10000, minNumPatients = 40,batchNormalize = c('none'));
saveRDS(output,file=paste0(saveDirMerged, "/mergedExprMatrix_minVar001_17_studies_no_norm.rds",compress=TRUE))

load("pam50Short_genes.RData")

#now our data matrix "list" is just one gigantic matrix
dataMatrixList <- list(mergedNoNorm=output$mergedExprMatrix)

clustFeaturesList <- list(pam50Short=pam50Short)

#we know these are strong clusters. have  minMeanClustConsensus=.8
kmeansConsensus <-  clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                             pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                             numSims=500,maxNumClusters=20,
                                             outputFile=outputFile,distMethod=c("euclidean"),
                                             hclustAlgorithm=c("average"),
                                             consensusHclustAlgorithm=c("average"),
                                             minClustConsensus=.7, minMeanClustConsensus=.85,
                                             corUse="everything",pItem=.9,maxPAC=.15)


saveRDS(kmeansConsensus,file=paste0(saveDirMerged,"/curatedBreastData_kmeansConsensus_mergedNoNorm_pam50Short_nstart1",".rds"),compress=TRUE)

#how do the centroids look?

# combat must be run on intersecting genes only; so limited pam50 gene set is a confounding factor.
#but still, if subtypes changes between no norm and combat with limited pam50 gene set, this is still an issue:
mergedNoNorm <- output$mergedExprMatrix
Pam50_subtypes_noNorm <- assignCentroidSubtype(t(mergedNoNorm),minNumGenes=30,centroidRData="pam50_centroids_updatedSymbols.RData");


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

#try with the full pam50 subtypes, too?
#we made this in the breast processing script:
load(paste0(saveDirGlobal, "/pam50Full_subtypeDF.RData.gzip"))
subtypes <- subtypeDF$subtype[na.omit(match(colnames(mergedNoNorm),subtypeDF$sampleName))]


table(subtypes[kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]]])
table(subtypes[kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]]])
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
ggplot(data=attrDF_pam50Full,aes(subtype))+geom_bar() + facet_grid(.~clustNum)

#########
#now with BMC normalization
dataMatrixList <- readRDS(paste0(saveDirGlobal,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.rds"))
output <- merge_datasetList(datasetList=dataMatrixList,minNumGenes = 10000, minNumPatients = 40,batchNormalize = c('BMC'));
saveRDS(output,file=paste0(saveDirMerged,"/mergedExprMatrix_minVar001_17_studies_BMC_norm.rds"),compress=TRUE)


load("pam50Short_genes.RData")

dataMatrixList <- list(mergedBMC=output$mergedExprMatrix)

clustFeaturesList <- list(pam50Short=pam50Short)

#we know these are strong clusters. have  minMeanClustConsensus=.8
kmeansConsensus <-  clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                           pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                           numSims=500,maxNumClusters=20,
                                           outputFile=outputFile,distMethod=c("euclidean"),
                                           hclustAlgorithm=c("average"),
                                           consensusHclustAlgorithm=c("average"),
                                           minClustConsensus=.7, minMeanClustConsensus=.85,
                                           corUse="everything",pItem=.9,maxPAC=.15)


saveRDS(kmeansConsensus,file=paste0(saveDirMerged, "/curatedBreastData_kmeansConsensus_mergedBMCNorm_pam50Short_nstart1",".rds"),compress=TRUE)

#########
#combat
###did not detect a batch effect!! p-value was .055....so try .1:
output <- merge_datasetList(datasetList=dataMatrixList,minNumGenes = 10000, minNumPatients = 40,batchNormalize = c('combat'),combatPvalueThresh=.1);
save(output,file=paste0(saveDirMerged, "/home/kplaney/breast_analysis//mergedExprMatrix_minVar001_17_studies_combat_norm.rds",compress=TRUE)
message("Done with Combat")

load("pam50Short_genes.RData")

dataMatrixList <- list(mergedCombat=output$mergedExprMatrix)

clustFeaturesList <- list(pam50Short=pam50Short)

#we know these are strong clusters. have  minMeanClustConsensus=.8
kmeansConsensus <-  clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                           pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                           numSims=500,maxNumClusters=20,
                                           outputFile=outputFile,distMethod=c("euclidean"),
                                           hclustAlgorithm=c("average"),
                                           consensusHclustAlgorithm=c("average"),
                                           minClustConsensus=.7, minMeanClustConsensus=.85,
                                           corUse="everything",pItem=.9,maxPAC=.15)


saveRDS(kmeansConsensus,file=paste0(saveDirMerged, "/curatedBreastData_kmeansConsensus_mergedCombatNorm_pam50Short_nstart1",".rds"),compress=TRUE)


####Analyses/inspections:
##first off, let's just see how the actual PAM50 subtypings differ when we applied 
#these batch normalizations:

# combat must be run on intersecting genes only; so limited pam50 gene set is a confounding factor.
#but still, if subtypes changes between no norm and combat with limited pam50 gene set, this is still an issue:
##also: look at effects of BMC normalization:
load(paste0(saveDirGlobal,"/mergedExprMatrix_minVar001_17_studies_BMC_norm.RData.gzip"))
mergedBMC <- output$mergedExprMatrix
load(paste0(saveDirGlobal,"/mergedExprMatrix_minVar001_17_studies_no_norm.RData.gzip"))
mergedNoNorm <- output$mergedExprMatrix

Pam50_subtypes_noNorm <- assignCentroidSubtype(t(mergedNoNorm),minNumGenes=30,centroidRData="pam50_centroids_updatedSymbols.RData");
Pam50_subtypes_BMCNorm <- assignCentroidSubtype(t(mergedBMC),minNumGenes=30,centroidRData="pam50_centroids_updatedSymbols.RData");
#wow...a lot changed!!! of course, if all pam50 genes were in all datasets, these numbers may slightly changed, but unlikely by much.
length(which(Pam50_subtypes_noNorm$subtype[,"subtypeLabels"] != Pam50_subtypes_BMCNorm$subtypes[,"subtypeLabels"]))

load(paste0(saveDirGlobal,"/mergedExprMatrix_minVar001_17_studies_combat_norm.RData.gzip"))
mergedCombat <- output$mergedExprMatrix
Pam50_subtypes_combatNorm <- assignCentroidSubtype(t(mergedCombat),minNumGenes=30,centroidRData="pam50_centroids_updatedSymbols.RData");

#how changed from no norm using the restricted gene set?
#indices matched up?
all(colnames(mergedCombat)==colnames(mergedNoNorm))
all(colnames(mergedBMC)==colnames(mergedNoNorm))
#starting with 2237 samples. this is what we're comparing to with our analysis.
#no norm here is the restricted PAM50 subtypes
length(Pam50_subtypes_noNorm$subtypes[,"subtypeLabels"])
#combat changed 673 samples
length(which(Pam50_subtypes_combatNorm$subtypes[,"subtypeLabels"]!=Pam50_subtypes_noNorm$subtypes[,"subtypeLabels"]))
#changed 821
length(which(Pam50_subtypes_BMCNorm$subtypes[,"subtypeLabels"]!=Pam50_subtypes_noNorm$subtypes[,"subtypeLabels"]))

#what about between ComBat and BMC?
#changed 848
length(which(Pam50_subtypes_BMCNorm$subtypes[,"subtypeLabels"]!=Pam50_subtypes_combatNorm$subtypes[,"subtypeLabels"]))


#load up original full subtypes with unrestricted gene set too.
load(paste0(saveDirGlobal, ("pam50FullAndShort_subtypeDF.RData.gzip"))
#short vs full difference? 309
length(which(subtypeDF_master$subtype!=Pam50_subtypes_noNorm$subtypes[,"subtypeLabels"]))
#vs combat? 723 have changed now
length(which(subtypeDF_master$subtype!=Pam50_subtypes_combatNorm$subtypes[,"subtypeLabels"]))
#vs BMC? 893 now
length(which(subtypeDF_master$subtype!=Pam50_subtypes_BMCNorm$subtypes[,"subtypeLabels"]))


#indices already match up
all(subtypeDF_master$sampleName==colnames(mergedNoNorm))
#and our current pam50 short subtypings match up with this older run:
all(subtypeDF_master$subtype_short==Pam50_subtypes_noNorm$subtypes[,"subtypeLabels"])
pam50ShortSubtypeDF <- cbind(colnames(mergedNoNorm),
                             subtypeDF_master$studyNum,
                             subtypeDF_master$subtype,
                             Pam50_subtypes_noNorm$subtypes[,"subtypeLabels"],
                             Pam50_subtypes_BMCNorm$subtypes[,"subtypeLabels"],
                             Pam50_subtypes_combatNorm$subtypes[,"subtypeLabels"])

colnames(pam50ShortSubtypeDF) <- c("sampleName",
                                   "datasetNum","noNorm_fullGeneSetAllowed",
                                   "noNorm_35genes","BMC_norm_35genes","combat_norm_35genes")
save(pam50ShortSubtypeDF,file=paste0(saveDirMerged,"/mergedPam50SubtypeMatrix.RData.gzip"),compress="gzip")
write.table(pam50ShortSubtypeDF,file=paste0(saveDirMerged,"/mergedPam50SubtypeMatrix.txt"),quote=FALSE,
                                            col.names=TRUE,row.names=FALSE)
   
textOut <- capture.output(table(subtypeDF_master$subtype))
cat(textOut,sep="\n",file=paste0(saveDirMerged,"/noNormPam50Full_subtypeTable.txt"))

textOut <- capture.output(table(subtypeDF_master$subtype_short))
cat(textOut,sep="\n",file=paste0(saveDirMerged,"/noNormPam50Short_subtypeTable.txt"))

textOut <- capture.output(table(Pam50_subtypes_BMCNorm$subtypes[,"subtypeLabels"]))
cat(textOut,sep="\n",file=paste0(saveDirMerged,"/mergedBMCNormPam50Short_subtypeTable.txt"))

textOut <- capture.output(table(Pam50_subtypes_combatNorm$subtypes[,"subtypeLabels"]))
cat(textOut,sep="\n",file=paste0(saveDirMerged,"/mergedcombatNormPam50Short_subtypeTable.txt"))
####

####now: semi-supervised centroid analyses for all 3 batch methods PLUS version that isn't shortened.
###i.e. if we treated the post-batch normalized PAM50 subtypes as clusters, 
#how well would this split patients by survival?
###not shortened.
#indices do match up
patientIndices <- data.frame(as.numeric(subtypeDF_master$subtypeNum),
                       subtypeDF_master$sampleName)


colnames(patientIndices) <- c("community","sampleName")
                       
esets <- readRDS(paste0(saveDirGlobal, "curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.rds"))


phenoMasterDF <- createPhenoMasterTableFromMatrixList(esetList=esets,sampleKeyColName="dbUniquePatientID")
library("plyr")
sampleClustCommPhenoData <- join(patientIndices,phenoMasterDF,by="sampleName",type="left",match="all")

groupingTerm="community"

groupings <- as.factor(sampleClustCommPhenoData[, groupingTerm])

expName <- "noNorm_PAM50centroidsFULL_survival"

ROC_output <- runBreastCancerBinarySurvModels(sampleClustCommPhenoData,expName,saveDirMerged)
saveRDS(ROC_output,compress=TRUE,file=paste0(saveDirMerged,"/",expName,"_ROC_stats.rds"))



###no norm/shortened.
patientIndices <- data.frame(as.numeric(Pam50_subtypes_noNorm$subtypes[,"subtypes"]),
                       colnames(mergedNoNorm))


colnames(patientIndices) <- c("community","sampleName")
                       
esets <- readRDS(paste0(saveDirGlobal, "curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.rds"))


phenoMasterDF <- createPhenoMasterTableFromMatrixList(esetList=esets,sampleKeyColName="dbUniquePatientID")
library("plyr")
sampleClustCommPhenoData <- join(patientIndices,phenoMasterDF,by="sampleName",type="left",match="all")

groupingTerm="community"

groupings <- as.factor(sampleClustCommPhenoData[, groupingTerm])

expName <- "noNorm_PAM50centroids_survival"

ROC_output <- runBreastCancerBinarySurvModels(sampleClustCommPhenoData,expName,saveDirMerged)
saveRDS(ROC_output,compress=TRUE,file=paste0(saveDirMerged,"/",expName,"_ROC_stats.rds"))


###BMC
patientIndices <- data.frame(as.numeric(Pam50_subtypes_BMCNorm$subtypes[,"subtypes"]),
                       colnames(mergedBMC))


colnames(patientIndices) <- c("community","sampleName")
                       
esets <- readRDS(paste0(saveDirGlobal, "curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.rds"))

phenoMasterDF <- createPhenoMasterTableFromMatrixList(esetList=esets,sampleKeyColName="dbUniquePatientID")
library("plyr")
sampleClustCommPhenoData <- join(patientIndices,phenoMasterDF,by="sampleName",type="left",match="all")

groupingTerm="community"

groupings <- as.factor(sampleClustCommPhenoData[, groupingTerm])

expName <- "BMC_PAM50centroids_survival"

ROC_output <- runBreastCancerBinarySurvModels(sampleClustCommPhenoData,expName,saveDirMerged)
saveRDS(ROC_output,compress=TRUE,file=paste0(saveDirMerged,"/",expName,"_ROC_stats.rds"))

##Combat
patientIndices <- data.frame(as.numeric(Pam50_subtypes_combatNorm$subtypes[,"subtypes"]),
                       colnames(mergedCombat))


colnames(patientIndices) <- c("community","sampleName")
                       
esets <- readRDS(paste0(saveDirGlobal, "curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.rds"))

phenoMasterDF <- createPhenoMasterTableFromMatrixList(esetList=esets,sampleKeyColName="dbUniquePatientID")

sampleClustCommPhenoData <- join(patientIndices,phenoMasterDF,by="sampleName",type="left",match="all")

groupingTerm="community"

groupings <- as.factor(sampleClustCommPhenoData[, groupingTerm])

expName <- "Combat_PAM50centroids_survival"

ROC_output <- runBreastCancerBinarySurvModels(sampleClustCommPhenoData,expName,saveDirMerged)
saveRDS(ROC_output,compress=TRUE,file=paste0(saveDirMerged,"/",expName,"_ROC_stats.rds"))

########noNorm plotting for merged consensus clustering results
kmeansConsensus <- readRDS(paste0(saveDirMerged, "/curatedBreastData_kmeansConsensus_mergedNoNorm_pam50Short_nstart1.rds"))

table(Pam50_subtypes_noNorm$subtypes[kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]] ,"subtypeLabels"])
table(Pam50_subtypes_noNorm$subtypes[kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]] ,"subtypeLabels"])
#and this cluster is super tiny!
table(Pam50_subtypes_noNorm$subtype[kmeansConsensus$clustSampleIndexList_PACR[[1]][[3]] ,"subtypeLabels"])

attrDF <- data.frame(c(rep.int(1,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]])),
                       rep.int(2,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]])),
                       rep.int(3,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[3]]))
),
c(Pam50_subtypes_noNorm$subtypes[kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]] ,"subtypeLabels"],
  Pam50_subtypes_noNorm$subtypes[kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]] ,"subtypeLabels"],
  Pam50_subtypes_noNorm$subtypes[kmeansConsensus$clustSampleIndexList_PACR[[1]][[3]] ,"subtypeLabels"]))
library("ggplot2")
colnames(attrDF) <- c("clustNum","subtype")

variableColorMatrix <- c("#970eec","#ec0e58","#0e58ec","#7d7d7d","#0eec97")
names(variableColorMatrix) <- c("LumB","LumA","Her2","Normal","Basal")

plotG <-    ggplot(attrDF,aes(factor(subtype),fill=factor(subtype)))+geom_bar() + facet_grid(.~clustNum,scales="free_x")+
  labs(y="Number of samples", fill=paste0("pam50 subtype"),
       title=paste0("Pam50 subtype with 35 gene list",
                    " by cluster for merged \nno dataset normalization and k-means"))+
  theme(axis.text.x = element_text(colour = "black",size=8,angle=45,vjust=1,hjust=1),
        axis.title=element_blank())+
  theme(legend.title=element_text(size=12))+ scale_fill_manual(values = variableColorMatrix)+  theme(strip.text.x = element_text(size = 17, colour = "black",face="bold"))+
  theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
  theme(axis.text.x = element_text(colour = "black",size=12,angle=45,vjust=1,hjust=1),axis.title.x= element_blank(),
        axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1))+
  theme(panel.grid.major = element_line(colour = 0),panel.grid.minor = element_line(colour = 0))+
  theme(plot.title=element_text(colour="black",size=12,vjust=1)) + coord_cartesian(ylim=c(0,2000))
#above: font size needs to be smaller for subtypes
png(filename=paste0(saveDirMerged,"/mergedNoNorm_subtype_breakdowns_",".png"),width=1000,height=1000,res=160)

plot( plotG)

dev.off()

#stacked
plotG <-    ggplot(attrDF,aes(factor(clustNum),fill=factor(subtype)),scales="free_x")+geom_bar() + facet_grid(.~clustNum,scales="free_x")+
  labs(y="Number of samples", fill=paste0("pam50 subtype"),
       title=paste0("Pam50 subtype with 35 gene list",
                    " by cluster for merged \nno dataset normalization and k-means"))+
  theme(axis.text.x = element_text(colour = "black",size=8,angle=45,vjust=1,hjust=1),
        axis.title=element_blank())+
  theme(legend.title=element_text(size=12))+ scale_fill_manual(values = variableColorMatrix)+  theme(strip.text.x = element_text(size = 17, colour = "black",face="bold"))+
  theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
  theme(axis.text.x = element_text(colour = "black",size=12,angle=45,vjust=1,hjust=1),axis.title.x= element_blank(),
        axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1))+
  theme(panel.grid.major = element_line(colour = 0),panel.grid.minor = element_line(colour = 0))+
  theme(plot.title=element_text(colour="black",size=12,vjust=1)) + coord_cartesian(ylim=c(0,2000))
#above: font size needs to be smaller for subtypes
png(filename=paste0(saveDirMerged,"/mergedNoNorm_subtype_breakdowns_stacked_",".png"),width=1000,height=1000,res=160)

plot( plotG)

dev.off()

##AUC data
#dataset indices still the same
patientIndices <- data.frame(c(rep.int(1,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]])),
                       rep.int(2,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]])),
                       rep.int(3,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[3]]))),
                       colnames(mergedNoNorm))


colnames(patientIndices) <- c("community","sampleName")
                       
esets <- readRDS(paste0(saveDirGlobal, "curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.rds"))

phenoMasterDF <- createPhenoMasterTableFromMatrixList(esetList=esets,sampleKeyColName="dbUniquePatientID")
library("plyr")
sampleClustCommPhenoData <- join(patientIndices,phenoMasterDF,by="sampleName",type="left",match="all")

groupingTerm="community"

groupings <- as.factor(sampleClustCommPhenoData[, groupingTerm])


expName <- "noNorm_survival"
ROC_output <- runBreastCancerBinarySurvModels(sampleClustCommPhenoData,expName,saveDirMerged)
saveRDS(ROC_output,compress=TRUE,file=paste0(saveDirMerged,"/",expName,"_ROC_stats.rds"))


########combat plotting
kmeansConsensus <- readRDS(paste0(saveDirMerged, "/curatedBreastData_kmeansConsensus_mergedCombatNorm_pam50Short_nstart1.rds"))


attrDF <- data.frame(c(rep.int(1,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]])),
                       rep.int(2,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]])),
                       rep.int(3,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[3]]))
),
c(Pam50_subtypes_combatNorm$subtypes[kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]] ,"subtypeLabels"],
  Pam50_subtypes_combatNorm$subtypes[kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]] ,"subtypeLabels"],
  Pam50_subtypes_combatNorm$subtypes[kmeansConsensus$clustSampleIndexList_PACR[[1]][[3]] ,"subtypeLabels"]))

colnames(attrDF) <- c("clustNum","subtype")
variableColorMatrix <- c("#970eec","#ec0e58","#0e58ec","#7d7d7d","#0eec97")
names(variableColorMatrix) <- c("LumB","LumA","Her2","Normal","Basal")

plotG <-    ggplot(attrDF,aes(factor(subtype),fill=factor(subtype)))+geom_bar() + facet_grid(.~clustNum,scales="free_x")+
  labs(y="Number of samples", fill=paste0("pam50 subtype"),
       title=paste0("Pam50 subtype with 35 gene list",
                    " by cluster for merged \ncombat dataset normalization and k-means"))+
  theme(axis.text.x = element_text(colour = "black",size=8,angle=45,vjust=1,hjust=1),
        axis.title=element_blank())+
  theme(legend.title=element_text(size=12))+ scale_fill_manual(values = variableColorMatrix)+  theme(strip.text.x = element_text(size = 17, colour = "black",face="bold"))+
  theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
  theme(axis.text.x = element_text(colour = "black",size=12,angle=45,vjust=1,hjust=1),axis.title.x= element_blank(),
        axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1))+
  theme(panel.grid.major = element_line(colour = 0),panel.grid.minor = element_line(colour = 0))+
  theme(plot.title=element_text(colour="black",size=12,vjust=1)) +coord_cartesian(ylim=c(0,2000))
#above: font size needs to be smaller for subtypes
png(filename=paste0(saveDirMerged,"/mergedCombat_subtype_breakdowns_",".png"),width=1000,height=1000,res=160)

plot( plotG)

dev.off()

plotG <-    ggplot(attrDF,aes(factor(clustNum),fill=factor(subtype)),scales="free_x")+geom_bar() + facet_grid(.~clustNum,scales="free_x")+
  labs(y="Number of samples", fill=paste0("pam50 subtype"),
       title=paste0("Pam50 subtype with 35 gene list",
                    " by cluster for merged \ncombat normalization and k-means"))+
  theme(axis.text.x = element_text(colour = "black",size=8,angle=45,vjust=1,hjust=1),
        axis.title=element_blank())+
  theme(legend.title=element_text(size=12))+ scale_fill_manual(values = variableColorMatrix)+  theme(strip.text.x = element_text(size = 17, colour = "black",face="bold"))+
  theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
  theme(axis.text.x = element_text(colour = "black",size=12,angle=45,vjust=1,hjust=1),axis.title.x= element_blank(),
        axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1))+
  theme(panel.grid.major = element_line(colour = 0),panel.grid.minor = element_line(colour = 0))+
  theme(plot.title=element_text(colour="black",size=12,vjust=1)) + coord_cartesian(ylim=c(0,2000))

png(filename=paste0(saveDirMerged,"/mergedCombat_subtype_breakdowns_stacked_",".png"),width=1000,height=1000,res=160)

plot( plotG)

dev.off()
###ALSO: with original non-transformed PAM50 statuses
kmeansConsensus <- readRDS(paste0(saveDirMerged, "/curatedBreastData_kmeansConsensus_mergedCombatNorm_pam50Short_nstart1.rds")
#do patient names all match up? yup
all(colnames(mergedCombat)==colnames(mergedNoNorm))

attrDF <- data.frame(c(rep.int(1,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]])),
                       rep.int(2,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]])),
                       rep.int(3,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[3]]))
),
c(Pam50_subtypes_noNorm$subtypes[kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]] ,"subtypeLabels"],
  Pam50_subtypes_noNorm$subtypes[kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]] ,"subtypeLabels"],
  Pam50_subtypes_noNorm$subtypes[kmeansConsensus$clustSampleIndexList_PACR[[1]][[3]] ,"subtypeLabels"]))

colnames(attrDF) <- c("clustNum","subtype")
variableColorMatrix <- c("#970eec","#ec0e58","#0e58ec","#7d7d7d","#0eec97")
names(variableColorMatrix) <- c("LumB","LumA","Her2","Normal","Basal")

plotG <-    ggplot(attrDF,aes(factor(subtype),fill=factor(subtype)))+geom_bar() + facet_grid(.~clustNum,scales="free_x")+
  labs(y="Number of samples", fill=paste0("pam50 subtype"),
       title=paste0("Pam50 subtype with 35 gene list",
                    " by cluster for merged \ncombat dataset normalization and k-means"))+
  theme(axis.text.x = element_text(colour = "black",size=8,angle=45,vjust=1,hjust=1),
        axis.title=element_blank())+
  theme(legend.title=element_text(size=12))+ scale_fill_manual(values = variableColorMatrix)+  theme(strip.text.x = element_text(size = 17, colour = "black",face="bold"))+
  theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
  theme(axis.text.x = element_text(colour = "black",size=12,angle=45,vjust=1,hjust=1),axis.title.x= element_blank(),
        axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1))+
  theme(panel.grid.major = element_line(colour = 0),panel.grid.minor = element_line(colour = 0))+
  theme(plot.title=element_text(colour="black",size=12,vjust=1)) +coord_cartesian(ylim=c(0,2000))
#above: font size needs to be smaller for subtypes
png(filename=paste0(saveDirMerged,"/mergedCombat_subtype_breakdowns_assignedBEFOREcombat_",".png"),width=1000,height=1000,res=160)

plot( plotG)

dev.off()

plotG <-    ggplot(attrDF,aes(factor(clustNum),fill=factor(subtype)),scales="free_x")+geom_bar() + facet_grid(.~clustNum,scales="free_x")+
  labs(y="Number of samples", fill=paste0("pam50 subtype"),
       title=paste0("Pam50 subtype with 35 gene list",
                    " by cluster for merged \ncombat normalization and k-means"))+
  theme(axis.text.x = element_text(colour = "black",size=8,angle=45,vjust=1,hjust=1),
        axis.title=element_blank())+
  theme(legend.title=element_text(size=12))+ scale_fill_manual(values = variableColorMatrix)+  theme(strip.text.x = element_text(size = 17, colour = "black",face="bold"))+
  theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
  theme(axis.text.x = element_text(colour = "black",size=12,angle=45,vjust=1,hjust=1),axis.title.x= element_blank(),
        axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1))+
  theme(panel.grid.major = element_line(colour = 0),panel.grid.minor = element_line(colour = 0))+
  theme(plot.title=element_text(colour="black",size=12,vjust=1)) + coord_cartesian(ylim=c(0,2000))

png(filename=paste0(saveDirMerged,"/mergedCombat_subtype_breakdowns_assignedBEFOREcombat_stacked_",".png"),width=1000,height=1000,res=160)

plot( plotG)

dev.off()


###survival code
#3 clusters
#dataset indices still the same
patientIndices <- data.frame(c(rep.int(1,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]])),
                       rep.int(2,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]])),
                       rep.int(3,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[3]]))),
                       colnames(mergedCombat))


colnames(patientIndices) <- c("community","sampleName")
                       
esets <- readRDS(paste0(saveDirGlobal, "curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.rds"))

phenoMasterDF <- createPhenoMasterTableFromMatrixList(esetList=esets,sampleKeyColName="dbUniquePatientID")
library("plyr")
sampleClustCommPhenoData <- join(patientIndices,phenoMasterDF,by="sampleName",type="left",match="all")

groupingTerm="community"

groupings <- as.factor(sampleClustCommPhenoData[, groupingTerm])

expName <- "combat_survival"

ROC_output <- runBreastCancerBinarySurvModels(sampleClustCommPhenoData,expName,saveDirMerged)
saveRDS(ROC_output,compress=TRUE,file=paste0(saveDirMerged,"/",expName,"_ROC_stats.rds"))

####BMC norm plotting
kmeansConsensus <- readRDS(paste0(saveDirMerged, "/curatedBreastData_kmeansConsensus_mergedBMCNorm_pam50Short_nstart1.rds")

#cluster indices match up with merged mattrix; this is what we fed in the clustering algorithm
#only 2 clusters!
attrDF <- data.frame(c(rep.int(1,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]])),
                       rep.int(2,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]]))
                       ),
c(Pam50_subtypes_BMCNorm$subtypes[kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]] ,"subtypeLabels"],
  Pam50_subtypes_BMCNorm$subtypes[kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]] ,"subtypeLabels"]
  ))

colnames(attrDF) <- c("clustNum","subtype")
variableColorMatrix <- c("#970eec","#ec0e58","#0e58ec","#7d7d7d","#0eec97")
names(variableColorMatrix) <- c("LumB","LumA","Her2","Normal","Basal")

plotG <-    ggplot(attrDF,aes(factor(subtype),fill=factor(subtype)))+geom_bar() + facet_grid(.~clustNum,scales="free_x")+
  labs(y="Number of samples", fill=paste0("pam50 subtype"),
       title=paste0("Pam50 subtype with 35 gene list",
                    " by cluster for merged \nBMC dataset normalization and k-means"))+
  theme(axis.text.x = element_text(colour = "black",size=8,angle=45,vjust=1,hjust=1),
        axis.title=element_blank())+
  theme(legend.title=element_text(size=12))+ scale_fill_manual(values = variableColorMatrix)+  theme(strip.text.x = element_text(size = 17, colour = "black",face="bold"))+
  theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
  theme(axis.text.x = element_text(colour = "black",size=12,angle=45,vjust=1,hjust=1),axis.title.x= element_blank(),
        axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1))+
  theme(panel.grid.major = element_line(colour = 0),panel.grid.minor = element_line(colour = 0))+
  theme(plot.title=element_text(colour="black",size=12,vjust=1)) +coord_cartesian(ylim=c(0,2000))
#above: font size needs to be smaller for subtypes
png(filename=paste0(saveDirMerged,"/mergedBMC_subtype_breakdowns_",".png"),width=1000,height=1000,res=160)

plot( plotG)

dev.off()

plotG <-    ggplot(attrDF,aes(factor(clustNum),fill=factor(subtype)),scales="free_x")+geom_bar() + facet_grid(.~clustNum,scales="free_x")+
  labs(y="Number of samples", fill=paste0("pam50 subtype"),
       title=paste0("Pam50 subtype with 35 gene list",
                    " by cluster for merged \nBMC normalization and k-means"))+
  theme(axis.text.x = element_text(colour = "black",size=8,angle=45,vjust=1,hjust=1),
        axis.title=element_blank())+
  theme(legend.title=element_text(size=12))+ scale_fill_manual(values = variableColorMatrix)+  theme(strip.text.x = element_text(size = 17, colour = "black",face="bold"))+
  theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
  theme(axis.text.x = element_text(colour = "black",size=12,angle=45,vjust=1,hjust=1),axis.title.x= element_blank(),
        axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1))+
  theme(panel.grid.major = element_line(colour = 0),panel.grid.minor = element_line(colour = 0))+
  theme(plot.title=element_text(colour="black",size=12,vjust=1)) + coord_cartesian(ylim=c(0,2000))

png(filename=paste0(saveDirMerged,"/mergedBMC_subtype_breakdowns_stacked_",".png"),width=1000,height=1000,res=160)

plot( plotG)

dev.off()
####ALSO: with original non-transformed PAM50 statuses
#do patient names name up? yup
all(colnames(mergedBMC)==colnames(mergedNoNorm))

kmeansConsensus <- readRDS(paste0(saveDirMerged, "/curatedBreastData_kmeansConsensus_mergedBMCNorm_pam50Short_nstart1.rds"))

#cluster indices match up with merged mattrix; this is what we fed in the clustering algorithm
#only 2 clusters!
attrDF <- data.frame(c(rep.int(1,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]])),
                       rep.int(2,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]]))
),
c(Pam50_subtypes_noNorm$subtypes[kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]] ,"subtypeLabels"],
  Pam50_subtypes_noNorm$subtypes[kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]] ,"subtypeLabels"]
))

colnames(attrDF) <- c("clustNum","subtype")
variableColorMatrix <- c("#970eec","#ec0e58","#0e58ec","#7d7d7d","#0eec97")
names(variableColorMatrix) <- c("LumB","LumA","Her2","Normal","Basal")

plotG <-    ggplot(attrDF,aes(factor(subtype),fill=factor(subtype)))+geom_bar() + facet_grid(.~clustNum,scales="free_x")+
  labs(y="Number of samples", fill=paste0("pam50 subtype"),
       title=paste0("Pam50 subtype with 35 gene list",
                    " by cluster for merged \nBMC dataset normalization and k-means"))+
  theme(axis.text.x = element_text(colour = "black",size=8,angle=45,vjust=1,hjust=1),
        axis.title=element_blank())+
  theme(legend.title=element_text(size=12))+ scale_fill_manual(values = variableColorMatrix)+  theme(strip.text.x = element_text(size = 17, colour = "black",face="bold"))+
  theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
  theme(axis.text.x = element_text(colour = "black",size=12,angle=45,vjust=1,hjust=1),axis.title.x= element_blank(),
        axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1))+
  theme(panel.grid.major = element_line(colour = 0),panel.grid.minor = element_line(colour = 0))+
  theme(plot.title=element_text(colour="black",size=12,vjust=1)) +coord_cartesian(ylim=c(0,2000))
#above: font size needs to be smaller for subtypes
png(filename=paste0(saveDirMerged,"/mergedBMC_subtype_breakdowns_assignedBEFOREBMC_",".png"),width=1000,height=1000,res=160)

plot( plotG)

dev.off()
#stacked
plotG <-    ggplot(attrDF,aes(factor(clustNum),fill=factor(subtype)),scales="free_x")+geom_bar() + facet_grid(.~clustNum,scales="free_x")+
  labs(y="Number of samples", fill=paste0("pam50 subtype"),
       title=paste0("Pam50 subtype with 35 gene list",
                    " by cluster for merged \nBMC normalization and k-means"))+
  theme(axis.text.x = element_text(colour = "black",size=8,angle=45,vjust=1,hjust=1),
        axis.title=element_blank())+
  theme(legend.title=element_text(size=12))+ scale_fill_manual(values = variableColorMatrix)+  theme(strip.text.x = element_text(size = 17, colour = "black",face="bold"))+
  theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
  theme(axis.text.x = element_text(colour = "black",size=12,angle=45,vjust=1,hjust=1),axis.title.x= element_blank(),
        axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1))+
  theme(panel.grid.major = element_line(colour = 0),panel.grid.minor = element_line(colour = 0))+
  theme(plot.title=element_text(colour="black",size=12,vjust=1)) + coord_cartesian(ylim=c(0,2000))

png(filename=paste0(saveDirMerged,"/mergedBMC_subtype_breakdowns_assignedBEFOREBMC_stacked_",".png"),width=1000,height=1000,res=160)

plot( plotG)

dev.off()



##gut check: got all patients? yes
length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[3]]) + length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]]) + length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]])
###how many of GSE32646 are in each combat cluster? they have an "809" tag
length(grep("809",names(kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]])))
#two are not from GSE32646 that contain 809:
names(kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]])[grep("809",names(kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]]))]
length(grep("809",names(kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]])))
names(kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]])[grep("809",names(kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]]))]
length(grep("809",names(kmeansConsensus$clustSampleIndexList_PACR[[1]][[3]])))
names(kmeansConsensus$clustSampleIndexList_PACR[[1]][[3]])[grep("809",names(kmeansConsensus$clustSampleIndexList_PACR[[1]][[3]]))]

#survival data
#2 clusters
#dataset indices still the same
patientIndices <- data.frame(c(rep.int(1,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]])),
                       rep.int(2,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]]))),
                       colnames(mergedBMC))


colnames(patientIndices) <- c("community","sampleName")

esets <- readRDS(paste0(saveDirGlobal, "curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.rds"))

phenoMasterDF <- createPhenoMasterTableFromMatrixList(esetList=esets,sampleKeyColName="dbUniquePatientID")
library("plyr")
sampleClustCommPhenoData <- join(patientIndices,phenoMasterDF,by="sampleName",type="left",match="all")

groupingTerm="community"

groupings <- as.factor(sampleClustCommPhenoData[, groupingTerm])


expName <- "BMC_survival"
ROC_output <- runBreastCancerBinarySurvModels(sampleClustCommPhenoData,expName,saveDirMerged)
saveRDS(ROC_output,compress=TRUE,file=paste0(saveDirMerged,"/",expName,"_ROC_stats.rds"))

