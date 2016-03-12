###look at Pam50 breakdown
library("Coincide")

#this script follows from: breastProcessAndGeneFeature_script.R

#other libraries called later down: plyr, ggplot2
saveDirMerged <- "/home/ywrfc09/breast_analysis/mergedMatrix/"
##assign subtypes
saveDirGlobal <- "/home/ywrfc09/breast_analysis"
# combat must be run on intersecting genes only; so limited pam50 gene set is a confounding factor.
#but still, if subtypes changes between no norm and combat with limited pam50 gene set, this is still an issue:
##also: look at effects of BMC normalization:
load(paste0(saveDirGlobal,"/mergedExprMatrix_minVar001_17_studies_BMC_norm.RData.gzip"))
mergedBMC <- output$mergedExprMatrix
load(paste0(saveDirGlobal,"/mergedExprMatrix_minVar001_17_studies_no_norm.RData.gzip"))
mergedNoNorm <- output$mergedExprMatrix

Pam50_subtypes_noNorm <- assignCentroidSubtype(t(mergedNoNorm),minNumGenes=30,centroidRData="/home/ywrfc09/breast_analysis/PAM50_analyses/pam50_centroids_updatedSymbols.RData");
Pam50_subtypes_BMCNorm <- assignCentroidSubtype(t(mergedBMC),minNumGenes=30,centroidRData="/home/ywrfc09/breast_analysis/PAM50_analyses/pam50_centroids_updatedSymbols.RData");
#wow...a lot changed!!! of course, if all pam50 genes were in all datasets, these numbers may slightly changed, but unlikely by much.
length(which(Pam50_subtypes_noNorm$subtype[,"subtypeLabels"] != Pam50_subtypes_BMCNorm$subtypes[,"subtypeLabels"]))

load(paste0(saveDirGlobal,"/mergedExprMatrix_minVar001_17_studies_combat_norm.RData.gzip"))
mergedCombat <- output$mergedExprMatrix
Pam50_subtypes_combatNorm <- assignCentroidSubtype(t(mergedCombat),minNumGenes=30,centroidRData="/home/ywrfc09/breast_analysis/PAM50_analyses/pam50_centroids_updatedSymbols.RData");

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
load("/home/ywrfc09/breast_analysis/PAM50_analyses/pam50FullAndShort_subtypeDF.RData.gzip")
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

####supervised centroid analyses for all 3 methods PLUS version that isn't shortened.

###not shortened.
#indices do match up
patientIndices <- data.frame(as.numeric(subtypeDF_master$subtypeNum),
                       subtypeDF_master$sampleName)


colnames(patientIndices) <- c("community","sampleName")
                       
load("/home/ywrfc09/breast_analysis/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
  esets=esets_minVar001_17_studies

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
                       
load("/home/ywrfc09/breast_analysis/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
  esets=esets_minVar001_17_studies

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
                       
load("/home/ywrfc09/breast_analysis/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
  esets=esets_minVar001_17_studies

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
                       
load("/home/ywrfc09/breast_analysis/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
  esets=esets_minVar001_17_studies

phenoMasterDF <- createPhenoMasterTableFromMatrixList(esetList=esets,sampleKeyColName="dbUniquePatientID")

sampleClustCommPhenoData <- join(patientIndices,phenoMasterDF,by="sampleName",type="left",match="all")

groupingTerm="community"

groupings <- as.factor(sampleClustCommPhenoData[, groupingTerm])

expName <- "Combat_PAM50centroids_survival"

ROC_output <- runBreastCancerBinarySurvModels(sampleClustCommPhenoData,expName,saveDirMerged)
saveRDS(ROC_output,compress=TRUE,file=paste0(saveDirMerged,"/",expName,"_ROC_stats.rds"))

########noNorm plotting
load("/home/ywrfc09/breast_analysis/mergedMatrix/curatedBreastData_kmeansConsensus_mergedNoNorm_pam50Short_nstart12015-05-04.RData.gzip")

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
png(filename=paste0(saveDirMerged,"/mergedNoNorm_subtype_breakdowns_",Sys.Date(),".png"),width=1000,height=1000,res=160)

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
png(filename=paste0(saveDirMerged,"/mergedNoNorm_subtype_breakdowns_stacked_",Sys.Date(),".png"),width=1000,height=1000,res=160)

plot( plotG)

dev.off()

##AUC data
#dataset indices still the same
patientIndices <- data.frame(c(rep.int(1,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]])),
                       rep.int(2,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]])),
                       rep.int(3,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[3]]))),
                       colnames(mergedNoNorm))


colnames(patientIndices) <- c("community","sampleName")
                       
load("/home/ywrfc09/breast_analysis/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
  esets=esets_minVar001_17_studies

phenoMasterDF <- createPhenoMasterTableFromMatrixList(esetList=esets,sampleKeyColName="dbUniquePatientID")
library("plyr")
sampleClustCommPhenoData <- join(patientIndices,phenoMasterDF,by="sampleName",type="left",match="all")

groupingTerm="community"

groupings <- as.factor(sampleClustCommPhenoData[, groupingTerm])


expName <- "noNorm_survival"
ROC_output <- runBreastCancerBinarySurvModels(sampleClustCommPhenoData,expName,saveDirMerged)
saveRDS(ROC_output,compress=TRUE,file=paste0(saveDirMerged,"/",expName,"_ROC_stats.rds"))


########combat plotting
load("/home/ywrfc09/breast_analysis/mergedMatrix/curatedBreastData_kmeansConsensus_mergedCombatNorm_pam50Short_nstart12015-05-04.RData.gzip")


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
png(filename=paste0(saveDirMerged,"/mergedCombat_subtype_breakdowns_",Sys.Date(),".png"),width=1000,height=1000,res=160)

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

png(filename=paste0(saveDirMerged,"/mergedCombat_subtype_breakdowns_stacked_",Sys.Date(),".png"),width=1000,height=1000,res=160)

plot( plotG)

dev.off()
###ALSO: with original non-transformed PAM50 statuses
load("/home/ywrfc09/breast_analysis/mergedMatrix/curatedBreastData_kmeansConsensus_mergedCombatNorm_pam50Short_nstart12015-05-04.RData.gzip")
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
png(filename=paste0(saveDirMerged,"/mergedCombat_subtype_breakdowns_assignedBEFOREcombat_",Sys.Date(),".png"),width=1000,height=1000,res=160)

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

png(filename=paste0(saveDirMerged,"/mergedCombat_subtype_breakdowns_assignedBEFOREcombat_stacked_",Sys.Date(),".png"),width=1000,height=1000,res=160)

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
                       
load("/home/ywrfc09/breast_analysis/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
  esets=esets_minVar001_17_studies

phenoMasterDF <- createPhenoMasterTableFromMatrixList(esetList=esets,sampleKeyColName="dbUniquePatientID")
library("plyr")
sampleClustCommPhenoData <- join(patientIndices,phenoMasterDF,by="sampleName",type="left",match="all")

groupingTerm="community"

groupings <- as.factor(sampleClustCommPhenoData[, groupingTerm])

expName <- "combat_survival"

ROC_output <- runBreastCancerBinarySurvModels(sampleClustCommPhenoData,expName,saveDirMerged)
saveRDS(ROC_output,compress=TRUE,file=paste0(saveDirMerged,"/",expName,"_ROC_stats.rds"))

####BMC norm plotting
load("/home/ywrfc09/breast_analysis/mergedMatrix//curatedBreastData_kmeansConsensus_mergedBMCNorm_pam50Short_nstart12015-05-04.RData.gzip")

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
png(filename=paste0(saveDirMerged,"/mergedBMC_subtype_breakdowns_",Sys.Date(),".png"),width=1000,height=1000,res=160)

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

png(filename=paste0(saveDirMerged,"/mergedBMC_subtype_breakdowns_stacked_",Sys.Date(),".png"),width=1000,height=1000,res=160)

plot( plotG)

dev.off()
####ALSO: with original non-transformed PAM50 statuses
#do patient names name up? yup
all(colnames(mergedBMC)==colnames(mergedNoNorm))


load("/home/kplaney/breast_analysis/curatedBreastData_kmeansConsensus_mergedBMCNorm_pam50Short_nstart12015-05-04.RData.gzip")

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
png(filename=paste0(saveDirMerged,"/mergedBMC_subtype_breakdowns_assignedBEFOREBMC_",Sys.Date(),".png"),width=1000,height=1000,res=160)

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

png(filename=paste0(saveDirMerged,"/mergedBMC_subtype_breakdowns_assignedBEFOREBMC_stacked_",Sys.Date(),".png"),width=1000,height=1000,res=160)

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
                       
load("/home/ywrfc09/breast_analysis/curatedBreastData_esets_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
  esets=esets_minVar001_17_studies

phenoMasterDF <- createPhenoMasterTableFromMatrixList(esetList=esets,sampleKeyColName="dbUniquePatientID")
library("plyr")
sampleClustCommPhenoData <- join(patientIndices,phenoMasterDF,by="sampleName",type="left",match="all")

groupingTerm="community"

groupings <- as.factor(sampleClustCommPhenoData[, groupingTerm])


expName <- "BMC_survival"
ROC_output <- runBreastCancerBinarySurvModels(sampleClustCommPhenoData,expName,saveDirMerged)
saveRDS(ROC_output,compress=TRUE,file=paste0(saveDirMerged,"/",expName,"_ROC_stats.rds"))

###survival/AUC
runBreastCancerBinarySurvModels <- function(sampleClustCommPhenoData,expName="test",saveDirMerged="./"){
  library("pROC")
  saveDirMerged <- paste0(saveDirMerged,"/",expName,"_",Sys.Date(),"/")
  system(paste0("mkdir ",saveDirMerged))
  options(bitmapType="cairo")
  plotToScreen <- FALSE
  summaryFile <- paste0(saveDirMerged,"/",expName, "_",Sys.Date(),"_summary.txt")
  experimentName <- expName
  linearModels <- list()
  linearModelsSumm <- list()
  ROC_list <- list()
   #message("running survival data on breast")
      message("Running RFS models")
      cat(gsub(" ","_",paste0("\nLinear analysis with only RFS:\n")),append=TRUE,file=summaryFile)
      
      dataMatrix <- sampleClustCommPhenoData[which(!is.na((sampleClustCommPhenoData$RFS))), ]
      cat(gsub(" ","_",paste0("Studies used:\n",unique(dataMatrix$studyNum),"\n")),append=TRUE,file=summaryFile)
      cat(gsub(" ","_",paste0("\nMeta-clusters used:\n",unique(dataMatrix$community),"\n")),append=TRUE,file=summaryFile)
      
      dataMatrix <- data.frame(dataMatrix)
      dataMatrix$community <- as.factor(dataMatrix$community)
      dataMatrix$RFS <- as.factor(dataMatrix$RFS)
      
      linearModels[["RFS_logitAlone"]] <- tryCatch(glm(RFS~community,data=dataMatrix,family=binomial(link="logit")),
                                                   error = function(e) {
                                                     return(NA)
                                                   }
      )
  
      
      if(any(!is.na(  linearModels[["RFS_logitAlone"]]))){
        ##ROC curves
        #above: font size needs to be smaller for subtypes
        predictor <- predict(linearModels[["RFS_logitAlone"]],type="response")
        
        png(filename=paste0(saveDirMerged,"/",experimentName,"_RFS_ROC_",Sys.Date(),".png"),width=1000,height=1000,res=200)
        roc_data <- roc(response=dataMatrix$RFS,predictor=predictor,plot=TRUE,na.rm=TRUE,auc=TRUE)
        dev.off()
        textOut <- capture.output(roc_data)
        cat("roc_info\n",textOut,sep="\n",append=TRUE,file=summaryFile)
        cat("\nAUC : ",roc_data$auc,"\n",append=TRUE,file=summaryFile)
        
        linearModelsSumm[["RFS_logitAlone"]] <- summary(linearModels[["RFS_logitAlone"]] ) 
        textOut <- capture.output( linearModelsSumm [["RFS_logitAlone"]]$coefficients)
        cat(textOut,sep="\n",append=TRUE,file=summaryFile)
        textOut <- capture.output( linearModelsSumm [["RFS_logitAlone"]])
        cat(textOut,sep="\n",append=TRUE,file=summaryFile)
        linearModelsSumm[["RFS_numPatients"]] <- nrow(dataMatrix)
        cat("\nNumber_of_patients_for_RFS__model: ",nrow(dataMatrix),"\n",append=TRUE,file=summaryFile)
        linearModelsSumm[["RFS_aloneMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
        
        ROC_list[["RFS_logitAlone"]] <- roc_data
        
      }else{
        cat("\nRFS model alone returned NA.\n")
        cat("\nRFS model alone returned NA.\n",append=TRUE,file=outputFile)
      }
      ##with therapies
      #interesting: community becomes insignificant if add in therapies now: this could be
      #because community is most likely correlated with treatments; it's certainly correlated with pam50 status.
      #low Rsquared.
      cat(gsub(" ","_",paste0("\nLinear analysis with RFS and treatment:\n")),append=TRUE,file=summaryFile)
      dataMatrix <- dataMatrix[which(!is.na(dataMatrix$chemotherapyClass)), ]
      dataMatrix <- dataMatrix[which(!is.na(dataMatrix$anti_estrogen)), ]
      dataMatrix <- dataMatrix[which(!is.na(dataMatrix$anti_HER2)), ]
      dataMatrix <- data.frame(dataMatrix)
      dataMatrix$community <- as.factor(dataMatrix$community)
      dataMatrix$RFS <- as.factor(dataMatrix$RFS)
      dataMatrix$hist_grade <- as.factor(dataMatrix$hist_grade)
      dataMatrix$chemotherapyClass <- as.factor(dataMatrix$chemotherapyClass)
      dataMatrix$anti_estrogen <- as.factor(dataMatrix$anti_estrogen)
      dataMatrix$anti_HER2 <- as.factor(dataMatrix$anti_HER2)
      cat(gsub(" ","_",paste0("\nMeta-clusters used:\n",unique(dataMatrix$community),"\n")),append=TRUE,file=summaryFile)
      
      linearModelsSumm [["RFS_withRx_numPatients"]] <- nrow(dataMatrix)
      cat("\nNumber_of_patients_for_RFS_with_Rx_model: ",nrow(dataMatrix),"\n",append=TRUE,file=summaryFile)
      
      if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))>1 
         && length(unique(dataMatrix$anti_HER2))>1 ){
        
        linearModels[["RFS_logitWithRx"]] <- tryCatch(glm(RFS~community+chemotherapyClass+anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
        
        #no community
        reducedModel <- tryCatch(glm(RFS~chemotherapyClass+anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
         linearModels[["RFS_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["RFS_logitWithRx"]],test="Chisq")
           
     
        #now try with  grade.    
        cat(gsub(" ","_",paste0("\nLinear analysis with RFS and histological grade :\n")),append=TRUE,file=summaryFile)
        dataMatrixHist <- dataMatrix[which(!is.na((dataMatrix$hist_grade))),]
 
        
        cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
        
        linearModelsSumm[["RFS_withGrade_numPatients"]] <- nrow(dataMatrixHist)
        cat("\nNumber_of_patients_for_RFS_with_grade_model: ",nrow(dataMatrixHist),"\n",append=TRUE,file=summaryFile)
        
        linearModels[["RFS_logitWithRxGrade"]] <- tryCatch( glm(RFS~community+chemotherapyClass+anti_estrogen+anti_HER2+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
                                                            error = function(e) {
                                                              return(NA)
                                                            }
        )
        
        
      }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))>1 
               && length(unique(dataMatrix$anti_HER2))==1 ){
        
        
        cat(gsub(" ","_",paste0("\nFor RFS analysis, all anti_HER2 values the same: ",unique(unique(dataMatrix$anti_HER2)),"\n")),
            append=TRUE,file=summaryFile)
        linearModels[["RFS_logitWithRx"]] <- tryCatch(glm(RFS~community+chemotherapyClass+anti_estrogen,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
                                                      
        )
        
               #no community
        reducedModel <- tryCatch(glm(RFS~chemotherapyClass+anti_estrogen,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
         linearModels[["RFS_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["RFS_logitWithRx"]],test="Chisq")
        
        #now try with  grade.    
        cat(gsub(" ","_",paste0("\nLinear analysis with RFS and histological grade :\n")),append=TRUE,file=summaryFile)
        dataMatrixHist <- dataMatrix[which(!is.na((dataMatrix$hist_grade))),]
 
        
        cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
        
        linearModelsSumm[["RFS_withGrade_numPatients"]] <- nrow(dataMatrixHist)
        cat("\nNumber_of_patients_for_RFS_with_grade_model: ",nrow(dataMatrixHist),"\n",append=TRUE,file=summaryFile)
        linearModels[["RFS_logitWithRxGrade"]] <- tryCatch( glm(RFS~community+chemotherapyClass+anti_estrogen+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
                                                            error = function(e) {
                                                              return(NA)
                                                            }
        )
        
        
        
      }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))==1 
               && length(unique(dataMatrix$anti_HER2))>1 ){
        
        cat(gsub(" ","_",paste0("\nFor RFS analysis, all anti_estrogen values the same: ",unique(unique(dataMatrix$anti_estrogen)),"\n")),
            append=TRUE,file=summaryFile)
        
        linearModels[["RFS_logitWithRx"]] <- tryCatch(glm(RFS~community+chemotherapyClass+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
                                                      
        )
        
                   #no community
        reducedModel <- tryCatch(glm(RFS~chemotherapyClass+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
         linearModels[["RFS_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["RFS_logitWithRx"]],test="Chisq")
        
 
        
        #now try with  grade.    
        cat(gsub(" ","_",paste0("\nLinear analysis with RFS and histological grade :\n")),append=TRUE,file=summaryFile)
        dataMatrixHist <- dataMatrix[which(!is.na((dataMatrix$hist_grade))),]
 
        
        cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
        
        linearModelsSumm[["RFS_withGrade_numPatients"]] <- nrow(dataMatrixHist)
        cat("\nNumber_of_patients_for_RFS_with_grade_model: ",nrow(dataMatrixHist),"\n",append=TRUE,file=summaryFile)
        linearModels[["RFS_logitWithRxGrade"]] <- tryCatch( glm(RFS~community+chemotherapyClass+anti_HER2+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
                                                            error = function(e) {
                                                              return(NA)
                                                            }
        )
        
        
      }else if(length(unique(dataMatrix$chemotherapyClass))==1 && length(unique(dataMatrix$anti_estrogen))>1 
               && length(unique(dataMatrix$anti_HER2))>1 ){
        
        
        cat(gsub(" ","_",paste0("\nFor RFS analysis, all chemo values the same: ",unique(unique(dataMatrix$chemotherapyClass)),"\n")),
            append=TRUE,file=summaryFile)
        
        linearModels[["RFS_logitWithRx"]] <- tryCatch(glm(RFS~community+anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
                                                      
        )
        
                           #no community
        reducedModel <- tryCatch(glm(RFS~anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
         linearModels[["RFS_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["RFS_logitWithRx"]],test="Chisq")
        #now try with  grade.    
        cat(gsub(" ","_",paste0("\nLinear analysis with RFS and histological grade :\n")),append=TRUE,file=summaryFile)
        dataMatrixHist <- dataMatrix[which(!is.na((dataMatrix$hist_grade))),]
 
        
        cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
        
        linearModelsSumm[["RFS_withGrade_numPatients"]] <- nrow(dataMatrixHist)
        cat("\nNumber_of_patients_for_RFS_with_grade_model: ",nrow(dataMatrixHist),"\n",append=TRUE,file=summaryFile)
        linearModels[["RFS_logitWithRxGrade"]] <- tryCatch( glm(RFS~community+anti_estrogen+anti_HER2+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
                                                            error = function(e) {
                                                              return(NA)
                                                            }
        )
        
        
        
      }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))==1 
               && length(unique(dataMatrix$anti_HER2))==1 ){
        
        
        cat(gsub(" ","_",paste0("\nFor RFS analysis, the anti HER2 and estrogen variables have only one unique value each.\n")),append=TRUE,
            file=summaryFile)
        
        
        linearModels[["RFS_logitWithRx"]] <- tryCatch(glm(RFS~community+chemotherapyClass,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
                                                      
        )
         reducedModel <- tryCatch(glm(RFS~chemotherapyClass,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
         linearModels[["RFS_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["RFS_logitWithRx"]],test="Chisq")
    
        
        #now try with  grade.    
        cat(gsub(" ","_",paste0("\nLinear analysis with RFS and histological grade :\n")),append=TRUE,file=summaryFile)
        dataMatrixHist <- dataMatrix[which(!is.na((dataMatrix$hist_grade))),]
 
        
        cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
        
        linearModelsSumm[["RFS_withGrade_numPatients"]] <- nrow(dataMatrixHist)
        cat("\nNumber_of_patients_for_RFS_with_grade_model: ",nrow(dataMatrixHist),"\n",append=TRUE,file=summaryFile)
        linearModels[["RFS_logitWithRxGrade"]] <- tryCatch( glm(RFS~community+chemotherapyClass+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
                                                            error = function(e) {
                                                              return(NA)
                                                            }
        )
        
        
        
        
      }else{
        
        linearModels[["RFS_logitWithRx"]] <- NA
        
      }
      
      cat(gsub(" ","_",paste0("\nLinear analysis with RFS and treatment:\n")),append=TRUE,file=summaryFile)
      if(any(!is.na(linearModels[["RFS_logitWithRx"]]))){
        
        linearModelsSumm[["RFS_logitWithRx"]] <- summary( linearModels[["RFS_logitWithRx"]])
        textOut <- capture.output( linearModelsSumm[["RFS_logitWithRx"]]$coefficients)
        cat(textOut,sep="\n",append=TRUE,file=summaryFile)
        textOut <- capture.output( linearModelsSumm [["RFS_logitWithRx"]])
        cat(textOut,sep="\n",append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nLinear analysis with RFS and treatment ANOVA for meta-clusters:\n")),append=TRUE,file=summaryFile)
        textOut <- capture.output(linearModels[["RFS_logitWithRx_ANOVA"]])
        cat(textOut,sep="\n",append=TRUE,file=summaryFile)
        
        linearModelsSumm[["RFS_WithRxMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",") 
        predictor <- predict(linearModels[["RFS_logitWithRx"]],type="response")
        png(filename=paste0(saveDirMerged,"/",experimentName,"_RFS_withRx_ROC_",Sys.Date(),".png"),width=1000,height=1000,res=200)
        roc_data <- roc(response=dataMatrix$RFS,predictor=predictor,plot=TRUE,na.rm=TRUE,auc=TRUE)
        dev.off()
        textOut <- capture.output(roc_data)   
        ROC_list[["RFS_logitWithRx"]] <- roc_data
        cat("\n roc_info:",textOut,sep="\n",append=TRUE,file=summaryFile)
        cat("\nAUC : ",roc_data$auc,"\n",append=TRUE,file=summaryFile)
        
        linearModelsSumm[["RFS_logitWithRx"]] <- summary( linearModels[["RFS_logitWithRx"]])
        textOut <- capture.output( linearModelsSumm[["RFS_logitWithRx"]]$coefficients)
        cat(textOut,sep="\n",append=TRUE,file=summaryFile)
        
        linearModelsSumm[["RFS_WithRxMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
        
      }else{
        
        cat(gsub(" ","_",paste0("\nRFS analysis with therapies returned NA.\n")),append=TRUE,file=summaryFile)
      }
      
      #now record grade.
      cat(gsub(" ","_",paste0("\nLinear analysis with RFS and histological grade :\n")),append=TRUE,file=summaryFile)
      cat(gsub(" ","_",paste0("\nUnique hist grade values :",unique(dataMatrixHist$hist_grade),"\n")),append=TRUE,file=summaryFile)
      
      
      if(any(!is.na(linearModels[["RFS_logitWithRxGrade"]]))){
        
        predictor <- predict(linearModels[["RFS_logitWithRxGrade"]],type="response")
        png(filename=paste0(saveDirMerged,"/",experimentName,"_RFS_withGrade_ROC_",Sys.Date(),".png"),width=1000,height=1000,res=200)
        roc_data <- roc(response=dataMatrixHist$RFS,predictor=predictor,plot=TRUE,na.rm=TRUE,auc=TRUE)
        dev.off()
          ROC_list[["RFS_logitWithRxGrade"]] <- roc_data
        textOut <- capture.output(roc_data)
        cat("roc_info:\n",textOut,sep="\n",append=TRUE,file=summaryFile)
        cat("\nAUC: ",roc_data$auc,"\n",append=TRUE,file=summaryFile)
        
        linearModelsSumm[["RFS_logitWithRxGrade"]] <- summary( linearModels[["RFS_logitWithRxGrade"]])
        textOut <- capture.output( linearModelsSumm[["RFS_logitWithRxGrade"]]$coefficients)
        cat(textOut,sep="\n",append=TRUE,file=summaryFile) 
        linearModelsSumm[["RFS_WithGradeMetaClusters_used"]] <- paste(unique(dataMatrixHist$community),collapse=",")
        
      }else{
        
        cat(gsub(" ","_",paste0("\nHist grade with RFS analysis returned NA.\n")),append=TRUE,file=summaryFile)  
      }
      
      
      
      message("Running DFS models")
      ###DFS
      cat(gsub(" ","_",paste0("\nLinear analysis with only DFS:\n")),append=TRUE,file=summaryFile)
      
      dataMatrix <- sampleClustCommPhenoData[which(!is.na((sampleClustCommPhenoData$DFS))), ]
      dataMatrix$community <- as.factor(dataMatrix$community)
      dataMatrix$DFS <- as.factor(dataMatrix$DFS)
      dataMatrix$hist_grade <- as.factor(dataMatrix$hist_grade)
      dataMatrix$chemotherapyClass <- as.factor(dataMatrix$chemotherapyClass)
      dataMatrix$anti_estrogen <- as.factor(dataMatrix$anti_estrogen)
      dataMatrix$anti_HER2 <- as.factor(dataMatrix$anti_HER2)
      cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrix$studyNum),"\n")),append=TRUE,file=summaryFile)
      cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrix$community),"\n")),append=TRUE,file=summaryFile)
      
      #strongly predicts DFS.
      linearModels[["DFS_logitAlone"]] <- glm(DFS~community,data=dataMatrix,family=binomial(link="logit"))
      png(filename=paste0(saveDirMerged,"/",experimentName,"_DFS_ROC_",Sys.Date(),".png"),width=1000,height=1000,res=200)
      predictor <- predict(linearModels[["DFS_logitAlone"]],type="response")
      roc_data <- roc(response=dataMatrix$DFS,predictor=predictor,plot=TRUE,na.rm=TRUE,auc=TRUE)
      dev.off()
     ROC_list[["DFS_logitAlone"]] <- roc_data
      textOut <- capture.output(roc_data)
      cat("roc_info\n",textOut,sep="\n",append=TRUE,file=summaryFile)
      cat("\nAUC: ",roc_data$auc,"\n",append=TRUE,file=summaryFile)
      linearModelsSumm[["DFS_logitAlone"]] <- summary( linearModels[["DFS_logitAlone"]] ) 
      textOut <- capture.output( linearModelsSumm[["DFS_logitAlone"]]$coefficients)
      cat(textOut,sep="\n",append=TRUE,file=summaryFile)
      textOut <- capture.output( linearModelsSumm[["DFS_logitAlone"]])
      cat(textOut,sep="\n",append=TRUE,file=summaryFile)
      linearModelsSumm[["DFS_numPatients"]] <- nrow(dataMatrix)
      cat("\nNumber of patients: ",nrow(dataMatrix),"\n",append=TRUE,file=summaryFile)
      linearModelsSumm[["DFS_aloneMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
      
      ##with therapies
      #interesting: community becomes insignificant if add in therapies now: this could be
      #because community is most likely correlated with treatments; it's certainly correlated with pam50 status.
      #low Rsquared.
      cat(gsub(" ","_",paste0("\nLinear analysis with DFS and treatment:\n")),append=TRUE,file=summaryFile)
      dataMatrix <- dataMatrix[which(!is.na(dataMatrix$chemotherapyClass)), ]
      dataMatrix <- dataMatrix[which(!is.na(dataMatrix$anti_estrogen)), ]
      dataMatrix <- dataMatrix[which(!is.na(dataMatrix$anti_HER2)), ]
      dataMatrix <- data.frame(dataMatrix)
      dataMatrix$community <- as.factor(dataMatrix$community)
      dataMatrix$DFS <- as.factor(dataMatrix$DFS)
      dataMatrix$hist_grade <- as.factor(dataMatrix$hist_grade)
      dataMatrix$chemotherapyClass <- as.factor(dataMatrix$chemotherapyClass)
      dataMatrix$anti_estrogen <- as.factor(dataMatrix$anti_estrogen)
      dataMatrix$anti_HER2 <- as.factor(dataMatrix$anti_HER2)
      cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrix$studyNum),"\n")),append=TRUE,file=summaryFile)
      cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrix$community),"\n")),append=TRUE,file=summaryFile)
      
      linearModelsSumm[["DFS_withRx_numPatients"]] <- nrow(dataMatrix)
      cat("\nNumber_of_patients_for_DFS_with_Rx_model: ",nrow(dataMatrix),"\n",append=TRUE,file=summaryFile)
      
      if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))>1 && 
         length(unique(dataMatrix$anti_HER2))>1 ){
        
        linearModels[["DFS_logitWithRx"]] <- tryCatch(glm(DFS~community+chemotherapyClass+anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
        
        #remove communities for ANOVA
           reducedModel <- tryCatch(glm(DFS~chemotherapyClass+anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
         linearModels[["DFS_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["DFS_logitWithRx"]],test="Chisq")
    
    
        #now try with grade.
        cat(gsub(" ","_",paste0("\nLinear analysis with DFS and histological grade :\n")),append=TRUE,file=summaryFile)
        dataMatrixHist <- dataMatrix[which(!is.na((dataMatrix$hist_grade))),]
 

        cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
        
        linearModelsSumm[["DFS_withGrade_numPatients"]] <- nrow(dataMatrixHist)
        cat("\nNumber_of_patients_for_DFS_withGrade_model: ",nrow(dataMatrixHist),"\n",append=TRUE,file=summaryFile)
        linearModels[["DFS_logitWithRxGrade"]] <- tryCatch( glm(DFS~community+chemotherapyClass+anti_estrogen+anti_HER2+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
                                                            error = function(e) {
                                                              return(NA)
                                                            }
        )
        
        
        
      }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))>1 &&
               length(unique(dataMatrix$anti_HER2))==1 ){
        
        
        cat(gsub(" ","_",paste0("\nFor DFS analysis, all anti_HER2 values the same: ",unique(unique(dataMatrix$anti_HER2)),"\n")),
            append=TRUE,file=summaryFile)
        linearModels[["DFS_logitWithRx"]] <- tryCatch(glm(DFS~community+chemotherapyClass+anti_estrogen,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
                #remove communities for ANOVA
           reducedModel <- tryCatch(glm(DFS~chemotherapyClass+anti_estrogen,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
         linearModels[["DFS_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["DFS_logitWithRx"]],test="Chisq")
    

        #now try with grade.
        cat(gsub(" ","_",paste0("\nLinear analysis with DFS and histological grade :\n")),append=TRUE,file=summaryFile)
        dataMatrixHist <- dataMatrix[which(!is.na((dataMatrix$hist_grade))),]
 
        
        cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
        
        linearModelsSumm[["DFS_withGrade_numPatients"]] <- nrow(dataMatrixHist)
        cat("\nNumber_of_patients_for_DFS_withGrade_model: ",nrow(dataMatrixHist),"\n",append=TRUE,file=summaryFile)
        linearModels[["DFS_logitWithRxGrade"]] <- tryCatch( glm(DFS~community+chemotherapyClass+anti_estrogen+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
                                                            error = function(e) {
                                                              return(NA)
                                                            }
        )
        
        
        
        
      }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))==1 
               && length(unique(dataMatrix$anti_HER2))>1 ){
        
        cat(gsub(" ","_",paste0("\nFor DFS analysis, all anti_estrogen values the same: ",unique(unique(dataMatrix$anti_estrogen)),"\n")),
            append=TRUE,file=summaryFile)
        
        linearModels[["DFS_logitWithRx"]] <- tryCatch(glm(DFS~community+chemotherapyClass+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
                    #remove communities for ANOVA
           reducedModel <- tryCatch(glm(DFS~chemotherapyClass+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
         linearModels[["DFS_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["DFS_logitWithRx"]],test="Chisq")
    

        #now try with grade.
        cat(gsub(" ","_",paste0("\nLinear analysis with DFS and histological grade :\n")),append=TRUE,file=summaryFile)
        dataMatrixHist <- dataMatrix[which(!is.na((dataMatrix$hist_grade))),]
 
        
        cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
        
        linearModelsSumm[["DFS_withGrade_numPatients"]] <- nrow(dataMatrixHist)
        cat("\nNumber_of_patients_for_DFS_withGrade_model: ",nrow(dataMatrixHist),"\n",append=TRUE,file=summaryFile)
        linearModels[["DFS_logitWithRxGrade"]] <- tryCatch( glm(DFS~community+chemotherapyClass+anti_HER2+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
                                                            error = function(e) {
                                                              return(NA)
                                                            }
        )
        
        
        
      }else if(length(unique(dataMatrix$chemotherapyClass))==1 && length(unique(dataMatrix$anti_estrogen))>1 
               && length(unique(dataMatrix$anti_HER2))>1 ){
        
        
        cat(gsub(" ","_",paste0("\nFor DFS analysis, all chemo values the same: ",unique(unique(dataMatrix$chemotherapyClass)),"\n")),
            append=TRUE,file=summaryFile)
        
        linearModels[["DFS_logitWithRx"]] <- tryCatch(glm(DFS~community+anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
                          #remove communities for ANOVA
           reducedModel <- tryCatch(glm(DFS~anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
         linearModels[["DFS_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["DFS_logitWithRx"]],test="Chisq")
    

        #now try with grade.
        cat(gsub(" ","_",paste0("\nLinear analysis with DFS and histological grade :\n")),append=TRUE,file=summaryFile)
        dataMatrixHist <- dataMatrix[which(!is.na((dataMatrix$hist_grade))),]
 
        
        cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
        
        linearModelsSumm[["DFS_withGrade_numPatients"]] <- nrow(dataMatrixHist)
        cat("\nNumber_of_patients_for_DFS_withGrade_model: ",nrow(dataMatrixHist),"\n",append=TRUE,file=summaryFile)
        linearModels[["DFS_logitWithRxGrade"]] <- tryCatch( glm(DFS~community+anti_estrogen+anti_HER2+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
                                                            error = function(e) {
                                                              return(NA)
                                                            }
        )
        
        
        
      }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))==1 
               && length(unique(dataMatrix$anti_HER2))==1 ){
        
        
        cat(gsub(" ","_",paste0("\nFor DFS analysis, the anti HER2 and estrogen variables have only one unique value each.\n")),append=TRUE,
            file=summaryFile)
        
        
        linearModels[["DFS_logitWithRx"]] <- tryCatch(glm(DFS~community+chemotherapyClass,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
                               #remove communities for ANOVA
           reducedModel <- tryCatch(glm(DFS~chemotherapyClass,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
         linearModels[["DFS_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["DFS_logitWithRx"]],test="Chisq")
    

        #now try with grade.
        cat(gsub(" ","_",paste0("\nLinear analysis with DFS and histological grade :\n")),append=TRUE,file=summaryFile)
        dataMatrixHist <- dataMatrix[which(!is.na((dataMatrix$hist_grade))),]
 
        
        cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
        
        linearModelsSumm[["DFS_withGrade_numPatients"]] <- nrow(dataMatrixHist)
        cat("\nNumber_of_patients_for_DFS_withGrade_model: ",nrow(dataMatrixHist),"\n",append=TRUE,file=summaryFile)
        linearModels[["DFS_logitWithRxGrade"]] <- tryCatch( glm(DFS~community+chemotherapyClass+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
                                                            error = function(e) {
                                                              return(NA)
                                                            }
        )
        
        
        
      }else{
        
        # don't run models
        
      }
      
      cat(gsub(" ","_",paste0("\nLinear analysis with DFS and treatment:\n")),append=TRUE,file=summaryFile)
      if(any(!is.na(linearModels[["DFS_logitWithRx"]]))){
        
        predictor <- predict(linearModels[["DFS_logitWithRx"]],type="response")
        png(filename=paste0(saveDirMerged,"/",experimentName,"_DFS_withRx_ROC_",Sys.Date(),".png"),width=1000,height=1000,res=200)
        roc_data <- roc(response=dataMatrix$DFS,predictor=predictor,plot=TRUE,na.rm=TRUE,auc=TRUE)
        dev.off()
          ROC_list[["DFS_logitWithRx"]] <- roc_data
        textOut <- capture.output(roc_data)
        cat("roc_info\n",textOut,sep="\n",append=TRUE,file=summaryFile)
        cat("\nAUC: ",roc_data$auc,"\n",append=TRUE,file=summaryFile)
        
        linearModelsSumm[["DFS_logitWithRx"]] <- summary( linearModels[["DFS_logitWithRx"]])
        textOut <- capture.output( linearModelsSumm[["DFS_logitWithRx"]]$coefficients)
        cat(textOut,sep="\n",append=TRUE,file=summaryFile)
                textOut <- capture.output( linearModelsSumm[["DFS_logitWithRx"]])
        cat(textOut,sep="\n",append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nLinear analysis with DFS and treatment ANOVA for meta-clusters:\n")),append=TRUE,file=summaryFile)
        textOut <- capture.output(linearModels[["DFS_logitWithRx_ANOVA"]])
        cat(textOut,sep="\n",append=TRUE,file=summaryFile)
        linearModelsSumm[["DFS_WithRxMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
        
      }else{
        
        cat(gsub(" ","_",paste0("\nRx with DFS analysis returned NA.\n")),append=TRUE,file=summaryFile)  
      }
      
      
      #now try with grade.
      cat(gsub(" ","_",paste0("\nLinear analysis with DFS and histological grade :\n")),append=TRUE,file=summaryFile)
      cat(gsub(" ","_",paste0("\nUnique hist grade values :",unique(dataMatrixHist$hist_grade),"\n")),append=TRUE,file=summaryFile)
      
      
      if(!is.na(linearModels[["DFS_logitWithRxGrade"]])){
        
        linearModelsSumm[["DFS_logitWithRxGrade"]] <- summary( linearModels[["DFS_logitWithRxGrade"]])
        textOut <- capture.output( linearModelsSumm[["DFS_logitWithRxGrade"]]$coefficients)
        cat(textOut,sep="\n",append=TRUE,file=summaryFile) 
        linearModelsSumm[["DFS_WithGradeMetaClusters_used"]] <- paste(unique(dataMatrixHist$community),collapse=",")
        
        predictor <- predict(linearModels[["DFS_logitWithRxGrade"]],type="response")
        png(filename=paste0(saveDirMerged,"/",experimentName,"_DFS_withGrade_ROC_",Sys.Date(),".png"),width=1000,height=1000,res=200)
        roc_data <- roc(response=dataMatrixHist$DFS,predictor=predictor,plot=TRUE,na.rm=TRUE,auc=TRUE)
        dev.off()
          ROC_list[["DFS_logitWithRxGrade"]] <- roc_data
        textOut <- capture.output(roc_data)
        cat("roc_info\n",textOut,sep="\n",append=TRUE,file=summaryFile)
        cat("\nAUC: ",roc_data$auc,"\n",append=TRUE,file=summaryFile)
        
      }else{
        
        cat(gsub(" ","_",paste0("\nHist grade with DFS analysis returned NA\n.")),append=TRUE,file=summaryFile)  
      }
      
      
      message("Running pCR models")
      #####pCR
      cat(gsub(" ","_","\nLinear analysis with only pCR:\n"),append=TRUE,file=summaryFile)
      
      dataMatrix <- sampleClustCommPhenoData[which(!is.na((sampleClustCommPhenoData$pCR))), ]
      dataMatrix$community <- as.factor(dataMatrix$community)
      dataMatrix$pCR <- as.factor(dataMatrix$pCR)
      dataMatrix$hist_grade <- as.factor(dataMatrix$hist_grade)
      dataMatrix$chemotherapyClass <- as.factor(dataMatrix$chemotherapyClass)
      dataMatrix$anti_estrogen <- as.factor(dataMatrix$anti_estrogen)
      dataMatrix$anti_HER2 <- as.factor(dataMatrix$anti_HER2)
      cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrix$studyNum),"\n")),append=TRUE,file=summaryFile)
      cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrix$community),"\n")),append=TRUE,file=summaryFile)
      
      #strongly predicts pCR.
      
  #hmm...turns out that sometimes only one meta-cluster has pCR values.
      linearModels[["pCR_logitAlone"]] <- glm(pCR~community,data=dataMatrix,family=binomial(link="logit"))
      if(!is.na(linearModelslinearModels[["pCR_logitAlone"]])){
              
      predictor <- predict(linearModels[["pCR_logitAlone"]],type="response")
      png(filename=paste0(saveDirMerged,"/",experimentName,"_pCR_ROC_",Sys.Date(),".png"),width=1000,height=1000,res=200)
      roc_data <- roc(response=dataMatrix$pCR,predictor=predictor,plot=TRUE,na.rm=TRUE,auc=TRUE)
      dev.off()
      ROC_list[["pCR_logitAlone"]] <- roc_data
      textOut <- capture.output(roc_data)
      cat("roc_info\n",textOut,sep="\n",append=TRUE,file=summaryFile)
      cat("\nAUC: ",roc_data$auc,"\n",append=TRUE,file=summaryFile)
      
      
      linearModelsSumm[["pCR_logitAlone"]] <- summary(  linearModels[["pCR_logitAlone"]] ) 
      textOut <- capture.output( linearModelsSumm[["pCR_logitAlone"]]$coefficients)
      cat(textOut,sep="\n",append=TRUE,file=summaryFile)
        textOut <- capture.output( linearModelsSumm[["pCR_logitAlone"]])
      cat(textOut,sep="\n",append=TRUE,file=summaryFile)
      linearModelsSumm[["pCR_numPatients"]] <- nrow(dataMatrix)
      cat("\nNumber_of_patients_for_pCR_model: ",nrow(dataMatrix),"\n",append=TRUE,file=summaryFile)
      linearModelsSumm[["pCR_aloneMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
      
      }else{
        
         cat(gsub(" ","_",paste0("\n pCR alone analysis returned NA\n.")),append=TRUE,file=summaryFile)  
      }
      ##with therapies
      #interesting: community becomes insignificant if add in therapies now: this could be
      #because community is most likely correlated with treatments; it's certainly correlated with pam50 status.
      #low Rsquared.
      cat(gsub(" ","_",paste0("\nLinear analysis with pCR and treatment:\n")),append=TRUE,file=summaryFile)
      dataMatrix <- dataMatrix[which(!is.na(dataMatrix$chemotherapyClass)), ]
      dataMatrix <- dataMatrix[which(!is.na(dataMatrix$anti_estrogen)), ]
      dataMatrix <- dataMatrix[which(!is.na(dataMatrix$anti_HER2)), ]
      dataMatrix <- data.frame(dataMatrix)
      dataMatrix$community <- as.factor(dataMatrix$community)
      dataMatrix$pCR <- as.factor(dataMatrix$pCR)
      dataMatrix$hist_grade <- as.factor(dataMatrix$hist_grade)
      dataMatrix$chemotherapyClass <- as.factor(dataMatrix$chemotherapyClass)
      dataMatrix$anti_estrogen <- as.factor(dataMatrix$anti_estrogen)
      dataMatrix$anti_HER2 <- as.factor(dataMatrix$anti_HER2)
      cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrix$studyNum),"\n")),append=TRUE,file=summaryFile)
      cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrix$community),"\n")),append=TRUE,file=summaryFile)
      
      linearModelsSumm[["pCR_withRx_numPatients"]] <- nrow(dataMatrix)
      cat(gsub(" ","_",paste0("\nNumber_of_patients_for_pCR_withRx_model: ",nrow(dataMatrix),"\n")),append=TRUE,file=summaryFile)
      
      if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))>1 &&
         length(unique(dataMatrix$anti_HER2))>1 ){
        
        linearModels[["pCR_logitWithRx"]] <- tryCatch(glm(pCR~community+chemotherapyClass+anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
        
                                    #remove communities for ANOVA
           reducedModel <- tryCatch(glm(pCR~chemotherapyClass+anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
         linearModels[["pCR_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["pCR_logitWithRx"]],test="Chisq")
    

        #now try with grade.
        cat(gsub(" ","_",paste0("\nLinear analysis with pCR and histological grade :\n")),append=TRUE,file=summaryFile)
        dataMatrixHist <- dataMatrix[which(!is.na(dataMatrix$hist_grade)),]
        cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
        
        linearModelsSumm[["pCR_withGrade_numPatients"]] <- nrow(dataMatrixHist)
        cat(gsub(" ","_",paste0("\nNumber_of_patients_for_pCR_withGrade_model: ",nrow(dataMatrixHist),"\n")),append=TRUE,file=summaryFile)
        linearModels[["pCR_logitWithRxGrade"]] <- tryCatch( glm(pCR~community+chemotherapyClass+anti_estrogen+anti_HER2+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
                                                            error = function(e) {
                                                              return(NA)
                                                            }
        )
        
        
      }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))>1 
               && length(unique(dataMatrix$anti_HER2))==1 ){
        
        
        cat(gsub(" ","_",paste0("\nFor pCR analysis, all anti_HER2 values the same: ",unique(unique(dataMatrix$anti_HER2)),"\n")),
            append=TRUE,file=summaryFile)
        linearModels[["pCR_logitWithRx"]] <- tryCatch(glm(pCR~community+chemotherapyClass+anti_estrogen,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
                                            #remove communities for ANOVA
           reducedModel <- tryCatch(glm(pCR~chemotherapyClass+anti_estrogen,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
         linearModels[["pCR_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["pCR_logitWithRx"]],test="Chisq")
        #now try with stage, grade.
        cat(gsub(" ","_",paste0("\nLinear analysis with pCR and histological grade :\n")),append=TRUE,file=summaryFile)
        dataMatrixHist <- dataMatrix[which(!is.na(dataMatrix$hist_grade)),]
 
        cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
        
        linearModelsSumm[["pCR_withGrade_numPatients"]] <- nrow(dataMatrixHist)
        cat(gsub(" ","_",paste0("\nNumber_of_patients_for_pCR_withGrade_model: ",nrow(dataMatrixHist),"\n")),append=TRUE,file=summaryFile)
        linearModels[["pCR_logitWithRxGrade"]] <- tryCatch( glm(pCR~community+chemotherapyClass+anti_estrogen+anti_HER2+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
                                                            error = function(e) {
                                                              return(NA)
                                                            }
        )
        
        
      }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))==1 
               && length(unique(dataMatrix$anti_HER2))>1 ){
        
        cat(gsub(" ","_",paste0("\nFor pCR analysis, all anti_estrogen values the same: ",unique(unique(dataMatrix$anti_estrogen)),"\n")),
            append=TRUE,file=summaryFile)
        
        linearModels[["pCR_logitWithRx"]] <- tryCatch(glm(pCR~community+chemotherapyClass+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
                                                  #remove communities for ANOVA
           reducedModel <- tryCatch(glm(pCR~chemotherapyClass+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
         linearModels[["pCR_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["pCR_logitWithRx"]],test="Chisq")

        #now try with stage, grade.
        cat(gsub(" ","_",paste0("\nLinear analysis with pCR and histological grade :\n")),append=TRUE,file=summaryFile)
        dataMatrixHist <- dataMatrix[which(!is.na(dataMatrix$hist_grade)),]
 
        cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
        
        linearModelsSumm[["pCR_withGrade_numPatients"]] <- nrow(dataMatrixHist)
        cat(gsub(" ","_",paste0("\nNumber_of_patients_for_pCR_withGrade_model: ",nrow(dataMatrixHist),"\n")),append=TRUE,file=summaryFile)
        linearModels[["pCR_logitWithRxGrade"]] <- tryCatch( glm(pCR~community+chemotherapyClass+anti_HER2+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
                                                            error = function(e) {
                                                              return(NA)
                                                            }
        )
        
        
        
        
        
      }else if(length(unique(dataMatrix$chemotherapyClass))==1 && length(unique(dataMatrix$anti_estrogen))>1 
               && length(unique(dataMatrix$anti_HER2))>1 ){
        
        
        cat(gsub(" ","_",paste0("\nFor pCR analysis, all chemo values the same: ",unique(dataMatrix$chemotherapyClass),"\n")),
            append=TRUE,file=summaryFile)
        
        linearModels[["pCR_logitWithRx"]] <- tryCatch(glm(pCR~community+anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
   #remove community for ANOVA
        reducedModel <- tryCatch(glm(pCR~anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
         linearModels[["pCR_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["pCR_logitWithRx"]],test="Chisq")

 
        #now try with stage, grade.
        cat(gsub(" ","_",paste0("\nLinear analysis with pCR and histological grade :\n")),append=TRUE,file=summaryFile)
        dataMatrixHist <- dataMatrix[which(!is.na(dataMatrix$hist_grade)),]
 
        cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
        
        linearModelsSumm[["pCR_withGrade_numPatients"]] <- nrow(dataMatrix)
        cat(gsub(" ","_",paste0("\nNumber_of_patients_for_pCR_withGrade_model: ",nrow(dataMatrixHist),"\n")),append=TRUE,file=summaryFile)
        linearModels[["pCR_logitWithRxGrade"]] <- tryCatch( glm(pCR~community+anti_estrogen+anti_HER2+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
                                                            error = function(e) {
                                                              return(NA)
                                                            }
        )
        
        
        
        
        
      }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))==1 
               && length(unique(dataMatrix$anti_HER2))==1 ){
        
        
        cat(gsub(" ","_",paste0("\nFor pCR analysis, the anti HER2 and estrogen variables have only one unique value each.\n")),append=TRUE,
            file=summaryFile)
        
        
        linearModels[["pCR_logitWithRx"]] <- tryCatch(glm(pCR~community+chemotherapyClass,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )

   #remove community for ANOVA
        reducedModel <- tryCatch(glm(pCR~chemotherapyClass,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
         linearModels[["pCR_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["pCR_logitWithRx"]],test="Chisq")


        #now try with stage, grade.
        cat(gsub(" ","_",paste0("\nLinear analysis with pCR and histological grade :\n")),append=TRUE,file=summaryFile)
        dataMatrixHist <- dataMatrix[which(!is.na(dataMatrix$hist_grade)),]
 
        cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
        
        linearModelsSumm[["pCR_withGrade_numPatients"]] <- nrow(dataMatrixHist)
        cat(gsub(" ","_",paste0("\nNumber_of_patients_for_pCR_withGrade_model: ",nrow(dataMatrixHist),"\n")),append=TRUE,file=summaryFile)
        linearModels[["pCR_logitWithRxGrade"]] <- tryCatch( glm(pCR~community+chemotherapyClass+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
                                                            error = function(e) {
                                                              return(NA)
                                                            }
        )
        
        
        
        
      }else{
        
        # don't run models
        
      }
      
      cat(gsub(" ","_",paste0("\nLinear analysis with pCR and treatment:\n")),append=TRUE,file=summaryFile)
      if(any(!is.na(linearModels[["pCR_logitWithRx"]]))){
        
        linearModelsSumm[["pCR_logitWithRx"]] <- summary( linearModels[["pCR_logitWithRx"]])
        textOut <- capture.output( linearModelsSumm[["pCR_logitWithRx"]]$coefficients)
        cat(textOut,sep="\n",append=TRUE,file=summaryFile) 
        linearModelsSumm[["pCR_WithRxMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
                textOut <- capture.output( linearModelsSumm [["pCR_logitWithRx"]])
        cat(textOut,sep="\n",append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nLinear analysis with pCR and treatment ANOVA for meta-clusters:\n")),append=TRUE,file=summaryFile)
        textOut <- capture.output(linearModels[["pCR_logitWithRx_ANOVA"]])
        cat(textOut,sep="\n",append=TRUE,file=summaryFile)
        
        predictor <- predict(linearModels[["pCR_logitWithRx"]],type="response")
        png(filename=paste0(saveDirMerged,"/",experimentName,"_pCR_withRx_ROC_",Sys.Date(),".png"),width=1000,height=1000,res=200)
        roc_data <- roc(response=dataMatrix$pCR,predictor=predictor,plot=TRUE,na.rm=TRUE,auc=TRUE)
        dev.off()
          ROC_list[["pCR_logitWithRx"]] <- roc_data
        textOut <- capture.output(roc_data)
        cat("roc_info\n",textOut,sep="\n",append=TRUE,file=summaryFile)
        cat("\nAUC: ",roc_data$auc,"\n",append=TRUE,file=summaryFile)
        
        
        
      }else{
        
        cat(gsub(" ","_",paste0("\nRx with pCR analysis returned NA\n.")),append=TRUE,file=summaryFile)  
      }
      
      
      #now try with stage, grade.
      cat(gsub(" ","_",paste0("\nLinear analysis with pCR and histological grade :\n")),append=TRUE,file=summaryFile)
      cat(gsub(" ","_",paste0("\nUnique hist grade values :",unique(dataMatrixHist$hist_grade),"\n")),append=TRUE,file=summaryFile)

      if(any(!is.na(linearModels[["pCR_logitWithRxGrade"]]))){
        
        linearModelsSumm[["pCR_logitWithRxGrade"]] <- summary( linearModels[["pCR_logitWithRxGrade"]])
        textOut <- capture.output( linearModelsSumm[["pCR_logitWithRxGrade"]]$coefficients)
        cat(textOut,sep="\n",append=TRUE,file=summaryFile) 
        linearModelsSumm[["pCR_logitWithRxGrade"]] <- paste(unique(dataMatrixHist$community),collapse=",")
        
        predictor <- predict(linearModels[["pCR_logitWithRxGrade"]],type="response")
        png(filename=paste0(saveDirMerged,"/",experimentName,"_pCR_withGrade_ROC_",Sys.Date(),".png"),width=1000,height=1000,res=200)
        roc_data <- roc(response=dataMatrixHist$pCR,predictor=predictor,plot=TRUE,na.rm=TRUE,auc=TRUE)
        dev.off()
             ROC_list[["pCR_logitWithRxGrade"]] <- roc_data
        textOut <- capture.output(roc_data)
        cat("roc_info\n",textOut,sep="\n",append=TRUE,file=summaryFile)
        cat("\nAUC: ",roc_data$auc,"\n",append=TRUE,file=summaryFile)
        
      }else{
        
        cat(gsub(" ","_",paste0("\nHist grade with pCR analysis returned NA.")),append=TRUE,file=summaryFile)  
      }

  output <- list(sampleClustCommPhenoData=sampleClustCommPhenoData,ROC_data=ROC_list,
                 linearModelSummaries=linearModelsSumm,linearSurvModels=linearModels)
  return(output)
}
