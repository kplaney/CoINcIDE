###look at Pam50 breakdown
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_visualization.R")
saveDir <- "/home/kplaney/breast_analysis/merged/"
##assign subtypes
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

#how changed from original??
#indices matched up?
all(colnames(mergedCombat)==colnames(mergedNoNorm))
all(colnames(mergedBMC)==colnames(mergedNoNorm))
#starting with 2237 samples
length(Pam50_subtypes_noNorm$subtypes[,"subtypeLabels"])
#combat changed 673 samples
length(which(Pam50_subtypes_combatNorm$subtypes[,"subtypeLabels"]!=Pam50_subtypes_noNorm$subtypes[,"subtypeLabels"]))
#changed 821
length(which(Pam50_subtypes_BMCNorm$subtypes[,"subtypeLabels"]!=Pam50_subtypes_noNorm$subtypes[,"subtypeLabels"]))



#load up original full subtypes too.
load("/home/kplaney/breast_analysis/pam50FullAndShort_subtypeDF.RData.gzip")
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
save(pam50ShortSubtypeDF,file=paste0(saveDir,"/mergedPam50SubtypeMatrix.RData.gzip"),compress="gzip")
write.table(pam50ShortSubtypeDF,file=paste0(saveDir,"/mergedPam50SubtypeMatrix.txt"),quote=FALSE,
                                            col.names=TRUE,row.names=FALSE)
   
textOut <- capture.output(table(subtypeDF_master$subtype))
cat(textOut,sep="\n",file=paste0(saveDir,"/noNormPam50Full_subtypeTable.txt"))

textOut <- capture.output(table(subtypeDF_master$subtype_short))
cat(textOut,sep="\n",file=paste0(saveDir,"/noNormPam50Short_subtypeTable.txt"))

textOut <- capture.output(table(Pam50_subtypes_BMCNorm$subtypes[,"subtypeLabels"]))
cat(textOut,sep="\n",file=paste0(saveDir,"/mergedBMCNormPam50Short_subtypeTable.txt"))

textOut <- capture.output(table(Pam50_subtypes_combatNorm$subtypes[,"subtypeLabels"]))
cat(textOut,sep="\n",file=paste0(saveDir,"/mergedcombatNormPam50Short_subtypeTable.txt"))
####

########noNorm plotting
load("/home/kplaney/breast_analysis/curatedBreastData_kmeansConsensus_mergedNoNorm_pam50Short_nstart12015-05-04.RData.gzip")

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
  theme(plot.title=element_text(colour="black",size=12,vjust=1))
#above: font size needs to be smaller for subtypes
png(filename=paste0(saveDir,"/mergedNoNorm_subtype_breakdowns_",Sys.Date(),".png"),width=1000,height=1000,res=160)

plot( plotG)

dev.off()



########combat plotting
load("/home/kplaney/breast_analysis/curatedBreastData_kmeansConsensus_mergedCombatNorm_pam50Short_nstart12015-05-04.RData.gzip")


attrDF <- data.frame(c(rep.int(1,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]])),
                       rep.int(2,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]])),
                       rep.int(3,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[3]]))
),
c(Pam50_subtypes_combatNorm$subtype[kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]] ,"subtypeLabels"],
  Pam50_subtypes_combatNorm$subtype[kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]] ,"subtypeLabels"],
  Pam50_subtypes_combatNorm$subtype[kmeansConsensus$clustSampleIndexList_PACR[[1]][[3]] ,"subtypeLabels"]))

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
  theme(plot.title=element_text(colour="black",size=12,vjust=1))
#above: font size needs to be smaller for subtypes
png(filename=paste0(saveDir,"/mergedCombat_subtype_breakdowns_",Sys.Date(),".png"),width=1000,height=1000,res=160)

plot( plotG)

dev.off()

####BMC norm plotting
load("/home/kplaney/breast_analysis/curatedBreastData_kmeansConsensus_mergedBMCNorm_pam50Short_nstart12015-05-04.RData.gzip")

#cluster indices match up with merged mattrix; this is what we fed in the clustering algorithm
#only 2 clusters!
attrDF <- data.frame(c(rep.int(1,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]])),
                       rep.int(2,length(kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]]))
                       ),
c(Pam50_subtypes_BMCNorm$subtype[kmeansConsensus$clustSampleIndexList_PACR[[1]][[1]] ,"subtypeLabels"],
  Pam50_subtypes_BMCNorm$subtype[kmeansConsensus$clustSampleIndexList_PACR[[1]][[2]] ,"subtypeLabels"]
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
  theme(plot.title=element_text(colour="black",size=12,vjust=1))
#above: font size needs to be smaller for subtypes
png(filename=paste0(saveDir,"/mergedBMC_subtype_breakdowns_",Sys.Date(),".png"),width=1000,height=1000,res=160)

plot( plotG)

dev.off()


########inspecting how normalization affects pam50 groups:


