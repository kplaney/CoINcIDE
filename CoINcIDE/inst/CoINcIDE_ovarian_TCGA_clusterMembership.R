
TCGA_clusterAssign <- read.table("/home/kplaney/ovarian_analysis//gdac.broadinstitute.org_OV-TP.Aggregate_Molecular_Subtype_Clusters.Level_4.2014101700.0.0/OV-TP.mergedcluster.txt",
                                 header=TRUE)


TCGA_clusterAssign_NMF_cHclust <- TCGA_clusterAssign[,c(1:3)]
#remove last number, and replace "-" with "."
library("limma")
TCGA_clusterAssign_NMF_cHclust[,1] <- substr(TCGA_clusterAssign_NMF_cHclust[,1],start=1,stop=12)
TCGA_clusterAssign_NMF_cHclust[,1] <- gsub("-",".",TCGA_clusterAssign_NMF_cHclust[,1] )
load("/home/kplaney/ovarian_analysis/curatedOvarianData_kmeansConsensus_nstart1_200Features_2015-04-28.RData.gzip")

tmpSampleName <- c()
tmpClustAssign <-c()

#TCGA is dataset 24 here.
for(c in 1:length(kmeansConsensus$clustSampleIndexList_PACR[[24]])){
  
  tmpSampleName <- append(tmpSampleName,names(kmeansConsensus$clustSampleIndexList_PACR[[24]][[c]]))
  tmpClustAssign <- append(tmpClustAssign,rep.int(c,times=length(kmeansConsensus$clustSampleIndexList_PACR[[24]][[c]])))
  
}

TCGA_200F_cKmeans <- cbind(tmpSampleName,tmpClustAssign)

##now do 2000 features
load("/home/kplaney/ovarian_analysis/curatedOvarianData_kmeansConsensus_nstart1_2000Features_2015-04-30.RData.gzip")

tmpSampleName <- c()
tmpClustAssign <-c()

#TCGA is dataset 23 here.
for(c in 1:length(kmeansConsensus$clustSampleIndexList_PACR[[23]])){
  
  tmpSampleName <- append(tmpSampleName,names(kmeansConsensus$clustSampleIndexList_PACR[[23]][[c]]))
  tmpClustAssign <- append(tmpClustAssign,rep.int(c,times=length(kmeansConsensus$clustSampleIndexList_PACR[[23]][[c]])))
  
}

TCGA_2000F_cKmeans <- cbind(tmpSampleName,tmpClustAssign)

#we should have the same # patients
nrow(TCGA_2000F_cKmeans)==nrow(TCGA_200F_cKmeans)

#hmmm...more patients in the original clustering - not sure why, but these were removed from
#the curatedOvarianData, so remove them here.
#596
#ha some patients are NA - perhaps this is because they are normal patients??
#I removed these anyways for CoINcIDE
nrow(TCGA_clusterAssign_NMF_cHclust)
#569
nrow(TCGA_2000F_cKmeans)

#in CoINcIDE clustering runs, all patients have a match but aren't in the same order
all(TCGA_2000F_cKmeans[,1]==TCGA_200F_cKmeans[,1])
#[1] FALSE
 any(is.na(match(TCGA_2000F_cKmeans[,1],TCGA_200F_cKmeans[,1])))
#[1] FALSE

TCGA_200F_cKmeans <- TCGA_200F_cKmeans[na.omit(match(TCGA_2000F_cKmeans[,1],TCGA_200F_cKmeans[,1])), ]
all(TCGA_2000F_cKmeans[,1]==TCGA_200F_cKmeans[,1])
#remove patients that are in TCGA's clustering but not ours. they may have also been removed by combat.

#569 patients left
dim(TCGA_clusterAssign_NMF_cHclust)

TCGA_clusterAssign_NMF_cHclust <- TCGA_clusterAssign_NMF_cHclust[na.omit(match(TCGA_200F_cKmeans[,1],
                                                                       TCGA_clusterAssign_NMF_cHclust[,1])), ]

masterTable <- cbind(TCGA_2000F_cKmeans,TCGA_200F_cKmeans[,2],TCGA_clusterAssign_NMF_cHclust[,c(2:3)])
colnames(masterTable) <- c("Patient","cKmeans_2000F","cKmeans_200F","cNMF","cHclust")


plotData <- data.matrix(masterTable[,c(2:5)])
#why do 2 patients have NA data?
which(is.na(plotData[,3]))
#501 483 
#168 502 
which(is.na(plotData[,4]))
#501 483 
#168 502 
#remove these two patients
plotData <- plotData[-which(is.na(plotData[,4])), ]
library("gplots")
library("RColorBrewer")
colorCodeF <- brewer.pal(8,"Dark2")
heatmap.2(t(plotData),Rowv=NULL,Colv=NULL,dendrogram='none',cexRow=1,srtRow=45, trace="none",
          scale="none",notecol="black",key=FALSE,col=colorCodeF)
#hmmm...we can see that clusters correspond well, but the actual cluster labels (1,2,3) are switched sometimes
#hand-tune this just for the plot so blocks of color are consistent

#cNMF: change cluster 2 to cluster 1.  change cluster 1 to cluster 2.
#cHclust: change cluster 2 to cluster 3. change cluster 3 to cluster 2.

plotDataOrig <- plotData
plotData[which(plotDataOrig[,3]==1) ,3] <- 2
plotData[which(plotDataOrig[,3]==2) ,3] <- 1
plotData[which(plotDataOrig[,4]==2) ,4] <- 3
plotData[which(plotDataOrig[,4]==3) ,4] <- 2

#now we're all good!

heatmap.2(t(plotData),Rowv=NULL,Colv=NULL,dendrogram='none',cexRow=1,srtRow=45, trace="none",
          scale="none",notecol="black",key=FALSE,col=colorCodeF)

png(filename=paste0("/home/kplaney/ovarian_analysis/TCGA_clusterMembership_",Sys.Date(),".png"),width=1000,height=800,res=160)
#,width = 700, height = 1000,res=160);

heatmap.2(t(plotData),Rowv=NULL,Colv=NULL,dendrogram='none',cexRow=1,srtRow=45, trace="none",
          scale="none",notecol="black",key=FALSE,col=colorCodeF)

dev.off()



