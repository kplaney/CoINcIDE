
#####NOTE: pam50 merged clustering is in merged script.
###pam50 k-means is in the breast test K script. We just add on centroid clustering and extra gap test runs below.

####make subtype dataframe first
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_visualization.R")
for(d in 1:length(dataMatrixList)){
  
  tmp <- assignCentroidSubtype(t(dataMatrixList[[d]]),minNumGenes=30,centroidRData="/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData");
 
  tmp <- data.frame(colnames(dataMatrixList[[d]]),rep.int(d,times=ncol(dataMatrixList[[d]])),as.numeric(as.character(tmp$subtypes[,1])),tmp$subtypes[,2],stringsAsFactors=FALSE)
  
  if(d>1){
    
   subtypeDF <- rbind(subtypeDF,tmp)
   
  }else{
    
    subtypeDF <- tmp
    
  }
}

colnames(subtypeDF) <- c("sampleName","studyNum","subtypeNum","subtype")
save(subtypeDF,file="/home/kplaney/breast_analysis/pam50Full_subtypeDF.RData.gzip",compress="gzip")

load("/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData")
load("/home/kplaney/breast_analysis/pam50Short_genes.RData")
centroidMatrix <- centroidMatrix[na.omit(match(pam50Short,centroidMatrix[,1])), ]
save(centroidMatrix,file="/home/kplaney/breast_analysis/pam50Short_updatedSymbols_centroidMatrix.RData")
for(d in 1:length(dataMatrixList)){
  
  tmp <- assignCentroidSubtype(t(dataMatrixList[[d]]),minNumGenes=30,centroidRData="/home/kplaney/breast_analysis/pam50Short_updatedSymbols_centroidMatrix.RData");
  
  tmp <- data.frame(colnames(dataMatrixList[[d]]),rep.int(d,times=ncol(dataMatrixList[[d]])),as.numeric(as.character(tmp$subtypes[,1])),tmp$subtypes[,2],stringsAsFactors=FALSE)
  
  if(d>1){
    
    subtypeDF_short <- rbind(subtypeDF_short,tmp)
    
  }else{
    
    subtypeDF_short <- tmp
    
  }
}

colnames(subtypeDF_short) <- c("sampleName","studyNum","subtypeNum_short","subtype_short")
#all rows ordered correctly?
all(subtypeDF$sampleName==subtypeDF_short$sampleName)
subtypeDF_master <- cbind(subtypeDF,subtypeDF_short[,c(3:4)])

save(subtypeDF_master,file="/home/kplaney/breast_analysis/pam50FullAndShort_subtypeDF.RData.gzip",compress="gzip")

#NOTE: around ~300 samples are re-categorized if use short gene list.
length(which(subtypeDF_master$subtype!=subtypeDF_master$subtype_short))
table(subtypeDF_master$subtype)
#basically: less basal, luminal A with short subtypings:
table(subtypeDF_master$subtype_short)


#########centroid subtyping (only use pam50 full):
load("/home/kplaney/breast_analysis/pam50FullAndShort_subtypeDF.RData.gzip")

subtype_studySplit <- split(subtypeDF_master[,c("subtype","sampleName")],f=subtypeDF_master[,"studyNum"])
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

clustSampleIndexList <- list()
clustFeatureIndexList <- list()

length(subtype_studySplit)==length(dataMatrixList)


load("/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData");
pam50Genes <- centroidMatrix[ ,1];

for(d in 1:length(dataMatrixList)){
  
  clustSampleIndexList[[d]] <- list()
  clustFeatureIndexList[[d]] <- list()
  studySubtypes <- unique(subtype_studySplit[[d]][,"subtype"])
  
  featureIndices <- na.omit(match(pam50Genes,rownames(dataMatrixList[[d]])))
  
  for(s in 1:length(studySubtypes)){
    
    clustFeatureIndexList[[d]][[s]] <- featureIndices
    clustSampleIndexList[[d]][[s]] <- na.omit(match(
      subtype_studySplit[[d]][which(subtype_studySplit[[d]][,"subtype"]==studySubtypes[s]),"sampleName"],
      colnames(dataMatrixList[[d]])))
    
  }
  
  
}

  
  clustFeaturesList[[d]] <- pam50GeneShort
  
}

saveDir <- "/home/kplaney/breast_analysis/"
kmeansConsensus<- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                         pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                         numSims=500,maxNumClusters=10,
                                         outputFile=paste0(saveDir,"/",numFeatures,"F_clust.txt"),distMethod=c("euclidean"),
                                         hclustAlgorithm=c("average"),
                                         consensusHclustAlgorithm=c("average"),
                                         minClustConsensus=.7, minMeanClustConsensus=.85,
                                         corUse="everything",pItem=.9,maxPAC=.15)

save(kmeansConsensus,file=paste0(saveDir,"/curatedbreastData_kmeansConsensus_nstart1pItem9",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")

#############################now try gap test for sh
