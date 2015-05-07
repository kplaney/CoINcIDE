
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

#no need to do k-means here; we created naive clustering assignments with the pam50 centroids
saveDir <- "/home/kplaney/breast_analysis/"
pam50Full_centroidCluster <- list()
pam50Full_centroidCluster$clustSampleIndexList <- clustSampleIndexList
pam50Full_centroidCluster$clustFeatureIndexList <- clustFeatureIndexList 
save(pam50Full_centroidCluster,file=paste0(saveDir,"/pam50Full_centroidCluster.RData.gzip"),compress="gzip")

#############################now try gap test not for actual k-means pam50 clustering (consensus is in selectK script.)

##with short, nstart=25
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
load("/home/kplaney/breast_analysis/pam50Short_genes.RData")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")
saveDir <- "/home/kplaney/breast_analysis/"


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- pam50Short
  
}

#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansGapTest_pam50_short_Nstart25 <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                                                    pickKMethod=c("gap"),iter.max=20,nstart=25,
                                                                    numSims=500,maxNumClusters=15,
                                                                    outputFile="/home/kplaney/breast_analysis/test.txt"
                                                                  )


save(kmeansGapTest_pam50_short_Nstart25,
     file=paste0(saveDir,"/gapTestKmeans_pam50Short_nstart25",Sys.Date(),".RData.gzip"),compress="gzip")


###with short, nstart =1
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
load("/home/kplaney/breast_analysis/pam50Short_genes.RData")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")
saveDir <- "/home/kplaney/breast_analysis/"


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- pam50Short
  
}

#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansGapTest_pam50_short_Nstart1 <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                                             pickKMethod=c("gap"),iter.max=20,nstart=1,
                                                             numSims=500,maxNumClusters=15,
                                                             outputFile="/home/kplaney/breast_analysis/test.txt"
)


save(kmeansGapTest_pam50_short_Nstart1,
     file=paste0(saveDir,"/gapTestKmeans_pam50Short_nstart1",Sys.Date(),".RData.gzip"),compress="gzip")


####with full, nstart=25
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
load("/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData")
pam50GenesFull <- centroidMatrix[,1]
saveDir <- "/home/kplaney/breast_analysis/"
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")
numFeatures <- "pam50Full"



clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- pam50GenesFull
  
}


kmeansGapTest_pam50_full_Nstart25 <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                                             pickKMethod=c("gap"),iter.max=20,nstart=25,
                                                             numSims=500,maxNumClusters=15,
                                                             outputFile="/home/kplaney/breast_analysis/test.txt"
)


save(kmeansGapTest_pam50_full_Nstart25,
     file=paste0(saveDir,"/gapTestKmeans_pam50Full_nstart25",Sys.Date(),".RData.gzip"),compress="gzip")





####with full, nstart=1
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
load("/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData")
pam50GenesFull <- centroidMatrix[,1]
saveDir <- "/home/kplaney/breast_analysis/"
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")
numFeatures <- "pam50Full"


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- pam50GenesFull
  
}

#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansGapTest_pam50_full_Nstart1 <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                                            pickKMethod=c("gap"),iter.max=20,nstart=1,
                                                            numSims=500,maxNumClusters=15,
                                                            outputFile="/home/kplaney/breast_analysis/test.txt"
)


save(kmeansGapTest_pam50_full_Nstart1,
     file=paste0(saveDir,"/gapTestKmeans_pam50Full_nstart1",Sys.Date(),".RData.gzip"),compress="gzip")
