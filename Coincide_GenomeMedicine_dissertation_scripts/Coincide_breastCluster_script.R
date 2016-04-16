
#clustering the beast cancer datasets using various feature sets
#clustering method was selected using the Coincide_selectK_cluster_script.R

library("Coincide")

##CHANGE these paths to your user directory
saveDirGlobal <- "/home/ywrfc09/breast_analysis/"
saveDir_PAM50 <- "/home/ywrfc09/breast_analysis/PAM50_analyses/"
saveDir_20 <- "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes"
saveDir_no20 <-  "/home/ywrfc09/breast_analysis/metaRankNoTop20Genes"
saveDirMerged <- "/home/ywrfc09/breast_analysis/mergedMatrix/"
outputFile <- paste0(saveDirMerged ,"/merged_outMessages.txt")

dataMatrixList <- readRDS(paste0(saveDirGlobal,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.rds"))

#below named "centroidMatrix" when load it:
load("pam50_centroids_updatedSymbols.RData")
load("pam50Short_genes.RData")

###NOTE: pam50 consensus clustering with PAM50 genes using Coincide
#was run in the "select K" script. Merged matrix clustering was run in the 
#breastMergedAnalysis_script

###so just do pam50 semi-supervised centroid clustering here:

####make subtype dataframe first - i.e. assign our samples to centroid groups
for(d in 1:length(dataMatrixList)){
  
  tmp <- assignCentroidSubtype(t(dataMatrixList[[d]]),minNumGenes=30,centroidRData="pam50_centroids_updatedSymbols.RData");
  
  tmp <- data.frame(colnames(dataMatrixList[[d]]),rep.int(d,times=ncol(dataMatrixList[[d]])),as.numeric(as.character(tmp$subtypes[,1])),tmp$subtypes[,2],stringsAsFactors=FALSE)
  
  if(d>1){
    
    subtypeDF <- rbind(subtypeDF,tmp)
    
  }else{
    
    subtypeDF <- tmp
    
  }
}


colnames(subtypeDF) <- c("sampleName","studyNum","subtypeNum","subtype")

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

save(subtypeDF_master,file=paste0(saveDirGlobal,"pam50FullAndShort_subtypeDF.RData.gzip"),compress="gzip")

#NOTE: around ~300 samples are re-categorized if use short gene list.
length(which(subtypeDF_master$subtype!=subtypeDF_master$subtype_short))
table(subtypeDF_master$subtype)
#basically: less basal, luminal A with short subtypings:
table(subtypeDF_master$subtype_short)


#########centroid subtyping (only use pam50 full):
subtype_studySplit <- split(subtypeDF_master[,c("subtype","sampleName")],f=subtypeDF_master[,"studyNum"])

clustSampleIndexList <- list()
clustFeatureIndexList <- list()

length(subtype_studySplit)==length(dataMatrixList)


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
  clustFeaturesList[[d]] <- pam50GeneShort
  
}



#no need to do k-means here; we created naive clustering assignments with the pam50 centroids

pam50Full_centroidCluster <- list()
pam50Full_centroidCluster$clustSampleIndexList <- clustSampleIndexList
pam50Full_centroidCluster$clustFeatureIndexList <- clustFeatureIndexList 
save(pam50Full_centroidCluster,file=paste0(saveDir_PAM50,"/pam50Full_centroidCluster.RData.gzip"),compress="gzip")



###now on to true clustering:
###########do for each number of features (just change numFeatures variable)
#because it gave much clearer clusters, used the meta-features with top 20
#genes by variance frome each dataset.
numFeatures <- 2000

metaFeatures <- readRDS(paste0(saveDir_20,"/metaFeatures_",numFeatures,".rds"))
##do for each 200,500,1000,2000 (load different metaFeatures RData object each time.)
#remove datasets with too many missing top gene features


if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- metaFeatures$finalFeatures
  
}

#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansConsensus<- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                                                    pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                                                    numSims=500,maxNumClusters=10,
                                                                    outputFile=paste0(saveDirGlobal,"/",numFeatures,"F_clust.txt"),distMethod=c("euclidean"),
                                                                    hclustAlgorithm=c("average"),
                                                                    consensusHclustAlgorithm=c("average"),
                                                                    minClustConsensus=.7, minMeanClustConsensus=.85,
                                                                    corUse="everything",pItem=.9,maxPAC=.15)

saveRDS(kmeansConsensus,file=paste0(saveDirGlobal,"/curatedbreastData_kmeansConsensus_nstart1pItem9",numFeatures,"Features_",Sys.Date(),".rds"),compress=TRUE)

#######1000
numFeatures <- 1000
metaFeatures <- readRDS(paste0(saveDir_20,"/metaFeatures_",numFeatures,".rds"))
##do for each 200,500,1000,2000 (load different metaFeatures RData object each time.)
#remove datasets with too many missing top gene features

if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- metaFeatures$finalFeatures
  
}


#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansConsensus<- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                         pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                         numSims=500,maxNumClusters=10,
                                         outputFile=paste0(saveDirGlobal,"/",numFeatures,"F_clust.txt"),distMethod=c("euclidean"),
                                         hclustAlgorithm=c("average"),
                                         consensusHclustAlgorithm=c("average"),
                                         minClustConsensus=.7, minMeanClustConsensus=.85,
                                         corUse="everything",pItem=.9,maxPAC=.15)

saveRDS(kmeansConsensus,file=paste0(saveDirGlobal,"/curatedbreastData_kmeansConsensus_nstart1pItem9",numFeatures,"Features_",Sys.Date(),".rds"),compress=TRUE)


####500
numFeatures <- 500
metaFeatures <- readRDS(paste0(saveDir_20,"/metaFeatures_",numFeatures,".rds"))
##do for each 200,500,1000,2000 (load different metaFeatures RData object each time.)
#remove datasets with too many missing top gene features


if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- metaFeatures$finalFeatures
  
}

#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansConsensus<- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                         pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                         numSims=500,maxNumClusters=10,
                                         outputFile=paste0(saveDirGlobal,"/",numFeatures,"F_clust.txt"),distMethod=c("euclidean"),
                                         hclustAlgorithm=c("average"),
                                         consensusHclustAlgorithm=c("average"),
                                         minClustConsensus=.7, minMeanClustConsensus=.85,
                                         corUse="everything",pItem=.9,maxPAC=.15)

saveRDS(kmeansConsensus,file=paste0(saveDirGlobal,"/curatedbreastData_kmeansConsensus_nstart1pItem9",numFeatures,"Features_",Sys.Date(),".rds"),compress=TRUE)


####200
numFeatures <- 200
metaFeatures <- readRDS(paste0(saveDir_20,"/metaFeatures_",numFeatures,".rds"))
##do for each 200,500,1000,2000 (load different metaFeatures RData object each time.)
#remove datasets with too many missing top gene features


if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- metaFeatures$finalFeatures
  
}


#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansConsensus<- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                         pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                         numSims=500,maxNumClusters=10,
                                         outputFile=paste0(saveDirGlobal,"/",numFeatures,"F_clust.txt"),distMethod=c("euclidean"),
                                         hclustAlgorithm=c("average"),
                                         consensusHclustAlgorithm=c("average"),
                                         minClustConsensus=.7, minMeanClustConsensus=.85,
                                         corUse="everything",pItem=.9,maxPAC=.15)

saveRDS(kmeansConsensus,file=paste0(saveDirGlobal,"/curatedbreastData_kmeansConsensus_nstart1pItem9",numFeatures,"Features_",Sys.Date(),".rds"),compress=TRUE)


###250,300
numFeatures <- 250
metaFeatures <- readRDS(paste0(saveDir_20,"/metaFeatures_",numFeatures,".rds"))
##do for each 200,500,1000,2000 (load different metaFeatures RData object each time.)
#remove datasets with too many missing top gene features


if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- metaFeatures$finalFeatures
  
}


#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansConsensus<- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                         pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                         numSims=500,maxNumClusters=10,
                                         outputFile=paste0(saveDirGlobal,"/",numFeatures,"F_clust.txt"),distMethod=c("euclidean"),
                                         hclustAlgorithm=c("average"),
                                         consensusHclustAlgorithm=c("average"),
                                         minClustConsensus=.7, minMeanClustConsensus=.85,
                                         corUse="everything",pItem=.9,maxPAC=.15)

saveRDS(kmeansConsensus,file=paste0(saveDirGlobal,"/curatedbreastData_kmeansConsensus_nstart1pItem9",numFeatures,"Features_",Sys.Date(),".rds"),compress=TRUE)

##300
numFeatures <- 300
metaFeatures <- readRDS(paste0(saveDir_20,"/metaFeatures_",numFeatures,".rds"))
##do for each 200,500,1000,2000 (load different metaFeatures RData object each time.)
#remove datasets with too many missing top gene features


if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- metaFeatures$finalFeatures
  
}


#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansConsensus<- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                         pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                         numSims=500,maxNumClusters=10,
                                         outputFile=paste0(saveDirGlobal,"/",numFeatures,"F_clust.txt"),distMethod=c("euclidean"),
                                         hclustAlgorithm=c("average"),
                                         consensusHclustAlgorithm=c("average"),
                                         minClustConsensus=.7, minMeanClustConsensus=.85,
                                         corUse="everything",pItem=.9,maxPAC=.15)

saveRDS(kmeansConsensus,file=paste0(saveDirGlobal,"/curatedbreastData_kmeansConsensus_nstart1pItem9",numFeatures,"Features_",Sys.Date(),".rds"),compress=TRUE)


#########kmeans again, but NOT including PAM50 genes
#why? want to see if removing these well-known gene still gives "good" clusters
##shoot for around 70 features, that should give us ~50 without PAM50?
#need to run meta-feature script again then.

#ran meta-feature analysis for 1000,500,1000,2000 features.

#OK: so really only 4 PAM50 genes in this top set anyways...
metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=54,
                                                   numTopFeatFromEachDataset=0,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile="/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/selectFeaturesMetaVarianceOut.txt")

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))


PAM50genes <- centroidMatrix[,"geneSymbol"]

#remove any 50 genes.
noPAM50_set <- setdiff(metaFeatures$finalFeatures,PAM50genes)
#out of curiosity: which ones DID intersect?
#[1] "NAT1"  "KRT14" "SFRP1" "ESR1" 
intersect(metaFeatures$finalFeatures,PAM50genes)
#OK back to non-PAM50 analysis:
metaFeatures$finalFeatures <- noPAM50_set
saveRDS(metaFeatures,file=paste0(globalSaveDir,"/metaFeatures_50_NOPAM50_NOTOP20GENES.rds"),compress=TRUE)
write.table(metaFeatures$finalFeatures,quote=FALSE,file=paste0(saveDirGlobal,"metaFeatures_50_NOPAM50_NOTOP20GENES.txt"))


if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- metaFeatures$finalFeatures
  
}

numFeatures <- 50
#we know these are strong clusters. have  minMeanClustConsensus around .85
#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansConsensus<- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                         pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                         numSims=500,maxNumClusters=10,
                                         outputFile=paste0(saveDirGlobal,"/",numFeatures,"F_clust.txt"),distMethod=c("euclidean"),
                                         hclustAlgorithm=c("average"),
                                         consensusHclustAlgorithm=c("average"),
                                         minClustConsensus=.7, minMeanClustConsensus=.85,
                                         corUse="everything",pItem=.9,maxPAC=.15)

saveRDS(kmeansConsensus,file=paste0(saveDirGlobal,"/curatedbreastData_kmeansConsensus_nstart1pItem9",numFeatures,"Features_NOPAM50_50genesNOTOP20_",Sys.Date(),".rds"),compress=TRUE)



#200 features
#ALSO: no PAM50
numFeatures <- 200
metaFeatures <- readRDS(paste0(saveDir_20,"/metaFeatures_",numFeatures,".rds"))

PAM50genes <- centroidMatrix[,"geneSymbol"]

noPAM50_200MetaF <- setdiff(metaFeatures$finalFeatures,PAM50genes)

if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

write.table(noPAM50_200MetaF,quote=FALSE,row.names=FALSE, file=paste0(saveDirGlobal,"metaFeatures_200F_NOPAM50.txt"))
clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- noPAM50_200MetaF
  
}


#we know these are strong clusters. have  minMeanClustConsensus around .85
#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansConsensus<- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                         pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                         numSims=500,maxNumClusters=10,
                                         outputFile=paste0(saveDirGlobal,"/",numFeatures,"F_clust.txt"),distMethod=c("euclidean"),
                                         hclustAlgorithm=c("average"),
                                         consensusHclustAlgorithm=c("average"),
                                         minClustConsensus=.7, minMeanClustConsensus=.85,
                                         corUse="everything",pItem=.9,maxPAC=.15)

saveRDS(kmeansConsensus,file=paste0(saveDirGlobal,"/curatedbreastData_kmeansConsensus_nstart1pItem9",numFeatures,"Features_NOPAM50_",Sys.Date(),".rds"),compress=TRUE)


##try 500, no PAM50
numFeatures <- 500
metaFeatures <- readRDS(paste0(saveDir_20,"/metaFeatures_",numFeatures,".rds"))

PAM50genes <- centroidMatrix[,"geneSymbol"]

noPAM50_500MetaF <- setdiff(metaFeatures$finalFeatures,PAM50genes)

if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- noPAM50_500MetaF
  
}

#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansConsensus<- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                         pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                         numSims=500,maxNumClusters=10,
                                         outputFile=paste0(saveDirGlobal,"/",numFeatures,"F_clust.txt"),distMethod=c("euclidean"),
                                         hclustAlgorithm=c("average"),
                                         consensusHclustAlgorithm=c("average"),
                                         minClustConsensus=.7, minMeanClustConsensus=.85,
                                         corUse="everything",pItem=.9,maxPAC=.15)

saveRDS(kmeansConsensus,file=paste0(saveDirGlobal,"/curatedbreastData_kmeansConsensus_nstart1pItem9",numFeatures,"Features_NOPAM50_",Sys.Date(),".rds"),compress=TRUE)


##try 1000, no PAM50
numFeatures <- 1000

metaFeatures <- readRDS(paste0(saveDir_20,"/metaFeatures_",numFeatures,".rds"))

PAM50genes <- centroidMatrix[,"geneSymbol"]

noPAM50_500MetaF <- setdiff(metaFeatures$finalFeatures,PAM50genes)

if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- noPAM50_500MetaF
  
}

#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansConsensus<- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                         pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                         numSims=500,maxNumClusters=10,
                                         outputFile=paste0(saveDirGlobal,"/",numFeatures,"F_clust.txt"),distMethod=c("euclidean"),
                                         hclustAlgorithm=c("average"),
                                         consensusHclustAlgorithm=c("average"),
                                         minClustConsensus=.7, minMeanClustConsensus=.85,
                                         corUse="everything",pItem=.9,maxPAC=.15)

saveRDS(kmeansConsensus,file=paste0(saveDirGlobal,"/curatedbreastData_kmeansConsensus_nstart1pItem9",numFeatures,"Features_NOPAM50_",Sys.Date(),".rds"),compress=TRUE)



##try 2000, no PAM50
numFeatures <- 2000

metaFeatures <- readRDS(paste0(saveDir_20,"/metaFeatures_",numFeatures,".rds"))

PAM50genes <- centroidMatrix[,"geneSymbol"]

noPAM50_2000MetaF <- setdiff(metaFeatures$finalFeatures,PAM50genes)


if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}
write.table(noPAM50_2000MetaF,quote=FALSE,row.names=FALSE, file=paste0(saveDirGlobal,"metaFeatures_2000F_NOPAM50.txt"))

clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- noPAM50_500MetaF
  
}

#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansConsensus<- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                         pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                         numSims=500,maxNumClusters=10,
                                         outputFile=paste0(saveDirGlobal,"/",numFeatures,"F_clust.txt"),distMethod=c("euclidean"),
                                         hclustAlgorithm=c("average"),
                                         consensusHclustAlgorithm=c("average"),
                                         minClustConsensus=.7, minMeanClustConsensus=.85,
                                         corUse="everything",pItem=.9,maxPAC=.15)

saveRDS(kmeansConsensus,file=paste0(saveDirGlobal,"/curatedbreastData_kmeansConsensus_nstart1pItem9",numFeatures,"Features_NOPAM50_",Sys.Date(),".rds"),compress=TRUE)

