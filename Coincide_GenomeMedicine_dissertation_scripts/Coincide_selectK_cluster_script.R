
##CHANGE these paths to your user directory
library("CoINcIDE")
saveDirGlobal <- "/home/ywrfc09/breast_analysis/"
saveDir_PAM50 <- "/home/ywrfc09/breast_analysis/PAM50_analyses/"
saveDir_20 <- "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes"
saveDir_no20 <-  "/home/ywrfc09/breast_analysis/metaRankNoTop20Genes"
outputFile <- "/home/kplaney/breast_analysis/clust_test_outMessages.txt"


####look at expected K for each dataset
##this script follows the breast breastProcessAndGeneFeatures_script.R
library("Biobase")
#load all processed eset data that has clinical data
#load("/home/ywrfc09/breast_analysis//curatedBreastData_esets_proc.RData.gzip")
library("curatedBreastData")
data("curatedBreastDataExprSetList")


###first: broadly look at the pam50 panel table stats and ER and HER2 IHC
#table stats for a study (when they were measured.) These are helpful proxies
#to guess on a broad level about how many clusters should be in the dataset.

#look at MD-Anderson dataset first.
table(pData(esets[["study_25055_GPL96_MDACC_M"]])$pam50)
table(pData(esets[["study_25055_GPL96_MDACC_M"]])$ER_preTrt)
table(pData(esets[["study_25055_GPL96_MDACC_M"]])$HER2_preTrt)


table(pData(esets[["study_22226_GPL1708_all"]])$pam50)
table(pData(esets[["study_22226_GPL1708_all"]])$ER_preTrt)
table(pData(esets[["study_22226_GPL1708_all"]])$HER2_preTrt)


table(pData(esets[["study_19615_GPL570_all"]])$pam50)
table(pData(esets[["study_19615_GPL570_all"]])$ER_preTrt)
table(pData(esets[["study_19615_GPL570_all"]])$HER2_preTrt)

table(pData(esets[["study_12093_GPL96_all"]])$pam50)
table(pData(esets[["study_12093_GPL96_all"]])$ER_preTrt)
table(pData(esets[["study_12093_GPL96_all"]])$HER2_preTrt)

table(pData(esets[["study_25065_GPL96_MDACC"]])$pam50)
table(pData(esets[["study_25065_GPL96_MDACC"]])$ER_preTrt)
table(pData(esets[["study_25065_GPL96_MDACC"]])$HER2_preTrt)

table(pData(esets[["study_20181_GPL96_all"]])$pam50)
table(pData(esets[["study_20181_GPL96_all"]])$ER_preTrt)
table(pData(esets[["study_20181_GPL96_all"]])$HER2_preTrt)

table(pData(esets[["study_16446_GPL570_all"]])$pam50)
table(pData(esets[["study_16446_GPL570_all"]])$ER_preTrt)
table(pData(esets[["study_16446_GPL570_all"]])$HER2_preTrt)

table(pData(esets[["study_2034_GPL96_all"]])$pam50)
table(pData(esets[["study_2034_GPL96_all"]])$ER_preTrt)
table(pData(esets[["study_2034_GPL96_all"]])$HER2_preTrt)


###consensus clustering on pam50 genes.

load(paste0(saveDirGlobal,
            "/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
load("pam50Short_genes.RData")

numFeatures <- "pam50Short"

clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- pam50Short
  
}


#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansConsensuspam50_short_Nstart15pItem9 <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                                                    pickKMethod=c("consensus"),iter.max=20,nstart=15,
                                                                    numSims=500,maxNumClusters=10,
                                                                    outputFile=outputFile,distMethod=c("euclidean"),
                                                                    hclustAlgorithm=c("average"),
                                                                    consensusHclustAlgorithm=c("average"),
                                                                    minClustConsensus=.7, minMeanClustConsensus=.85,
                                                                    corUse="everything",pItem=.9,maxPAC=.15)


save(kmeansConsensuspam50_short_Nstart15pItem9,
     file=paste0(saveDir_PAM50,"/kmeansConsensuspam50_short_Nstart15pItem9_",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")


kmeansConsensuspam50_short_Nstart1pItem9 <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                                                   pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                                                   numSims=500,maxNumClusters=10,
                                                                   outputFile=outputFile,distMethod=c("euclidean"),
                                                                   hclustAlgorithm=c("average"),
                                                                   consensusHclustAlgorithm=c("average"),
                                                                   minClustConsensus=.7, minMeanClustConsensus=.85,
                                                                   corUse="everything",pItem=.9,maxPAC=.15)


save(kmeansConsensuspam50_short_Nstart1pItem9,
     file=paste0(saveDir_PAM50,"/kmeansConsensuspam50_short_Nstart1pItem9_",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")


kmeansConsensuspam50_short_Nstart15pItem8 <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                                                    pickKMethod=c("consensus"),iter.max=20,nstart=15,
                                                                    numSims=500,maxNumClusters=10,
                                                                    outputFile=outputFile,distMethod=c("euclidean"),
                                                                    hclustAlgorithm=c("average"),
                                                                    consensusHclustAlgorithm=c("average"),
                                                                    minClustConsensus=.7, minMeanClustConsensus=.85,
                                                                    corUse="everything",pItem=.8,maxPAC=.15)


save(kmeansConsensuspam50_short_Nstart15pItem8 ,
     file=paste0(saveDir_PAM50,"/kmeansConsensuspam50_short_Nstart15pItem8_",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")


##hierarchical
hclustConsensuspam50_short_pItem9 <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("hc"),
                                                            pickKMethod=c("consensus"),iter.max=20,nstart=15,
                                                            numSims=500,maxNumClusters=10,
                                                            outputFile=outputFile,distMethod=c("euclidean"),
                                                            hclustAlgorithm=c("ward.D"),
                                                            consensusHclustAlgorithm=c("average"),
                                                            minClustConsensus=.7, minMeanClustConsensus=.85,
                                                            corUse="everything",pItem=.9,maxPAC=.15)


save(hclustConsensuspam50_short_pItem9,file=paste0(saveDir_PAM50,"/hclustConsensuspam50_short_pItem9_",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")

hclustConsensuspam50_short_pItem8 <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("hc"),
                                                            pickKMethod=c("consensus"),iter.max=20,nstart=15,
                                                            numSims=500,maxNumClusters=10,
                                                            outputFile=outputFile,distMethod=c("euclidean"),
                                                            hclustAlgorithm=c("ward.D"),
                                                            consensusHclustAlgorithm=c("average"),
                                                            minClustConsensus=.7, minMeanClustConsensus=.85,
                                                            corUse="everything",pItem=.8,maxPAC=.15)


save(hclustConsensuspam50_short_pItem8,file=paste0(saveDir_PAM50,"/hclustConsensuspam50_short_pItem8_",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")

#############
####pam50 full
load(paste0(saveDirGlobal,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
load("pam50_centroids_updatedSymbols.RData")
pam50GenesFull <- centroidMatrix[,1]

numFeatures <- "pam50Full"


#COME BACK: remove clusters with too few pam50 genes?

clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- pam50GenesFull
  
}

kmeansConsensuspam50_full_Nstart15pItem9 <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                                                   pickKMethod=c("consensus"),iter.max=20,nstart=15,
                                                                   numSims=500,maxNumClusters=10,
                                                                   outputFile=outputFile,distMethod=c("euclidean"),
                                                                   hclustAlgorithm=c("average"),
                                                                   consensusHclustAlgorithm=c("average"),
                                                                   minClustConsensus=.7, minMeanClustConsensus=.85,
                                                                   corUse="everything",pItem=.9,maxPAC=.15)


save(kmeansConsensuspam50_full_Nstart15pItem9,
     file=paste0(saveDir_PAM50,"/kmeansConsensuspam50_full_Nstart15pItem9_",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")


kmeansConsensuspam50_full_Nstart1pItem9 <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                                                  pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                                                  numSims=500,maxNumClusters=10,
                                                                  outputFile=outputFile,distMethod=c("euclidean"),
                                                                  hclustAlgorithm=c("average"),
                                                                  consensusHclustAlgorithm=c("average"),
                                                                  minClustConsensus=.7, minMeanClustConsensus=.85,
                                                                  corUse="everything",pItem=.9,maxPAC=.15)


save(kmeansConsensuspam50_full_Nstart1pItem9,
     file=paste0(saveDir_PAM50,"/kmeansConsensuspam50_full_Nstart1pItem9_",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")




hclustConsensuspam50_full_pItem9 <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("hc"),
                                                           pickKMethod=c("consensus"),iter.max=20,nstart=15,
                                                           numSims=500,maxNumClusters=10,
                                                           outputFile=outputFile,distMethod=c("euclidean"),
                                                           hclustAlgorithm=c("ward.D"),
                                                           consensusHclustAlgorithm=c("average"),
                                                           minClustConsensus=.7, minMeanClustConsensus=.85,
                                                           corUse="everything",pItem=.9,maxPAC=.15)


save(hclustConsensuspam50_full_pItem9 ,file=paste0(saveDir_PAM50,"/hclustConsensuspam50_full_pItem9_",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")



######gap test
#############
####pam50 full
load(paste0(saveDirGlobal,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
load("pam50_centroids_updatedSymbols.RData")
pam50GenesFull <- centroidMatrix[,1]

numFeatures <- "pam50Full"


#COME BACK: remove clusters with too few pam50 genes?

clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- pam50GenesFull
  
}

kmeansConsensuspam50_full_Nstart15pItem9 <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                                                   pickKMethod=c("consensus"),iter.max=20,nstart=15,
                                                                   numSims=500,maxNumClusters=10,
                                                                   outputFile=outputFile,distMethod=c("euclidean"),
                                                                   hclustAlgorithm=c("average"),
                                                                   consensusHclustAlgorithm=c("average"),
                                                                   minClustConsensus=.7, minMeanClustConsensus=.85,
                                                                   corUse="everything",pItem=.9,maxPAC=.15)


save(kmeansConsensuspam50_full_Nstart15pItem9,
     file=paste0(saveDir_PAM50,"/kmeansConsensuspam50_full_Nstart15pItem9_",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")


kmeansConsensuspam50_full_Nstart1pItem9 <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                                                  pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                                                  numSims=500,maxNumClusters=10,
                                                                  outputFile=outputFile,distMethod=c("euclidean"),
                                                                  hclustAlgorithm=c("average"),
                                                                  consensusHclustAlgorithm=c("average"),
                                                                  minClustConsensus=.7, minMeanClustConsensus=.85,
                                                                  corUse="everything",pItem=.9,maxPAC=.15)


save(kmeansConsensuspam50_full_Nstart1pItem9,
     file=paste0(saveDir_PAM50,"/kmeansConsensuspam50_full_Nstart1pItem9_",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")




hclustConsensuspam50_full_pItem9 <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("hc"),
                                                           pickKMethod=c("consensus"),iter.max=20,nstart=15,
                                                           numSims=500,maxNumClusters=10,
                                                           outputFile=outputFile,distMethod=c("euclidean"),
                                                           hclustAlgorithm=c("ward.D"),
                                                           consensusHclustAlgorithm=c("average"),
                                                           minClustConsensus=.7, minMeanClustConsensus=.85,
                                                           corUse="everything",pItem=.9,maxPAC=.15)


save(hclustConsensuspam50_full_pItem9 ,file=paste0(saveDir_PAM50,"/hclustConsensuspam50_full_pItem9_",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")


####gap test clustering
#############################

##with short, nstart=25
load(paste0(saveDirGlobal,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
load("pam50Short_genes.RData")


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- pam50Short
  
}

#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansGapTest_pam50_short_Nstart25 <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                                             pickKMethod=c("gap"),iter.max=20,nstart=25,
                                                             numSims=500,maxNumClusters=15,
                                                             outputFile=outputFile
)


save(kmeansGapTest_pam50_short_Nstart25,
     file=paste0(saveDir_PAM50,"/gapTestKmeans_pam50Short_nstart25",Sys.Date(),".RData.gzip"),compress="gzip")


###with short, nstart =1
load(paste0(saveDirGlobal,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
load("pam50Short_genes.RData")

saveDir_PAM50 <- "/home/kplaney/breast_analysis/"


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- pam50Short
  
}

#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansGapTest_pam50_short_Nstart1 <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                                            pickKMethod=c("gap"),iter.max=20,nstart=1,
                                                            numSims=500,maxNumClusters=15,
                                                            outputFile=outputFile
)


save(kmeansGapTest_pam50_short_Nstart1,
     file=paste0(saveDir_PAM50,"/gapTestKmeans_pam50Short_nstart1",Sys.Date(),".RData.gzip"),compress="gzip")

####hierarchical
load(paste0(saveDirGlobal,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
load("pam50Short_genes.RData")



clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- pam50Short
  
}

hclustGapTest_pam50_short <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("hc"),
                                                    pickKMethod=c("gap"),iter.max=20,nstart=25,distMethod=c("euclidean"),hclustAlgorithm=c("average"),
                                                    numSims=500,maxNumClusters=15,
                                                    outputFile=outputFile
)


save(hclustGapTest_pam50_short,
     file=paste0(saveDir_PAM50,"/hclustGapTest_pam50_short_",Sys.Date(),".RData.gzip"),compress="gzip")



####with full, nstart=25
load(paste0(saveDirGlobal,"/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
load("pam50_centroids_updatedSymbols.RData")
pam50GenesFull <- centroidMatrix[,1]

numFeatures <- "pam50Full"



clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- pam50GenesFull
  
}


kmeansGapTest_pam50_full_Nstart25 <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                                            pickKMethod=c("gap"),iter.max=20,nstart=25,
                                                            numSims=500,maxNumClusters=15,
                                                            outputFile=outputFile
)


save(kmeansGapTest_pam50_full_Nstart25,
     file=paste0(saveDir_PAM50,"/gapTestKmeans_pam50Full_nstart25",Sys.Date(),".RData.gzip"),compress="gzip")





####with full, nstart=1
load(paste0(saveDirGlobal,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
load("pam50_centroids_updatedSymbols.RData")
pam50GenesFull <- centroidMatrix[,1]
numFeatures <- "pam50Full"


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- pam50GenesFull
  
}

#we know these are strong clusters. have  minMeanClustConsensus around .85
kmeansGapTest_pam50_full_Nstart1 <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                                           pickKMethod=c("gap"),iter.max=20,nstart=1,
                                                           numSims=500,maxNumClusters=15,
                                                           outputFile=outputFile
)


save(kmeansGapTest_pam50_full_Nstart1,
     file=paste0(saveDir_PAM50,"/gapTestKmeans_pam50Full_nstart1",Sys.Date(),".RData.gzip"),compress="gzip")


####with hierachical

load(paste0(saveDirGlobal,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip"))
load("pam50_centroids_updatedSymbols.RData")
pam50GenesFull <- centroidMatrix[,1]


numFeatures <- "pam50Full"



clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- pam50GenesFull
  
}


hclustGapTest_pam50_full_Nstart25 <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("hc"),
                                                            pickKMethod=c("gap"),iter.max=20,nstart=25,distMethod=c("euclidean"),hclustAlgorithm=c("average"),
                                                            numSims=500,maxNumClusters=15,
                                                            outputFile=outputFile
)


save(hclustGapTest_pam50_full_Nstart25,
     file=paste0(saveDir_PAM50,"/gapTesthclust_pam50Full_",Sys.Date(),".RData.gzip"),compress="gzip")


#################################
###compare/look at results from all clustering methods
load(saveDirGlobal,"/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
names(dataMatrixList)
studyNames <- c("study_2034_GPL96_all","study_25055_GPL96_MDACC_M","study_22226_GPL1708_all",
                "study_20181_GPL96_all","study_19615_GPL570_all" ,"study_16446_GPL570_all",
                "study_12093_GPL96_all","study_25065_GPL96_MDACC")

datasetIndices <- na.omit(match(studyNames,names(dataMatrixList)))
load(paste0(saveDir_PAM50,"/gapTestKmeans_pam50Full_nstart12015-05-04.RData.gzip"))

unlist(kmeansGapTest_pam50_full_Nstart1$bestK)[datasetIndices]

load(paste0(saveDir_PAM50,"/gapTestKmeans_pam50Full_nstart252015-05-04.RData.gzip"))
unlist(kmeansGapTest_pam50_full_Nstart25$bestK)[datasetIndices]

load(paste0(saveDir_PAM50,"/gapTestKmeans_pam50Short_nstart12015-05-04.RData.gzip"))

unlist(kmeansGapTest_pam50_short_Nstart1$bestK)[datasetIndices]

load(paste0(saveDir_PAM50,"/gapTestKmeans_pam50Short_nstart252015-05-04.RData.gzip"))

unlist(kmeansGapTest_pam50_short_Nstart25$bestK)[datasetIndices]


##now hclust

#(there's no actual "nstart" for hclust; it's all the same.)

unlist(kmeansGapTest_pam50_full_Nstart1$bestK)[datasetIndices]

load(paste0(saveDir_PAM50,"/gapTesthclust_pam50Full_2015-05-15.RData.gzip"))
unlist(hclustGapTest_pam50_full_Nstart25$bestK)[datasetIndices]

load(paste0(saveDir_PAM50,"/hclustGapTest_pam50_short_2015-05-15.RData.gzip"))

unlist(hclustGapTest_pam50_short$bestK)[datasetIndices]




