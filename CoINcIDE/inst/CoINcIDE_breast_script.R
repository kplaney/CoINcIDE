library("curatedBreastData")
#had to download source: wget http://www.bioconductor.org/packages/release/data/experiment/src/contrib/curatedBreastData_1.0.0.tar.gz

load("/home/kplaney/R/x86_64-redhat-linux-gnu-library/curatedBreastData/data/curatedBreastDataExprSetList.rda")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")

#NOTE: function may conflict with other package
#this takes a while for this large of a database! (few hours)
esets <- processExpressionSetList(exprSetList=curatedBreastDataExprSetList,outputFileDirectory="/home/kplaney/breast_analysis/",
                                  minVar=.001,featureDataFieldName="gene_symbol",uniquePDataID="patient_ID")


save(esets,file="/home/kplaney/breast_analysis/curatedBreastData_esets_proc.RData.gzip",compress="gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_featureSelection.R")

#now format just as a list of data matrices.
dataMatrixList <- exprSetListToMatrixList(esets,featureDataFieldName="gene_symbol")

names(dataMatrixList) <- names(esets)


##ALSO: merge matrices first to help decide which studies to keep overall

#also merge this one
output <- merge_datasetList(datasetList=dataMatrixList,minNumGenes = 10000, minNumPatients = 40,batchNormalize = c('none'));
save(output,file="/home/kplaney/breast_analysis//mergedExprMatrix_minVar001_17_studies_no_norm.RData.gzip",compress="gzip")

#have NOT run these last two yet:
output <- merge_datasetList(datasetList=dataMatrixList,minNumGenes = 10000, minNumPatients = 40,batchNormalize = c('BMC'));
save(output,file="/home/kplaney/breast_analysis//mergedExprMatrix_minVar001_17_studies_BMC_norm.RData.gzip",compress="gzip")

output <- merge_datasetList(datasetList=dataMatrixList,minNumGenes = 10000, minNumPatients = 40,batchNormalize = c('combat'));
save(output,file="/home/kplaney/breast_analysis//mergedExprMatrix_minVar001_17_studies_combat_norm.RData.gzip",compress="gzip")

#NOW: also remove these smaller datasets from the esets list before save

if(length(output$removeDatasetIndices>0)){
  
  dataMatrixList <- dataMatrixList[-output$removeDatasetIndices]
  
}
save(dataMatrixList,file="/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip",compress="gzip")




load("/home/kplaney/breast_analysis/curatedBreastData_esets_proc.RData.gzip")



###save this dataMatrixList now.

#ran meta-feature analysis for 1000,500,1000,2000 features.
metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=200,
                                                   numTopFeatFromEachDataset=20,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile="/home/kplaney/breast_analysis//selectFeaturesMetaVarianceOut.txt")

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

#(saved each one after ran)
save(metaFeatures,file="/home/kplaney/breast_analysis/metaFeatures_200.RData.gzip",compress="gzip")
metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=500,
                                                   numTopFeatFromEachDataset=20,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile="/home/kplaney/breast_analysis//selectFeaturesMetaVarianceOut.txt")

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

save(metaFeatures,file="/home/kplaney/breast_analysis/metaFeatures_500.RData.gzip",compress="gzip")

metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=1000,
                                                   numTopFeatFromEachDataset=20,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile="/home/kplaney/breast_analysis//selectFeaturesMetaVarianceOut.txt")

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

save(metaFeatures,file="/home/kplaney/breast_analysis/metaFeatures_1000.RData.gzip",compress="gzip")

metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=2000,
                                                   numTopFeatFromEachDataset=20,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile="/home/kplaney/breast_analysis//selectFeaturesMetaVarianceOut.txt")

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

save(metaFeatures,file="/home/kplaney/breast_analysis/metaFeatures_2000.RData.gzip",compress="gzip")


###########clustering
#just change the numFeatures for each run:
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")
numFeatures <- 200
load(paste0("/home/kplaney/breast_analysis/metaFeatures_",numFeatures,".RData.gzip")
load("/home/kplaney/breast_analysis/curatedBreastData_esets_proc.RData.gzip")
#now format just as a list of data matrices.
dataMatrixList <- exprSetListToMatrixList(esets,featureDataFieldName="gene")

names(dataMatrixList) <- names(esets)
##do for each 200,500,1000,2000 (load different metaFeatures RData object each time.)
#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- metaFeatures$finalFeatures
  
}

kmeansConsensus <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList=clustFeaturesList,clustMethod=c("km"),
                                          pickKMethod=c("consensus"),
                                          numSims=500,maxNumClusters=15,
                                          outputFile="/home/kplaney/breast_analysis/test.txt",iter.max=30,nstart=25,distMethod=c("euclidean"),
                                          hclustAlgorithm=c("average"),
                                          consensusHclustAlgorithm=c("average"),
                                          minClustConsensus=.7, minMeanClustConsensus=.7,corUse="everything",pItem=.9)


save(kmeansConsensus,file=paste0("/home/kplaney/breast_analysis/curatedBreastData_kmeansConsensus_",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")



######pam50 clustering: 
##"dummy" clustering
load("/home/data/breast_microarrayDB//bc_expr_baseMinVar01_studiesWithAtLeast10kGenes40Patients.RData.gzip");

fullMatrixList <- list();
for(f in 1:length(bc_expr_baseMinVar01)){
  
  fullMatrixList[[f]] <- t(bc_expr_baseMinVar01[[f]]$expr);
  
}
names(fullMatrixList) <- names(bc_expr_baseMinVar01);
minNumGenes <- 25;
centroidData <- "/home/data/breast_microarrayDB/pam50_centroids.RData";
source("/home/kplaney/gitRepos/IGP_network/igp_network/clust_robust.R");
######cluster samples.
clustMatrixList_pam50 <- lapply(fullMatrixList,FUN=function(origDataMatrix,minNumGenes,centroidData){
  
  
  centroidAssignments <- assignCentroidSubtype(origDataMatrix,minNumGenes=minNumGenes,centroidRData=centroidData);
  cat("\n",unique(centroidAssignments$subtypes[,2]),"\n")
  clustMatrixList <- list();
  
  if(all(!is.na((centroidAssignments)))){
    
    for(c in 1:length(unique(centroidAssignments$subtypes[,2]))){
      
      
      clustMatrixList[[c]] <- data.matrix(origDataMatrix[which(centroidAssignments$subtypes[,2]==unique(centroidAssignments$subtypes[,2])[c]),centroidAssignments$centroidMatrix_genes,drop=FALSE]);
      
      if(all(is.null(dim(clustMatrixList[[c]])))){
        
        stop("\nNot getting a 2D clust matrix.")
      }
    }
    
    
    names(clustMatrixList) <- unique(centroidAssignments$subtypes[,2]);
    
  }
  
  return(clustMatrixList);
  
} ,minNumGenes=minNumGenes,centroidData=centroidData);

names(clustMatrixList_pam50) <- names(fullMatrixList);

output <- list(clustMatrixList_pam50=clustMatrixList_pam50,dataMatrixList=fullMatrixList);
save(output,file="/home/data/breast_microarrayDB/output/pam50_subtypes/pam50_clustMatrixList.RData.gzip");


##########
#########druggable 
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")
numFeatures <- "druggable"
#NOTE: you'll need to save this data object in the R package
load("/home/data/genomeReferences/annotationHuman/druggableGenomeList_updatedGeneSymbols.RData.gzip");
load("/home/kplaney/breast_analysis/curatedBreastData_esets_proc.RData.gzip")

#now format just as a list of data matrices.
dataMatrixList <- exprSetListToMatrixList(esets,featureDataFieldName="gene")

names(dataMatrixList) <- names(esets)

#use only the datasets that have enough features (looked at datasets ahead of time.)
dataMatrixList <- dataMatrixList[c(1,2,c(4:16))]
clustFeaturesList <- list()

#should be ~93 intersecting genes across all these datasets.
#looked like filtered out lowly varying genes - like 50th percentile?
for(d in 1:length(dataMatrixList)){
  
  druggableGenes <- intersect(druggableGenes,rownames(dataMatrixList[[d]]))
  
}
#check: genes should be in here...perhaps filter genes further with filter by variancea again. so that this code can run on its own.
#load("/home/data/breast_microarrayDB/druggableGeneShort_min50th_percVar_noStudy3_17_18.RData");

for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- druggableGenes
  
}


kmeansConsensus <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList=clustFeaturesList,clustMethod=c("km"),
                                          pickKMethod=c("consensus"),
                                          numSims=500,maxNumClusters=15,
                                          outputFile="/home/kplaney/breast_analysis/test.txt",iter.max=30,nstart=25,distMethod=c("euclidean"),
                                          hclustAlgorithm=c("average"),
                                          consensusHclustAlgorithm=c("average"),
                                          minClustConsensus=.7, minMeanClustConsensus=.7,corUse="everything",pItem=.9)


save(kmeansConsensus,file=paste0("/home/kplaney/breast_analysis/curatedBreastData_kmeansConsensus_",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")



###compare with pam50 short on merged matrix
#first, not normalized:
load("/home/data/breast_microarrayDB//mergedExprMatrix_minVar01_18_studies_no_norm.RData.gzip");

Pam50_subtypes <- assignCentroidSubtype(t(output$mergedExprMatrix),minNumGenes=30,centroidRData="/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData");
table(Pam50_subtypes$subtypes[,2]);
mergeData <- output;
patientData <- data.frame(mergeData$study,factor(as.numeric(factor(mergeData$study))),mergeData$GSMID);
colnames(patientData) <- c("fullStudyName","studyID","GSMID");
load("/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData");
#not a ton of genes that overlap! only around 190...
pam50Genes <- centroidMatrix[ ,1];

dataList <- list(output$mergedExprMatrix);
#which of the pam50 genes are in ALL studies?
pam50GeneShort <- pam50Genes[na.omit(match(rownames(output$mergedExprMatrix),pam50Genes))];
cat("\n",proc.time(),"\n");

#genes in rows?
dataMatrixList <- list(output$mergedExprMatrix)
clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- pam50GeneShort
  
}

kmeansConsensus <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList=clustFeaturesList,clustMethod=c("km"),
                                          pickKMethod=c("consensus"),
                                          numSims=500,maxNumClusters=15,
                                          outputFile="/home/kplaney/breast_analysis/test.txt",iter.max=30,nstart=25,distMethod=c("euclidean"),
                                          hclustAlgorithm=c("average"),
                                          consensusHclustAlgorithm=c("average"),
                                          minClustConsensus=.7, minMeanClustConsensus=.7,corUse="everything",pItem=.9)

##then, normalized:
#also try with combat?
load("/home/data/breast_microarrayDB//mergedExprMatrix_minVar01_18_studies_BMC_norm.RData.gzip");

mergeData <- output;
patientData <- data.frame(mergeData$study,factor(as.numeric(factor(mergeData$study))),mergeData$GSMID);
colnames(patientData) <- c("fullStudyName","studyID","GSMID");
load("/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData");
#not a ton of genes that overlap! only around 190...
pam50Genes <- centroidMatrix[ ,1];

#which of the pam50 genes are in ALL studies?
pam50GeneShort <- pam50Genes[na.omit(match(rownames(output$mergedExprMatrix),pam50Genes))];
cat("\n",proc.time(),"\n");

dataMatrixList <- list(output$mergedExprMatrix)
clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- pam50GeneShort
  
}


kmeansConsensus <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList=clustFeaturesList,clustMethod=c("km"),
                                          pickKMethod=c("consensus"),
                                          numSims=500,maxNumClusters=15,
                                          outputFile="/home/kplaney/breast_analysis/test.txt",iter.max=30,nstart=25,distMethod=c("euclidean"),
                                          hclustAlgorithm=c("average"),
                                          consensusHclustAlgorithm=c("average"),
                                          minClustConsensus=.7, minMeanClustConsensus=.7,corUse="everything",pItem=.9)

#try combat too
load("/home/data/breast_microarrayDB/mergedExprMatrix_minVar01_18_studies_combat_norm.RData.gzip")

mergeData <- output;
patientData <- data.frame(mergeData$study,factor(as.numeric(factor(mergeData$study))),mergeData$GSMID);
colnames(patientData) <- c("fullStudyName","studyID","GSMID");
load("/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData");
#not a ton of genes that overlap! only around 190...
pam50Genes <- centroidMatrix[ ,1];

#which of the pam50 genes are in ALL studies?
pam50GeneShort <- pam50Genes[na.omit(match(rownames(output$mergedExprMatrix),pam50Genes))];
cat("\n",proc.time(),"\n");

dataMatrixList <- list(output$mergedExprMatrix)
clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- pam50GeneShort
  
}


kmeansConsensus <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList=clustFeaturesList,clustMethod=c("km"),
                                          pickKMethod=c("consensus"),iter.max=15,nstart=15,
                                          numSims=500,maxNumClusters=15,
                                          outputFile="/home/kplaney/breast_analysis/test.txt",iter.max=30,nstart=25,distMethod=c("euclidean"),
                                          hclustAlgorithm=c("average"),
                                          consensusHclustAlgorithm=c("average"),
                                          minClustConsensus=.7, minMeanClustConsensus=.7,corUse="everything",pItem=.9)




save(kmeansConsensus,file=paste0("/home/kplaney/breast_analysis/curatedBreastData_kmeansConsensus_",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")

#############################
####CoINcIDE (not merged) with pam50 short for comparison


load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
load("/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData")
pam50GenesFull <- centroidMatrix[,1]

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")
numFeatures <- "pam50Full"



#COME BACK: remove clusters with too few pam50 genes?

clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- pam50GenesFull
  
}

#we know these are strong clusters. have  minMeanClustConsensus=.8
kmeansConsensus <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                          pickKMethod=c("consensus"),iter.max=20,nstart=15,
                                          numSims=500,maxNumClusters=10,
                                          outputFile="/home/kplaney/breast_analysis/test.txt",distMethod=c("euclidean"),
                                          hclustAlgorithm=c("average"),
                                          consensusHclustAlgorithm=c("average"),
                                          minClustConsensus=.7, minMeanClustConsensus=.85,corUse="everything",pItem=.9,maxPAC=.1)


save(kmeansConsensus,file=paste0("/home/kplaney/breast_analysis/curatedBreastData_kmeansConsensus_",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")


######
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
load("/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData")
pam50GenesFull <- centroidMatrix[,1]

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")
numFeatures <- "pam50Full"



#COME BACK: remove clusters with too few pam50 genes?

clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- pam50GenesFull
  
}

#we know these are strong clusters. have  minMeanClustConsensus=.8
hclustConsensus <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("hc"),
                                          pickKMethod=c("consensus"),iter.max=20,nstart=15,
                                          numSims=500,maxNumClusters=10,
                                          outputFile="/home/kplaney/breast_analysis/test.txt",distMethod=c("euclidean"),
                                          hclustAlgorithm=c("ward.D"),
                                          consensusHclustAlgorithm=c("average"),
                                          minClustConsensus=.7, minMeanClustConsensus=.85,corUse="everything",pItem=.9,maxPAC=.1)


save(hclustConsensus,file=paste0("/home/kplaney/breast_analysis/curatedBreastData_hclustConsensus_",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")


##########
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
#load("/home/kplaney/breast_analysis/mergedExprMatrix_minVar001_17_studies_no_norm.RData.gzip")
#load("/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData")
#pam50FullGenes <- centroidMatrix[,1]
#pam50Short <- pam50FullGenes[na.omit(match(rownames(output$mergedExprMatrix),pam50FullGenes))]
#save(pam50Short,file="/home/kplaney/breast_analysis/pam50Short_genes.RData")
#load("/home/data/breast_microarrayDB/output/ISMB/kmeans_merged_pam50Short.RData.gzip");

#pam50Short <- output$pam50GeneShort;
load("/home/kplaney/breast_analysis/pam50Short_genes.RData")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")
numFeatures <- "pam50Short_nstart15_pItem8"



#COME BACK: remove clusters with too few pam50 genes?

clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- pam50Short
  
}

#we know these are strong clusters. have  minMeanClustConsensus=.8
kmeansConsensus <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("km"),
                                          pickKMethod=c("consensus"),iter.max=20,nstart=1,
                                          numSims=500,maxNumClusters=10,
                                          outputFile="/home/kplaney/breast_analysis/test.txt",distMethod=c("euclidean"),
                                          hclustAlgorithm=c("average"),
                                          consensusHclustAlgorithm=c("average"),
                                          minClustConsensus=.7, minMeanClustConsensus=.85,corUse="everything",pItem=.8,maxPAC=.1)


save(kmeansConsensus,file=paste0("/home/kplaney/breast_analysis/curatedBreastData_kmeansConsensus_",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")


####
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
load("/home/kplaney/breast_analysis/pam50Short_genes.RData")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")
numFeatures <- "pam50Short_hclustpItem8"



#COME BACK: remove clusters with too few pam50 genes?

clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- pam50Short
  
}

#we know these are strong clusters. have  minMeanClustConsensus=.8
hclustConsensus <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("hc"),
                                          pickKMethod=c("consensus"),iter.max=20,nstart=15,
                                          numSims=500,maxNumClusters=10,
                                          outputFile="/home/kplaney/breast_analysis/test.txt",distMethod=c("euclidean"),
                                          hclustAlgorithm=c("ward.D"),
                                          consensusHclustAlgorithm=c("average"),
                                          minClustConsensus=.7, minMeanClustConsensus=.85,corUse="everything",pItem=.8,maxPAC=.1)


save(hclustConsensus,file=paste0("/home/kplaney/breast_analysis/curatedBreastData_hclustConsensus_",numFeatures,"Features_",Sys.Date(),".RData.gzip"),compress="gzip")



#original run:
load("/home/data/breast_microarrayDB/output/pam50_subtypes/pam50Short_kmeans_allStudies.RData.gzip")


hclustConsensus <- clustMatrixListWrapper(dataMatrixList,clustFeaturesList,clustMethod=c("hc"),
                                          pickKMethod=c("consensus"),
                                          numSims=100,maxNumClusters=15,
                                          outputFile="/home/kplaney/breast_analysis/hc_test.txt",iter.max=30,nstart=25,distMethod=c("euclidean"),
                                          hclustAlgorithm=c("ward.D"),
                                          consensusHclustAlgorithm=c("average"),
                                          minClustConsensus=.7, minMeanClustConsensus=.8,corUse="everything",pItem=.9)



hclustConsensus <- clustMatrixListWrapper(dataMatrix,clustFeaturesList,clustMethod=c("hc"),
                                          pickKMethod=c("consensus"),
                                          numSims=100,maxNumClusters=15,
                                          outputFile="/home/kplaney/breast_analysis/hc_test.txt",iter.max=30,nstart=25,distMethod=c("euclidean"),
                                          hclustAlgorithm=c("ward.D"),
                                          consensusHclustAlgorithm=c("average"),
                                          minClustConsensus=.7, minMeanClustConsensus=.8,corUse="everything",pItem=.9)

hclustConsensus <- clustMatrixListWrapper(dataMatrix,clustFeaturesList,clustMethod=c("kmeans"),
                                          pickKMethod=c("gap"),
                                          numSims=50,maxNumClusters=15,
                                          outputFile="/home/kplaney/breast_analysis/hc_test.txt",iter.max=30,nstart=25,distMethod=c("euclidean"),
                                          hclustAlgorithm=c("ward.D"),
                                          consensusHclustAlgorithm=c("average"),
                                          minClustConsensus=.7, minMeanClustConsensus=.8,corUse="everything",pItem=.9)



clustSamplesIndexList <- list()
dataMatrixList <- output$dataMatrixList
for(d in 1:length(dataMatrixList)){
  
  for(c in 1:length(output$clustMatrixListpam50Short[[d]])){
    
    clustSamplesIndexList[[d]] <- na.omit(match(colnames(output$clustMatrixListpam50Short[[d]][[c]])),rownames(dataMatrixList[[d]]))
    
  }
  
}

kmeansConsensus <- clustMatrixListWrapper(dataMatrix,clustFeaturesList,clustMethod=c("km"),
                                          pickKMethod=c("consensus"),
                                          numSims=500,maxNumClusters=15,
                                          outputFile="/home/kplaney/breast_analysis/test.txt",iter.max=30,nstart=25,distMethod=c("euclidean"),
                                          hclustAlgorithm=c("average"),
                                          consensusHclustAlgorithm=c("average"),
                                          minClustConsensus=.7, minMeanClustConsensus=.8,corUse="everything",pItem=.9)




