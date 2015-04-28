

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

saveDir <- "/home/kplaney/breast_analysis/"
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

load(paste0(saveDir,"/kmeansConsensuspam50_full_Nstart15pItem9_pam50ShortFeatures_04-28-2014.RData.gzip")
     
#all features should intersect here
fractFeatIntersectThresh=.8
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=8
minTrueSimilThresh=.3
#maxNullFractSize=.2
maxTrueSimilThresh=Inf
minTrueSimilThresh=.25
includeRefClustInNull=TRUE
numSims=500

sigMethod <- "meanMatrix"
edgeMethod <- "pearson"







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
