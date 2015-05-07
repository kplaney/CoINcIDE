

fractFeatIntersectThresh=.8
#we know these all share at least 35 genes.
numFeatIntersectThresh=150
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=3
#this is "softer" than breast cancer.
maxNullFractSize=.2
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"

#200 features, pearson:
edgeMethod <- "pearson"
load("/home/kplaney/breast_analysis//curatedbreastData_kmeansConsensus_nstart1pItem9200Features_2015-05-04.RData.gzip")
load("/home/kplaney/breast_analysis/metaFeatures_200.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
#load data matrix list
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

clustSampleIndexList <-  kmeansConsensus$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensus$clustFeatureIndexList_PACR

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast200F_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                      edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                      sigMethod=sigMethod,maxNullFractSize=maxNullFractSize,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                      outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                      numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                      clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


#had a bug in paste:
save(breast200F_pearson_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_200F_pearson_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")

#now do  spearman
edgeMethod <- "spearman"
load("/home/kplaney/breast_analysis//curatedbreastData_kmeansConsensus_nstart1_200Features_2015-04-28.RData.gzip")
load("/home/kplaney/breast_analysis/metaFeatures_200.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")


load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

clustSampleIndexList <-  kmeansConsensus$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensus$clustFeatureIndexList_PACR



source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast200F_spearman_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                       edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                       sigMethod=sigMethod,maxNullFractSize=maxNullFractSize,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                       outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                       numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                       clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


save(breast200F_spearman_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_200F_spearman_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")


########500:
fractFeatIntersectThresh=.8
#we know these all share at least 35 genes.
numFeatIntersectThresh=425
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=3
#this is "softer" than breast cancer.
maxNullFractSize=.2
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"

# 500features, pearson:
edgeMethod <- "pearson"
load("/home/kplaney/breast_analysis//curatedbreastData_kmeansConsensus_nstart1pItem9500Features_2015-05-04.RData.gzip")
load("/home/kplaney/breast_analysis/metaFeatures_500.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

clustSampleIndexList <-  kmeansConsensus$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensus$clustFeatureIndexList_PACR

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast500F_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                      edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                      sigMethod=sigMethod,maxNullFractSize=maxNullFractSize,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                      outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                      numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                      clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


#had a bug in paste:
save(breast500F_pearson_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_500F_pearson_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")


##############
#now do  spearman
edgeMethod <- "spearman"
load("/home/kplaney/breast_analysis//curatedbreastData_kmeansConsensus_nstart1_500Features_2015-04-29.RData.gzip")
load("/home/kplaney/breast_analysis/metaFeatures_500.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")

load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}
clustSampleIndexList <-  kmeansConsensus$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensus$clustFeatureIndexList_PACR



source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast500F_spearman_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                       edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                       sigMethod=sigMethod,maxNullFractSize=maxNullFractSize,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                       outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                       numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                       clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


save(breast500F_spearman_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_500F_spearman_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")

########1000
fractFeatIntersectThresh=.8
#we know these all share at least 35 genes.
numFeatIntersectThresh=800
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=3
#this is "softer" than breast cancer.
maxNullFractSize=.2
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"

# 500features, pearson:
edgeMethod <- "pearson"
load("/home/kplaney/breast_analysis//curatedbreastData_kmeansConsensus_nstart1pItem91000Features_2015-05-05.RData.gzip")
load("/home/kplaney/breast_analysis/metaFeatures_1000.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

clustSampleIndexList <-  kmeansConsensus$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensus$clustFeatureIndexList_PACR

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast1000F_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                       edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                       sigMethod=sigMethod,maxNullFractSize=maxNullFractSize,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                       outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                       numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                       clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)



save(breast1000F_pearson_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_1000F_pearson_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")


####2000:
fractFeatIntersectThresh=.85
#we know these all share at least 35 genes.
numFeatIntersectThresh=1700
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=3
#this is "softer" than breast cancer.
maxNullFractSize=.2
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"

# 2000features, pearson:
edgeMethod <- "pearson"
load("/home/kplaney/breast_analysis/curatedbreastData_kmeansConsensus_nstart1pItem92000Features_2015-05-05.RData.gzip")
load("/home/kplaney/breast_analysis/metaFeatures_2000.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
#load("/home/kplaney/breast_analysis/esets_proc_TCGAcombat.RData.gzip")

load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}
clustSampleIndexList <-  kmeansConsensus$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensus$clustFeatureIndexList_PACR

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast2000F_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                       edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                       sigMethod=sigMethod,maxNullFractSize=maxNullFractSize,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                       outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                       numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                       clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


#had a bug in paste:
save(breast2000F_pearson_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_2000F_pearson_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")

