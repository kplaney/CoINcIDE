
#CHANGE: /home/kplaney/breast_analysis to /home/kplaney/breast_analysis_withTop20Genes/
fractFeatIntersectThresh=.8
#we know these all share at least 35 genes.
numFeatIntersectThresh=150
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=3
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/breast_analysis_withTop20Genes//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"

#200 features, pearson:
edgeMethod <- "pearson"
load("/home/kplaney/breast_analysis_withTop20Genes///curatedbreastData_kmeansConsensus_nstart1pItem9200Features_2015-05-04.RData.gzip")
load("/home/kplaney/breast_analysis_withTop20Genes//metaFeatures_200.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
#load data matrix list
load("/home/kplaney/breast_analysis_withTop20Genes//curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

clustSampleIndexList <-  kmeansConsensus$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensus$clustFeatureIndexList_PACR

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast200F_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                      edgeMethod=edgeMethod,maxNullFractSize=.2,numSims=500,centroidMethod=c("mean"))


#intersect: p-value, fract, meanMetric

save(breast200F_pearson_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_200F_pearson_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")

#####now do  spearman
edgeMethod <- "spearman"
load("/home/kplaney/breast_analysis//curatedbreastData_kmeansConsensus_nstart1pItem9200Features_2015-05-04.RData.gzip")
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
                                                       sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                       outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                       numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                       clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


save(breast200F_spearman_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_200F_spearman_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")

###now: pearson, but with centroid method
fractFeatIntersectThresh=.8
#we know these all share at least 35 genes.
numFeatIntersectThresh=150
clustSizeThresh=5
clustSizeFractThresh=.05
#more cores: centroid methods run slower!
numParallelCores=4
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "centroid"

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


breast200F_pearson_centroid <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                         edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                         sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                         outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                         numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                         clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)



save(breast200F_pearson_centroid ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_200F_pearson_centroid_",Sys.Date(),".RData.gzip"),compress="gzip")

########500:
fractFeatIntersectThresh=.8
#we know these all share at least 35 genes.
numFeatIntersectThresh=425
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=3
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
                                                      sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                      outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                      numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                      clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)



save(breast500F_pearson_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_500F_pearson_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")


##############
#now do  spearman
edgeMethod <- "spearman"
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


breast500F_spearman_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                       edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                       sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                       outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                       numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                       clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


save(breast500F_spearman_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_500F_spearman_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")

###now pearson, but with centroid:
fractFeatIntersectThresh=.8
#we know these all share at least 35 genes.
numFeatIntersectThresh=425
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=4
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "centroid"

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


breast500F_pearson_centroid <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                         edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                         sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                         outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                         numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                         clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


#had a bug in paste:
save(breast500F_pearson_centroid,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_500F_pearson_centroid_",Sys.Date(),".RData.gzip"),compress="gzip")


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
                                                       sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                       outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                       numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                       clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)



save(breast1000F_pearson_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_1000F_pearson_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")



###1000 spearman:

edgeMethod <- "spearman"
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


breast1000F_spearman_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                          edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                          sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                          outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                          numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                          clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)



save(breast1000F_spearman_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_1000F_spearman_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")



###now pearson, but with centroid:
fractFeatIntersectThresh=.8
#we know these all share at least 35 genes.
numFeatIntersectThresh=800
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=4
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "centroid"

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


breast1000F_pearson_centroid <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                          edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                          sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                          outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                          numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                          clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)



save(breast1000F_pearson_centroid ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_1000F_pearson_centroid_",Sys.Date(),".RData.gzip"),compress="gzip")


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
load("/home/kplaney/breast_analysis//curatedbreastData_kmeansConsensus_nstart1pItem92000Features_2015-05-05.RData.gzip")
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
                                                       sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                       outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                       numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                       clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


#had a bug in paste:
save(breast2000F_pearson_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_2000F_pearson_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")


###spearman 2000

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
edgeMethod <- "spearman"
load("/home/kplaney/breast_analysis//curatedbreastData_kmeansConsensus_nstart1pItem92000Features_2015-05-05.RData.gzip")
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


breast2000F_spearman_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                          edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                          sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                          outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                          numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                          clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)



save(breast2000F_spearman_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_2000F_spearman_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")

###now with pearson, but centroid
fractFeatIntersectThresh=.85
numFeatIntersectThresh=1700
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=4
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "centroid"

# 2000features, pearson:
edgeMethod <- "pearson"
load("/home/kplaney/breast_analysis//curatedbreastData_kmeansConsensus_nstart1pItem92000Features_2015-05-05.RData.gzip")
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


breast2000F_pearson_centroid <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                          edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                          sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                          outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                          numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                          clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)



save(breast2000F_pearson_centroid,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_2000F_pearson_meanMatrix_centroid_",Sys.Date(),".RData.gzip"),compress="gzip")

#########
####gap test results as opposed to k-means consensus.
##gap test with k-means

###200 features
fractFeatIntersectThresh=.8
#we know these all share at least 35 genes.
numFeatIntersectThresh=150
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=3
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"

#200 features, pearson:
edgeMethod <- "pearson"
load("/home/kplaney/breast_analysis/curatedbreastData_kmeansGap_nstart25_200_features_2015-05-18.RData.gzip")
load("/home/kplaney/breast_analysis/metaFeatures_200.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
#load data matrix list
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

clustSampleIndexList <-  kmeansGap$clustSampleIndexList
clustFeatureIndexList <- kmeansGap$clustFeatureIndexList

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast200F_kmeansGap_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                                   edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                                   sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                                   outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                                   numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                                   clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)



save(breast200F_kmeansGap_pearson_meanMatrix,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_breast200F_kmeansGap_pearson_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")


###500 features
fractFeatIntersectThresh=.8
#we know these all share at least 35 genes.
numFeatIntersectThresh=425
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=3
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"

# 500features, pearson:
edgeMethod <- "pearson"
load("/home/kplaney/breast_analysis/curatedbreastData_kmeansGap_nstart25_500_features_2015-05-18.RData.gzip")
load("/home/kplaney/breast_analysis/metaFeatures_500.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

clustSampleIndexList <-  kmeansGap$clustSampleIndexList
clustFeatureIndexList <- kmeansGap$clustFeatureIndexList

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast500F_kmeansGap_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                                   edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                                   sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                                   outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                                   numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                                   clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)



save(breast500F_kmeansGap_pearson_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_breast500F_kmeansGap_pearson_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")


#TO DO: 1000,2000 have not finished running yet.
#1000

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


edgeMethod <- "pearson"
#TO DO: this hasn't finished running yet:
load("/home/kplaney/breast_analysis//curatedbreastData_kmeansGap_nstart25_1000_features_2015-05-18.RData.gzip")
load("/home/kplaney/breast_analysis/metaFeatures_1000.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

clustSampleIndexList <-  kmeansGap$clustSampleIndexList
clustFeatureIndexList <- kmeansGap$clustFeatureIndexList

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast1000F_kmeansGap_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                                    edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                                    sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                                    outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                                    numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                                    clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)



save(breast1000F_kmeansGap_pearson_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_breast1000F_kmeansGap_pearson_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")


###2000
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
load("/home/kplaney/breast_analysis/curatedbreastData_kmeansGap_nstart25_2000_features_")
load("/home/kplaney/breast_analysis/metaFeatures_2000.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
#load("/home/kplaney/breast_analysis/esets_proc_TCGAcombat.RData.gzip")

load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}
clustSampleIndexList <-  hclustOut$clustSampleIndexList
clustFeatureIndexList <- hclustOut$clustFeatureIndexList

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast2000F_kmeansGap_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                                    edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                                    sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                                    outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                                    numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                                    clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)



save(breast2000F_kmeansGap_pearson_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_2000F_kmeansGap_pearson_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")


##########
##gap test using hierarchical clustering:
###200 features
fractFeatIntersectThresh=.8
#we know these all share at least 35 genes.
numFeatIntersectThresh=150
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=3
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"

#200 features, pearson:
edgeMethod <- "pearson"
load("/home/kplaney/breast_analysis/curatedbreastData_hclust_200Features_2015-05-15.RData.gzip")
load("/home/kplaney/breast_analysis/metaFeatures_200.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
#load data matrix list
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){

dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]

}

clustSampleIndexList <-  hclustOut$clustSampleIndexList
clustFeatureIndexList <- hclustOut$clustFeatureIndexList

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast200F_hclustGap_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)



save(breast200F_hclustGap_pearson_meanMatrix,file=
paste0("/home/kplaney/breast_analysis/adjMatrices_breast200F_hclustGap_pearson_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")


###500 features
fractFeatIntersectThresh=.8
#we know these all share at least 35 genes.
numFeatIntersectThresh=425
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=3
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"

# 500features, pearson:
edgeMethod <- "pearson"
load("/home/kplaney/breast_analysis/curatedbreastData_hclust_500Features_2015-05-15.RData.gzip")
load("/home/kplaney/breast_analysis/metaFeatures_500.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

clustSampleIndexList <-  hclustOut$clustSampleIndexList
clustFeatureIndexList <- hclustOut$clustFeatureIndexList

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast500F_hclustGap_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                         edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                         sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                         outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                         numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                         clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)



save(breast500F_hclustGap_pearson_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_breast500F_hclustGap_pearson_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")


#
#10000

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


edgeMethod <- "pearson"
load("/home/kplaney/breast_analysis/curatedbreastData_hclust_1000Features_2015-05-15.RData.gzip")
load("/home/kplaney/breast_analysis/metaFeatures_1000.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

clustSampleIndexList <-  hclustOut$clustSampleIndexList
clustFeatureIndexList <- hclustOut$clustFeatureIndexList

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast1000F_hclustGap_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                          edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                          sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                          outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                          numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                          clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)



save(breast1000F_hclustGap_pearson_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_breast1000F_hclustGap_pearson_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")


###2000
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
load("/home/kplaney/breast_analysis/curatedbreastData_hclust_2000Features_2015-05-16.RData.gzip")
load("/home/kplaney/breast_analysis/metaFeatures_2000.RData.gzip")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
#load("/home/kplaney/breast_analysis/esets_proc_TCGAcombat.RData.gzip")

load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}
clustSampleIndexList <-  hclustOut$clustSampleIndexList
clustFeatureIndexList <- hclustOut$clustFeatureIndexList

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")


breast2000F_hclustGap_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                          edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                          sigMethod=sigMethod,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                          outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                          numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                          clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)



save(breast2000F_hclustGap_pearson_meanMatrix ,file=
       paste0("/home/kplaney/breast_analysis/adjMatrices_2000F_hclustGap_pearson_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")



