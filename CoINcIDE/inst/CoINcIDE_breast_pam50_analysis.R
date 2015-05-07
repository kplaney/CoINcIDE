
####pam50 centroids clustering
load("/home/kplaney/breast_analysis/pam50Full_centroidCluster.RData.gzip")
#load data matrix list
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")

clustSampleIndexList <- pam50Full_centroidCluster$clustSampleIndexList
clustFeatureIndexList <- pam50Full_centroidCluster$clustFeatureIndexList
#35/50 is 70%
fractFeatIntersectThresh=.69
#we know these all share at least 35 genes.
numFeatIntersectThresh=34
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
edgeMethod <- "pearson"
#not used here
maxNullFractSize=.2

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")
pam50Full_centroidCluster_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                        edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                        sigMethod=sigMethod,maxNullFractSize=maxNullFractSize,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                        outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                        numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                        clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


save(pam50Full_centroidCluster_pearson_meanMatrix,file=paste0("/home/kplaney/breast_analysis/adjMatrices_pam50Full_centroidCluster_pearson_meanMatrix_",Sys.Date(),".RData.gzip"),compress="gzip")

#######
####now do k-means pam50 results:
#load data matrix list
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
#35/50 is 70%
fractFeatIntersectThresh=.69
#we know these all share at least 35 genes.
numFeatIntersectThresh=34
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=3
#not used here
maxNullFractSize=.2
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"


#full features, pearson
edgeMethod <- "pearson"
load("/home/kplaney/breast_analysis/kmeansConsensuspam50_full_Nstart1pItem9_pam50FullFeatures_2015-04-28.RData.gzip")

clustSampleIndexList <- kmeansConsensuspam50_full_Nstart1pItem9$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensuspam50_full_Nstart1pItem9$clustFeatureIndexList_PACR

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")

pam50Full_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                        edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                        sigMethod=sigMethod,maxNullFractSize=maxNullFractSize,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                        outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                        numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                        clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)

save(pam50Full_spearman_meanMatrix="/home/kplaney/breast_analysis/adjMatrices_pam50Full_pearson_meanMatrix_",Sys.Date(),".RData.gzip",compress="gzip")


#now spearman:
edgeMethod <- "spearman"
load("/home/kplaney/breast_analysis/hclustConsensuspam50_full_pItem9_pam50FullFeatures_2015-04-28.RData.gzip")

clustSampleIndexList <- kmeansConsensuspam50_full_Nstart1pItem9$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensuspam50_full_Nstart1pItem9$clustFeatureIndexList_PACR

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")

pam50Full_spearman_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                        edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                        sigMethod=sigMethod,maxNullFractSize=maxNullFractSize,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                        outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                        numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                        clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


save(pam50Full_spearman_meanMatrix="/home/kplaney/breast_analysis/adjMatrices_pam50Full_spearman_meanMatrix_",Sys.Date(),".RData.gzip",compress="gzip")

#now pam50 short:

#load data matrix list
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
#35/50 is 70%
fractFeatIntersectThresh=.69
#we know these all share at least 35 genes.
numFeatIntersectThresh=34
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


#full features, pearson
edgeMethod <- "pearson"
load("/home/kplaney/breast_analysis/kmeansConsensuspam50_short_Nstart1pItem9_pam50ShortFeatures_2015-04-28.RData.gzip")

clustSampleIndexList <- kmeansConsensuspam50_short_Nstart1pItem9$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensuspam50_short_Nstart1pItem9$clustFeatureIndexList_PACR

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")

pam50Short_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                         edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                         sigMethod=sigMethod,maxNullFractSize=maxNullFractSize,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                         outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                         numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                         clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


save(pam50Short_pearson_meanMatrix="/home/kplaney/breast_analysis/adjMatrices_pam50Short_pearson_meanMatrix_",Sys.Date(),".RData.gzip",compress="gzip")

######
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
#35/50 is 70%
fractFeatIntersectThresh=.69
#we know these all share at least 35 genes.
numFeatIntersectThresh=34
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


#full features, pearson
edgeMethod <- "spearman"
load("/home/kplaney/breast_analysis/kmeansConsensuspam50_short_Nstart1pItem9_pam50ShortFeatures_2015-04-28.RData.gzip")

clustSampleIndexList <- kmeansConsensuspam50_short_Nstart1pItem9$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensuspam50_short_Nstart1pItem9$clustFeatureIndexList_PACR

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")

pam50Short_spearman_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                         edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                         sigMethod=sigMethod,maxNullFractSize=maxNullFractSize,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                         outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                         numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                         clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)


save(pam50Short_spearman_meanMatrix="/home/kplaney/breast_analysis/adjMatrices_pam50Short_spearman_meanMatrix_",Sys.Date(),".RData.gzip",compress="gzip")


##########gap test versions

load("/home/kplaney/breast_analysis/gapTestKmeans_pam50Short_nstart252015-05-04.RData.gzip")
load("/home/kplaney/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
#35/50 is 70%
fractFeatIntersectThresh=.69
#we know these all share at least 35 genes.
numFeatIntersectThresh=34
clustSizeThresh=5
clustSizeFractThresh=.05
numParallelCores=3
#not used here
maxNullFractSize=.2
maxTrueSimilThresh=Inf
minTrueSimilThresh=.3
includeRefClustInNull=TRUE
numSims=500
checkNA=FALSE
outputFile="/home/kplaney/breast_analysis//CoINcIDE_messages.txt"
sigMethod <- "meanMatrix"


#full features, pearson
edgeMethod <- "pearson"
load("/home/kplaney/breast_analysis/kmeansConsensuspam50_full_Nstart1pItem9_pam50FullFeatures_2015-04-28.RData.gzip")

clustSampleIndexList <- kmeansConsensuspam50_full_Nstart1pItem9$clustSampleIndexList_PACR
clustFeatureIndexList <- kmeansConsensuspam50_full_Nstart1pItem9$clustFeatureIndexList_PACR

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")

pam50Full_pearson_meanMatrix <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                                                        edgeMethod=edgeMethod,numParallelCores=numParallelCores,minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                                        sigMethod=sigMethod,maxNullFractSize=maxNullFractSize,numSims=numSims,includeRefClustInNull=includeRefClustInNull,
                                                        outputFile=outputFile,fractFeatIntersectThresh=fractFeatIntersectThresh,
                                                        numFeatIntersectThresh=numFeatIntersectThresh,clustSizeThresh=clustSizeThresh, 
                                                        clustSizeFractThresh=clustSizeFractThresh,checkNA=FALSE)



