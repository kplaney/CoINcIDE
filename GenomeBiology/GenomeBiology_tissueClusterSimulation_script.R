
source("/home/ywrfc09/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_simulation.R")
source("/home/ywrfc09/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")
source("/home/ywrfc09/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_communityDetection.R")

##test out FDR.
highQualityDataset <- createTissueSimDatasets(numSimDatasets=2,
                                   eigenValueMin = -.001,simType=c("highQualityClust"),
                                   numPerClust = c(50,50,50,50),
                                   stddevNoise=0,numRows=200)

dataMatrixList=highQualityDataset$dataMatrixList
clustSampleIndexList=highQualityDataset$clustSampleIndexList
clustFeatureIndexList=highQualityDataset$clustFeatureIndexList

real_adjMatrices <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=highQualityDataset$clustSampleIndexList,clustFeatureIndexList=highQualityDataset$clustFeatureIndexList,
                                    edgeMethod=c("pearson"),
                                    numSims=500,
                                    outputFile="~/CoINcIDE_messages.txt",
                                    checkNA=FALSE,centroidMethod=c("mean"))

real_adjMatrices <- CoINcIDE_getAdjMatrices(dataMatrixList=dataMatrixList,clustSampleIndexList=highQualityDataset$clustSampleIndexList,clustFeatureIndexList=highQualityDataset$clustFeatureIndexList,
                                            edgeMethod=c("spearman"),
                                            numSims=500,
                                            outputFile="~/CoINcIDE_messages.txt",
                                            checkNA=FALSE,centroidMethod=c("mean"))



highQuality <- runTissueClusterSimROC(saveDir = "/home/kplaney/lungSims/",numSimDatasets=10,
                          eigenValueMin = -.001,simType=c("highQualityClust"),
                          noiseVector = seq(from=0,to=2.5,by=.2),numPerClust = c(50,50,50,50),
                          numWrapperSims=100,numSims=500,fractFeatIntersectThresh=.7,numFeatIntersectThresh=190,
                          clustSizeThresh=5,clustSizeFractThresh=.05,numParallelCores=3,minTrueSimilThresh=.3,
                          maxTrueSimilThresh=Inf,includeRefClustInNull=TRUE, edgeMethod=c("pearson"),
                          indEdgePvalueThresh=.1,meanEdgePairPvalueThresh=.05,restrictEdges=FALSE
)
  
mixQuality <- runTissueClusterSimROC(saveDir = "/home/kplaney/lungSims/",numSimDatasets=10,
                             eigenValueMin = -.001,simType=c("mixedClustQualityClust"),
                             noiseVector = seq(from=0,to=2.5,by=.2),numPerClust = c(50,50,50,50),
                             numWrapperSims=100,numSims=500,fractFeatIntersectThresh=.7,numFeatIntersectThresh=190,
                             clustSizeThresh=5,clustSizeFractThresh=.05,numParallelCores=3,minTrueSimilThresh=.3,
                             maxTrueSimilThresh=Inf,includeRefClustInNull=TRUE, edgeMethod=c("pearson"),
                             indEdgePvalueThresh=.1,meanEdgePairPvalueThresh=.05,restrictEdges=FALSE
)

unevenSize <- runTissueClusterSimROC(saveDir = "/home/kplaney/lungSims/",numSimDatasets=10,
                             eigenValueMin = -.001,simType=c("unevenSizeClust"),
                             noiseVector = seq(from=0,to=2.5,by=.2),numPerClust = c(50,50,50,50),
                             numWrapperSims=100,numSims=500,fractFeatIntersectThresh=.7,numFeatIntersectThresh=190,
                             clustSizeThresh=5,clustSizeFractThresh=.05,numParallelCores=3,minTrueSimilThresh=.3,
                             maxTrueSimilThresh=Inf,includeRefClustInNull=TRUE, edgeMethod=c("pearson"),
                             indEdgePvalueThresh=.1,meanEdgePairPvalueThresh=.05,restrictEdges=FALSE
)

