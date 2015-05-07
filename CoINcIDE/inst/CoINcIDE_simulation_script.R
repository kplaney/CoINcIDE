
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_simulation.R")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_computeEdges.R")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_communityDetection.R")

highQuality <- runLungSimROC(saveDir = "/home/kplaney/lungSims/",numSimDatasets=10,
                          eigenValueMin = -.001,simType=c("highQualityClust"),
                          noiseVector = seq(from=0,to=2.5,by=.2),numPerClust = c(50,50,50,50),
                          numWrapperSims=100,numSims=500,fractFeatIntersectThresh=.7,numFeatIntersectThresh=190,
                          clustSizeThresh=5,clustSizeFractThresh=.05,numParallelCores=3,minTrueSimilThresh=.3,
                          maxTrueSimilThresh=Inf,includeRefClustInNull=TRUE, edgeMethod=c("pearson"),
                          indEdgePvalueThresh=.1,meanEdgePairPvalueThresh=.05,restrictEdges=FALSE
)
  
mixQuality <- runLungSimROC(saveDir = "/home/kplaney/lungSims/",numSimDatasets=10,
                             eigenValueMin = -.001,simType=c("mixedClustQualityClust"),
                             noiseVector = seq(from=0,to=2.5,by=.2),numPerClust = c(50,50,50,50),
                             numWrapperSims=100,numSims=500,fractFeatIntersectThresh=.7,numFeatIntersectThresh=190,
                             clustSizeThresh=5,clustSizeFractThresh=.05,numParallelCores=3,minTrueSimilThresh=.3,
                             maxTrueSimilThresh=Inf,includeRefClustInNull=TRUE, edgeMethod=c("pearson"),
                             indEdgePvalueThresh=.1,meanEdgePairPvalueThresh=.05,restrictEdges=FALSE
)

unevenQuality <- runLungSimROC(saveDir = "/home/kplaney/lungSims/",numSimDatasets=10,
                             eigenValueMin = -.001,simType=c("unevenSizeClust"),
                             noiseVector = seq(from=0,to=2.5,by=.2),numPerClust = c(50,50,50,50),
                             numWrapperSims=100,numSims=500,fractFeatIntersectThresh=.7,numFeatIntersectThresh=190,
                             clustSizeThresh=5,clustSizeFractThresh=.05,numParallelCores=3,minTrueSimilThresh=.3,
                             maxTrueSimilThresh=Inf,includeRefClustInNull=TRUE, edgeMethod=c("pearson"),
                             indEdgePvalueThresh=.1,meanEdgePairPvalueThresh=.05,restrictEdges=FALSE
)

