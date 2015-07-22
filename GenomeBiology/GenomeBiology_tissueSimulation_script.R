
library("CoINcIDE")
source("/home/ywrfc09//CoINcIDE/coincide/CoINcIDE_packageVersion//CoINcIDE/R/CoINcIDE_simulation.R")

minTrueSimilVector <- seq(from=0,to=1,by=.1)
saveDir <- "/home/ywrfc09/simulations"

for(m in 1:length(minTrueSimilVector)){
#this is NOT for communities - just for edges.
highQuality <- runTissueClusterSimROC(saveDir="./",numSimDatasets=10,
                                   eigenValueMin = -.001,simType=c("highQualityClust"),
                                   noiseVector = seq(from=0,to=2.5,by=.2),numPerClust = c(50,50,50,50),
                                   numWrapperSims=50,numSims=500,fractFeatIntersectThresh=.7,numFeatIntersectThresh=190,
                                   clustSizeThresh=0,clustSizeFractThresh=0,minTrueSimilThresh=minTrueSimilVector[m],
                                   maxTrueSimilThresh=Inf,edgeMethod=c("pearson"),centroidMethod=c("mean"),
                                   indEdgePvalueThresh=.01,meanEdgePairPvalueThresh=.01,
                                   minFractNN=.8,minRandNumClust=2,randNumClust=FALSE,minRandSize=1,maxRandSize=100
)

saveRDS(highQuality,file=paste0(saveDir,"/highQuality_minSimil_",minTrueSimilVector[m],".rds"),compress=TRUE)

}

for(m in 1:length(minTrueSimilVector)){
highQualityUnevenSize <- runTissueClusterSimROC(saveDir="./",numSimDatasets=10,
                                              eigenValueMin = -.001,simType=c("unevenSizeClust"),
                                              noiseVector = seq(from=0,to=2.5,by=.2),numPerClust = c(50,50,50,50),
                                              numWrapperSims=50,numSims=500,fractFeatIntersectThresh=.7,numFeatIntersectThresh=190,
                                              clustSizeThresh=0,clustSizeFractThresh=0,minTrueSimilThresh=minTrueSimilVector[m],
                                              maxTrueSimilThresh=Inf,edgeMethod=c("pearson"),centroidMethod=c("mean"),
                                              indEdgePvalueThresh=.01,meanEdgePairPvalueThresh=.01,
                                              minFractNN=.8,minRandNumClust=2,randNumClust=FALSE,minRandSize=1,maxRandSize=100
)

saveRDS(highQualityUnevenSize,file=paste0(saveDir,"/highQualityUnevenSize_minSimil_",minTrueSimilVector[m],".rds"),compress=TRUE)

}


for(m in 1:length(minTrueSimilVector)){
  
highQualityUnevenNumClustMin2 <- runTissueClusterSimROC(saveDir="./",numSimDatasets=10,
                                      eigenValueMin = -.001,simType=c("highQualityClust"),
                                      noiseVector = seq(from=0,to=2.5,by=.2),numPerClust = c(50,50,50,50),
                                      numWrapperSims=50,numSims=500,fractFeatIntersectThresh=.7,numFeatIntersectThresh=190,
                                      clustSizeThresh=0,clustSizeFractThresh=0,minTrueSimilThresh=minTrueSimilVector[m],
                                      maxTrueSimilThresh=Inf,edgeMethod=c("pearson"),centroidMethod=c("mean"),
                                      indEdgePvalueThresh=.01,meanEdgePairPvalueThresh=.01,
                                      minFractNN=.8,minRandNumClust=2,randNumClust=TRUE,minRandSize=1,maxRandSize=100
)

saveRDS(highQualityUnevenNumClustMin2,file=paste0(saveDir,"/highQualityUnevenNumClustMin2_minSimil_",minTrueSimilVector[m],".rds"),compress=TRUE)

}


for(m in 1:length(minTrueSimilVector)){
  
highQualityUnevenNumClustMin1 <- runTissueClusterSimROC(saveDir="./",numSimDatasets=10,
                                                        eigenValueMin = -.001,simType=c("highQualityClust"),
                                                        noiseVector = seq(from=0,to=2.5,by=.2),numPerClust = c(50,50,50,50),
                                                        numWrapperSims=50,numSims=500,fractFeatIntersectThresh=.7,numFeatIntersectThresh=190,
                                                        clustSizeThresh=0,clustSizeFractThresh=0,minTrueSimilThresh=minTrueSimilVector[m],
                                                        maxTrueSimilThresh=Inf,edgeMethod=c("pearson"),centroidMethod=c("mean"),
                                                        indEdgePvalueThresh=.01,meanEdgePairPvalueThresh=.01,
                                                        minFractNN=.8,minRandNumClust=1,randNumClust=TRUE,minRandSize=1,maxRandSize=100
)

saveRDS(highQualityUnevenNumClustMin1,file=paste0(saveDir,"/highQualityUnevenNumClustMin1_minSimil_",minTrueSimilVector[m],".rds"),compress=TRUE)

}


for(m in 1:length(minTrueSimilVector)){

highQualityUnevenSizeUnevenNumClustMin2 <- runTissueClusterSimROC(saveDir="./",numSimDatasets=10,
                                                        eigenValueMin = -.001,simType=c("unevenSizeClust"),
                                                        noiseVector = seq(from=0,to=2.5,by=.2),numPerClust = c(50,50,50,50),
                                                        numWrapperSims=50,numSims=500,fractFeatIntersectThresh=.7,numFeatIntersectThresh=190,
                                                        clustSizeThresh=0,clustSizeFractThresh=0,minTrueSimilThresh=minTrueSimilVector[m],
                                                        maxTrueSimilThresh=Inf,edgeMethod=c("pearson"),centroidMethod=c("mean"),
                                                        indEdgePvalueThresh=.01,meanEdgePairPvalueThresh=.01,
                                                        minFractNN=.8,minRandNumClust=2,randNumClust=TRUE,minRandSize=1,maxRandSize=100
)


saveRDS(highQualityUnevenSizeUnevenNumClustMin2,file=paste0(saveDir,"/highQualityUnevenSizeUnevenNumClustMin2_minSimil_",minTrueSimilVector[m],".rds"),compress=TRUE)

}

for(m in 1:length(minTrueSimilVector)){

highQualityUnevenSizeUnevenNumClustMin1 <- runTissueClusterSimROC(saveDir="./",numSimDatasets=10,
                                                                          eigenValueMin = -.001,simType=c("unevenSizeClust"),
                                                                          noiseVector = seq(from=0,to=2.5,by=.2),numPerClust = c(50,50,50,50),
                                                                          numWrapperSims=50,numSims=500,fractFeatIntersectThresh=.7,numFeatIntersectThresh=190,
                                                                          clustSizeThresh=0,clustSizeFractThresh=0,minTrueSimilThresh=minTrueSimilVector[m],
                                                                          maxTrueSimilThresh=Inf,edgeMethod=c("pearson"),centroidMethod=c("mean"),
                                                                          indEdgePvalueThresh=.01,meanEdgePairPvalueThresh=.01,
                                                                          minFractNN=.8,minRandNumClust=1,randNumClust=TRUE,minRandSize=1,maxRandSize=100
)

saveRDS(highQualityUnevenSizeUnevenNumClustMin1,file=paste0(saveDir,"/highQualityUnevenSizeUnevenNumClustMin1_minSimil_",minTrueSimilVector[m],".rds"),compress=TRUE)

}


for(m in 1:length(minTrueSimilVector)){
  
mixedQualityClust <- runTissueClusterSimROC(saveDir="./",numSimDatasets=10,
                                            eigenValueMin = -.001,simType=c("mixedClustQualityClust"),
                                            noiseVector = seq(from=0,to=2.5,by=.2),numPerClust = c(50,50,50,50),
                                            numWrapperSims=50,numSims=500,fractFeatIntersectThresh=.7,numFeatIntersectThresh=190,
                                            clustSizeThresh=0,clustSizeFractThresh=0,minTrueSimilThresh=minTrueSimilVector[m],
                                            maxTrueSimilThresh=Inf,edgeMethod=c("pearson"),centroidMethod=c("mean"),
                                            indEdgePvalueThresh=.01,meanEdgePairPvalueThresh=.01,
                                            minFractNN=.8,minRandNumClust=1,randNumClust=TRUE,minRandSize=1,maxRandSize=100
)

saveRDS(mixedQualityClust,file=paste0(saveDir,"/mixedQualityClust_minSimil_",minTrueSimilVector[m],".rds"),compress=TRUE)

}