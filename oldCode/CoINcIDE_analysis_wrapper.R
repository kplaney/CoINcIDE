

source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_communityDetection.R")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_visualization.R")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_metaClusterAnalysis.R")

#grab data matrix list, clust features list
CoINcIDE_rankFeatures = "/home/kplaney/ovarian_analysis/metaFeatures_1000.RData.gzip"
CoINcIDE_clusterOutput =  "/home/kplaney/ovarian_analysis/curatedOvarianData_kmeansConsensus_nstart1_1000Features_2015-04-29.RData.gzip"
esets <- "/home/kplaney/ovarian_analysis/esets_proc_TCGAcombat.RData.gzip"
CoINcIDE_computeEdgesObject = "/home/kplaney/ovarian_analysis/adjMatrices_1000F_pearson_meanMatrix_2015-04-29RData.gzip"
meanEdgePairPvalueThresh <- .01
indEdgePvalueThresh <- .05
minTrueSimilThresh <- .5
maxTrueSimilThresh <- Inf
#NOTE: I have found that removing clusters below size 10 does clear up signal.
clustSizeThresh <- 5
saveDir <- "/home/kplaney/ovarian_analysis/"
experimentTag <- 
load(CoINcIDE_rankFeatures)
metaFeatures <- metaFeatures
load(CoINcIDE_clusterOutput)
clusterOutput <- kmeansConsensus
load(esets)
networkColors <- "set3"
commMethod = "edgeBetween"
minNumUniqueStudiesPerCommunity=4
experimentName="",nodePlotSize=10,nodeFontSize=.7
ES_thresh <- 1
#load("/home/kplaney/ovarian_analysis/esets_proc_TCGAcombat.RData.gzip")

#make a unique directory for all files.
dir.create(paste0(saveDir,"/",experimentTag,"_",Sys.Date()))
#CoINcIDE_computeEdgesObject <- "/home/kplaney/ovarian_analysis/adjMatrices_1000F_pearson_meanMatrix_2015-04-29RData.gzip"
load(CoINcIDE_computeEdgesObject)
output <- ov_1000F_pearson_meanMatrix
inputVariablesDF <- output$inputVariablesDF
computeTrueSimilOutput <- output$computeTrueSimilOutput
pvalueMatrix <- output$pvalueMatrix
clustIndexMatrix <- output$clustIndexMatrix

#now format just as a list of data matrices.
dataMatrixList <- exprSetListToMatrixList(esets,featureDataFieldName="gene")

names(dataMatrixList) <- names(esets)
##do for each 200,500,1000,2000 (load different metaFeatures RData object each time.)
#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}

origToNewIndexMap <- cbind(1:length(dataMatrixList),na.omit(match(names(dataMatrixList),names(esets))))
message("Assuming we're using k-means consensus PACR to choose K.")
clustSampleIndexList <-  clusterOutput$clustSampleIndexList_PACR
clustFeatureIndexList <- clusterOutput$clustFeatureIndexList_PACR


fractFeatIntersectThresh <- inputVariablesDF$fractFeatIntersectThresh
numFeatIntersectThresh <- inputVariablesDF$numFeatIntersectThresh
clustSizeFractThresh <- inputVariablesDF$clustSizeFractThresh

cat("\nInput variables used to derive the adjacency matrix:\n",append=TRUE,file=outputFile)
textOut <- capture.output(inputVariablesDF)
cat(textOut,sep="\n",append=TRUE,file=outputFile)
cat("\nThere were ",nrow(clustIndexMatrix), " total input clusters from ",length(unique(clustIndexMatrix[,2])), " studies",append=TRUE,file=outputFile)
cat("\nThe total number of input features was ",length(metaFeatures$finalFeatures),append=TRUE,file=outputFile)
cat("\nAcross the entire square (nonsymmetric) p-value matrix, there are ",length(which(pvalueMatrix<=.1))," pvalues less than or equal to .1",append=TRUE,file=outputFile)
cat("\nAcross the entire square (nonsymmetric) p-value matrix, there are ",length(which(pvalueMatrix<=.05))," pvalues less than or equal to .05",append=TRUE,file=outputFile)
cat("\nAcross the entire square (symmetric) similarity matrix, there are ",length(which(computeTrueSimilOutput$similValueMatrix>=minTrueSimilThresh))," similarities greater than or equal to ",minTrueSimilThresh,append=TRUE,file=outputFile)

finalEdgeInfo <- assignFinalEdges(computeTrueSimilOutput=computeTrueSimilOutput,pvalueMatrix=pvalueMatrix,indEdgePvalueThresh=indEdgePvalueThresh,
                                  meanEdgePairPvalueThresh=meanEdgePairPvalueThresh,
                                  minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                  fractFeatIntersectThresh=fractFeatIntersectThresh,numFeatIntersectThresh=numFeatIntersectThresh ,
                                  clustSizeThresh=clustSizeThresh, clustSizeFractThresh= clustSizeFractThresh,saveDir=saveDir,fileTag="CoINcIDE_edges_",
                                  restrictEdges=FALSE
)

commInfo <- findCommunities(edgeMatrix=finalEdgeInfo$filterEdgeOutput$edgeMatrix,edgeWeightMatrix=finalEdgeInfo$filterEdgeOutput$edgeWeightMatrix,
                            clustIndexMatrix=output$clustIndexMatrix,fileTag=experimentTag,
                            saveDir=saveDir,minNumUniqueStudiesPerCommunity=minNumUniqueStudiesPerCommunity,clustMethodName="",
                            commMethod=commMethod,
                            makePlots=TRUE,saveGraphData=TRUE,plotToScreen=FALSE)

aggregateData <- returnSampleMemberMatrix(clustSampleIndexList=clustSampleIndexList,
                                          dataMatrixList=dataMatrixList,communityInfo=commInfo)
clustSizes <- table(aggregateData$sampleClustCommKey$globalClustNum)
#add clusterSizes to node attributes
commInfo$attrDF$size <- clustSizes[na.omit(match(commInfo$attrDF$clust,names(clustSizes)))]

networkStats <- advancedNetworkPlots(communityMembership=commInfo,
                                     brewPal = networkColors,
                                     saveDir=saveDir,
                                     plotToScreen=FALSE,experimentName=experimentTag)$network_stats

message("Overall network stats (I believe the origNumClusters variable may be off-need to debug:")
networkStats
#save variable
sampleClustCommKey<-aggregateData$sampleClustCommKey
binInfo <- binarizeMetaclustStudyStatus(aggregateData$sampleClustCommKey)

message("Running effect size analysis:")
ES_out <- computeMetaclustEffectSizes(metaClustSampleNames=binInfo$metaClustSampleNames,dataMatrixList=dataMatrixList,featureNames=metaFeatures$finalFeatures,minOtherClass=5,
                                      computeWilcoxon=FALSE)
sigGenes <- selectMetaclustSigGenes(computeMetaclustEffectSizesOutput=ES_out,qvalueThresh=1,
                                    ESthresh=ES_thresh,qvalueThresh=NA)

summaryESPos <- summarizePosESMetaclustGenes(selectMetaclustSigGenesOut=sigGenes,computeMetaclustEffectSizesOutput=ES_out)
  #write to table:
write.table(summaryESPos,file=paste0(saveDir,"/",experimentName,"_summaryGenes_ESpos_thresh_",ES_thresh,"_",Sys.Date(),".txt"),row.names=TRUE,quote=FALSE,col.names=TRUE)

save(summaryESPos,file=paste0(saveDir,"/",experimentName,"_ES_genesWithThresh_ES_",ES_thresh,".RData.gzip"),compress="gzip")

GSEA_out <- list()
for(c in 1:length(sigGenes$sigMetaclustGenes_pos)){
  
  testGeneSet <- sigGenes$sigMetaclustGenes_pos[[c]]
  cat("\n",length(testGeneSet)," genes for community ", c, " with >= ES of ",ES_thresh,"\n")
  GSEA_out[[c]] <- GSEA(testGeneVector=testGeneSet,refGeneLists=NULL,method=c("hypergeometric"),
                        refGeneListDir="/home/data/MSigDB/GSEA_base_MSigDB_lists_merged.RData.gzip")
  
}

save(GSEA_out,file=paste0(saveDir,"/",experimentName,"_GSEA_out.RData.gzip"),compress="gzip")
GSEA_out_unique <- list()
for(g in 1:length(GSEA_out)){
  
  tmp <- GSEA_out[[g]]$refGeneListNames[which(GSEA_outOrig[[g]]$qvalue<=.05)] 
  for(e in 1:length(GSEA_out)){
    
    if(e != g){
      
      tmp <- setdiff(tmp, GSEA_out[[e]]$refGeneListNames[which(GSEA_outOrig[[e]]$qvalue<=.05)] )
    }
  }
  GSEA_out_unique[[g]] <- tmp
  
}
save(GSEA_out_unique,file=paste0(saveDir,"/",experimentName,"_GSEA_out_uniqueSigLists_forEachMetaCluster.RData.gzip"),compress="gzip")

