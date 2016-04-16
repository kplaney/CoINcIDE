
#per reviewer requests:  an analysis that looks at Coincide ovarian samples 
#using the second mean correlation metric density "peak" at 0.7, and 
#only including Coincide subtypes with greater than 30 samples in the survival 
#analysis:

library("Coincide")
library("survival")
outputFile <- "~/CoINcIDE_messages.txt"
saveDir <- "/home/ywrfc09/ovarian_analysis/metaRankWithTop20Genes"
globalSaveDir <- "/home/ywrfc09/ovarian_analysis/"
#200 features, pearson:
load(paste0(saveDir,"/metaFeatures_2000.RData.gzip"))
clusterCoINcIDE_output =  readRDS("/home/ywrfc09/ovarian_analysis/metaRankWithTop20Genes/curatedOvarianData_kmeansConsensus_nstart1_2000Features_2015-04-30.rds")


clustSampleIndexList <-  clusterCoINcIDE_output$clustSampleIndexList_PACR
clustFeatureIndexList <- clusterCoINcIDE_output$clustFeatureIndexList_PACR

load(paste0(globalSaveDir,"/esets_proc_TCGAcombat.RData.gzip"))
esets = esets

CoINcIDE_output = readRDS("/home/ywrfc09/ovarian_analysis/metaRankWithTop20Genes/CoINcIDE_results_ovarian2014F_pearson_edgeMethod_mean_centroidMethod2015-07-08.rds")
experimentName <- "2014F_pear_meanCent_MM7_minNumPatientOutcomesPerMetaCluster29"
minTrueSimilThresh <- 0.7
eset_featureDataFieldName="gene"
networkColors = "Set2"
outcomesVarBinary="vital_status"
outcomesVarCont = "days_to_death"
ovarian <- TRUE
fisherTestVariables <- c("histological_type","tumorstage","recurrence_status","grade","age_at_initial_pathologic_diagnosis")
fisherTestVariableLegendNames <- c("hist\ntype","tumor\nstage","recurrence","tumor\ngrade","age")
fisherTestVariableTitleNames <- c("histological type","tumor stage", "recurrence status","tumor grade","age at diagnosis")
eset_uniquePatientID= "unique_patient_ID"

metaFeatures=metaFeatures;esets=esets;CoINcIDE_output=CoINcIDE_output ; clusterCoINcIDE_output=clusterCoINcIDE_output;
                                                                               meanEdgePairPvalueThresh = .01;indEdgePvalueThresh = .01; maxTrueSimilThresh = Inf;minFractNN=.8;
                                                                               clustSizeThresh = 0;saveDir =saveDir;experimentName = experimentName;networkColors = networkColors;
                                                                               commMethod = "edgeBetween"; minNumUniqueStudiesPerCommunity=3; nodePlotSize=10;nodeFontSize=.7;ES_thresh = 0.5;eset_featureDataFieldName=eset_featureDataFieldName;
                                                                               survivalAnalysis=TRUE;outcomesVarBinary=outcomesVarBinary;outcomesVarCont = outcomesVarCont;
                                                                               CutoffPointYears=5; eset_uniquePatientID=eset_uniquePatientID; fisherTestVariables = fisherTestVariables;
                                                                               ovarian=ovarian;fisherTestVariableLegendNames=fisherTestVariableLegendNames;fisherTestVariableTitleNames=fisherTestVariableTitleNames;
                                                                               GSEAanalysis=FALSE;clinVarPlots=FALSE; fractFeatIntersectThresh=.6;numFeatIntersectThresh =0;clustSizeFractThresh =0;
                                                                               findCommWithWeights=TRUE; plotSimilEdgeWeight = TRUE;plotToScreen=FALSE;fractEdgesInVsOutEdge=0; fractEdgesInVsOutComm=0;minNumEdgesForCluster=1;nodeSizeScaleFactor=1




##copied most of metaFeatures analysis wrapper function.

#but added this variable:
  minNumPatients <- 30
if(!plotToScreen){
  #otherwise png() won't work
  options(bitmapType="cairo")
  
}

dir.create(paste0(saveDir,"/",experimentName,"_",Sys.Date()),showWarnings=TRUE);
saveDir <- paste0(saveDir,"/",experimentName,"_",Sys.Date())

outputFile <- paste0(saveDir,"/",experimentName,"_outMessages_",Sys.Date(),".txt")

summaryFile <- paste0(saveDir,"/",experimentName,"_summary_",Sys.Date(),".txt")

inputVariablesDF = list()
inputVariablesDF$meanEdgePairPvalueThresh = meanEdgePairPvalueThresh
inputVariablesDF$indEdgePvalueThresh = indEdgePvalueThresh
inputVariablesDF$minTrueSimilThresh = minTrueSimilThresh
inputVariablesDF$maxTrueSimilThresh = maxTrueSimilThresh
inputVariablesDF$outcomesVarBinary=outcomesVarBinary
inputVariablesDF$outcomesVarCont=outcomesVarCont
inputVariablesDF$CutoffPointYears=CutoffPointYears
inputVariablesDF$ES_thresh=ES_thresh
inputVariablesDF$minNumUniqueStudiesPerCommunity=minNumUniqueStudiesPerCommunity
inputVariablesDF$commMethod = commMethod
inputVariablesDF$ovarian = ovarian
inputVariablesDF$fractFeatIntersectThresh <- fractFeatIntersectThresh
inputVariablesDF$numFeatIntersectThresh <- numFeatIntersectThresh 
inputVariablesDF$clustSizeFractThresh <- clustSizeFractThresh 
inputVariablesDF$clustSizeThresh <- clustSizeThresh 
cat(gsub(" ","_","CoINcIDE input variables used for this analysis:\n"),
    append=TRUE,file=outputFile)
textOut <-capture.output(inputVariablesDF)
cat("\n",textOut,sep="\n",append=TRUE,file=outputFile)


computeTrueSimilOutput = CoINcIDE_output$computeTrueSimilOutput
pvalueMatrix = CoINcIDE_output$pvalueMatrix
clustIndexMatrix = CoINcIDE_output$clustIndexMatrix

if(!is.null(esets)){
  
  #now format just as a list of data matrices.
  dataMatrixList <- exprSetListToMatrixList(esets,featureDataFieldName=eset_featureDataFieldName)
  
  names(dataMatrixList) <- names(esets)
  ##do for each 200,500,1000,2000 (load different metaFeatures RData object each time.)
  #remove datasets with too many missing top gene features
  if(length(metaFeatures$datasetListIndicesToRemove)>0){
    
    dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
    
  }
  
  
  cat(gsub(" ","_",paste0("\nStarted with  ",length(dataMatrixList), " datasets\n:")),
      append=TRUE,file=outputFile)
  textOut <- capture.output(names(dataMatrixList))
  cat("\n",textOut,sep="\n",append=TRUE,file=outputFile)
  
  #need mapping later for pheno data.
  if(length(dataMatrixList) != length(na.omit(match(names(dataMatrixList),names(esets))))){
    
    stop("old to new indexing for esets and data matrixlist not matching up.")
  }
  
  origToNewIndexMap <- cbind(1:length(dataMatrixList),na.omit(match(names(dataMatrixList),names(esets))),names(esets)[na.omit(match(names(dataMatrixList),names(esets)))])
  tmp <- origToNewIndexMap[ ,c(1,3)]
  
}else{
  
  tmp <- cbind(c(1:length(dataMatrixList)),names(dataMatrixList))
  
}
colnames(tmp) <- c("datasetNumber","datasetName")
write.table(tmp,file=paste0(saveDir,"/",experimentName,"_datasetNameNumberKey.txt"),quote=FALSE,row.names=TRUE,col.names=FALSE)

cat(gsub(" ","_",paste0("\nAssuming we're using k-means consensus PACR to choose K\n.")),file=outputFile,append=TRUE)
textOut <- capture.output(tmp)
cat(textOut,sep="\n",append=TRUE,file=summaryFile)
clustSampleIndexList <-  clusterCoINcIDE_output$clustSampleIndexList_PACR
clustFeatureIndexList <- clusterCoINcIDE_output$clustFeatureIndexList_PACR
#cat("\nK for each dataset: \n",append=TRUE,file=outputFile)
#capture.output(names(clusterCoINcIDE_output$bestK_PACR),file=outputFile,append=TRUE)
# capture.output(unlist(clusterCoINcIDE_output$bestK_PACR),file=outputFile,append=TRUE)

#write this to a file
numClusters <- clusterCoINcIDE_output$bestK_PACR
numClusters <- as.matrix(unlist(numClusters))
colnames(numClusters) <- "PACR_bestK"
write.table(numClusters,row.names=TRUE,col.name=TRUE,quote=FALSE,file=paste0(saveDir,"/",experimentName,"_numClustersForDatasets_",Sys.Date(),".txt"))
textOut <- capture.output(numClusters)
cat(gsub(" ","_","\nDataset number of clusters:\n"),textOut,sep="\n",
    append=TRUE,file=summaryFile)


cat(gsub(" ","_",paste0("\nThere were ",nrow(clustIndexMatrix), " total input clusters from ",length(unique(clustIndexMatrix[,2])), " studies")),append=TRUE,file=summaryFile)
cat(gsub(" ","_",paste0("\nThe total number of input features was ",length(metaFeatures$finalFeatures))),append=TRUE,file=summaryFile)
cat(gsub(" ","_",paste0("\nAcross the entire square (nonsymmetric) p-value matrix, there are ",length(which(pvalueMatrix<=.1))," pvalues less than or equal to .1")),append=TRUE,file=summaryFile)
cat(gsub(" ","_",paste0("\nAcross the entire square (nonsymmetric) p-value matrix, there are ",length(which(pvalueMatrix<=.05))," pvalues less than or equal to .05")),append=TRUE,file=summaryFile)
cat(gsub(" ","_",paste0("\nAcross the entire square (nonsymmetric) p-value matrix, there are ",length(which(pvalueMatrix<=.01))," pvalues less than or equal to .01")),append=TRUE,file=summaryFile)
cat(gsub(" ","_",paste0("\nAcross the entire square (symmetric) similarity matrix, there are ",length(which(computeTrueSimilOutput$similValueMatrix>=minTrueSimilThresh))/2,
                        " similarities greater than or equal to ",minTrueSimilThresh,"not double-counting.")),append=TRUE,file=summaryFile)

finalEdgeInfo <- assignFinalEdges(computeTrueSimilOutput=computeTrueSimilOutput,pvalueMatrix=pvalueMatrix,indEdgePvalueThresh=indEdgePvalueThresh,
                                  meanEdgePairPvalueThresh=meanEdgePairPvalueThresh,
                                  minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                  fractFeatIntersectThresh=fractFeatIntersectThresh,numFeatIntersectThresh=numFeatIntersectThresh ,
                                  clustSizeThresh=clustSizeThresh, clustSizeFractThresh= clustSizeFractThresh,saveDir=saveDir,fileTag="CoINcIDE_edges_",
                                  minFractNN =.7,minNumEdgesForCluster=minNumEdgesForCluster
)

commInfo <- findCommunities(edgeMatrix=finalEdgeInfo$filterEdgeOutput$edgeMatrix,edgeWeightMatrix=finalEdgeInfo$filterEdgeOutput$edgeWeightMatrix,
                            clustIndexMatrix=CoINcIDE_output$clustIndexMatrix,fileTag=experimentName,
                            saveDir=saveDir,minNumUniqueStudiesPerCommunity=minNumUniqueStudiesPerCommunity,experimentName=experimentName,
                            commMethod=commMethod,
                            makePlots=TRUE,saveGraphData=TRUE,plotToScreen=plotToScreen, findCommWithWeights=findCommWithWeights, plotSimilEdgeWeight = plotSimilEdgeWeight,
                            fractEdgesInVsOutComm=fractEdgesInVsOutComm, fractEdgesInVsOutEdge=fractEdgesInVsOutEdge)

textOut <- capture.output(commInfo$finalCommEdgeInfo)
cat(gsub(" ","_","\nFinal community stats:\n"),textOut,sep="\n",append=TRUE,file=paste0(saveDir,"/",experimentName,"_finalCommStats.txt"))
cat(gsub(" ","_","\nFinal community stats:\n"),textOut,sep="\n",
    append=TRUE,file=summaryFile)

aggregateData <- returnSampleMemberMatrix(clustSampleIndexList=clustSampleIndexList,
                                          dataMatrixList=dataMatrixList,communityInfo=commInfo)
clustSizes <- table(aggregateData$sampleClustCommKey$globalClustNum)
textOut <- capture.output(table(aggregateData$sampleClustCommKey$community,exclude=FALSE))
cat(gsub(" ","_","\nMeta-cluster sizes:\n"),textOut,sep="\n",
    append=TRUE,file=summaryFile)
#add clusterSizes to node attributes
commInfo$attrDF$size <- clustSizes[na.omit(match(as.numeric(as.character(commInfo$attrDF$clust)),as.numeric(names(clustSizes))))]
textOut <- capture.output(table(aggregateData$sampleClustCommKey$community,aggregateData$sampleClustCommKey$studyNum,exclude=FALSE))
cat(gsub(" ","_","\nMeta-cluster number of samples by studies:\n"),textOut,sep="\n",append=TRUE,file=paste0(saveDir,"/",experimentName,"_metaclusterNumSamplesByStudies.txt"))
cat(gsub(" ","_","\nMeta-cluster number of samples by studies:\n"),textOut,sep="\n",
    append=TRUE,file=summaryFile)
networkStats <- advancedNetworkPlots(communityMembership=commInfo,
                                     brewPal = networkColors,
                                     saveDir=saveDir,clustIndexMatrix=clustIndexMatrix,
                                     plotToScreen=plotToScreen,experimentName=experimentName,
                                     plotEdgeWeight=plotSimilEdgeWeight,nodeSizeScaleFactor=nodeSizeScaleFactor)$network_stats


write.table(networkStats,quote=FALSE,file=paste0(saveDir,"/",experimentName,"_networkStats_",Sys.Date(),".txt"))

textOut <- capture.output(networkStats)
cat(gsub(" ","_","\nNetwork stats:\n"),textOut,sep="\n",
    file=summaryFile,append=TRUE)
#save variable
sampleClustCommKey<-aggregateData$sampleClustCommKey
binInfo <- binarizeMetaclustStudyStatus(aggregateData$sampleClustCommKey)


if(clinVarPlots || survivalAnalysis){
  
  clinicalTables <- list()
  
  #already have esets loaded:
  #load("/home/kplaney/ovarian_analysis/esets_proc_TCGAcombat.RData.gzip")
  phenoMasterDF <- createPhenoMasterTableFromMatrixList(esetList=esets)
  
  #load("/home/kplaney/ovarian_analysis/curatedOvarian_phenoMasterDF.RData.gzip")
  #study numbers won't align here because some filtered out (had no robust clusters); need to "translate"
  
  origToNewIndexMap <- data.frame(origToNewIndexMap,stringsAsFactors=FALSE)
  colnames(origToNewIndexMap) <- c("studyNum","origStudyNum","esetName")
  #want to intersect on "new" study numbers, as this is what sampleClustCommKey has.
  #NOTE: graph will override this join function....it also has a "join"!
  #RDavidWebservice was the culprit...
  #detach("package:graph",unload=TRUE,force=TRUE)
  #left: only take studies in the sampleClustCommKey.
  sampleClustCommKey <- join(sampleClustCommKey,origToNewIndexMap,by="studyNum",type="left",match="all")
  sampleClustCommPhenoData <- addClinicalVarToNodeAttributes(sampleClustCommKey,phenoMasterDF=phenoMasterDF)
  
  
  sampleClustCommPhenoDataOrig <- sampleClustCommPhenoData

}

expName <- gsub("_"," ",experimentName)

result <- list()

linearModels <- list()
linearModelsSumm <- list()

if(survivalAnalysis){
  cat(gsub(" ","_","\nSurvival analysis:\n"),append=TRUE,file=summaryFile)
  if(ovarian){
    
    message("running survival analysis on ovarian")
    #groupings <- groupings[which(!is.na(sampleClustCommPhenoData[,outcomesVarBinary]))];
    #outcomesData <- sampleClustCommPhenoData[which(!is.na(sampleClustCommPhenoData[,outcomesVarBinary])),];
    
    #keep samples with NA days to event for now?
    #hmm...one meta-cluster is left out if use "days_to_death"...
    #groupings <- as.numeric(as.factor(groupings[which(!is.na(outcomesData[,outcomesVarCont]))]))
    #outcomesData <- outcomesData[which(!is.na(outcomesData[,outcomesVarCont])),];
    
    #if binary is character string categories: make it a factor first, then numeric,
    #otherwise coxph function will throw errors.
    #  nonCensoredTerm=1
    sampleClustCommPhenoData[which(sampleClustCommPhenoData[ ,outcomesVarBinary]=="deceased"),outcomesVarBinary] <- 1
    sampleClustCommPhenoData[which(sampleClustCommPhenoData[ ,outcomesVarBinary]=="living"),outcomesVarBinary] <- 0
    
    ##ALTER HERE: only take patients from meta-clusters with lots of samples
    numSurv <- table( sampleClustCommPhenoData[,"community"], sampleClustCommPhenoData[,"vital_status"])
    #minNumPatients variable denotes this:
    commToKeep <- rownames(numSurv)[which(rowSums(numSurv)>minNumPatients)]
    
    for(c in 1:length(commToKeep)){
    
      if(c==1){
        
        tmp <- sampleClustCommPhenoData[which(sampleClustCommPhenoData[,"community"]==commToKeep[c]), ]
        
      }else{
        
      tmp <- rbind(tmp,sampleClustCommPhenoData[which(sampleClustCommPhenoData[,"community"]==commToKeep[c]), ])
                   
      }
                                      
    }
    
    sampleClustCommPhenoData <- tmp
    
    groupingTerm="community"
    
    #only take samples with the groupingTerm you're looking at.
    sampleClustCommPhenoData <- sampleClustCommPhenoData[which(!is.na(sampleClustCommPhenoData[, groupingTerm])), ]
    #remove samples with NA values.
    groupings <- as.factor(sampleClustCommPhenoData[, groupingTerm])
    
    cat(gsub(" ", "_",paste0("\n",length(groupings), " patients included across all final meta-clusters.\n")),
        append=TRUE,file=summaryFile)
    
    outcomesDataShort <- data.frame(as.numeric(sampleClustCommPhenoData[,outcomesVarBinary]),as.numeric(sampleClustCommPhenoData[,outcomesVarCont])
    );
    
    #sometimes the names are duplicated across studies - remove this line
    #rownames(outcomesDataShort ) <- outcomesData[,eset_uniquePatientID];
    colnames(outcomesDataShort) <- c("Censoring","TimeToLastContactOrEvent")
    
    nonCensoredTerm=1
    censoredTerm=0
    Survival <- outcomesDataShort
    if(!is.factor(groupings)){
      
      groupings <- as.factor(groupings)
      
    }
    nonNAs <- intersect(which(!is.na(sampleClustCommPhenoData[,outcomesVarBinary])), which(!is.na(sampleClustCommPhenoData[,outcomesVarCont])))
    numOutcomesPerMetacluster <- table(sampleClustCommPhenoData[nonNAs,"studyNum"],groupings[nonNAs])
    
    textOut <- capture.output(numOutcomesPerMetacluster)
    cat("numNonNAs_info\n",textOut,sep="\n",file=paste0(saveDir,"/",experimentName,"_numNonNAOutcomesByMetacluster_",Sys.Date(),".txt"))
    
    cat(gsub(" ","_","\nnumNonNAs_info\n:\n"),textOut,sep="\n",
        append=TRUE,file=summaryFile)
    
    #creating the survival objects with the time and censoring variables
    #EVENT is nonCensoredTerm
    #For counting process data, event indicates whether an event occurred at the end of the interval
    #if provide a true/false like this: Surv knows that this is the event variable
    #From Surv documentation:
    #When the type argument is missing the code assumes a type based on the following rules:
    #If there are two unnamed arguments, they will match time and event in that order. 
    OverallSurvival <- Surv(Survival$TimeToLastContactOrEvent,Survival$Censoring==nonCensoredTerm);
    #creating a survival object cutoff at a certain point
    CutoffPoint <- CutoffPointYears*365;
    CutoffSamples=Survival$TimeToLastContactOrEvent>CutoffPoint & !is.na(Survival$TimeToLastContactOrEvent)
    SurvivalCutoff=Survival
    SurvivalCutoff$TimeToLastContactOrEvent[CutoffSamples]=CutoffPoint
    SurvivalCutoff$Censoring[CutoffSamples]=censoredTerm
    #   #"Surv" creates a survival object. really for binary outcomes data.
    OverallSurvivalCutoff=Surv(SurvivalCutoff$TimeToLastContactOrEvent,SurvivalCutoff$Censoring==nonCensoredTerm)
    
    coxfit=coxph(OverallSurvival~groupings, data=Survival)
    #message("coxfit summary for overall survival.")
    sumCoxfit_overall <- summary(coxfit)
    result$cox_overall_likelihoodPvalue <- sumCoxfit_overall$logtest["pvalue"]
    result$cox_overall_nevent <- sumCoxfit_overall$nevent
    #see page 7 - helps clarify that $logtest is NOT the logrank test - sctest is:
    #https://cran.r-project.org/web/packages/HSAUR/vignettes/Ch_survival_analysis.pdf
    result$cox_overall_logRankTestPvalue <- sumCoxfit_overall$sctest["pvalue"]
    result$cox_overall_rsq <- sumCoxfit_overall$rsq["rsq"]
    result$cox_overall_waldPvalue <- sumCoxfit_overall$waldtest["pvalue"]
    #plot(cox.zph(coxfit))
    kmfit=survdiff(OverallSurvival ~ groupings)
    #message("kaplan meier p-value for overall survival:")
    result$km_overall_pvalue <- 1 - pchisq(kmfit$chisq, length(kmfit$n) - 1)
    
    #  message("calculating the sign of the survival relationship")
    
    colorCodes <- rev(brewer.pal(length(unique(commInfo$attrDF[,"community"])),networkColors));
    #want color codes to be same as in networ plot.
    colorCodes <- colorCodes[na.omit(match(unique(groupings),unique(commInfo$attrDF[,"community"])))]
    
    mfit=survfit(OverallSurvival ~ groupings)
    
    if(!plotToScreen){
      
      png(filename=paste0(saveDir,"/",experimentName,"_kaplanMeier_OS_",Sys.Date(),".png"),width=1000,height=1000,res=160)
      
      plot(mfit,main=paste0("Overall survival for ",expName), lty=c(1:length(unique(groupings))),col= colorCodes)
      
      legend(60,.35,unique(groupings),lty=c(1:length(unique(groupings))),col= colorCodes)
      
      dev.off()
      
      
    }else{
      
      legend(60,.35,unique(groupings),lty=c(1:length(unique(groupings))),col= colorCodes)
      plot(mfit,main=paste0("Overall survival for ",expName), lty=c(1:length(unique(groupings))),col= colorCodes)
      
      
      
    }
    
    
    mfitCut <- survfit(OverallSurvivalCutoff ~ groupings)
    
    
    if(!plotToScreen){
      
      png(filename=paste0(saveDir,"/",experimentName,"_kaplanMeier_OS",CutoffPointYears,"years_",Sys.Date(),".png"),width=1000,height=1000,res=160)
      
      plot(mfitCut,main=paste0("Overall survival for ",expName, "\nwith cutoff at ",
                               CutoffPointYears," years"),lty=c(1:length(unique(groupings))),col= colorCodes)
      
      legend(25,.35,unique(groupings),lty=c(1:length(unique(groupings))),col= colorCodes)
      
      dev.off()
      
    }else{
      
      
      legend(25,.35,unique(groupings),lty=c(1:length(unique(groupings))),col= colorCodes)
      plot(mfitCut,main=paste0("Overall survival for ",expName, "\nwith cutoff at ",
                               CutoffPointYears," years"),lty=c(1:length(unique(groupings))),col= colorCodes)
      
      
    }
    coxfit=coxph(OverallSurvivalCutoff~groupings, data= SurvivalCutoff)
    # message("coxfit summary for survival cutoff at ",CutoffPointYears,"years:")
    
    #plot(cox.zph(coxfit))
    sumCoxfit_cut  <- summary(coxfit)
    result$cox_cut_likelihoodPvalue <- sumCoxfit_cut$logtest["pvalue"]
    result$cox_cut_nevent <- sumCoxfit_cut$nevent
    result$cox_cut_logRankTestPvalue <- sumCoxfit_cut$sctest["pvalue"]
    result$cox_cut_rsq <- sumCoxfit_cut$rsq["rsq"]
    result$cox_cut_waldPvalue <- sumCoxfit_cut$waldtest["pvalue"]
    
    kmfit=survdiff(OverallSurvivalCutoff ~ groupings)
    
    result$km_pv_cutoff <- 1 - pchisq(kmfit$chisq, length(kmfit$n) - 1) 
    
    
    tmp <- as.matrix(unlist(result))
    write.table(tmp,file=paste0(saveDir,"/",experimentName,"_survivalAnalysis_",Sys.Date(),".txt"),col.names=FALSE)
    textOut <- capture.output(tmp)
    cat(gsub(" ","_","\nSuvival models:\n"),textOut,sep="\n",
        append=TRUE,file=summaryFile)
    
    #end of survival analysis.
  }else{
    #breast analysis
  }
}
output <- list(sampleClustCommPhenoData=sampleClustCommPhenoData,aggregateData=aggregateData,commInfo=commInfo,finalEdgeInfo=finalEdgeInfo,CoINcIDE_computeEdgesObject=CoINcIDE_output,
               networkStats=networkStats,survivalResults=result,
               linearModelSummaries=linearModelsSumm,linearSurvModels=linearModels,
               clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
               clinicalTables=clinicalTables)

