library("pROC")
library("survival")
library("Coincide")

##wrapper function that analyzes a wide range of stats after run Coincide on the 
#curatedBreast or Ovarian data.  This is very tailored to producing outputs for the 
#Genome Medicine publication, and will not immediately translate to analyzing any 
#type of dataset, as the survival analyses are very data-specific..

#also: note that an path to reference gene lists for GSEA is needed, unless you set 
#GSEAanalysis=FALSE
metaFeaturesAnalysisWrapper <- function(metaFeatures,esets=NULL,CoINcIDE_output , clusterCoINcIDE_output,
                                        meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,nodeSizeScaleFactor=1,
                                        clustSizeThresh = 5,saveDir = "/home/kplaney/ovarian_analysis/",experimentName = "ovarian_2000F",networkColors = "Set3",
                                        commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=3, minMedianNumEdgesPerNodeInCommunity=3,nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName="gene",
                                        survivalAnalysis=TRUE,GSEAanalysis=FALSE,clinVarPlots=TRUE,outcomesVarBinary="vital_status",outcomesVarCont = "days_to_death",
                                        CutoffPointYears=5, eset_uniquePatientID="unique_patient_ID", ovarian=TRUE, fisherTestVariables = c("histological_type","tumorstage","recurrence_status","grade"),
                                        fisherTestVariableLegendNames=fisherTestVariables,fisherTestVariableTitleNames=fisherTestVariables,plotStackedYLimit=NA,minFractNN =.7,
                                        fractFeatIntersectThresh=.5,numFeatIntersectThresh =0,clustSizeFractThresh =0, findCommWithWeights=FALSE, plotSimilEdgeWeight = TRUE,plotToScreen=FALSE,
                                        refGeneListDir="~/GSEA_base_MSigDB_lists_merged.RData.gzip",minNumEdgesForCluster=1,fractEdgesInVsOutEdge=0, fractEdgesInVsOutComm=0
){
  
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
                                    minFractNN =minFractNN,minNumEdgesForCluster=minNumEdgesForCluster
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
  
  if(GSEAanalysis){
    message("Running effect size analysis")
    ES_out <- computeMetaclustEffectSizes(metaClustSampleNames=binInfo$metaClustSampleNames,dataMatrixList=dataMatrixList,featureNames=metaFeatures$finalFeatures,minOtherClass=5,
                                          computeWilcoxon=FALSE)
    sigGenes <- selectMetaclustSigGenes(computeMetaclustEffectSizesOutput=ES_out,qvalueThresh=1,
                                        ESthresh=ES_thresh,includeWilcoxon=FALSE)
    
    summaryESPos <- summarizePosESMetaclustGenes(selectMetaclustSigGenesOut=sigGenes,computeMetaclustEffectSizesOutput=ES_out)
    
    summaryESNeg <- summarizeNegESMetaclustGenes(selectMetaclustSigGenesOut=sigGenes,computeMetaclustEffectSizesOutput=ES_out)
    #write to table:
    write.table(summaryESPos,file=paste0(saveDir,"/",experimentName,"_summaryGenes_ESpos_thresh_",ES_thresh,"_",Sys.Date(),".txt"),row.names=TRUE,quote=FALSE,col.names=TRUE)
    
    write.table(summaryESNeg,file=paste0(saveDir,"/",experimentName,"_summaryGenes_ESNeg_thresh_-",ES_thresh,"_",Sys.Date(),".txt"),row.names=TRUE,quote=FALSE,col.names=TRUE)
   
   # for(g in 1:length(sigGenes$sigMetaclustGenesES_pos)){
      #what can't remove indices??
    #  tmp <- as.matrix(unlist(sigGenes$sigMetaclustGenesES_pos[[g]]))
    #  write.table(tmp,file=paste0(saveDir,"/",experimentName,"_genesESpos_thresh_",ES_thresh,"_community_",g,".txt"),row.names=FALSE,col.names=FALSE,quote=FALSE)
      
    #}
    
    
    save(summaryESPos,file=paste0(saveDir,"/",experimentName,"_ES_genesWithThresh_ES_",ES_thresh,".RData.gzip"),compress="gzip")
    message("running GSEA")
    GSEA_out <- list()
    for(c in 1:length(sigGenes$sigMetaclustGenesES_pos)){
      
      testGeneSet <- sigGenes$sigMetaclustGenesES_pos[[c]]
      cat(gsub(" ","_",paste0(length(testGeneSet)," genes for community ", c, " with >= ES of ",ES_thresh,"\n")),
          append=TRUE,file=outputFile)
      GSEA_out[[c]] <- GSEA(testGeneVector=testGeneSet,refGeneLists=NULL,method=c("hypergeometric"),
                            refGeneListDir=refGeneListDir)
      
    }
    
    save(GSEA_out,file=paste0(saveDir,"/",experimentName,"_GSEA_out.RData.gzip"),compress="gzip")
    GSEA_out_unique <- list()
    for(g in 1:length(GSEA_out)){
      
      tmp <- GSEA_out[[g]]$refGeneListNames[which(GSEA_out[[g]]$qvalue<=.05)] 
      for(e in 1:length(GSEA_out)){
        
        if(e != g){
          
          tmp <- setdiff(tmp, GSEA_out[[e]]$refGeneListNames[which(GSEA_out[[e]]$qvalue<=.05)] )
        }
      }
      GSEA_out_unique[[g]] <- tmp
      cat(gsub(" ","_",paste0("\n",length(tmp)," unique signficant GSEA ref gene lists community ", g, " with ES of >=",ES_thresh,"\n")),
          append=TRUE,file=outputFile)
      
    }
    save(GSEA_out_unique,file=paste0(saveDir,"/",experimentName,"_GSEA_out_uniqueSigLists_forEachMetaCluster.RData.gzip"),compress="gzip")
    
    for(g in 1:length(GSEA_out_unique)){
      #what can't remove indices??
      write.table(as.character(unlist(GSEA_out_unique[[g]])),file=paste0(saveDir,"/",experimentName,"_GSEA_uniqueSig_community_",g,".txt"),row.names=TRUE,quote=FALSE,col.names=FALSE)
      
    }
    
    GSEA_out_sig <- list()
    for(g in 1:length(GSEA_out)){
      
      GSEA_out_sig[[g]] <- GSEA_out[[g]]$refGeneListNames[which(GSEA_out[[g]]$qvalue<=.05)] 
      cat(gsub(" ","_",paste0("\n",length(tmp)," overall (not necessarily unique) signficant GSEA ref gene lists community ", g, " with >= ES of ",ES_thresh,"\n")),
          append=TRUE,file=outputFile)
      
      
    }
    save(GSEA_out_sig,file=paste0(saveDir,"/",experimentName,"_GSEA_out_allSigLists_forEachMetaCluster.RData.gzip"),compress="gzip")
    
    
    for(g in 1:length(GSEA_out_sig)){
      #what can't remove indices??
      write.table(as.character(unlist(GSEA_out_sig[[g]])),file=paste0(saveDir,"/",experimentName,"_GSEA_allSig_community_",g,".txt"),row.names=TRUE,quote=FALSE,col.names=FALSE)
      
    }
    #cat(GSEA_out_sig,sep="\n",file=paste0(saveDir,"/",experimentName,"_GSEA_allSig.txt"),append=FALSE)
    
  }else{
    
    GSEA_out_unique <- NA
    sigGenes <- NA
    
  }
  ###now on to plotting
   
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
  
  groupingTerm="community"
  
  #only take samples with the groupingTerm you're looking at.
  sampleClustCommPhenoData <- sampleClustCommPhenoData[which(!is.na(sampleClustCommPhenoData[, groupingTerm])), ]
  #remove samples with NA values.
  groupings <- as.factor(sampleClustCommPhenoData[, groupingTerm])
  
  cat(gsub(" ", "_",paste0("\n",length(groupings), " patients included across all final meta-clusters.\n")),
      append=TRUE,file=summaryFile)
  
  }else{
    
    sampleClustCommPhenoDataOrig <- NULL
  }
  
  expName <- gsub("_"," ",experimentName)
  
  if(clinVarPlots){
    
    if(!ovarian){
      #need to add pam50 centroids
      
      load("/home/ywrfc09/breast_analysis/PAM50_analyses/pam50FullAndShort_subtypeDF.RData.gzip")
      #hmmm...pam50 short seems to work better? perhaps because they all share the same genes?
      #careful: any duplicated sample names (across different studies)
      #any(duplicated(subtypeDF_master$sampleName))
      #any(duplicated(aggregateData$sampleClustCommKey$sampleName))
      colnames(subtypeDF_master)[2] <- "origStudyNum"
      #now: this data frame only has the "original" study number.
      #type=left: only want samples that are in the sampleClustCommPhenoData
      sampleClustCommPhenoData <- join(sampleClustCommPhenoData,subtypeDF_master,by=c("origStudyNum","sampleName"),type="left",
                                       match="all") 
      
      #add on subtype variables
      fisherTestVariables <- append(fisherTestVariables,c("subtype","subtype_short"))
      fisherTestVariableTitleNames <- append(fisherTestVariableTitleNames,c("pam50 subtype","35-gene pam50 subtype"))
      fisherTestVariableLegendNames <- append( fisherTestVariableLegendNames,c("pam50\nsubtype","pam50\nsubtype"))     
      
    }   
    
    
    
    for(f in 1:length(fisherTestVariables)){
      
      if(!fisherTestVariables[f] %in% colnames(sampleClustCommPhenoData)){
        #move on to next loop: don't stop code way at bottom.
        message("Clinical variable ", fisherTestVariables[f] , " not found in pheno data.")
        next;
      }
      
      #if any are "", make NA.
      sampleClustCommPhenoData[which(sampleClustCommPhenoData[,fisherTestVariables[f]]==""),
                               fisherTestVariables[f]] <- NA
      
      clinicalTables[[fisherTestVariables[f]]] <- table(sampleClustCommPhenoData[,groupingTerm],factor(sampleClustCommPhenoData[,fisherTestVariables[f]],exclude=NULL))
      #want NA as a factor, so set exclude=NULL, not NA.
      textOut <- capture.output(clinicalTables[[fisherTestVariables[f]]])
      cat(textOut,sep="\n",file=paste0(saveDir,"/",experimentName,"_tableStats_",fisherTestVariables[f],"_",Sys.Date(),".txt"),append=FALSE)
      
      cat(gsub(" ","_",paste0("\nTableStats for ",fisherTestVariables[f],":\n")),textOut,sep="\n",
          file=summaryFile,append=TRUE)
      
      dataF <- data.frame(sampleClustCommPhenoData[,groupingTerm],
                          sampleClustCommPhenoData[,fisherTestVariables[f]],
                          stringsAsFactors=FALSE)
      colnames(dataF) <- c("meta_cluster","clinVar")
      
      if(any(!is.na(dataF[,2]))){
        
        
        #pam50 colors:  t(data.frame(c("#970eec","#ec0e58","#0e58ec","#7d7d7d","#0eec97")))
        #(and name for subtypes-look up order in old clust robust code.)
        
        if(fisherTestVariables[f] != "subtype" &&  fisherTestVariables[f] != "subtype_short"){
          
          #no scale_fill_manual(values = variableColorMatrix)+
          #fill must be a factor.
          plotG <-    ggplot(dataF,aes(factor(clinVar),fill=factor(clinVar)))+geom_bar() + facet_grid(.~meta_cluster,scales="free_x")+
            labs(y="Number of samples",fill=paste0("",fisherTestVariableLegendNames[f],""),
                 title=paste0(fisherTestVariableTitleNames[f],
                              " by meta-cluster \nfor ",expName,"\n"))+
            theme(axis.text.x = element_text(colour = "black",size=8,angle=45,vjust=1,hjust=1),
                  axis.title=element_blank())+
            theme(legend.title=element_text(size=12))+   theme(strip.text.x = element_text(size = 17, colour = "black",face="bold"))+
            theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
            theme(axis.text.x = element_text(colour = "black",size=12,angle=45,vjust=1,hjust=1),axis.title.x= element_blank(),
                  axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1))+
            theme(panel.grid.major = element_line(colour = 0),panel.grid.minor = element_line(colour = 0))+
            theme(plot.title=element_text(colour="black",size=18,vjust=1))
          
          if(!plotToScreen){
            
            png(filename=paste0(saveDir,"/",experimentName,"_",fisherTestVariables[f],"_breakdowns_",Sys.Date(),".png"),width=1000,height=1000,res=160)
            
            plot(plotG)
            
            dev.off()
            
          }else{
            #plot at least the subtype figures to screen.
            plot(plotG)
          }
          
          if(is.na(plotStackedYLimit)){
            
            plotGStacked <-  ggplot(dataF,aes(factor(meta_cluster),fill=factor(clinVar)),scales="free_x")+geom_bar() + facet_grid(.~meta_cluster,scales="free_x")+
              labs(y="Number of samples", x=paste0("",fisherTestVariables[f],""),fill=paste0("",fisherTestVariableLegendNames[f],""),
                   title=paste0(fisherTestVariableTitleNames[f],
                                " by meta-cluster \nfor ",expName,"\n"))+
              theme(axis.text.x = element_text(colour = "black",size=8,angle=45,vjust=1,hjust=1),
                    axis.title=element_blank())+
              theme(legend.title=element_text(size=12))+  theme(strip.text.x = element_text(size = 17, colour = "black",face="bold"))+
              theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
              theme(axis.text.x = element_text(colour = "black",size=12,angle=45,vjust=1,hjust=1),axis.title.x= element_blank(),
                    axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1))+
              theme(panel.grid.major = element_line(colour = 0),panel.grid.minor = element_line(colour = 0))+
              theme(plot.title=element_text(colour="black",size=20,vjust=1))
            
          }else{
            
            plotGStacked <-  ggplot(dataF,aes(factor(meta_cluster),fill=factor(clinVar)),scales="free_x")+geom_bar() + facet_grid(.~meta_cluster,scales="free_x")+
              labs(y="Number of samples", x=paste0("",fisherTestVariables[f],""),fill=paste0("",fisherTestVariableLegendNames[f],""),
                   title=paste0(fisherTestVariableTitleNames[f],
                                " by meta-cluster \nfor ",expName,"\n"))+
              theme(axis.text.x = element_text(colour = "black",size=8,angle=45,vjust=1,hjust=1),
                    axis.title=element_blank())+
              theme(legend.title=element_text(size=12))+  theme(strip.text.x = element_text(size = 17, colour = "black",face="bold"))+
              theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
              theme(axis.text.x = element_text(colour = "black",size=12,angle=45,vjust=1,hjust=1),axis.title.x= element_blank(),
                    axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1))+
              theme(panel.grid.major = element_line(colour = 0),panel.grid.minor = element_line(colour = 0))+
              theme(plot.title=element_text(colour="black",size=20,vjust=1))+coord_cartesian(ylim=c(0,plotStackedYLimit))
            
            
          }
          
          if(!plotToScreen){
            
            png(filename=paste0(saveDir,"/",experimentName,"_",fisherTestVariables[f],"_breakdowns_stacked_",Sys.Date(),".png"),width=1000,height=1000,res=160)
            
            plot( plotGStacked)
            
            dev.off()
            
          }else{
            
            plot(plotGStacked)
          }
          
          
        }else{
          
          variableColorMatrix <- c("#970eec","#ec0e58","#0e58ec","#7d7d7d","#0eec97")
          names(variableColorMatrix) <- c("LumB","LumA","Her2","Normal","Basal")
          
          plotG <-    ggplot(dataF,aes(factor(clinVar),fill=factor(clinVar)))+geom_bar() + facet_grid(.~meta_cluster,scales="free_x")+
            labs(y="Number of samples", fill=paste0("",fisherTestVariableLegendNames[f],""),
                 title=paste0(fisherTestVariableTitleNames[f],
                              " by meta-cluster \nfor ",expName))+
            theme(axis.text.x = element_text(colour = "black",size=8,angle=45,vjust=1,hjust=1),
                  axis.title=element_blank())+
            theme(legend.title=element_text(size=12))+ scale_fill_manual(values = variableColorMatrix)+  theme(strip.text.x = element_text(size = 17, colour = "black",face="bold"))+
            theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
            theme(axis.text.x = element_text(colour = "black",size=12,angle=45,vjust=1,hjust=1),axis.title.x= element_blank(),
                  axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1))+
            theme(panel.grid.major = element_line(colour = 0),panel.grid.minor = element_line(colour = 0))+
            theme(plot.title=element_text(colour="black",size=20,vjust=1))+coord_cartesian(ylim=c(0,1000))
          #+coord_cartesian(ylim=c(0,900)) #if want axes to all be the same.
          #above: font size needs to be smaller for subtypes
          
          if(!plotToScreen){
            
            png(filename=paste0(saveDir,"/",experimentName,"_",fisherTestVariables[f],"_breakdowns_",Sys.Date(),".png"),width=1000,height=1000,res=160)
            
            plot(plotG)
            
            dev.off()
            
          }else{
            
            plot(plotG)
            
          }
          
          if(is.na(plotStackedYLimit)){
            
            plotGStacked <- ggplot(dataF,aes(factor(meta_cluster),fill=factor(clinVar)),scales="free_x")+geom_bar() + facet_grid(.~meta_cluster,scales="free_x")+
              labs(y="Number of samples", fill=paste0("",fisherTestVariableLegendNames[f],""),
                   title=paste0(fisherTestVariableTitleNames[f],
                                " by meta-cluster \nfor ",expName,"\n"))+
              theme(axis.text.x = element_text(colour = "black",size=8,angle=45,vjust=1,hjust=1),
                    axis.title=element_blank())+
              theme(legend.title=element_text(size=12))+ scale_fill_manual(values = variableColorMatrix)+  theme(strip.text.x = element_text(size = 17, colour = "black",face="bold"))+
              theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
              theme(axis.text.x = element_text(colour = "black",size=12,angle=45,vjust=1,hjust=1),axis.title.x= element_blank(),
                    axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1))+
              theme(panel.grid.major = element_line(colour = 0),panel.grid.minor = element_line(colour = 0))+
              theme(plot.title=element_text(colour="black",size=20,vjust=1))
            
          }else{
            plotGStacked <- ggplot(dataF,aes(factor(meta_cluster),fill=factor(clinVar)),scales="free_x")+geom_bar() + facet_grid(.~meta_cluster,scales="free_x")+
              labs(y="Number of samples", fill=paste0("",fisherTestVariableLegendNames[f],""),
                   title=paste0(fisherTestVariableTitleNames[f],
                                " by meta-cluster \nfor ",expName,"\n"))+
              theme(axis.text.x = element_text(colour = "black",size=8,angle=45,vjust=1,hjust=1),
                    axis.title=element_blank())+
              theme(legend.title=element_text(size=12))+ scale_fill_manual(values = variableColorMatrix)+  theme(strip.text.x = element_text(size = 17, colour = "black",face="bold"))+
              theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
              theme(axis.text.x = element_text(colour = "black",size=12,angle=45,vjust=1,hjust=1),axis.title.x= element_blank(),
                    axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1))+
              theme(panel.grid.major = element_line(colour = 0),panel.grid.minor = element_line(colour = 0))+
              theme(plot.title=element_text(colour="black",size=20,vjust=1))+coord_cartesian(ylim=c(0,plotStackedYLimit))
            
            
            
          }
          
          if(!plotToScreen){
            
            png(filename=paste0(saveDir,"/",experimentName,"_",fisherTestVariables[f],"_breakdowns_stacked_",Sys.Date(),".png"),width=1000,height=1000,res=160)
            
            plot( plotGStacked)
            
            dev.off()
            
          }else{
            
            plot(plotGStacked)
            
          }
          
        }
        
        #chi-squared test. fisher's test doesn't work here...it looks like bar charts, plots will be more informative.
        tmp <- factor(sampleClustCommPhenoData[,groupingTerm])[which(!is.na(sampleClustCommPhenoData[,fisherTestVariables[f]]))]
        tmp1 <- factor(sampleClustCommPhenoData[,fisherTestVariables[f]][which(!is.na(sampleClustCommPhenoData[,fisherTestVariables[f]]))])
        if(length(tmp1)==0 || length(unique(tmp1)==1)){
          
          cat(gsub(" ","_",paste0("\nClinical variable ", fisherTestVariables[f], " had only NAs or all one value so not analyzing it.\n")),
              append=TRUE,file=outputFile)
          next;
        }
        #need to simulate p-value other gets wonky.
        result[[paste0("fisher_",fisherTestVariables[f],"_pvalue")]] <- chisq.test(tmp,tmp1,simulate.p.value=TRUE)$p.value
        
        
        #end of if is NA
      }else{
        cat(gsub(" ","_",paste0("\nVariable ",fisherTestVariables[f], " is all NAs.\n")),append=TRUE,file=outputFile)
        result[[paste0("fisher_",fisherTestVariables[f],"_pvalue")]]  <- NA
      }
      #end of looping over fisher variables.
    }
    
  }
  
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
     breast_survivalData <- NA
    }else{
      #breast analysis
      #message("running survival data on breast")
      message("Running RFS models")
#       cat(gsub(" ","_",paste0("\nLinear analysis with only RFS:\n")),append=TRUE,file=summaryFile)
#       
#       dataMatrix <- sampleClustCommPhenoData[which(!is.na((sampleClustCommPhenoData$RFS))), ]
#       cat(gsub(" ","_",paste0("Studies used:\n",unique(dataMatrix$studyNum),"\n")),append=TRUE,file=summaryFile)
#       cat(gsub(" ","_",paste0("\nMeta-clusters used:\n",unique(dataMatrix$community),"\n")),append=TRUE,file=summaryFile)
#       
#       dataMatrix <- data.frame(dataMatrix)
#       dataMatrix$community <- as.factor(dataMatrix$community)
#       dataMatrix$RFS <- as.factor(dataMatrix$RFS)
#       
#       linearModels[["RFS_logitAlone"]] <- tryCatch(glm(RFS~community,data=dataMatrix,family=binomial(link="logit")),
#                                                    error = function(e) {
#                                                      return(NA)
#                                                    }
#       )
#       
#       if(any(!is.na(  linearModels[["RFS_logitAlone"]]))){
#         ##ROC curves
#         #above: font size needs to be smaller for subtypes
#         predictor <- predict(linearModels[["RFS_logitAlone"]],type="response")
#         
#         png(filename=paste0(saveDir,"/",experimentName,"_RFS_ROC_",Sys.Date(),".png"),width=1000,height=1000,res=200)
#         roc_data <- roc(response=dataMatrix$RFS,predictor=predictor,plot=TRUE,na.rm=TRUE,auc=TRUE)
#         dev.off()
#         textOut <- capture.output(roc_data)
#         cat("roc_info\n",textOut,sep="\n",append=TRUE,file=summaryFile)
#         cat("\nAUC : ",roc_data$auc,"\n",append=TRUE,file=summaryFile)
#         
#         linearModelsSumm[["RFS_logitAlone"]] <- summary(linearModels[["RFS_logitAlone"]] ) 
#         textOut <- capture.output( linearModelsSumm [["RFS_logitAlone"]]$coefficients)
#         cat(textOut,sep="\n",append=TRUE,file=summaryFile)
#         textOut <- capture.output( linearModelsSumm [["RFS_logitAlone"]])
#         cat(textOut,sep="\n",append=TRUE,file=summaryFile)
#         linearModelsSumm[["RFS_numPatients"]] <- nrow(dataMatrix)
#         cat("\nNumber_of_patients_for_RFS__model: ",nrow(dataMatrix),"\n",append=TRUE,file=summaryFile)
#         linearModelsSumm[["RFS_aloneMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
#         
#       }else{
#         cat("\nRFS model alone returned NA.\n")
#         cat("\nRFS model alone returned NA.\n",append=TRUE,file=outputFile)
#       }
#       ##with therapies
#       #interesting: community becomes insignificant if add in therapies now: this could be
#       #because community is most likely correlated with treatments; it's certainly correlated with pam50 status.
#       #low Rsquared.
#       cat(gsub(" ","_",paste0("\nLinear analysis with RFS and treatment:\n")),append=TRUE,file=summaryFile)
#       dataMatrix <- dataMatrix[which(!is.na(dataMatrix$chemotherapyClass)), ]
#       dataMatrix <- dataMatrix[which(!is.na(dataMatrix$anti_estrogen)), ]
#       dataMatrix <- dataMatrix[which(!is.na(dataMatrix$anti_HER2)), ]
#       dataMatrix <- data.frame(dataMatrix)
#       dataMatrix$community <- as.factor(dataMatrix$community)
#       dataMatrix$RFS <- as.factor(dataMatrix$RFS)
#       dataMatrix$hist_grade <- as.factor(dataMatrix$hist_grade)
#       dataMatrix$chemotherapyClass <- as.factor(dataMatrix$chemotherapyClass)
#       dataMatrix$anti_estrogen <- as.factor(dataMatrix$anti_estrogen)
#       dataMatrix$anti_HER2 <- as.factor(dataMatrix$anti_HER2)
#       cat(gsub(" ","_",paste0("\nMeta-clusters used:\n",unique(dataMatrix$community),"\n")),append=TRUE,file=summaryFile)
#       
#       linearModelsSumm [["RFS_withRx_numPatients"]] <- nrow(dataMatrix)
#       cat("\nNumber_of_patients_for_RFS_with_Rx_model: ",nrow(dataMatrix),"\n",append=TRUE,file=summaryFile)
#       
#       if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))>1 
#          && length(unique(dataMatrix$anti_HER2))>1 ){
#         
#         linearModels[["RFS_logitWithRx"]] <- tryCatch(glm(RFS~community+chemotherapyClass+anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
#                                                       error = function(e) {
#                                                         return(NA)
#                                                       }
#         )
#         
#         #no community
#         reducedModel <- tryCatch(glm(RFS~chemotherapyClass+anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
#                                                       error = function(e) {
#                                                         return(NA)
#                                                       }
#         )
#          linearModels[["RFS_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["RFS_logitWithRx"]],test="Chisq")
#            
#      
#         #now try with  grade.    
#         cat(gsub(" ","_",paste0("\nLinear analysis with RFS and histological grade :\n")),append=TRUE,file=summaryFile)
#         dataMatrixHist <- dataMatrix[which(!is.na((dataMatrix$hist_grade))),]
#  
#         
#         cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
#         cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
#         
#         linearModelsSumm[["RFS_withGrade_numPatients"]] <- nrow(dataMatrixHist)
#         cat("\nNumber_of_patients_for_RFS_with_grade_model: ",nrow(dataMatrixHist),"\n",append=TRUE,file=summaryFile)
#         
#         linearModels[["RFS_logitWithRxGrade"]] <- tryCatch( glm(RFS~community+chemotherapyClass+anti_estrogen+anti_HER2+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
#                                                             error = function(e) {
#                                                               return(NA)
#                                                             }
#         )
#         
#         
#       }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))>1 
#                && length(unique(dataMatrix$anti_HER2))==1 ){
#         
#         
#         cat(gsub(" ","_",paste0("\nFor RFS analysis, all anti_HER2 values the same: ",unique(unique(dataMatrix$anti_HER2)),"\n")),
#             append=TRUE,file=summaryFile)
#         linearModels[["RFS_logitWithRx"]] <- tryCatch(glm(RFS~community+chemotherapyClass+anti_estrogen,data=dataMatrix,family=binomial(link="logit")),
#                                                       error = function(e) {
#                                                         return(NA)
#                                                       }
#                                                       
#         )
#         
#                #no community
#         reducedModel <- tryCatch(glm(RFS~chemotherapyClass+anti_estrogen,data=dataMatrix,family=binomial(link="logit")),
#                                                       error = function(e) {
#                                                         return(NA)
#                                                       }
#         )
#          linearModels[["RFS_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["RFS_logitWithRx"]],test="Chisq")
#         
#         #now try with  grade.    
#         cat(gsub(" ","_",paste0("\nLinear analysis with RFS and histological grade :\n")),append=TRUE,file=summaryFile)
#         dataMatrixHist <- dataMatrix[which(!is.na((dataMatrix$hist_grade))),]
#  
#         
#         cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
#         cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
#         
#         linearModelsSumm[["RFS_withGrade_numPatients"]] <- nrow(dataMatrixHist)
#         cat("\nNumber_of_patients_for_RFS_with_grade_model: ",nrow(dataMatrixHist),"\n",append=TRUE,file=summaryFile)
#         linearModels[["RFS_logitWithRxGrade"]] <- tryCatch( glm(RFS~community+chemotherapyClass+anti_estrogen+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
#                                                             error = function(e) {
#                                                               return(NA)
#                                                             }
#         )
#         
#         
#         
#       }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))==1 
#                && length(unique(dataMatrix$anti_HER2))>1 ){
#         
#         cat(gsub(" ","_",paste0("\nFor RFS analysis, all anti_estrogen values the same: ",unique(unique(dataMatrix$anti_estrogen)),"\n")),
#             append=TRUE,file=summaryFile)
#         
#         linearModels[["RFS_logitWithRx"]] <- tryCatch(glm(RFS~community+chemotherapyClass+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
#                                                       error = function(e) {
#                                                         return(NA)
#                                                       }
#                                                       
#         )
#         
#                    #no community
#         reducedModel <- tryCatch(glm(RFS~chemotherapyClass+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
#                                                       error = function(e) {
#                                                         return(NA)
#                                                       }
#         )
#          linearModels[["RFS_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["RFS_logitWithRx"]],test="Chisq")
#         
#  
#         
#         #now try with  grade.    
#         cat(gsub(" ","_",paste0("\nLinear analysis with RFS and histological grade :\n")),append=TRUE,file=summaryFile)
#         dataMatrixHist <- dataMatrix[which(!is.na((dataMatrix$hist_grade))),]
#  
#         
#         cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
#         cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
#         
#         linearModelsSumm[["RFS_withGrade_numPatients"]] <- nrow(dataMatrixHist)
#         cat("\nNumber_of_patients_for_RFS_with_grade_model: ",nrow(dataMatrixHist),"\n",append=TRUE,file=summaryFile)
#         linearModels[["RFS_logitWithRxGrade"]] <- tryCatch( glm(RFS~community+chemotherapyClass+anti_HER2+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
#                                                             error = function(e) {
#                                                               return(NA)
#                                                             }
#         )
#         
#         
#       }else if(length(unique(dataMatrix$chemotherapyClass))==1 && length(unique(dataMatrix$anti_estrogen))>1 
#                && length(unique(dataMatrix$anti_HER2))>1 ){
#         
#         
#         cat(gsub(" ","_",paste0("\nFor RFS analysis, all chemo values the same: ",unique(unique(dataMatrix$chemotherapyClass)),"\n")),
#             append=TRUE,file=summaryFile)
#         
#         linearModels[["RFS_logitWithRx"]] <- tryCatch(glm(RFS~community+anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
#                                                       error = function(e) {
#                                                         return(NA)
#                                                       }
#                                                       
#         )
#         
#                            #no community
#         reducedModel <- tryCatch(glm(RFS~anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
#                                                       error = function(e) {
#                                                         return(NA)
#                                                       }
#         )
#          linearModels[["RFS_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["RFS_logitWithRx"]],test="Chisq")
#         #now try with  grade.    
#         cat(gsub(" ","_",paste0("\nLinear analysis with RFS and histological grade :\n")),append=TRUE,file=summaryFile)
#         dataMatrixHist <- dataMatrix[which(!is.na((dataMatrix$hist_grade))),]
#  
#         
#         cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
#         cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
#         
#         linearModelsSumm[["RFS_withGrade_numPatients"]] <- nrow(dataMatrixHist)
#         cat("\nNumber_of_patients_for_RFS_with_grade_model: ",nrow(dataMatrixHist),"\n",append=TRUE,file=summaryFile)
#         linearModels[["RFS_logitWithRxGrade"]] <- tryCatch( glm(RFS~community+anti_estrogen+anti_HER2+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
#                                                             error = function(e) {
#                                                               return(NA)
#                                                             }
#         )
#         
#         
#         
#       }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))==1 
#                && length(unique(dataMatrix$anti_HER2))==1 ){
#         
#         
#         cat(gsub(" ","_",paste0("\nFor RFS analysis, the anti HER2 and estrogen variables have only one unique value each.\n")),append=TRUE,
#             file=summaryFile)
#         
#         
#         linearModels[["RFS_logitWithRx"]] <- tryCatch(glm(RFS~community+chemotherapyClass,data=dataMatrix,family=binomial(link="logit")),
#                                                       error = function(e) {
#                                                         return(NA)
#                                                       }
#                                                       
#         )
#          reducedModel <- tryCatch(glm(RFS~chemotherapyClass,data=dataMatrix,family=binomial(link="logit")),
#                                                       error = function(e) {
#                                                         return(NA)
#                                                       }
#         )
#          linearModels[["RFS_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["RFS_logitWithRx"]],test="Chisq")
#     
#         
#         #now try with  grade.    
#         cat(gsub(" ","_",paste0("\nLinear analysis with RFS and histological grade :\n")),append=TRUE,file=summaryFile)
#         dataMatrixHist <- dataMatrix[which(!is.na((dataMatrix$hist_grade))),]
#  
#         
#         cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
#         cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
#         
#         linearModelsSumm[["RFS_withGrade_numPatients"]] <- nrow(dataMatrixHist)
#         cat("\nNumber_of_patients_for_RFS_with_grade_model: ",nrow(dataMatrixHist),"\n",append=TRUE,file=summaryFile)
#         linearModels[["RFS_logitWithRxGrade"]] <- tryCatch( glm(RFS~community+chemotherapyClass+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
#                                                             error = function(e) {
#                                                               return(NA)
#                                                             }
#         )
#         
#         
#         
#         
#       }else{
#         
#         linearModels[["RFS_logitWithRx"]] <- NA
#         
#       }
#       
#       cat(gsub(" ","_",paste0("\nLinear analysis with RFS and treatment:\n")),append=TRUE,file=summaryFile)
#       if(any(!is.na(linearModels[["RFS_logitWithRx"]]))){
#         
#         linearModelsSumm[["RFS_logitWithRx"]] <- summary( linearModels[["RFS_logitWithRx"]])
#         textOut <- capture.output( linearModelsSumm[["RFS_logitWithRx"]]$coefficients)
#         cat(textOut,sep="\n",append=TRUE,file=summaryFile)
#         textOut <- capture.output( linearModelsSumm [["RFS_logitWithRx"]])
#         cat(textOut,sep="\n",append=TRUE,file=summaryFile)
#         cat(gsub(" ","_",paste0("\nLinear analysis with RFS and treatment ANOVA for meta-clusters:\n")),append=TRUE,file=summaryFile)
#         textOut <- capture.output(linearModels[["RFS_logitWithRx_ANOVA"]])
#         cat(textOut,sep="\n",append=TRUE,file=summaryFile)
#         
#         linearModelsSumm[["RFS_WithRxMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",") 
#         predictor <- predict(linearModels[["RFS_logitWithRx"]],type="response")
#         png(filename=paste0(saveDir,"/",experimentName,"_RFS_withRx_ROC_",Sys.Date(),".png"),width=1000,height=1000,res=200)
#         roc_data <- roc(response=dataMatrix$RFS,predictor=predictor,plot=TRUE,na.rm=TRUE,auc=TRUE)
#         dev.off()
#         textOut <- capture.output(roc_data)
#         cat("\n roc_info:",textOut,sep="\n",append=TRUE,file=summaryFile)
#         cat("\nAUC : ",roc_data$auc,"\n",append=TRUE,file=summaryFile)
#         
#         linearModelsSumm[["RFS_logitWithRx"]] <- summary( linearModels[["RFS_logitWithRx"]])
#         textOut <- capture.output( linearModelsSumm[["RFS_logitWithRx"]]$coefficients)
#         cat(textOut,sep="\n",append=TRUE,file=summaryFile)
#         
#         linearModelsSumm[["RFS_WithRxMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
#         
#       }else{
#         
#         cat(gsub(" ","_",paste0("\nRFS analysis with therapies returned NA.\n")),append=TRUE,file=summaryFile)
#       }
#       
#       #now record grade.
#       cat(gsub(" ","_",paste0("\nLinear analysis with RFS and histological grade :\n")),append=TRUE,file=summaryFile)
#       cat(gsub(" ","_",paste0("\nUnique hist grade values :",unique(dataMatrixHist$hist_grade),"\n")),append=TRUE,file=summaryFile)
#       
#       
#       if(any(!is.na(linearModels[["RFS_logitWithRxGrade"]]))){
#         
#         predictor <- predict(linearModels[["RFS_logitWithRxGrade"]],type="response")
#         png(filename=paste0(saveDir,"/",experimentName,"_RFS_withGrade_ROC_",Sys.Date(),".png"),width=1000,height=1000,res=200)
#         roc_data <- roc(response=dataMatrixHist$RFS,predictor=predictor,plot=TRUE,na.rm=TRUE,auc=TRUE)
#         dev.off()
#         textOut <- capture.output(roc_data)
#         cat("roc_info:\n",textOut,sep="\n",append=TRUE,file=summaryFile)
#         cat("\nAUC: ",roc_data$auc,"\n",append=TRUE,file=summaryFile)
#         
#         linearModelsSumm[["RFS_logitWithRxGrade"]] <- summary( linearModels[["RFS_logitWithRxGrade"]])
#         textOut <- capture.output( linearModelsSumm[["RFS_logitWithRxGrade"]]$coefficients)
#         cat(textOut,sep="\n",append=TRUE,file=summaryFile) 
#         linearModelsSumm[["RFS_WithGradeMetaClusters_used"]] <- paste(unique(dataMatrixHist$community),collapse=",")
#         
#       }else{
#         
#         cat(gsub(" ","_",paste0("\nHist grade with RFS analysis returned NA.\n")),append=TRUE,file=summaryFile)  
#       }
#       
#       
#       
#       message("Running DFS models")
#       ###DFS
#       cat(gsub(" ","_",paste0("\nLinear analysis with only DFS:\n")),append=TRUE,file=summaryFile)
#       
#       dataMatrix <- sampleClustCommPhenoData[which(!is.na((sampleClustCommPhenoData$DFS))), ]
#       dataMatrix$community <- as.factor(dataMatrix$community)
#       dataMatrix$DFS <- as.factor(dataMatrix$DFS)
#       dataMatrix$hist_grade <- as.factor(dataMatrix$hist_grade)
#       dataMatrix$chemotherapyClass <- as.factor(dataMatrix$chemotherapyClass)
#       dataMatrix$anti_estrogen <- as.factor(dataMatrix$anti_estrogen)
#       dataMatrix$anti_HER2 <- as.factor(dataMatrix$anti_HER2)
#       cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrix$studyNum),"\n")),append=TRUE,file=summaryFile)
#       cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrix$community),"\n")),append=TRUE,file=summaryFile)
#       
#       #strongly predicts DFS.
#       linearModels[["DFS_logitAlone"]] <- glm(DFS~community,data=dataMatrix,family=binomial(link="logit"))
#       png(filename=paste0(saveDir,"/",experimentName,"_DFS_ROC_",Sys.Date(),".png"),width=1000,height=1000,res=200)
#       predictor <- predict(linearModels[["DFS_logitAlone"]],type="response")
#       roc_data <- roc(response=dataMatrix$DFS,predictor=predictor,plot=TRUE,na.rm=TRUE,auc=TRUE)
#       dev.off()
#       textOut <- capture.output(roc_data)
#       cat("roc_info\n",textOut,sep="\n",append=TRUE,file=summaryFile)
#       cat("\nAUC: ",roc_data$auc,"\n",append=TRUE,file=summaryFile)
#       linearModelsSumm[["DFS_logitAlone"]] <- summary( linearModels[["DFS_logitAlone"]] ) 
#       textOut <- capture.output( linearModelsSumm[["DFS_logitAlone"]]$coefficients)
#       cat(textOut,sep="\n",append=TRUE,file=summaryFile)
#       textOut <- capture.output( linearModelsSumm[["DFS_logitAlone"]])
#       cat(textOut,sep="\n",append=TRUE,file=summaryFile)
#       linearModelsSumm[["DFS_numPatients"]] <- nrow(dataMatrix)
#       cat("\nNumber of patients: ",nrow(dataMatrix),"\n",append=TRUE,file=summaryFile)
#       linearModelsSumm[["DFS_aloneMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
#       
#       ##with therapies
#       #interesting: community becomes insignificant if add in therapies now: this could be
#       #because community is most likely correlated with treatments; it's certainly correlated with pam50 status.
#       #low Rsquared.
#       cat(gsub(" ","_",paste0("\nLinear analysis with DFS and treatment:\n")),append=TRUE,file=summaryFile)
#       dataMatrix <- dataMatrix[which(!is.na(dataMatrix$chemotherapyClass)), ]
#       dataMatrix <- dataMatrix[which(!is.na(dataMatrix$anti_estrogen)), ]
#       dataMatrix <- dataMatrix[which(!is.na(dataMatrix$anti_HER2)), ]
#       dataMatrix <- data.frame(dataMatrix)
#       dataMatrix$community <- as.factor(dataMatrix$community)
#       dataMatrix$DFS <- as.factor(dataMatrix$DFS)
#       dataMatrix$hist_grade <- as.factor(dataMatrix$hist_grade)
#       dataMatrix$chemotherapyClass <- as.factor(dataMatrix$chemotherapyClass)
#       dataMatrix$anti_estrogen <- as.factor(dataMatrix$anti_estrogen)
#       dataMatrix$anti_HER2 <- as.factor(dataMatrix$anti_HER2)
#       cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrix$studyNum),"\n")),append=TRUE,file=summaryFile)
#       cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrix$community),"\n")),append=TRUE,file=summaryFile)
#       
#       linearModelsSumm[["DFS_withRx_numPatients"]] <- nrow(dataMatrix)
#       cat("\nNumber_of_patients_for_DFS_with_Rx_model: ",nrow(dataMatrix),"\n",append=TRUE,file=summaryFile)
#       
#       if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))>1 && 
#          length(unique(dataMatrix$anti_HER2))>1 ){
#         
#         linearModels[["DFS_logitWithRx"]] <- tryCatch(glm(DFS~community+chemotherapyClass+anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
#                                                       error = function(e) {
#                                                         return(NA)
#                                                       }
#         )
#         
#         #remove communities for ANOVA
#            reducedModel <- tryCatch(glm(DFS~chemotherapyClass+anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
#                                                       error = function(e) {
#                                                         return(NA)
#                                                       }
#         )
#          linearModels[["DFS_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["DFS_logitWithRx"]],test="Chisq")
#     
#     
#         #now try with grade.
#         cat(gsub(" ","_",paste0("\nLinear analysis with DFS and histological grade :\n")),append=TRUE,file=summaryFile)
#         dataMatrixHist <- dataMatrix[which(!is.na((dataMatrix$hist_grade))),]
#  
# 
#         cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
#         cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
#         
#         linearModelsSumm[["DFS_withGrade_numPatients"]] <- nrow(dataMatrixHist)
#         cat("\nNumber_of_patients_for_DFS_withGrade_model: ",nrow(dataMatrixHist),"\n",append=TRUE,file=summaryFile)
#         linearModels[["DFS_logitWithRxGrade"]] <- tryCatch( glm(DFS~community+chemotherapyClass+anti_estrogen+anti_HER2+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
#                                                             error = function(e) {
#                                                               return(NA)
#                                                             }
#         )
#         
#         
#         
#       }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))>1 &&
#                length(unique(dataMatrix$anti_HER2))==1 ){
#         
#         
#         cat(gsub(" ","_",paste0("\nFor DFS analysis, all anti_HER2 values the same: ",unique(unique(dataMatrix$anti_HER2)),"\n")),
#             append=TRUE,file=summaryFile)
#         linearModels[["DFS_logitWithRx"]] <- tryCatch(glm(DFS~community+chemotherapyClass+anti_estrogen,data=dataMatrix,family=binomial(link="logit")),
#                                                       error = function(e) {
#                                                         return(NA)
#                                                       }
#         )
#                 #remove communities for ANOVA
#            reducedModel <- tryCatch(glm(DFS~chemotherapyClass+anti_estrogen,data=dataMatrix,family=binomial(link="logit")),
#                                                       error = function(e) {
#                                                         return(NA)
#                                                       }
#         )
#          linearModels[["DFS_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["DFS_logitWithRx"]],test="Chisq")
#     
# 
#         #now try with grade.
#         cat(gsub(" ","_",paste0("\nLinear analysis with DFS and histological grade :\n")),append=TRUE,file=summaryFile)
#         dataMatrixHist <- dataMatrix[which(!is.na((dataMatrix$hist_grade))),]
#  
#         
#         cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
#         cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
#         
#         linearModelsSumm[["DFS_withGrade_numPatients"]] <- nrow(dataMatrixHist)
#         cat("\nNumber_of_patients_for_DFS_withGrade_model: ",nrow(dataMatrixHist),"\n",append=TRUE,file=summaryFile)
#         linearModels[["DFS_logitWithRxGrade"]] <- tryCatch( glm(DFS~community+chemotherapyClass+anti_estrogen+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
#                                                             error = function(e) {
#                                                               return(NA)
#                                                             }
#         )
#         
#         
#         
#         
#       }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))==1 
#                && length(unique(dataMatrix$anti_HER2))>1 ){
#         
#         cat(gsub(" ","_",paste0("\nFor DFS analysis, all anti_estrogen values the same: ",unique(unique(dataMatrix$anti_estrogen)),"\n")),
#             append=TRUE,file=summaryFile)
#         
#         linearModels[["DFS_logitWithRx"]] <- tryCatch(glm(DFS~community+chemotherapyClass+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
#                                                       error = function(e) {
#                                                         return(NA)
#                                                       }
#         )
#                     #remove communities for ANOVA
#            reducedModel <- tryCatch(glm(DFS~chemotherapyClass+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
#                                                       error = function(e) {
#                                                         return(NA)
#                                                       }
#         )
#          linearModels[["DFS_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["DFS_logitWithRx"]],test="Chisq")
#     
# 
#         #now try with grade.
#         cat(gsub(" ","_",paste0("\nLinear analysis with DFS and histological grade :\n")),append=TRUE,file=summaryFile)
#         dataMatrixHist <- dataMatrix[which(!is.na((dataMatrix$hist_grade))),]
#  
#         
#         cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
#         cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
#         
#         linearModelsSumm[["DFS_withGrade_numPatients"]] <- nrow(dataMatrixHist)
#         cat("\nNumber_of_patients_for_DFS_withGrade_model: ",nrow(dataMatrixHist),"\n",append=TRUE,file=summaryFile)
#         linearModels[["DFS_logitWithRxGrade"]] <- tryCatch( glm(DFS~community+chemotherapyClass+anti_HER2+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
#                                                             error = function(e) {
#                                                               return(NA)
#                                                             }
#         )
#         
#         
#         
#       }else if(length(unique(dataMatrix$chemotherapyClass))==1 && length(unique(dataMatrix$anti_estrogen))>1 
#                && length(unique(dataMatrix$anti_HER2))>1 ){
#         
#         
#         cat(gsub(" ","_",paste0("\nFor DFS analysis, all chemo values the same: ",unique(unique(dataMatrix$chemotherapyClass)),"\n")),
#             append=TRUE,file=summaryFile)
#         
#         linearModels[["DFS_logitWithRx"]] <- tryCatch(glm(DFS~community+anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
#                                                       error = function(e) {
#                                                         return(NA)
#                                                       }
#         )
#                           #remove communities for ANOVA
#            reducedModel <- tryCatch(glm(DFS~anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
#                                                       error = function(e) {
#                                                         return(NA)
#                                                       }
#         )
#          linearModels[["DFS_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["DFS_logitWithRx"]],test="Chisq")
#     
# 
#         #now try with grade.
#         cat(gsub(" ","_",paste0("\nLinear analysis with DFS and histological grade :\n")),append=TRUE,file=summaryFile)
#         dataMatrixHist <- dataMatrix[which(!is.na((dataMatrix$hist_grade))),]
#  
#         
#         cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
#         cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
#         
#         linearModelsSumm[["DFS_withGrade_numPatients"]] <- nrow(dataMatrixHist)
#         cat("\nNumber_of_patients_for_DFS_withGrade_model: ",nrow(dataMatrixHist),"\n",append=TRUE,file=summaryFile)
#         linearModels[["DFS_logitWithRxGrade"]] <- tryCatch( glm(DFS~community+anti_estrogen+anti_HER2+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
#                                                             error = function(e) {
#                                                               return(NA)
#                                                             }
#         )
#         
#         
#         
#       }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))==1 
#                && length(unique(dataMatrix$anti_HER2))==1 ){
#         
#         
#         cat(gsub(" ","_",paste0("\nFor DFS analysis, the anti HER2 and estrogen variables have only one unique value each.\n")),append=TRUE,
#             file=summaryFile)
#         
#         
#         linearModels[["DFS_logitWithRx"]] <- tryCatch(glm(DFS~community+chemotherapyClass,data=dataMatrix,family=binomial(link="logit")),
#                                                       error = function(e) {
#                                                         return(NA)
#                                                       }
#         )
#                                #remove communities for ANOVA
#            reducedModel <- tryCatch(glm(DFS~chemotherapyClass,data=dataMatrix,family=binomial(link="logit")),
#                                                       error = function(e) {
#                                                         return(NA)
#                                                       }
#         )
#          linearModels[["DFS_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["DFS_logitWithRx"]],test="Chisq")
#     
# 
#         #now try with grade.
#         cat(gsub(" ","_",paste0("\nLinear analysis with DFS and histological grade :\n")),append=TRUE,file=summaryFile)
#         dataMatrixHist <- dataMatrix[which(!is.na((dataMatrix$hist_grade))),]
#  
#         
#         cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
#         cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
#         
#         linearModelsSumm[["DFS_withGrade_numPatients"]] <- nrow(dataMatrixHist)
#         cat("\nNumber_of_patients_for_DFS_withGrade_model: ",nrow(dataMatrixHist),"\n",append=TRUE,file=summaryFile)
#         linearModels[["DFS_logitWithRxGrade"]] <- tryCatch( glm(DFS~community+chemotherapyClass+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
#                                                             error = function(e) {
#                                                               return(NA)
#                                                             }
#         )
#         
#         
#         
#       }else{
#         
#         # don't run models
#         
#       }
#       
#       cat(gsub(" ","_",paste0("\nLinear analysis with DFS and treatment:\n")),append=TRUE,file=summaryFile)
#       if(any(!is.na(linearModels[["DFS_logitWithRx"]]))){
#         
#         predictor <- predict(linearModels[["DFS_logitWithRx"]],type="response")
#         png(filename=paste0(saveDir,"/",experimentName,"_DFS_withRx_ROC_",Sys.Date(),".png"),width=1000,height=1000,res=200)
#         roc_data <- roc(response=dataMatrix$DFS,predictor=predictor,plot=TRUE,na.rm=TRUE,auc=TRUE)
#         dev.off()
#         textOut <- capture.output(roc_data)
#         cat("roc_info\n",textOut,sep="\n",append=TRUE,file=summaryFile)
#         cat("\nAUC: ",roc_data$auc,"\n",append=TRUE,file=summaryFile)
#         
#         linearModelsSumm[["DFS_logitWithRx"]] <- summary( linearModels[["DFS_logitWithRx"]])
#         textOut <- capture.output( linearModelsSumm[["DFS_logitWithRx"]]$coefficients)
#         cat(textOut,sep="\n",append=TRUE,file=summaryFile)
#                 textOut <- capture.output( linearModelsSumm[["DFS_logitWithRx"]])
#         cat(textOut,sep="\n",append=TRUE,file=summaryFile)
#         cat(gsub(" ","_",paste0("\nLinear analysis with DFS and treatment ANOVA for meta-clusters:\n")),append=TRUE,file=summaryFile)
#         textOut <- capture.output(linearModels[["DFS_logitWithRx_ANOVA"]])
#         cat(textOut,sep="\n",append=TRUE,file=summaryFile)
#         linearModelsSumm[["DFS_WithRxMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
#         
#       }else{
#         
#         cat(gsub(" ","_",paste0("\nRx with DFS analysis returned NA.\n")),append=TRUE,file=summaryFile)  
#       }
#       
#       
#       #now try with grade.
#       cat(gsub(" ","_",paste0("\nLinear analysis with DFS and histological grade :\n")),append=TRUE,file=summaryFile)
#       cat(gsub(" ","_",paste0("\nUnique hist grade values :",unique(dataMatrixHist$hist_grade),"\n")),append=TRUE,file=summaryFile)
#       
#       
#       if(!is.na(linearModels[["DFS_logitWithRxGrade"]])){
#         
#         linearModelsSumm[["DFS_logitWithRxGrade"]] <- summary( linearModels[["DFS_logitWithRxGrade"]])
#         textOut <- capture.output( linearModelsSumm[["DFS_logitWithRxGrade"]]$coefficients)
#         cat(textOut,sep="\n",append=TRUE,file=summaryFile) 
#         linearModelsSumm[["DFS_WithGradeMetaClusters_used"]] <- paste(unique(dataMatrixHist$community),collapse=",")
#         
#         predictor <- predict(linearModels[["DFS_logitWithRxGrade"]],type="response")
#         png(filename=paste0(saveDir,"/",experimentName,"_DFS_withGrade_ROC_",Sys.Date(),".png"),width=1000,height=1000,res=200)
#         roc_data <- roc(response=dataMatrixHist$DFS,predictor=predictor,plot=TRUE,na.rm=TRUE,auc=TRUE)
#         dev.off()
#         textOut <- capture.output(roc_data)
#         cat("roc_info\n",textOut,sep="\n",append=TRUE,file=summaryFile)
#         cat("\nAUC: ",roc_data$auc,"\n",append=TRUE,file=summaryFile)
#         
#       }else{
#         
#         cat(gsub(" ","_",paste0("\nHist grade with DFS analysis returned NA\n.")),append=TRUE,file=summaryFile)  
#       }
#       
#       
#       message("Running pCR models")
#       #####pCR
#       cat(gsub(" ","_","\nLinear analysis with only pCR:\n"),append=TRUE,file=summaryFile)
#       
#       dataMatrix <- sampleClustCommPhenoData[which(!is.na((sampleClustCommPhenoData$pCR))), ]
#       dataMatrix$community <- as.factor(dataMatrix$community)
#       dataMatrix$pCR <- as.factor(dataMatrix$pCR)
#       dataMatrix$hist_grade <- as.factor(dataMatrix$hist_grade)
#       dataMatrix$chemotherapyClass <- as.factor(dataMatrix$chemotherapyClass)
#       dataMatrix$anti_estrogen <- as.factor(dataMatrix$anti_estrogen)
#       dataMatrix$anti_HER2 <- as.factor(dataMatrix$anti_HER2)
#       cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrix$studyNum),"\n")),append=TRUE,file=summaryFile)
#       cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrix$community),"\n")),append=TRUE,file=summaryFile)
#       
#       #strongly predicts pCR.
#       linearModels[["pCR_logitAlone"]] <- glm(pCR~community,data=dataMatrix,family=binomial(link="logit"))
#       
#       predictor <- predict(linearModels[["pCR_logitAlone"]],type="response")
#       png(filename=paste0(saveDir,"/",experimentName,"_pCR_ROC_",Sys.Date(),".png"),width=1000,height=1000,res=200)
#       roc_data <- roc(response=dataMatrix$pCR,predictor=predictor,plot=TRUE,na.rm=TRUE,auc=TRUE)
#       dev.off()
#       textOut <- capture.output(roc_data)
#       cat("roc_info\n",textOut,sep="\n",append=TRUE,file=summaryFile)
#       cat("\nAUC: ",roc_data$auc,"\n",append=TRUE,file=summaryFile)
#       
#       
#       linearModelsSumm[["pCR_logitAlone"]] <- summary(  linearModels[["pCR_logitAlone"]] ) 
#       textOut <- capture.output( linearModelsSumm[["pCR_logitAlone"]]$coefficients)
#       cat(textOut,sep="\n",append=TRUE,file=summaryFile)
#         textOut <- capture.output( linearModelsSumm[["pCR_logitAlone"]])
#       cat(textOut,sep="\n",append=TRUE,file=summaryFile)
#       linearModelsSumm[["pCR_numPatients"]] <- nrow(dataMatrix)
#       cat("\nNumber_of_patients_for_pCR_model: ",nrow(dataMatrix),"\n",append=TRUE,file=summaryFile)
#       linearModelsSumm[["pCR_aloneMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
#       ##with therapies
#       #interesting: community becomes insignificant if add in therapies now: this could be
#       #because community is most likely correlated with treatments; it's certainly correlated with pam50 status.
#       #low Rsquared.
#       cat(gsub(" ","_",paste0("\nLinear analysis with pCR and treatment:\n")),append=TRUE,file=summaryFile)
#       dataMatrix <- dataMatrix[which(!is.na(dataMatrix$chemotherapyClass)), ]
#       dataMatrix <- dataMatrix[which(!is.na(dataMatrix$anti_estrogen)), ]
#       dataMatrix <- dataMatrix[which(!is.na(dataMatrix$anti_HER2)), ]
#       dataMatrix <- data.frame(dataMatrix)
#       dataMatrix$community <- as.factor(dataMatrix$community)
#       dataMatrix$pCR <- as.factor(dataMatrix$pCR)
#       dataMatrix$hist_grade <- as.factor(dataMatrix$hist_grade)
#       dataMatrix$chemotherapyClass <- as.factor(dataMatrix$chemotherapyClass)
#       dataMatrix$anti_estrogen <- as.factor(dataMatrix$anti_estrogen)
#       dataMatrix$anti_HER2 <- as.factor(dataMatrix$anti_HER2)
#       cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrix$studyNum),"\n")),append=TRUE,file=summaryFile)
#       cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrix$community),"\n")),append=TRUE,file=summaryFile)
#       
#       linearModelsSumm[["pCR_withRx_numPatients"]] <- nrow(dataMatrix)
#       cat(gsub(" ","_",paste0("\nNumber_of_patients_for_pCR_withRx_model: ",nrow(dataMatrix),"\n")),append=TRUE,file=summaryFile)
#       
#       if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))>1 &&
#          length(unique(dataMatrix$anti_HER2))>1 ){
#         
#         linearModels[["pCR_logitWithRx"]] <- tryCatch(glm(pCR~community+chemotherapyClass+anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
#                                                       error = function(e) {
#                                                         return(NA)
#                                                       }
#         )
#         
#                                     #remove communities for ANOVA
#            reducedModel <- tryCatch(glm(pCR~chemotherapyClass+anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
#                                                       error = function(e) {
#                                                         return(NA)
#                                                       }
#         )
#          linearModels[["pCR_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["pCR_logitWithRx"]],test="Chisq")
#     
# 
#         #now try with grade.
#         cat(gsub(" ","_",paste0("\nLinear analysis with pCR and histological grade :\n")),append=TRUE,file=summaryFile)
#         dataMatrixHist <- dataMatrix[which(!is.na(dataMatrix$hist_grade)),]
#         cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
#         cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
#         
#         linearModelsSumm[["pCR_withGrade_numPatients"]] <- nrow(dataMatrixHist)
#         cat(gsub(" ","_",paste0("\nNumber_of_patients_for_pCR_withGrade_model: ",nrow(dataMatrixHist),"\n")),append=TRUE,file=summaryFile)
#         linearModels[["pCR_logitWithRxGrade"]] <- tryCatch( glm(pCR~community+chemotherapyClass+anti_estrogen+anti_HER2+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
#                                                             error = function(e) {
#                                                               return(NA)
#                                                             }
#         )
#         
#         
#       }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))>1 
#                && length(unique(dataMatrix$anti_HER2))==1 ){
#         
#         
#         cat(gsub(" ","_",paste0("\nFor pCR analysis, all anti_HER2 values the same: ",unique(unique(dataMatrix$anti_HER2)),"\n")),
#             append=TRUE,file=summaryFile)
#         linearModels[["pCR_logitWithRx"]] <- tryCatch(glm(pCR~community+chemotherapyClass+anti_estrogen,data=dataMatrix,family=binomial(link="logit")),
#                                                       error = function(e) {
#                                                         return(NA)
#                                                       }
#         )
#                                             #remove communities for ANOVA
#            reducedModel <- tryCatch(glm(pCR~chemotherapyClass+anti_estrogen,data=dataMatrix,family=binomial(link="logit")),
#                                                       error = function(e) {
#                                                         return(NA)
#                                                       }
#         )
#          linearModels[["pCR_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["pCR_logitWithRx"]],test="Chisq")
#         #now try with stage, grade.
#         cat(gsub(" ","_",paste0("\nLinear analysis with pCR and histological grade :\n")),append=TRUE,file=summaryFile)
#         dataMatrixHist <- dataMatrix[which(!is.na(dataMatrix$hist_grade)),]
#  
#         cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
#         cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
#         
#         linearModelsSumm[["pCR_withGrade_numPatients"]] <- nrow(dataMatrixHist)
#         cat(gsub(" ","_",paste0("\nNumber_of_patients_for_pCR_withGrade_model: ",nrow(dataMatrixHist),"\n")),append=TRUE,file=summaryFile)
#         linearModels[["pCR_logitWithRxGrade"]] <- tryCatch( glm(pCR~community+chemotherapyClass+anti_estrogen+anti_HER2+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
#                                                             error = function(e) {
#                                                               return(NA)
#                                                             }
#         )
#         
#         
#       }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))==1 
#                && length(unique(dataMatrix$anti_HER2))>1 ){
#         
#         cat(gsub(" ","_",paste0("\nFor pCR analysis, all anti_estrogen values the same: ",unique(unique(dataMatrix$anti_estrogen)),"\n")),
#             append=TRUE,file=summaryFile)
#         
#         linearModels[["pCR_logitWithRx"]] <- tryCatch(glm(pCR~community+chemotherapyClass+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
#                                                       error = function(e) {
#                                                         return(NA)
#                                                       }
#         )
#                                                   #remove communities for ANOVA
#            reducedModel <- tryCatch(glm(pCR~chemotherapyClass+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
#                                                       error = function(e) {
#                                                         return(NA)
#                                                       }
#         )
#          linearModels[["pCR_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["pCR_logitWithRx"]],test="Chisq")
# 
#         #now try with stage, grade.
#         cat(gsub(" ","_",paste0("\nLinear analysis with pCR and histological grade :\n")),append=TRUE,file=summaryFile)
#         dataMatrixHist <- dataMatrix[which(!is.na(dataMatrix$hist_grade)),]
#  
#         cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
#         cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
#         
#         linearModelsSumm[["pCR_withGrade_numPatients"]] <- nrow(dataMatrixHist)
#         cat(gsub(" ","_",paste0("\nNumber_of_patients_for_pCR_withGrade_model: ",nrow(dataMatrixHist),"\n")),append=TRUE,file=summaryFile)
#         linearModels[["pCR_logitWithRxGrade"]] <- tryCatch( glm(pCR~community+chemotherapyClass+anti_HER2+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
#                                                             error = function(e) {
#                                                               return(NA)
#                                                             }
#         )
#         
#         
#         
#         
#         
#       }else if(length(unique(dataMatrix$chemotherapyClass))==1 && length(unique(dataMatrix$anti_estrogen))>1 
#                && length(unique(dataMatrix$anti_HER2))>1 ){
#         
#         
#         cat(gsub(" ","_",paste0("\nFor pCR analysis, all chemo values the same: ",unique(dataMatrix$chemotherapyClass),"\n")),
#             append=TRUE,file=summaryFile)
#         
#         linearModels[["pCR_logitWithRx"]] <- tryCatch(glm(pCR~community+anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
#                                                       error = function(e) {
#                                                         return(NA)
#                                                       }
#         )
#    #remove community for ANOVA
#         reducedModel <- tryCatch(glm(pCR~anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
#                                                       error = function(e) {
#                                                         return(NA)
#                                                       }
#         )
#          linearModels[["pCR_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["pCR_logitWithRx"]],test="Chisq")
# 
#  
#         #now try with stage, grade.
#         cat(gsub(" ","_",paste0("\nLinear analysis with pCR and histological grade :\n")),append=TRUE,file=summaryFile)
#         dataMatrixHist <- dataMatrix[which(!is.na(dataMatrix$hist_grade)),]
#  
#         cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
#         cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
#         
#         linearModelsSumm[["pCR_withGrade_numPatients"]] <- nrow(dataMatrix)
#         cat(gsub(" ","_",paste0("\nNumber_of_patients_for_pCR_withGrade_model: ",nrow(dataMatrixHist),"\n")),append=TRUE,file=summaryFile)
#         linearModels[["pCR_logitWithRxGrade"]] <- tryCatch( glm(pCR~community+anti_estrogen+anti_HER2+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
#                                                             error = function(e) {
#                                                               return(NA)
#                                                             }
#         )
#         
#         
#         
#         
#         
#       }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))==1 
#                && length(unique(dataMatrix$anti_HER2))==1 ){
#         
#         
#         cat(gsub(" ","_",paste0("\nFor pCR analysis, the anti HER2 and estrogen variables have only one unique value each.\n")),append=TRUE,
#             file=summaryFile)
#         
#         
#         linearModels[["pCR_logitWithRx"]] <- tryCatch(glm(pCR~community+chemotherapyClass,data=dataMatrix,family=binomial(link="logit")),
#                                                       error = function(e) {
#                                                         return(NA)
#                                                       }
#         )
# 
#    #remove community for ANOVA
#         reducedModel <- tryCatch(glm(pCR~chemotherapyClass,data=dataMatrix,family=binomial(link="logit")),
#                                                       error = function(e) {
#                                                         return(NA)
#                                                       }
#         )
#          linearModels[["pCR_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["pCR_logitWithRx"]],test="Chisq")
# 
# 
#         #now try with stage, grade.
#         cat(gsub(" ","_",paste0("\nLinear analysis with pCR and histological grade :\n")),append=TRUE,file=summaryFile)
#         dataMatrixHist <- dataMatrix[which(!is.na(dataMatrix$hist_grade)),]
#  
#         cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
#         cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
#         
#         linearModelsSumm[["pCR_withGrade_numPatients"]] <- nrow(dataMatrixHist)
#         cat(gsub(" ","_",paste0("\nNumber_of_patients_for_pCR_withGrade_model: ",nrow(dataMatrixHist),"\n")),append=TRUE,file=summaryFile)
#         linearModels[["pCR_logitWithRxGrade"]] <- tryCatch( glm(pCR~community+chemotherapyClass+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
#                                                             error = function(e) {
#                                                               return(NA)
#                                                             }
#         )
#         
#         
#         
#         
#       }else{
#         
#         # don't run models
#         
#       }
#       
#       cat(gsub(" ","_",paste0("\nLinear analysis with pCR and treatment:\n")),append=TRUE,file=summaryFile)
#       if(any(!is.na(linearModels[["pCR_logitWithRx"]]))){
#         
#         linearModelsSumm[["pCR_logitWithRx"]] <- summary( linearModels[["pCR_logitWithRx"]])
#         textOut <- capture.output( linearModelsSumm[["pCR_logitWithRx"]]$coefficients)
#         cat(textOut,sep="\n",append=TRUE,file=summaryFile) 
#         linearModelsSumm[["pCR_WithRxMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
#                 textOut <- capture.output( linearModelsSumm [["pCR_logitWithRx"]])
#         cat(textOut,sep="\n",append=TRUE,file=summaryFile)
#         cat(gsub(" ","_",paste0("\nLinear analysis with pCR and treatment ANOVA for meta-clusters:\n")),append=TRUE,file=summaryFile)
#         textOut <- capture.output(linearModels[["pCR_logitWithRx_ANOVA"]])
#         cat(textOut,sep="\n",append=TRUE,file=summaryFile)
#         
#         predictor <- predict(linearModels[["pCR_logitWithRx"]],type="response")
#         png(filename=paste0(saveDir,"/",experimentName,"_pCR_withRx_ROC_",Sys.Date(),".png"),width=1000,height=1000,res=200)
#         roc_data <- roc(response=dataMatrix$pCR,predictor=predictor,plot=TRUE,na.rm=TRUE,auc=TRUE)
#         dev.off()
#         textOut <- capture.output(roc_data)
#         cat("roc_info\n",textOut,sep="\n",append=TRUE,file=summaryFile)
#         cat("\nAUC: ",roc_data$auc,"\n",append=TRUE,file=summaryFile)
#         
#         
#         
#       }else{
#         
#         cat(gsub(" ","_",paste0("\nRx with pCR analysis returned NA\n.")),append=TRUE,file=summaryFile)  
#       }
#       
#       
#       #now try with stage, grade.
#       cat(gsub(" ","_",paste0("\nLinear analysis with pCR and histological grade :\n")),append=TRUE,file=summaryFile)
#       cat(gsub(" ","_",paste0("\nUnique hist grade values :",unique(dataMatrixHist$hist_grade),"\n")),append=TRUE,file=summaryFile)
# 
#       if(any(!is.na(linearModels[["pCR_logitWithRxGrade"]]))){
#         
#         linearModelsSumm[["pCR_logitWithRxGrade"]] <- summary( linearModels[["pCR_logitWithRxGrade"]])
#         textOut <- capture.output( linearModelsSumm[["pCR_logitWithRxGrade"]]$coefficients)
#         cat(textOut,sep="\n",append=TRUE,file=summaryFile) 
#         linearModelsSumm[["pCR_logitWithRxGrade"]] <- paste(unique(dataMatrixHist$community),collapse=",")
#         
#         predictor <- predict(linearModels[["pCR_logitWithRxGrade"]],type="response")
#         png(filename=paste0(saveDir,"/",experimentName,"_pCR_withGrade_ROC_",Sys.Date(),".png"),width=1000,height=1000,res=200)
#         roc_data <- roc(response=dataMatrixHist$pCR,predictor=predictor,plot=TRUE,na.rm=TRUE,auc=TRUE)
#         dev.off()
#         textOut <- capture.output(roc_data)
#         cat("roc_info\n",textOut,sep="\n",append=TRUE,file=summaryFile)
#         cat("\nAUC: ",roc_data$auc,"\n",append=TRUE,file=summaryFile)
#         
#       }else{
#         
#         cat(gsub(" ","_",paste0("\nHist grade with pCR analysis returned NA.")),append=TRUE,file=summaryFile)  
#       }
#       
breast_survivalData=runBreastCancerBinarySurvModels(sampleClustCommPhenoData,expName,saveDir,summaryFile)
#       #end of if not ovarian.
    }
    #end of if survival analysis.
    
  }else{
    
    clinicalTables <- NULL
    
  }
  
  output <- list(sampleClustCommPhenoData=sampleClustCommPhenoDataOrig,aggregateData=aggregateData,commInfo=commInfo,finalEdgeInfo=finalEdgeInfo,CoINcIDE_computeEdgesObject=CoINcIDE_output,
                 networkStats=networkStats,sigGenes=sigGenes,GSEA_out_unique=GSEA_out_unique,survivalResults=result,
                 linearModelSummaries=linearModelsSumm,linearSurvModels=linearModels,
                 clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                 clinicalTables=clinicalTables,breast_survivalData=breast_survivalData)
  
  return(output)
}

runBreastCancerBinarySurvModels <- function(sampleClustCommPhenoData,expName="test",saveDir="./",summaryFile=paste0(saveDir,"/survivalOutput.txt")){
  library("pROC")
  options(bitmapType="cairo")
  plotToScreen <- FALSE
  experimentName <- expName
  linearModels <- list()
  linearModelsSumm <- list()
  ROC_list <- list()
   #message("running survival data on breast")
      message("Running RFS models")
      cat(gsub(" ","_",paste0("\nLinear analysis with only RFS:\n")),append=TRUE,file=summaryFile)
      
      dataMatrix <- sampleClustCommPhenoData[which(!is.na((sampleClustCommPhenoData$RFS))), ]
      cat(gsub(" ","_",paste0("Studies used:\n",unique(dataMatrix$studyNum),"\n")),append=TRUE,file=summaryFile)
      cat(gsub(" ","_",paste0("\nMeta-clusters used:\n",unique(dataMatrix$community),"\n")),append=TRUE,file=summaryFile)
      
      dataMatrix <- data.frame(dataMatrix)
      dataMatrix$community <- as.factor(dataMatrix$community)
      dataMatrix$RFS <- as.factor(dataMatrix$RFS)
      
      linearModels[["RFS_logitAlone"]] <- tryCatch(glm(RFS~community,data=dataMatrix,family=binomial(link="logit")),
                                                   error = function(e) {
                                                     return(NA)
                                                   }
      )
  
      
      if(any(!is.na(  linearModels[["RFS_logitAlone"]]))){
        ##ROC curves
        #above: font size needs to be smaller for subtypes
        predictor <- predict(linearModels[["RFS_logitAlone"]],type="response")
        
        png(filename=paste0(saveDir,"/",experimentName,"_RFS_ROC_",Sys.Date(),".png"),width=1000,height=1000,res=200)
        roc_data <- roc(response=dataMatrix$RFS,predictor=predictor,plot=TRUE,na.rm=TRUE,auc=TRUE)
        dev.off()
        textOut <- capture.output(roc_data)
        cat("roc_info\n",textOut,sep="\n",append=TRUE,file=summaryFile)
        cat("\nAUC : ",roc_data$auc,"\n",append=TRUE,file=summaryFile)
        
        linearModelsSumm[["RFS_logitAlone"]] <- summary(linearModels[["RFS_logitAlone"]] ) 
        textOut <- capture.output( linearModelsSumm [["RFS_logitAlone"]]$coefficients)
        cat(textOut,sep="\n",append=TRUE,file=summaryFile)
        textOut <- capture.output( linearModelsSumm [["RFS_logitAlone"]])
        cat(textOut,sep="\n",append=TRUE,file=summaryFile)
        linearModelsSumm[["RFS_numPatients"]] <- nrow(dataMatrix)
        cat("\nNumber_of_patients_for_RFS__model: ",nrow(dataMatrix),"\n",append=TRUE,file=summaryFile)
        linearModelsSumm[["RFS_aloneMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
        
        ROC_list[["RFS_logitAlone"]] <- roc_data
        
      }else{
        cat("\nRFS model alone returned NA.\n")
        cat("\nRFS model alone returned NA.\n",append=TRUE,file=outputFile)
      }
      ##with therapies
      #interesting: community becomes insignificant if add in therapies now: this could be
      #because community is most likely correlated with treatments; it's certainly correlated with pam50 status.
      #low Rsquared.
      cat(gsub(" ","_",paste0("\nLinear analysis with RFS and treatment:\n")),append=TRUE,file=summaryFile)
      dataMatrix <- dataMatrix[which(!is.na(dataMatrix$chemotherapyClass)), ]
      dataMatrix <- dataMatrix[which(!is.na(dataMatrix$anti_estrogen)), ]
      dataMatrix <- dataMatrix[which(!is.na(dataMatrix$anti_HER2)), ]
      dataMatrix <- data.frame(dataMatrix)
      dataMatrix$community <- as.factor(dataMatrix$community)
      dataMatrix$RFS <- as.factor(dataMatrix$RFS)
      dataMatrix$hist_grade <- as.factor(dataMatrix$hist_grade)
      dataMatrix$chemotherapyClass <- as.factor(dataMatrix$chemotherapyClass)
      dataMatrix$anti_estrogen <- as.factor(dataMatrix$anti_estrogen)
      dataMatrix$anti_HER2 <- as.factor(dataMatrix$anti_HER2)
      cat(gsub(" ","_",paste0("\nMeta-clusters used:\n",unique(dataMatrix$community),"\n")),append=TRUE,file=summaryFile)
      
      linearModelsSumm [["RFS_withRx_numPatients"]] <- nrow(dataMatrix)
      cat("\nNumber_of_patients_for_RFS_with_Rx_model: ",nrow(dataMatrix),"\n",append=TRUE,file=summaryFile)
      
      if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))>1 
         && length(unique(dataMatrix$anti_HER2))>1 ){
        
        linearModels[["RFS_logitWithRx"]] <- tryCatch(glm(RFS~community+chemotherapyClass+anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
        
        #no community
        reducedModel <- tryCatch(glm(RFS~chemotherapyClass+anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
         linearModels[["RFS_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["RFS_logitWithRx"]],test="Chisq")
           
     
        #now try with  grade.    
        cat(gsub(" ","_",paste0("\nLinear analysis with RFS and histological grade :\n")),append=TRUE,file=summaryFile)
        dataMatrixHist <- dataMatrix[which(!is.na((dataMatrix$hist_grade))),]
 
        
        cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
        
        linearModelsSumm[["RFS_withGrade_numPatients"]] <- nrow(dataMatrixHist)
        cat("\nNumber_of_patients_for_RFS_with_grade_model: ",nrow(dataMatrixHist),"\n",append=TRUE,file=summaryFile)
        
        linearModels[["RFS_logitWithRxGrade"]] <- tryCatch( glm(RFS~community+chemotherapyClass+anti_estrogen+anti_HER2+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
                                                            error = function(e) {
                                                              return(NA)
                                                            }
        )
        
        
      }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))>1 
               && length(unique(dataMatrix$anti_HER2))==1 ){
        
        
        cat(gsub(" ","_",paste0("\nFor RFS analysis, all anti_HER2 values the same: ",unique(unique(dataMatrix$anti_HER2)),"\n")),
            append=TRUE,file=summaryFile)
        linearModels[["RFS_logitWithRx"]] <- tryCatch(glm(RFS~community+chemotherapyClass+anti_estrogen,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
                                                      
        )
        
               #no community
        reducedModel <- tryCatch(glm(RFS~chemotherapyClass+anti_estrogen,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
         linearModels[["RFS_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["RFS_logitWithRx"]],test="Chisq")
        
        #now try with  grade.    
        cat(gsub(" ","_",paste0("\nLinear analysis with RFS and histological grade :\n")),append=TRUE,file=summaryFile)
        dataMatrixHist <- dataMatrix[which(!is.na((dataMatrix$hist_grade))),]
 
        
        cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
        
        linearModelsSumm[["RFS_withGrade_numPatients"]] <- nrow(dataMatrixHist)
        cat("\nNumber_of_patients_for_RFS_with_grade_model: ",nrow(dataMatrixHist),"\n",append=TRUE,file=summaryFile)
        linearModels[["RFS_logitWithRxGrade"]] <- tryCatch( glm(RFS~community+chemotherapyClass+anti_estrogen+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
                                                            error = function(e) {
                                                              return(NA)
                                                            }
        )
        
        
        
      }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))==1 
               && length(unique(dataMatrix$anti_HER2))>1 ){
        
        cat(gsub(" ","_",paste0("\nFor RFS analysis, all anti_estrogen values the same: ",unique(unique(dataMatrix$anti_estrogen)),"\n")),
            append=TRUE,file=summaryFile)
        
        linearModels[["RFS_logitWithRx"]] <- tryCatch(glm(RFS~community+chemotherapyClass+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
                                                      
        )
        
                   #no community
        reducedModel <- tryCatch(glm(RFS~chemotherapyClass+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
         linearModels[["RFS_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["RFS_logitWithRx"]],test="Chisq")
        
 
        
        #now try with  grade.    
        cat(gsub(" ","_",paste0("\nLinear analysis with RFS and histological grade :\n")),append=TRUE,file=summaryFile)
        dataMatrixHist <- dataMatrix[which(!is.na((dataMatrix$hist_grade))),]
 
        
        cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
        
        linearModelsSumm[["RFS_withGrade_numPatients"]] <- nrow(dataMatrixHist)
        cat("\nNumber_of_patients_for_RFS_with_grade_model: ",nrow(dataMatrixHist),"\n",append=TRUE,file=summaryFile)
        linearModels[["RFS_logitWithRxGrade"]] <- tryCatch( glm(RFS~community+chemotherapyClass+anti_HER2+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
                                                            error = function(e) {
                                                              return(NA)
                                                            }
        )
        
        
      }else if(length(unique(dataMatrix$chemotherapyClass))==1 && length(unique(dataMatrix$anti_estrogen))>1 
               && length(unique(dataMatrix$anti_HER2))>1 ){
        
        
        cat(gsub(" ","_",paste0("\nFor RFS analysis, all chemo values the same: ",unique(unique(dataMatrix$chemotherapyClass)),"\n")),
            append=TRUE,file=summaryFile)
        
        linearModels[["RFS_logitWithRx"]] <- tryCatch(glm(RFS~community+anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
                                                      
        )
        
                           #no community
        reducedModel <- tryCatch(glm(RFS~anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
         linearModels[["RFS_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["RFS_logitWithRx"]],test="Chisq")
        #now try with  grade.    
        cat(gsub(" ","_",paste0("\nLinear analysis with RFS and histological grade :\n")),append=TRUE,file=summaryFile)
        dataMatrixHist <- dataMatrix[which(!is.na((dataMatrix$hist_grade))),]
 
        
        cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
        
        linearModelsSumm[["RFS_withGrade_numPatients"]] <- nrow(dataMatrixHist)
        cat("\nNumber_of_patients_for_RFS_with_grade_model: ",nrow(dataMatrixHist),"\n",append=TRUE,file=summaryFile)
        linearModels[["RFS_logitWithRxGrade"]] <- tryCatch( glm(RFS~community+anti_estrogen+anti_HER2+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
                                                            error = function(e) {
                                                              return(NA)
                                                            }
        )
        
        
        
      }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))==1 
               && length(unique(dataMatrix$anti_HER2))==1 ){
        
        
        cat(gsub(" ","_",paste0("\nFor RFS analysis, the anti HER2 and estrogen variables have only one unique value each.\n")),append=TRUE,
            file=summaryFile)
        
        
        linearModels[["RFS_logitWithRx"]] <- tryCatch(glm(RFS~community+chemotherapyClass,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
                                                      
        )
         reducedModel <- tryCatch(glm(RFS~chemotherapyClass,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
         linearModels[["RFS_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["RFS_logitWithRx"]],test="Chisq")
    
        
        #now try with  grade.    
        cat(gsub(" ","_",paste0("\nLinear analysis with RFS and histological grade :\n")),append=TRUE,file=summaryFile)
        dataMatrixHist <- dataMatrix[which(!is.na((dataMatrix$hist_grade))),]
 
        
        cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
        
        linearModelsSumm[["RFS_withGrade_numPatients"]] <- nrow(dataMatrixHist)
        cat("\nNumber_of_patients_for_RFS_with_grade_model: ",nrow(dataMatrixHist),"\n",append=TRUE,file=summaryFile)
        linearModels[["RFS_logitWithRxGrade"]] <- tryCatch( glm(RFS~community+chemotherapyClass+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
                                                            error = function(e) {
                                                              return(NA)
                                                            }
        )
        
        
        
        
      }else{
        
        linearModels[["RFS_logitWithRx"]] <- NA
        
      }
      
      cat(gsub(" ","_",paste0("\nLinear analysis with RFS and treatment:\n")),append=TRUE,file=summaryFile)
      if(any(!is.na(linearModels[["RFS_logitWithRx"]]))){
        
        linearModelsSumm[["RFS_logitWithRx"]] <- summary( linearModels[["RFS_logitWithRx"]])
        textOut <- capture.output( linearModelsSumm[["RFS_logitWithRx"]]$coefficients)
        cat(textOut,sep="\n",append=TRUE,file=summaryFile)
        textOut <- capture.output( linearModelsSumm [["RFS_logitWithRx"]])
        cat(textOut,sep="\n",append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nLinear analysis with RFS and treatment ANOVA for meta-clusters:\n")),append=TRUE,file=summaryFile)
        textOut <- capture.output(linearModels[["RFS_logitWithRx_ANOVA"]])
        cat(textOut,sep="\n",append=TRUE,file=summaryFile)
        
        linearModelsSumm[["RFS_WithRxMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",") 
        predictor <- predict(linearModels[["RFS_logitWithRx"]],type="response")
        png(filename=paste0(saveDir,"/",experimentName,"_RFS_withRx_ROC_",Sys.Date(),".png"),width=1000,height=1000,res=200)
        roc_data <- roc(response=dataMatrix$RFS,predictor=predictor,plot=TRUE,na.rm=TRUE,auc=TRUE)
        dev.off()
        textOut <- capture.output(roc_data)   
        ROC_list[["RFS_logitWithRx"]] <- roc_data
        cat("\n roc_info:",textOut,sep="\n",append=TRUE,file=summaryFile)
        cat("\nAUC : ",roc_data$auc,"\n",append=TRUE,file=summaryFile)
        
        linearModelsSumm[["RFS_logitWithRx"]] <- summary( linearModels[["RFS_logitWithRx"]])
        textOut <- capture.output( linearModelsSumm[["RFS_logitWithRx"]]$coefficients)
        cat(textOut,sep="\n",append=TRUE,file=summaryFile)
        
        linearModelsSumm[["RFS_WithRxMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
        
      }else{
        
        cat(gsub(" ","_",paste0("\nRFS analysis with therapies returned NA.\n")),append=TRUE,file=summaryFile)
      }
      
      #now record grade.
      cat(gsub(" ","_",paste0("\nLinear analysis with RFS and histological grade :\n")),append=TRUE,file=summaryFile)
      cat(gsub(" ","_",paste0("\nUnique hist grade values :",unique(dataMatrixHist$hist_grade),"\n")),append=TRUE,file=summaryFile)
      
      
      if(any(!is.na(linearModels[["RFS_logitWithRxGrade"]]))){
        
        predictor <- predict(linearModels[["RFS_logitWithRxGrade"]],type="response")
        png(filename=paste0(saveDir,"/",experimentName,"_RFS_withGrade_ROC_",Sys.Date(),".png"),width=1000,height=1000,res=200)
        roc_data <- roc(response=dataMatrixHist$RFS,predictor=predictor,plot=TRUE,na.rm=TRUE,auc=TRUE)
        dev.off()
          ROC_list[["RFS_logitWithRxGrade"]] <- roc_data
        textOut <- capture.output(roc_data)
        cat("roc_info:\n",textOut,sep="\n",append=TRUE,file=summaryFile)
        cat("\nAUC: ",roc_data$auc,"\n",append=TRUE,file=summaryFile)
        
        linearModelsSumm[["RFS_logitWithRxGrade"]] <- summary( linearModels[["RFS_logitWithRxGrade"]])
        textOut <- capture.output( linearModelsSumm[["RFS_logitWithRxGrade"]]$coefficients)
        cat(textOut,sep="\n",append=TRUE,file=summaryFile) 
        linearModelsSumm[["RFS_WithGradeMetaClusters_used"]] <- paste(unique(dataMatrixHist$community),collapse=",")
        
      }else{
        
        cat(gsub(" ","_",paste0("\nHist grade with RFS analysis returned NA.\n")),append=TRUE,file=summaryFile)  
      }
      
      
      
      message("Running DFS models")
      ###DFS
      cat(gsub(" ","_",paste0("\nLinear analysis with only DFS:\n")),append=TRUE,file=summaryFile)
      
      dataMatrix <- sampleClustCommPhenoData[which(!is.na((sampleClustCommPhenoData$DFS))), ]
      dataMatrix$community <- as.factor(dataMatrix$community)
      dataMatrix$DFS <- as.factor(dataMatrix$DFS)
      dataMatrix$hist_grade <- as.factor(dataMatrix$hist_grade)
      dataMatrix$chemotherapyClass <- as.factor(dataMatrix$chemotherapyClass)
      dataMatrix$anti_estrogen <- as.factor(dataMatrix$anti_estrogen)
      dataMatrix$anti_HER2 <- as.factor(dataMatrix$anti_HER2)
      cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrix$studyNum),"\n")),append=TRUE,file=summaryFile)
      cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrix$community),"\n")),append=TRUE,file=summaryFile)
      
      #strongly predicts DFS.
      linearModels[["DFS_logitAlone"]] <- glm(DFS~community,data=dataMatrix,family=binomial(link="logit"))
      png(filename=paste0(saveDir,"/",experimentName,"_DFS_ROC_",Sys.Date(),".png"),width=1000,height=1000,res=200)
      predictor <- predict(linearModels[["DFS_logitAlone"]],type="response")
      roc_data <- roc(response=dataMatrix$DFS,predictor=predictor,plot=TRUE,na.rm=TRUE,auc=TRUE)
      dev.off()
     ROC_list[["DFS_logitAlone"]] <- roc_data
      textOut <- capture.output(roc_data)
      cat("roc_info\n",textOut,sep="\n",append=TRUE,file=summaryFile)
      cat("\nAUC: ",roc_data$auc,"\n",append=TRUE,file=summaryFile)
      linearModelsSumm[["DFS_logitAlone"]] <- summary( linearModels[["DFS_logitAlone"]] ) 
      textOut <- capture.output( linearModelsSumm[["DFS_logitAlone"]]$coefficients)
      cat(textOut,sep="\n",append=TRUE,file=summaryFile)
      textOut <- capture.output( linearModelsSumm[["DFS_logitAlone"]])
      cat(textOut,sep="\n",append=TRUE,file=summaryFile)
      linearModelsSumm[["DFS_numPatients"]] <- nrow(dataMatrix)
      cat("\nNumber of patients: ",nrow(dataMatrix),"\n",append=TRUE,file=summaryFile)
      linearModelsSumm[["DFS_aloneMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
      
      ##with therapies
      #interesting: community becomes insignificant if add in therapies now: this could be
      #because community is most likely correlated with treatments; it's certainly correlated with pam50 status.
      #low Rsquared.
      cat(gsub(" ","_",paste0("\nLinear analysis with DFS and treatment:\n")),append=TRUE,file=summaryFile)
      dataMatrix <- dataMatrix[which(!is.na(dataMatrix$chemotherapyClass)), ]
      dataMatrix <- dataMatrix[which(!is.na(dataMatrix$anti_estrogen)), ]
      dataMatrix <- dataMatrix[which(!is.na(dataMatrix$anti_HER2)), ]
      dataMatrix <- data.frame(dataMatrix)
      dataMatrix$community <- as.factor(dataMatrix$community)
      dataMatrix$DFS <- as.factor(dataMatrix$DFS)
      dataMatrix$hist_grade <- as.factor(dataMatrix$hist_grade)
      dataMatrix$chemotherapyClass <- as.factor(dataMatrix$chemotherapyClass)
      dataMatrix$anti_estrogen <- as.factor(dataMatrix$anti_estrogen)
      dataMatrix$anti_HER2 <- as.factor(dataMatrix$anti_HER2)
      cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrix$studyNum),"\n")),append=TRUE,file=summaryFile)
      cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrix$community),"\n")),append=TRUE,file=summaryFile)
      
      linearModelsSumm[["DFS_withRx_numPatients"]] <- nrow(dataMatrix)
      cat("\nNumber_of_patients_for_DFS_with_Rx_model: ",nrow(dataMatrix),"\n",append=TRUE,file=summaryFile)
      
      if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))>1 && 
         length(unique(dataMatrix$anti_HER2))>1 ){
        
        linearModels[["DFS_logitWithRx"]] <- tryCatch(glm(DFS~community+chemotherapyClass+anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
        
        #remove communities for ANOVA
           reducedModel <- tryCatch(glm(DFS~chemotherapyClass+anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
         linearModels[["DFS_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["DFS_logitWithRx"]],test="Chisq")
    
    
        #now try with grade.
        cat(gsub(" ","_",paste0("\nLinear analysis with DFS and histological grade :\n")),append=TRUE,file=summaryFile)
        dataMatrixHist <- dataMatrix[which(!is.na((dataMatrix$hist_grade))),]
 

        cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
        
        linearModelsSumm[["DFS_withGrade_numPatients"]] <- nrow(dataMatrixHist)
        cat("\nNumber_of_patients_for_DFS_withGrade_model: ",nrow(dataMatrixHist),"\n",append=TRUE,file=summaryFile)
        linearModels[["DFS_logitWithRxGrade"]] <- tryCatch( glm(DFS~community+chemotherapyClass+anti_estrogen+anti_HER2+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
                                                            error = function(e) {
                                                              return(NA)
                                                            }
        )
        
        
        
      }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))>1 &&
               length(unique(dataMatrix$anti_HER2))==1 ){
        
        
        cat(gsub(" ","_",paste0("\nFor DFS analysis, all anti_HER2 values the same: ",unique(unique(dataMatrix$anti_HER2)),"\n")),
            append=TRUE,file=summaryFile)
        linearModels[["DFS_logitWithRx"]] <- tryCatch(glm(DFS~community+chemotherapyClass+anti_estrogen,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
                #remove communities for ANOVA
           reducedModel <- tryCatch(glm(DFS~chemotherapyClass+anti_estrogen,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
         linearModels[["DFS_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["DFS_logitWithRx"]],test="Chisq")
    

        #now try with grade.
        cat(gsub(" ","_",paste0("\nLinear analysis with DFS and histological grade :\n")),append=TRUE,file=summaryFile)
        dataMatrixHist <- dataMatrix[which(!is.na((dataMatrix$hist_grade))),]
 
        
        cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
        
        linearModelsSumm[["DFS_withGrade_numPatients"]] <- nrow(dataMatrixHist)
        cat("\nNumber_of_patients_for_DFS_withGrade_model: ",nrow(dataMatrixHist),"\n",append=TRUE,file=summaryFile)
        linearModels[["DFS_logitWithRxGrade"]] <- tryCatch( glm(DFS~community+chemotherapyClass+anti_estrogen+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
                                                            error = function(e) {
                                                              return(NA)
                                                            }
        )
        
        
        
        
      }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))==1 
               && length(unique(dataMatrix$anti_HER2))>1 ){
        
        cat(gsub(" ","_",paste0("\nFor DFS analysis, all anti_estrogen values the same: ",unique(unique(dataMatrix$anti_estrogen)),"\n")),
            append=TRUE,file=summaryFile)
        
        linearModels[["DFS_logitWithRx"]] <- tryCatch(glm(DFS~community+chemotherapyClass+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
                    #remove communities for ANOVA
           reducedModel <- tryCatch(glm(DFS~chemotherapyClass+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
         linearModels[["DFS_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["DFS_logitWithRx"]],test="Chisq")
    

        #now try with grade.
        cat(gsub(" ","_",paste0("\nLinear analysis with DFS and histological grade :\n")),append=TRUE,file=summaryFile)
        dataMatrixHist <- dataMatrix[which(!is.na((dataMatrix$hist_grade))),]
 
        
        cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
        
        linearModelsSumm[["DFS_withGrade_numPatients"]] <- nrow(dataMatrixHist)
        cat("\nNumber_of_patients_for_DFS_withGrade_model: ",nrow(dataMatrixHist),"\n",append=TRUE,file=summaryFile)
        linearModels[["DFS_logitWithRxGrade"]] <- tryCatch( glm(DFS~community+chemotherapyClass+anti_HER2+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
                                                            error = function(e) {
                                                              return(NA)
                                                            }
        )
        
        
        
      }else if(length(unique(dataMatrix$chemotherapyClass))==1 && length(unique(dataMatrix$anti_estrogen))>1 
               && length(unique(dataMatrix$anti_HER2))>1 ){
        
        
        cat(gsub(" ","_",paste0("\nFor DFS analysis, all chemo values the same: ",unique(unique(dataMatrix$chemotherapyClass)),"\n")),
            append=TRUE,file=summaryFile)
        
        linearModels[["DFS_logitWithRx"]] <- tryCatch(glm(DFS~community+anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
                          #remove communities for ANOVA
           reducedModel <- tryCatch(glm(DFS~anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
         linearModels[["DFS_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["DFS_logitWithRx"]],test="Chisq")
    

        #now try with grade.
        cat(gsub(" ","_",paste0("\nLinear analysis with DFS and histological grade :\n")),append=TRUE,file=summaryFile)
        dataMatrixHist <- dataMatrix[which(!is.na((dataMatrix$hist_grade))),]
 
        
        cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
        
        linearModelsSumm[["DFS_withGrade_numPatients"]] <- nrow(dataMatrixHist)
        cat("\nNumber_of_patients_for_DFS_withGrade_model: ",nrow(dataMatrixHist),"\n",append=TRUE,file=summaryFile)
        linearModels[["DFS_logitWithRxGrade"]] <- tryCatch( glm(DFS~community+anti_estrogen+anti_HER2+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
                                                            error = function(e) {
                                                              return(NA)
                                                            }
        )
        
        
        
      }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))==1 
               && length(unique(dataMatrix$anti_HER2))==1 ){
        
        
        cat(gsub(" ","_",paste0("\nFor DFS analysis, the anti HER2 and estrogen variables have only one unique value each.\n")),append=TRUE,
            file=summaryFile)
        
        
        linearModels[["DFS_logitWithRx"]] <- tryCatch(glm(DFS~community+chemotherapyClass,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
                               #remove communities for ANOVA
           reducedModel <- tryCatch(glm(DFS~chemotherapyClass,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
         linearModels[["DFS_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["DFS_logitWithRx"]],test="Chisq")
    

        #now try with grade.
        cat(gsub(" ","_",paste0("\nLinear analysis with DFS and histological grade :\n")),append=TRUE,file=summaryFile)
        dataMatrixHist <- dataMatrix[which(!is.na((dataMatrix$hist_grade))),]
 
        
        cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
        
        linearModelsSumm[["DFS_withGrade_numPatients"]] <- nrow(dataMatrixHist)
        cat("\nNumber_of_patients_for_DFS_withGrade_model: ",nrow(dataMatrixHist),"\n",append=TRUE,file=summaryFile)
        linearModels[["DFS_logitWithRxGrade"]] <- tryCatch( glm(DFS~community+chemotherapyClass+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
                                                            error = function(e) {
                                                              return(NA)
                                                            }
        )
        
        
        
      }else{
        
        # don't run models
        
      }
      
      cat(gsub(" ","_",paste0("\nLinear analysis with DFS and treatment:\n")),append=TRUE,file=summaryFile)
      if(any(!is.na(linearModels[["DFS_logitWithRx"]]))){
        
        predictor <- predict(linearModels[["DFS_logitWithRx"]],type="response")
        png(filename=paste0(saveDir,"/",experimentName,"_DFS_withRx_ROC_",Sys.Date(),".png"),width=1000,height=1000,res=200)
        roc_data <- roc(response=dataMatrix$DFS,predictor=predictor,plot=TRUE,na.rm=TRUE,auc=TRUE)
        dev.off()
          ROC_list[["DFS_logitWithRx"]] <- roc_data
        textOut <- capture.output(roc_data)
        cat("roc_info\n",textOut,sep="\n",append=TRUE,file=summaryFile)
        cat("\nAUC: ",roc_data$auc,"\n",append=TRUE,file=summaryFile)
        
        linearModelsSumm[["DFS_logitWithRx"]] <- summary( linearModels[["DFS_logitWithRx"]])
        textOut <- capture.output( linearModelsSumm[["DFS_logitWithRx"]]$coefficients)
        cat(textOut,sep="\n",append=TRUE,file=summaryFile)
                textOut <- capture.output( linearModelsSumm[["DFS_logitWithRx"]])
        cat(textOut,sep="\n",append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nLinear analysis with DFS and treatment ANOVA for meta-clusters:\n")),append=TRUE,file=summaryFile)
        textOut <- capture.output(linearModels[["DFS_logitWithRx_ANOVA"]])
        cat(textOut,sep="\n",append=TRUE,file=summaryFile)
        linearModelsSumm[["DFS_WithRxMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
        
      }else{
        
        cat(gsub(" ","_",paste0("\nRx with DFS analysis returned NA.\n")),append=TRUE,file=summaryFile)  
      }
      
      
      #now try with grade.
      cat(gsub(" ","_",paste0("\nLinear analysis with DFS and histological grade :\n")),append=TRUE,file=summaryFile)
      cat(gsub(" ","_",paste0("\nUnique hist grade values :",unique(dataMatrixHist$hist_grade),"\n")),append=TRUE,file=summaryFile)
      
      
      if(!is.na(linearModels[["DFS_logitWithRxGrade"]])){
        
        linearModelsSumm[["DFS_logitWithRxGrade"]] <- summary( linearModels[["DFS_logitWithRxGrade"]])
        textOut <- capture.output( linearModelsSumm[["DFS_logitWithRxGrade"]]$coefficients)
        cat(textOut,sep="\n",append=TRUE,file=summaryFile) 
        linearModelsSumm[["DFS_WithGradeMetaClusters_used"]] <- paste(unique(dataMatrixHist$community),collapse=",")
        
        predictor <- predict(linearModels[["DFS_logitWithRxGrade"]],type="response")
        png(filename=paste0(saveDir,"/",experimentName,"_DFS_withGrade_ROC_",Sys.Date(),".png"),width=1000,height=1000,res=200)
        roc_data <- roc(response=dataMatrixHist$DFS,predictor=predictor,plot=TRUE,na.rm=TRUE,auc=TRUE)
        dev.off()
          ROC_list[["DFS_logitWithRxGrade"]] <- roc_data
        textOut <- capture.output(roc_data)
        cat("roc_info\n",textOut,sep="\n",append=TRUE,file=summaryFile)
        cat("\nAUC: ",roc_data$auc,"\n",append=TRUE,file=summaryFile)
        
      }else{
        
        cat(gsub(" ","_",paste0("\nHist grade with DFS analysis returned NA\n.")),append=TRUE,file=summaryFile)  
      }
      
      
      message("Running pCR models")
      #####pCR
      cat(gsub(" ","_","\nLinear analysis with only pCR:\n"),append=TRUE,file=summaryFile)
      
      dataMatrix <- sampleClustCommPhenoData[which(!is.na((sampleClustCommPhenoData$pCR))), ]
      dataMatrix$community <- as.factor(dataMatrix$community)
      dataMatrix$pCR <- as.factor(dataMatrix$pCR)
      dataMatrix$hist_grade <- as.factor(dataMatrix$hist_grade)
      dataMatrix$chemotherapyClass <- as.factor(dataMatrix$chemotherapyClass)
      dataMatrix$anti_estrogen <- as.factor(dataMatrix$anti_estrogen)
      dataMatrix$anti_HER2 <- as.factor(dataMatrix$anti_HER2)
      cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrix$studyNum),"\n")),append=TRUE,file=summaryFile)
      cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrix$community),"\n")),append=TRUE,file=summaryFile)
      
      #strongly predicts pCR.
      
  #hmm...turns out that sometimes only one meta-cluster has pCR values.
      linearModels[["pCR_logitAlone"]] <- glm(pCR~community,data=dataMatrix,family=binomial(link="logit"))
      if(!is.na(linearModels[["pCR_logitAlone"]])){
              
      predictor <- predict(linearModels[["pCR_logitAlone"]],type="response")
      png(filename=paste0(saveDir,"/",experimentName,"_pCR_ROC_",Sys.Date(),".png"),width=1000,height=1000,res=200)
      roc_data <- roc(response=dataMatrix$pCR,predictor=predictor,plot=TRUE,na.rm=TRUE,auc=TRUE)
      dev.off()
      ROC_list[["pCR_logitAlone"]] <- roc_data
      textOut <- capture.output(roc_data)
      cat("roc_info\n",textOut,sep="\n",append=TRUE,file=summaryFile)
      cat("\nAUC: ",roc_data$auc,"\n",append=TRUE,file=summaryFile)
      
      
      linearModelsSumm[["pCR_logitAlone"]] <- summary(  linearModels[["pCR_logitAlone"]] ) 
      textOut <- capture.output( linearModelsSumm[["pCR_logitAlone"]]$coefficients)
      cat(textOut,sep="\n",append=TRUE,file=summaryFile)
        textOut <- capture.output( linearModelsSumm[["pCR_logitAlone"]])
      cat(textOut,sep="\n",append=TRUE,file=summaryFile)
      linearModelsSumm[["pCR_numPatients"]] <- nrow(dataMatrix)
      cat("\nNumber_of_patients_for_pCR_model: ",nrow(dataMatrix),"\n",append=TRUE,file=summaryFile)
      linearModelsSumm[["pCR_aloneMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
      
      }else{
        
         cat(gsub(" ","_",paste0("\n pCR alone analysis returned NA\n.")),append=TRUE,file=summaryFile)  
      }
      ##with therapies
      #interesting: community becomes insignificant if add in therapies now: this could be
      #because community is most likely correlated with treatments; it's certainly correlated with pam50 status.
      #low Rsquared.
      cat(gsub(" ","_",paste0("\nLinear analysis with pCR and treatment:\n")),append=TRUE,file=summaryFile)
      dataMatrix <- dataMatrix[which(!is.na(dataMatrix$chemotherapyClass)), ]
      dataMatrix <- dataMatrix[which(!is.na(dataMatrix$anti_estrogen)), ]
      dataMatrix <- dataMatrix[which(!is.na(dataMatrix$anti_HER2)), ]
      dataMatrix <- data.frame(dataMatrix)
      dataMatrix$community <- as.factor(dataMatrix$community)
      dataMatrix$pCR <- as.factor(dataMatrix$pCR)
      dataMatrix$hist_grade <- as.factor(dataMatrix$hist_grade)
      dataMatrix$chemotherapyClass <- as.factor(dataMatrix$chemotherapyClass)
      dataMatrix$anti_estrogen <- as.factor(dataMatrix$anti_estrogen)
      dataMatrix$anti_HER2 <- as.factor(dataMatrix$anti_HER2)
      cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrix$studyNum),"\n")),append=TRUE,file=summaryFile)
      cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrix$community),"\n")),append=TRUE,file=summaryFile)
      
      linearModelsSumm[["pCR_withRx_numPatients"]] <- nrow(dataMatrix)
      cat(gsub(" ","_",paste0("\nNumber_of_patients_for_pCR_withRx_model: ",nrow(dataMatrix),"\n")),append=TRUE,file=summaryFile)
      
      if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))>1 &&
         length(unique(dataMatrix$anti_HER2))>1 ){
        
        linearModels[["pCR_logitWithRx"]] <- tryCatch(glm(pCR~community+chemotherapyClass+anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
        
                                    #remove communities for ANOVA
           reducedModel <- tryCatch(glm(pCR~chemotherapyClass+anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
         linearModels[["pCR_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["pCR_logitWithRx"]],test="Chisq")
    

        #now try with grade.
        cat(gsub(" ","_",paste0("\nLinear analysis with pCR and histological grade :\n")),append=TRUE,file=summaryFile)
        dataMatrixHist <- dataMatrix[which(!is.na(dataMatrix$hist_grade)),]
        cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
        
        linearModelsSumm[["pCR_withGrade_numPatients"]] <- nrow(dataMatrixHist)
        cat(gsub(" ","_",paste0("\nNumber_of_patients_for_pCR_withGrade_model: ",nrow(dataMatrixHist),"\n")),append=TRUE,file=summaryFile)
        linearModels[["pCR_logitWithRxGrade"]] <- tryCatch( glm(pCR~community+chemotherapyClass+anti_estrogen+anti_HER2+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
                                                            error = function(e) {
                                                              return(NA)
                                                            }
        )
        
        
      }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))>1 
               && length(unique(dataMatrix$anti_HER2))==1 ){
        
        
        cat(gsub(" ","_",paste0("\nFor pCR analysis, all anti_HER2 values the same: ",unique(unique(dataMatrix$anti_HER2)),"\n")),
            append=TRUE,file=summaryFile)
        linearModels[["pCR_logitWithRx"]] <- tryCatch(glm(pCR~community+chemotherapyClass+anti_estrogen,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
                                            #remove communities for ANOVA
           reducedModel <- tryCatch(glm(pCR~chemotherapyClass+anti_estrogen,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
         linearModels[["pCR_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["pCR_logitWithRx"]],test="Chisq")
        #now try with stage, grade.
        cat(gsub(" ","_",paste0("\nLinear analysis with pCR and histological grade :\n")),append=TRUE,file=summaryFile)
        dataMatrixHist <- dataMatrix[which(!is.na(dataMatrix$hist_grade)),]
 
        cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
        
        linearModelsSumm[["pCR_withGrade_numPatients"]] <- nrow(dataMatrixHist)
        cat(gsub(" ","_",paste0("\nNumber_of_patients_for_pCR_withGrade_model: ",nrow(dataMatrixHist),"\n")),append=TRUE,file=summaryFile)
        linearModels[["pCR_logitWithRxGrade"]] <- tryCatch( glm(pCR~community+chemotherapyClass+anti_estrogen+anti_HER2+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
                                                            error = function(e) {
                                                              return(NA)
                                                            }
        )
        
        
      }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))==1 
               && length(unique(dataMatrix$anti_HER2))>1 ){
        
        cat(gsub(" ","_",paste0("\nFor pCR analysis, all anti_estrogen values the same: ",unique(unique(dataMatrix$anti_estrogen)),"\n")),
            append=TRUE,file=summaryFile)
        
        linearModels[["pCR_logitWithRx"]] <- tryCatch(glm(pCR~community+chemotherapyClass+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
                                                  #remove communities for ANOVA
           reducedModel <- tryCatch(glm(pCR~chemotherapyClass+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
         linearModels[["pCR_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["pCR_logitWithRx"]],test="Chisq")

        #now try with stage, grade.
        cat(gsub(" ","_",paste0("\nLinear analysis with pCR and histological grade :\n")),append=TRUE,file=summaryFile)
        dataMatrixHist <- dataMatrix[which(!is.na(dataMatrix$hist_grade)),]
 
        cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
        
        linearModelsSumm[["pCR_withGrade_numPatients"]] <- nrow(dataMatrixHist)
        cat(gsub(" ","_",paste0("\nNumber_of_patients_for_pCR_withGrade_model: ",nrow(dataMatrixHist),"\n")),append=TRUE,file=summaryFile)
        linearModels[["pCR_logitWithRxGrade"]] <- tryCatch( glm(pCR~community+chemotherapyClass+anti_HER2+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
                                                            error = function(e) {
                                                              return(NA)
                                                            }
        )
        
        
        
        
        
      }else if(length(unique(dataMatrix$chemotherapyClass))==1 && length(unique(dataMatrix$anti_estrogen))>1 
               && length(unique(dataMatrix$anti_HER2))>1 ){
        
        
        cat(gsub(" ","_",paste0("\nFor pCR analysis, all chemo values the same: ",unique(dataMatrix$chemotherapyClass),"\n")),
            append=TRUE,file=summaryFile)
        
        linearModels[["pCR_logitWithRx"]] <- tryCatch(glm(pCR~community+anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
   #remove community for ANOVA
        reducedModel <- tryCatch(glm(pCR~anti_estrogen+anti_HER2,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
         linearModels[["pCR_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["pCR_logitWithRx"]],test="Chisq")

 
        #now try with stage, grade.
        cat(gsub(" ","_",paste0("\nLinear analysis with pCR and histological grade :\n")),append=TRUE,file=summaryFile)
        dataMatrixHist <- dataMatrix[which(!is.na(dataMatrix$hist_grade)),]
 
        cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
        
        linearModelsSumm[["pCR_withGrade_numPatients"]] <- nrow(dataMatrix)
        cat(gsub(" ","_",paste0("\nNumber_of_patients_for_pCR_withGrade_model: ",nrow(dataMatrixHist),"\n")),append=TRUE,file=summaryFile)
        linearModels[["pCR_logitWithRxGrade"]] <- tryCatch( glm(pCR~community+anti_estrogen+anti_HER2+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
                                                            error = function(e) {
                                                              return(NA)
                                                            }
        )
        
        
        
        
        
      }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))==1 
               && length(unique(dataMatrix$anti_HER2))==1 ){
        
        
        cat(gsub(" ","_",paste0("\nFor pCR analysis, the anti HER2 and estrogen variables have only one unique value each.\n")),append=TRUE,
            file=summaryFile)
        
        
        linearModels[["pCR_logitWithRx"]] <- tryCatch(glm(pCR~community+chemotherapyClass,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )

   #remove community for ANOVA
        reducedModel <- tryCatch(glm(pCR~chemotherapyClass,data=dataMatrix,family=binomial(link="logit")),
                                                      error = function(e) {
                                                        return(NA)
                                                      }
        )
         linearModels[["pCR_logitWithRx_ANOVA"]] <- anova(reducedModel, linearModels[["pCR_logitWithRx"]],test="Chisq")


        #now try with stage, grade.
        cat(gsub(" ","_",paste0("\nLinear analysis with pCR and histological grade :\n")),append=TRUE,file=summaryFile)
        dataMatrixHist <- dataMatrix[which(!is.na(dataMatrix$hist_grade)),]
 
        cat(gsub(" ","_",paste0("\nStudies used: ",unique(dataMatrixHist$studyNum),"\n")),append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nMeta-clusters used: ",unique(dataMatrixHist$community),"\n")),append=TRUE,file=summaryFile)
        
        linearModelsSumm[["pCR_withGrade_numPatients"]] <- nrow(dataMatrixHist)
        cat(gsub(" ","_",paste0("\nNumber_of_patients_for_pCR_withGrade_model: ",nrow(dataMatrixHist),"\n")),append=TRUE,file=summaryFile)
        linearModels[["pCR_logitWithRxGrade"]] <- tryCatch( glm(pCR~community+chemotherapyClass+hist_grade,data=dataMatrixHist,family=binomial(link="logit")),
                                                            error = function(e) {
                                                              return(NA)
                                                            }
        )
        
        
        
        
      }else{
        
        # don't run models
        
      }
      
      cat(gsub(" ","_",paste0("\nLinear analysis with pCR and treatment:\n")),append=TRUE,file=summaryFile)
      if(any(!is.na(linearModels[["pCR_logitWithRx"]]))){
        
        linearModelsSumm[["pCR_logitWithRx"]] <- summary( linearModels[["pCR_logitWithRx"]])
        textOut <- capture.output( linearModelsSumm[["pCR_logitWithRx"]]$coefficients)
        cat(textOut,sep="\n",append=TRUE,file=summaryFile) 
        linearModelsSumm[["pCR_WithRxMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
                textOut <- capture.output( linearModelsSumm [["pCR_logitWithRx"]])
        cat(textOut,sep="\n",append=TRUE,file=summaryFile)
        cat(gsub(" ","_",paste0("\nLinear analysis with pCR and treatment ANOVA for meta-clusters:\n")),append=TRUE,file=summaryFile)
        textOut <- capture.output(linearModels[["pCR_logitWithRx_ANOVA"]])
        cat(textOut,sep="\n",append=TRUE,file=summaryFile)
        
        predictor <- predict(linearModels[["pCR_logitWithRx"]],type="response")
        png(filename=paste0(saveDir,"/",experimentName,"_pCR_withRx_ROC_",Sys.Date(),".png"),width=1000,height=1000,res=200)
        roc_data <- roc(response=dataMatrix$pCR,predictor=predictor,plot=TRUE,na.rm=TRUE,auc=TRUE)
        dev.off()
          ROC_list[["pCR_logitWithRx"]] <- roc_data
        textOut <- capture.output(roc_data)
        cat("roc_info\n",textOut,sep="\n",append=TRUE,file=summaryFile)
        cat("\nAUC: ",roc_data$auc,"\n",append=TRUE,file=summaryFile)
        
        
        
      }else{
        
        cat(gsub(" ","_",paste0("\nRx with pCR analysis returned NA\n.")),append=TRUE,file=summaryFile)  
      }
      
      
      #now try with stage, grade.
      cat(gsub(" ","_",paste0("\nLinear analysis with pCR and histological grade :\n")),append=TRUE,file=summaryFile)
      cat(gsub(" ","_",paste0("\nUnique hist grade values :",unique(dataMatrixHist$hist_grade),"\n")),append=TRUE,file=summaryFile)

      if(any(!is.na(linearModels[["pCR_logitWithRxGrade"]]))){
        
        linearModelsSumm[["pCR_logitWithRxGrade"]] <- summary( linearModels[["pCR_logitWithRxGrade"]])
        textOut <- capture.output( linearModelsSumm[["pCR_logitWithRxGrade"]]$coefficients)
        cat(textOut,sep="\n",append=TRUE,file=summaryFile) 
        linearModelsSumm[["pCR_logitWithRxGrade"]] <- paste(unique(dataMatrixHist$community),collapse=",")
        
        predictor <- predict(linearModels[["pCR_logitWithRxGrade"]],type="response")
        png(filename=paste0(saveDir,"/",experimentName,"_pCR_withGrade_ROC_",Sys.Date(),".png"),width=1000,height=1000,res=200)
        roc_data <- roc(response=dataMatrixHist$pCR,predictor=predictor,plot=TRUE,na.rm=TRUE,auc=TRUE)
        dev.off()
             ROC_list[["pCR_logitWithRxGrade"]] <- roc_data
        textOut <- capture.output(roc_data)
        cat("roc_info\n",textOut,sep="\n",append=TRUE,file=summaryFile)
        cat("\nAUC: ",roc_data$auc,"\n",append=TRUE,file=summaryFile)
        
      }else{
        
        cat(gsub(" ","_",paste0("\nHist grade with pCR analysis returned NA.")),append=TRUE,file=summaryFile)  
      }

  output <- list(sampleClustCommPhenoData=sampleClustCommPhenoData,ROC_data=ROC_list,
                 linearModelSummaries=linearModelsSumm,linearSurvModels=linearModels)
  return(output)
}


