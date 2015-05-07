
library("ggplot2")
library("RColorBrewer")
library("plyr")
library("survival")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_geneExprProcess.R")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_communityDetection.R")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_visualization.R")
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_metaClusterAnalysis.R")


metaFeaturesAnalysisWrapper <- function(metaFeatures,esets,CoINcIDE_output , clusterCoINcIDE_output,
    meanEdgePairPvalueThresh = .01,indEdgePvalueThresh = .05, minTrueSimilThresh = .4, maxTrueSimilThresh = Inf,
    clustSizeThresh = 5,saveDir = "/home/kplaney/ovarian_analysis/",experimentName = "ovarian_2000F",networkColors = "Set3",
    commMethod = "edgeBetween", minNumUniqueStudiesPerCommunity=4, nodePlotSize=10,nodeFontSize=.7,ES_thresh = .5,eset_featureDataFieldName="gene",
    survivalAnalysis=TRUE,outcomesVarBinary="vital_status",outcomesVarCont = "days_to_death",
    CutoffPointYears=5, eset_uniquePatientID="unique_patient_ID", ovarian=TRUE, fisherTestVariables = c("histological_type","tumorstage","recurrence_status","grade"),
    fisherTestVariableLegendNames=fisherTestVariables,fisherTestVariableTitleNames=fisherTestVariables
    ){
  
  
 dir.create(paste0(saveDir,"/",experimentName,"_",Sys.Date()),showWarnings=TRUE);
  saveDir <- paste0(saveDir,"/",experimentName,"_",Sys.Date())
  
  outputFile <- paste0(saveDir,"/",experimentName,"_outMessages_",Sys.Date(),".txt")
  
  inputVariablesDF = CoINcIDE_output$inputVariablesDF
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
  cat(paste("CoINcIDE input variables used for this analysis: ",collapse="_"),
      append=TRUE,file=outputFile)
  textOut <-capture.output(inputVariablesDF)
  cat("\n",textOut,sep="\n",append=TRUE,file=outputFile)
  
  
  computeTrueSimilOutput = CoINcIDE_output$computeTrueSimilOutput
  pvalueMatrix = CoINcIDE_output$pvalueMatrix
  clustIndexMatrix = CoINcIDE_output$clustIndexMatrix

    #now format just as a list of data matrices.
    dataMatrixList <- exprSetListToMatrixList(esets,featureDataFieldName=eset_featureDataFieldName)
    
    names(dataMatrixList) <- names(esets)
    ##do for each 200,500,1000,2000 (load different metaFeatures RData object each time.)
    #remove datasets with too many missing top gene features
    if(length(metaFeatures$datasetListIndicesToRemove)>0){
      
      dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
      
    }

  
  cat(paste("\nStarted with  ",length(dataMatrixList), " datasets\n:",collapse="_"),
      append=TRUE,file=outputFile)
  textOut <- capture.output(names(dataMatrixList))
  cat("\n",textOut,sep="\n",append=TRUE,file=outputFile)
  
  #need mapping later for pheno data.
  if(length(dataMatrixList) != length(na.omit(match(names(dataMatrixList),names(esets))))){
    
    stop("old to new indexing for esets and data matrixlist not matching up.")
  }
  origToNewIndexMap <- cbind(1:length(dataMatrixList),na.omit(match(names(dataMatrixList),names(esets))),names(esets)[na.omit(match(names(dataMatrixList),names(esets)))])
  tmp <- origToNewIndexMap[ ,c(1,3)]
  colnames(tmp) <- c("datasetNumber","datasetName")
  write.table(tmp,file=paste0(saveDir,"/",experimentName,"_datasetNameNumberKey.txt"),quote=FALSE,row.names=TRUE,col.names=FALSE)
  cat(paste("\nAssuming we're using k-means consensus PACR to choose K\n.",collapse="_"),file=outputFile,append=TRUE)
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
  
  fractFeatIntersectThresh <- inputVariablesDF$fractFeatIntersectThresh
  numFeatIntersectThresh <- inputVariablesDF$numFeatIntersectThresh
  clustSizeFractThresh <- inputVariablesDF$clustSizeFractThresh
  
  cat(paste("\nThere were ",nrow(clustIndexMatrix), " total input clusters from ",length(unique(clustIndexMatrix[,2])), " studies",collapse="_"),append=TRUE,file=outputFile)
  cat(paste("\nThe total number of input features was ",length(metaFeatures$finalFeatures),collapse="_"),append=TRUE,file=outputFile)
  cat(paste("\nAcross the entire square (nonsymmetric) p-value matrix, there are ",length(which(pvalueMatrix<=.1))," pvalues less than or equal to .1",collapse="_"),append=TRUE,file=outputFile)
  cat(paste("\nAcross the entire square (nonsymmetric) p-value matrix, there are ",length(which(pvalueMatrix<=.05))," pvalues less than or equal to .05",collapse="_"),append=TRUE,file=outputFile)
  cat(paste("\nAcross the entire square (nonsymmetric) p-value matrix, there are ",length(which(pvalueMatrix<=.01))," pvalues less than or equal to .01",collapse="_"),append=TRUE,file=outputFile)
  cat(paste("\nAcross the entire square (symmetric) similarity matrix, there are ",length(which(computeTrueSimilOutput$similValueMatrix>=minTrueSimilThresh))/2,
      " similarities greater than or equal to ",minTrueSimilThresh,"not double-counting.",collapse="_"),append=TRUE,file=outputFile)
  
  finalEdgeInfo <- assignFinalEdges(computeTrueSimilOutput=computeTrueSimilOutput,pvalueMatrix=pvalueMatrix,indEdgePvalueThresh=indEdgePvalueThresh,
                                    meanEdgePairPvalueThresh=meanEdgePairPvalueThresh,
                                    minTrueSimilThresh=minTrueSimilThresh,maxTrueSimilThresh=maxTrueSimilThresh,
                                    fractFeatIntersectThresh=fractFeatIntersectThresh,numFeatIntersectThresh=numFeatIntersectThresh ,
                                    clustSizeThresh=clustSizeThresh, clustSizeFractThresh= clustSizeFractThresh,saveDir=saveDir,fileTag="CoINcIDE_edges_",
                                    restrictEdges=FALSE
  )
  
  commInfo <- findCommunities(edgeMatrix=finalEdgeInfo$filterEdgeOutput$edgeMatrix,edgeWeightMatrix=finalEdgeInfo$filterEdgeOutput$edgeWeightMatrix,
                              clustIndexMatrix=CoINcIDE_output$clustIndexMatrix,fileTag=experimentName,
                              saveDir=saveDir,minNumUniqueStudiesPerCommunity=minNumUniqueStudiesPerCommunity,experimentName=experimentName,
                              commMethod=commMethod,
                              makePlots=TRUE,saveGraphData=TRUE,plotToScreen=FALSE)
  
  aggregateData <- returnSampleMemberMatrix(clustSampleIndexList=clustSampleIndexList,
                                            dataMatrixList=dataMatrixList,communityInfo=commInfo)
  clustSizes <- table(aggregateData$sampleClustCommKey$globalClustNum)
  textOut <- capture.output(table(aggregateData$sampleClustCommKey$community,exclude=FALSE))
  cat(paste("\nMeta-cluster sizes:\n",collapse="_"),textOut,sep="\n",
      append=TRUE,file=paste0(saveDir,"/",experimentName,"_metaclusterSizes.txt"))
  #add clusterSizes to node attributes
  commInfo$attrDF$size <- clustSizes[na.omit(match(commInfo$attrDF$clust,names(clustSizes)))]
 textOut <- capture.output(table(aggregateData$sampleClustCommKey$community,aggregateData$sampleClustCommKey$studyNum,exclude=FALSE))
 cat(paste("\nMeta-cluster number of samples by studies:\n",collapse="_"),textOut,sep="\n",append=TRUE,file=paste0(saveDir,"/",experimentName,"_metaclusterNumSamplesByStudies.txt"))
 
  networkStats <- advancedNetworkPlots(communityMembership=commInfo,
                                       brewPal = networkColors,
                                       saveDir=saveDir,clustIndexMatrix=clustIndexMatrix,
                                       plotToScreen=FALSE,experimentName=experimentName)$network_stats
  
  
  write.table(networkStats,quote=FALSE,file=paste0(saveDir,"/",experimentName,"_networkStats_",Sys.Date(),".txt"))
  #save variable
  sampleClustCommKey<-aggregateData$sampleClustCommKey
  binInfo <- binarizeMetaclustStudyStatus(aggregateData$sampleClustCommKey)
  
  message("Running effect size analysis")
  ES_out <- computeMetaclustEffectSizes(metaClustSampleNames=binInfo$metaClustSampleNames,dataMatrixList=dataMatrixList,featureNames=metaFeatures$finalFeatures,minOtherClass=5,
                                        computeWilcoxon=FALSE)
  sigGenes <- selectMetaclustSigGenes(computeMetaclustEffectSizesOutput=ES_out,qvalueThresh=1,
                                      ESthresh=ES_thresh,includeWilcoxon=FALSE)
  
  summaryESPos <- summarizePosESMetaclustGenes(selectMetaclustSigGenesOut=sigGenes,computeMetaclustEffectSizesOutput=ES_out)
    #write to table:
  write.table(summaryESPos,file=paste0(saveDir,"/",experimentName,"_summaryGenes_ESpos_thresh_",ES_thresh,"_",Sys.Date(),".txt"),row.names=TRUE,quote=FALSE,col.names=TRUE)
  for(g in 1:length(sigGenes$sigMetaclustGenesES_pos)){
    #what can't remove indices??
    tmp <- as.matrix(unlist(sigGenes$sigMetaclustGenesES_pos[[g]]))
    write.table(tmp,file=paste0(saveDir,"/",experimentName,"_genesESpos_thresh_",ES_thresh,"_community_",g,".txt"),row.names=FALSE,col.names=FALSE,quote=FALSE)
    
  }
  
  save(summaryESPos,file=paste0(saveDir,"/",experimentName,"_ES_genesWithThresh_ES_",ES_thresh,".RData.gzip"),compress="gzip")
  message("running GSEA")
  GSEA_out <- list()
  for(c in 1:length(sigGenes$sigMetaclustGenesES_pos)){
    
    testGeneSet <- sigGenes$sigMetaclustGenesES_pos[[c]]
    cat(paste("\n",length(testGeneSet)," genes for community ", c, " with >= ES of ",ES_thresh,"\n",collapse="_"),
              append=TRUE,file=outputFile)
    GSEA_out[[c]] <- GSEA(testGeneVector=testGeneSet,refGeneLists=NULL,method=c("hypergeometric"),
                          refGeneListDir="/home/data/MSigDB/GSEA_base_MSigDB_lists_merged.RData.gzip")
    
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
    cat(paste("\n",length(tmp)," unique signficant GSEA ref gene lists community ", g, " with ES of >=",ES_thresh,"\n",collapse="_"),
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
    cat(paste("\n",length(tmp)," overall (not necessarily unique) signficant GSEA ref gene lists community ", g, " with >= ES of ",ES_thresh,"\n",collapse="_"),
              append=TRUE,file=outputFile)
    
    
  }
  save(GSEA_out_sig,file=paste0(saveDir,"/",experimentName,"_GSEA_out_allSigLists_forEachMetaCluster.RData.gzip"),compress="gzip")
  
  
  for(g in 1:length(GSEA_out_sig)){
    #what can't remove indices??
    write.table(as.character(unlist(GSEA_out_sig[[g]])),file=paste0(saveDir,"/",experimentName,"_GSEA_allSig_community_",g,".txt"),row.names=TRUE,quote=FALSE,col.names=FALSE)
    
  }
  #cat(GSEA_out_sig,sep="\n",file=paste0(saveDir,"/",experimentName,"_GSEA_allSig.txt"),append=FALSE)
  

  ###now on to survival analysis
  result <- list()
  clinicalTables <- list()
  linearModels <- list()
 linearModelsSumm <- list()
  if(survivalAnalysis){
    

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

    
    if(!ovarian){
      #need to add pam50 centroids
      
      load("/home/kplaney/breast_analysis/pam50FullAndShort_subtypeDF.RData.gzip")
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
    sampleClustCommPhenoDataOrig <- sampleClustCommPhenoData

    groupingTerm="community"
    
    #only take samples with the groupingTerm you're looking at.
    sampleClustCommPhenoData <- sampleClustCommPhenoData[which(!is.na(sampleClustCommPhenoData[, groupingTerm])), ]
    #remove samples with NA values.
    groupings <- sampleClustCommPhenoData[, groupingTerm]
    
    cat(paste("\n",length(groupings), " patients included across all final meta-clusters.\n",collapse="_"),
        append=TRUE,file=outputFile)
    
    expName <- gsub("_"," ",experimentName)
    for(f in 1:length(fisherTestVariables)){
      
      if(!fisherTestVariables[f] %in% colnames(sampleClustCommPhenoData)){
        #move on to next loop: don't stop code way at bottom.
        message("Clinical variable ", fisherTestVariables[f] , " not found in pheno data.")
        next;
      }
      
      clinicalTables[[fisherTestVariables[f]]] <- table(sampleClustCommPhenoData[,groupingTerm],sampleClustCommPhenoData[,fisherTestVariables[f]],exclude=FALSE)
      textOut <- capture.output(table(sampleClustCommPhenoData[,groupingTerm],sampleClustCommPhenoData[,fisherTestVariables[f]],exclude=FALSE))
      cat(textOut,sep="\n",file=paste0(saveDir,"/",experimentName,"_tableStats_",fisherTestVariables[f],"_",Sys.Date(),".txt"),append=FALSE)
      
      dataF <- data.frame(sampleClustCommPhenoData[,groupingTerm],
                          sampleClustCommPhenoData[,fisherTestVariables[f]],
                          stringsAsFactors=FALSE)
      colnames(dataF) <- c("meta_cluster","clinVar")
      
    if(any(!is.na(dataF[,2]))){
      #if any are ""
      sampleClustCommPhenoData[which(sampleClustCommPhenoData[,fisherTestVariables[f]]==""),
                               fisherTestVariables[f]] <- NA
      
      #pam50 colors:  t(data.frame(c("#970eec","#ec0e58","#0e58ec","#7d7d7d","#0eec97")))
      #(and name for subtypes-look up order in old clust robust code.)
      
      if(fisherTestVariables[f] != "subtype" &&  fisherTestVariables[f] != "subtype_short"){
        
      #no scale_fill_manual(values = variableColorMatrix)+
        #fill must be a factor.
        plotG <-    ggplot(dataF,aes(clinVar,fill=as.factor(clinVar)))+geom_bar() + facet_grid(.~meta_cluster,scales="free_x")+
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
        
        png(filename=paste0(saveDir,"/",experimentName,"_",fisherTestVariables[f],"_breakdowns_",Sys.Date(),".png"),width=1000,height=1000,res=160)
        
        plot(plotG)
        
        dev.off()
        
        plotGStacked <- ggplot(dataF,aes(x=meta_cluster,fill=as.factor(clinVar)))+geom_bar() + 
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
        
        png(filename=paste0(saveDir,"/",experimentName,"_",fisherTestVariables[f],"_breakdowns_stacked_",Sys.Date(),".png"),width=1000,height=1000,res=160)
        
        plot( plotGStacked)
        
        dev.off()
        
        
      }else{
   
        variableColorMatrix <- c("#970eec","#ec0e58","#0e58ec","#7d7d7d","#0eec97")
        names(variableColorMatrix) <- c("LumB","LumA","Her2","Normal","Basal")
      
        plotG <-    ggplot(dataF,aes(clinVar,fill=as.factor(clinVar)))+geom_bar() + facet_grid(.~meta_cluster,scales="free_x")+
          labs(y="Number of samples", fill=paste0("",fisherTestVariableLegendNames[f],""),
               title=paste0(fisherTestVariableTitleNames[f],
                            " by meta-cluster \nfor ",expName))+
          theme(axis.text.x = element_text(colour = "black",size=8,angle=45,vjust=1,hjust=1),
                axis.title=element_blank())+
          theme(legend.title=element_text(size=12))+ scale_fill_manual(values = variableColorMatrix)+  theme(strip.text.x = element_text(size = 17, colour = "black",face="bold"))+
          theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
          theme(axis.text.x = element_text(colour = "black",size=18,angle=45,vjust=1,hjust=1),axis.title.x= element_blank(),
                axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1))+
          theme(panel.grid.major = element_line(colour = 0),panel.grid.minor = element_line(colour = 0))+
          theme(plot.title=element_text(colour="black",size=20,vjust=1))
          
        png(filename=paste0(saveDir,"/",experimentName,"_",fisherTestVariables[f],"_breakdowns_",Sys.Date(),".png"),width=1000,height=1000,res=160)
        
        plot( plotG)
        
        dev.off()
        
        plotGStacked <- ggplot(dataF,aes(x=meta_cluster,fill=as.factor(clinVar)))+geom_bar() + 
          labs(y="Number of samples", fill=paste0("",fisherTestVariableLegendNames[f],""),
               title=paste0(fisherTestVariableTitleNames[f],
                            " by meta-cluster \nfor ",expName,"\n"))+
          theme(axis.text.x = element_text(colour = "black",size=8,angle=45,vjust=1,hjust=1),
                axis.title=element_blank())+
          theme(legend.title=element_text(size=12))+ scale_fill_manual(values = variableColorMatrix)+  theme(strip.text.x = element_text(size = 17, colour = "black",face="bold"))+
          theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
          theme(axis.text.x = element_text(colour = "black",size=18,angle=45,vjust=1,hjust=1),axis.title.x= element_blank(),
                axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1))+
          theme(panel.grid.major = element_line(colour = 0),panel.grid.minor = element_line(colour = 0))+
          theme(plot.title=element_text(colour="black",size=20,vjust=1))
        
        png(filename=paste0(saveDir,"/",experimentName,"_",fisherTestVariables[f],"_breakdowns_stacked_",Sys.Date(),".png"),width=1000,height=1000,res=160)
        
        plot( plotGStacked)
        
        dev.off()
      
      }

      #chi-squared test. fisher's test doesn't work here...it looks like bar charts, plots will be more informative.
      tmp <- factor(sampleClustCommPhenoData[,groupingTerm])[which(!is.na(sampleClustCommPhenoData[,fisherTestVariables[f]]))]
      tmp1 <- factor(sampleClustCommPhenoData[,fisherTestVariables[f]][which(!is.na(sampleClustCommPhenoData[,fisherTestVariables[f]]))])
         if(length(tmp1)==0 || length(unique(tmp1)==1)){
        
        cat(paste("\nClinical variable ", fisherTestVariables[f], " had only NAs or all one value so not analyzing it.\n",collapse="_"),
            append=TRUE,file=outputFile)
        next;
      }
      #need to simulate p-value other gets wonky.
      result[[paste0("fisher_",fisherTestVariables[f],"_pvalue")]] <- chisq.test(tmp,tmp1,simulate.p.value=TRUE)$p.value
    
  
  #end of if is NA
      }else{
        cat(paste("\nVariable ",fisherTestVariables[f], " is all NAs.\n",collapse="_"),append=TRUE,file=outputFile)
        result[[paste0("fisher_",fisherTestVariables[f],"_pvalue")]]  <- NA
      }
  #end of looping over fisher variables.
  }
    
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
      #creating the survival objects with the time and censoring variables
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
      result$cox_overall_logTestPvalue <- sumCoxfit_overall$logtest["pvalue"]
      result$cox_overall_nevent <- sumCoxfit_overall$nevent
      result$cox_overall_sctestPvalue <- sumCoxfit_overall$sctest["pvalue"]
      result$cox_overall_rsq <- sumCoxfit_overall$rsq["rsq"]
      result$cox_overall_waldPvalue <- sumCoxfit_overall$waldtest["pvalue"]
      #plot(cox.zph(coxfit))
      kmfit=survdiff(OverallSurvival ~ groupings)
      #message("kaplan meier p-value for overall survival:")
      result$km_overall_pvalue <- 1 - pchisq(kmfit$chisq, length(kmfit$n) - 1)
      
    #  message("calculating the sign of the survival relationship")
      mfit=survfit(OverallSurvival ~ groupings)
      
      png(filename=paste0(saveDir,"/",experimentName,"_kaplanMeier_OS_",Sys.Date(),".png"),width=1000,height=1000,res=160)
      
      plot(mfit,main=paste0("Overall survival for ",experimentName))
      
      dev.off()
      
      plot(mfit,main="overall survival")
      
      mfitCut <- survfit(OverallSurvivalCutoff ~ groupings)
      
      png(filename=paste0(saveDir,"/",experimentName,"_kaplanMeier_OS",CutoffPointYears,"years_",Sys.Date(),".png"),width=1000,height=1000,res=160)
      
      plot(mfitCut,main=paste0("Overall survival for ",experimentName, "\nwith cutoff at ",
                               CutoffPointYears," years"))
      
      dev.off()
      
      coxfit=coxph(OverallSurvivalCutoff~groupings, data= SurvivalCutoff)
      message("coxfit summary for survival cutoff at ",CutoffPointYears,"years:")
    
      #plot(cox.zph(coxfit))
    sumCoxfit_cut  <- summary(coxfit)
    result$cox_cut_logTestPvalue <- sumCoxfit_cut$logtest["pvalue"]
    result$cox_cut_nevent <- sumCoxfit_cut$nevent
    result$cox_cut_sctestPvalue <- sumCoxfit_cut$sctest["pvalue"]
    result$cox_cut_rsq <- sumCoxfit_cut$rsq["rsq"]
    result$cox_cut_waldPvalue <- sumCoxfit_cut$waldtest["pvalue"]
    
      kmfit=survdiff(OverallSurvivalCutoff ~ groupings)
    
      result$km_pv_cutoff <- 1 - pchisq(kmfit$chisq, length(kmfit$n) - 1) 
    
    
    tmp <- as.matrix(unlist(result))
    write.table(tmp,file=paste0(saveDir,"/",experimentName,"_survivalAnalysis_",Sys.Date(),".txt"),col.names=FALSE)
    
    
   #end of survival analysis.
  }else{
 
    #message("running survival data on breast")
    #OS: not enough variables.
 
    survFile <- paste0(saveDir,"/",experimentName,"_linearSurvivalStats_",Sys.Date(),".txt")
    cat(paste("\nLinear analysis with only RFS:\n",collapse="_"),append=TRUE,file=survFile)
    
    dataMatrix <- sampleClustCommPhenoData[which(!is.na((sampleClustCommPhenoData$RFS))), ]
    cat(paste("\nMeta-clusters used: ",paste(unique(dataMatrix$community),collapse=","),"\n",collapse="_"),append=TRUE,file=survFile)
    
    #strongly predicts RFS.
    linearModels[["RFS_linearAlone"]] <- lm(RFS~community,data=dataMatrix)
     linearModelsSumm[["RFS_linearAlone"]] <- summary(  linearModels[["RFS_linearAlone"]] ) 
    textOut <- capture.output( linearModelsSumm [["RFS_linearAlone"]]$coefficients)
    cat(textOut,sep="\n",append=TRUE,file=survFile)
     linearModelsSumm[["RFS_numPatients"]] <- nrow(dataMatrix)
     linearModelsSumm[["RFS_aloneMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
    ##with therapies
    #interesting: community becomes insignificant if add in therapies now: this could be
    #because community is most likely correlated with treatments; it's certainly correlated with pam50 status.
    #low Rsquared.
    cat(paste("\nLinear analysis with RFS and treatment:\n",collapse="_"),append=TRUE,file=survFile)
    dataMatrix <- dataMatrix[which(!is.na(dataMatrix$chemotherapyClass)), ]
    dataMatrix <- dataMatrix[which(!is.na(dataMatrix$anti_estrogen)), ]
    dataMatrix <- dataMatrix[which(!is.na(dataMatrix$anti_HER2)), ]
    dataMatrix <- data.frame(dataMatrix)
    cat(paste("\nMeta-clusters used: ",paste(unique(dataMatrix$community),collapse=","),"\n",collapse="_"),append=TRUE,file=survFile)
    
     linearModelsSumm [["RFS_withRx_numPatients"]] <- nrow(dataMatrix)
    
    if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))>1 
       && length(unique(dataMatrix$anti_HER2))>1 ){
         
         linearModels[["RFS_linearWithRx"]] <- tryCatch(lm(RFS~community+chemotherapyClass+anti_estrogen+anti_HER2,data=dataMatrix),
                                                        error = function(e) {
                                                          return(NA)
                                                        }
         )

        if(!is.na(linearModels[["RFS_linearWithRx"]])){
          linearModelsSumm[["RFS_linearWithRx"]] <- summary( linearModels[["RFS_linearWithRx"]])
         textOut <- capture.output( linearModelsSumm[["RFS_linearWithRx"]]$coefficients)
         cat(textOut,sep="\n",append=TRUE,file=survFile)
         
          linearModelsSumm[["RFS_WithRxMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
         
        }else{
          
          cat(paste("\nRFS analysis with therapies returned NA.\n",collapse="_"),append=TRUE,file=survFile)
        }
         
       }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))>1 
                && length(unique(dataMatrix$anti_HER2))==1 ){
         
         
         cat(paste("\nFor RFS analysis, all anti_HER2 values the same: ",unique(unique(dataMatrix$anti_HER2)),"\n",collapse="_"),
             append=TRUE,file=survFile)
         linearModels[["RFS_linearWithRx"]] <- lm(RFS~community+chemotherapyClass+anti_estrogen,data=dataMatrix)
          linearModelsSumm[["RFS_linearWithRx"]] <- summary( linearModels[["RFS_linearWithRx"]])
         textOut <- capture.output( linearModelsSumm[["RFS_linearWithRx"]]$coefficients)
         cat(textOut,sep="\n",append=TRUE,file=survFile)
         
          linearModelsSumm[["RFS_WithRxMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
         
         
       }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))==1 
                && length(unique(dataMatrix$anti_HER2))>1 ){
         
         cat(paste("\nFor RFS analysis, all anti_estrogen values the same: ",unique(unique(dataMatrix$anti_estrogen)),"\n",collapse="_"),
             append=TRUE,file=survFile)
         
         linearModels[["RFS_linearWithRx"]] <- lm(RFS~community+chemotherapyClass+anti_HER2,data=dataMatrix)
          linearModelsSumm[["RFS_linearWithRx"]] <- summary( linearModels[["RFS_linearWithRx"]])
         textOut <- capture.output( linearModelsSumm[["RFS_linearWithRx"]]$coefficients)
         cat(textOut,sep="\n",append=TRUE,file=survFile)
         
          linearModelsSumm[["RFS_WithRxMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
         
         
         
       }else if(length(unique(dataMatrix$chemotherapyClass))==1 && length(unique(dataMatrix$anti_estrogen))>1 
                && length(unique(dataMatrix$anti_HER2))>1 ){
         
         
         cat(paste("\nFor RFS analysis, all chemo values the same: ",unique(unique(dataMatrix$chemotherapyClass)),"\n",collapse="_"),
             append=TRUE,file=survFile)
         
         linearModels[["RFS_linearWithRx"]] <- lm(RFS~community+anti_estrogen+anti_HER2,data=dataMatrix)
          linearModelsSumm[["RFS_linearWithRx"]] <- summary( linearModels[["RFS_linearWithRx"]])
         textOut <- capture.output( linearModelsSumm[["RFS_linearWithRx"]]$coefficients)
         cat(textOut,sep="\n",append=TRUE,file=survFile)
         
          linearModelsSumm[["RFS_WithRxMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
         
         
       }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))==1 
                && length(unique(dataMatrix$anti_HER2))==1 ){
         
         
         cat(paste("\nFor RFS analysis, the anti HER2 and estrogen variables have only one unique value each.\n",collapse="_"),append=TRUE,
             file=survFile)
         
         
         linearModels[["RFS_linearWithRx"]] <- lm(RFS~community+chemotherapyClass,data=dataMatrix)
          linearModelsSumm[["RFS_linearWithRx"]] <- summary( linearModels[["RFS_linearWithRx"]])
         textOut <- capture.output( linearModelsSumm[["RFS_linearWithRx"]]$coefficients)
         cat(textOut,sep="\n",append=TRUE,file=survFile)
         
          linearModelsSumm[["RFS_WithRxMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
         
         
       }else{
         
         # don't run models
         
       }



  #now try with stage, grade.
  cat(paste("\nLinear analysis with RFS and histological grade :\n",collapse="_"),append=TRUE,file=survFile)
  dataMatrix <- sampleClustCommPhenoData[which(!is.na((sampleClustCommPhenoData$RFS))), ]
  dataMatrix <- sampleClustCommPhenoData[which(!is.na((sampleClustCommPhenoData$hist_grade))),]
  dataMatrix <- data.frame(dataMatrix)
  cat(paste("\nMeta-clusters used: ",paste(unique(dataMatrix$community),collapse=","),"\n",collapse="_"),append=TRUE,file=survFile)
  
   linearModelsSumm[["RFS_withGrade_numPatients"]] <- nrow(dataMatrix)
  linearModels[["RFS_linearWithGrade"]] <- tryCatch( lm(RFS~community+hist_grade,data=dataMatrix),
                                                     error = function(e) {
                                                       return(NA)
                                                     }
  )
  if(!is.na(linearModels[["RFS_linearWithGrade"]])){
    linearModelsSumm[["RFS_linearWithGrade"]] <- summary( linearModels[["RFS_linearWithGrade"]])
    textOut <- capture.output( linearModelsSumm[["RFS_linearWithGrade"]]$coefficients)
    cat(textOut,sep="\n",append=TRUE,file=survFile) 
    linearModelsSumm[["RFS_WithGradeMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
    
  }else{
    
    cat(paste("\nHist grade with RFS analysis returned NA.\n",collapse="_"),append=TRUE,file=survFile)  
  }

###DFS
   cat(paste("\nLinear analysis with only DFS:\n",collapse="_"),append=TRUE,file=survFile)
   
   dataMatrix <- sampleClustCommPhenoData[which(!is.na((sampleClustCommPhenoData$DFS))), ]
   cat(paste("\nMeta-clusters used: ",paste(unique(dataMatrix$community),collapse=","),"\n",collapse="_"),append=TRUE,file=survFile)
   
   #strongly predicts DFS.
   linearModels[["DFS_linearAlone"]] <- lm(DFS~community,data=dataMatrix)
    linearModelsSumm[["DFS_linearAlone"]] <- summary(  linearModels[["DFS_linearAlone"]] ) 
   textOut <- capture.output( linearModelsSumm[["DFS_linearAlone"]]$coefficients)
   cat(textOut,sep="\n",append=TRUE,file=survFile)
    linearModelsSumm[["DFS_numPatients"]] <- nrow(dataMatrix)
    linearModelsSumm[["DFS_aloneMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
   ##with therapies
   #interesting: community becomes insignificant if add in therapies now: this could be
   #because community is most likely correlated with treatments; it's certainly correlated with pam50 status.
   #low Rsquared.
   cat(paste("\nLinear analysis with DFS and treatment:\n",collapse="_"),append=TRUE,file=survFile)
   dataMatrix <- dataMatrix[which(!is.na(dataMatrix$chemotherapyClass)), ]
   dataMatrix <- dataMatrix[which(!is.na(dataMatrix$anti_estrogen)), ]
   dataMatrix <- dataMatrix[which(!is.na(dataMatrix$anti_HER2)), ]
   dataMatrix <- data.frame(dataMatrix)
   cat(paste("\nMeta-clusters used: ",paste(unique(dataMatrix$community),collapse=","),"\n",collapse="_"),append=TRUE,file=survFile)
   
    linearModelsSumm[["DFS_withRx_numPatients"]] <- nrow(dataMatrix)
   
   if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))>1 && 
        length(unique(dataMatrix$anti_HER2))>1 ){
        
   linearModels[["DFS_linearWithRx"]] <- tryCatch(lm(DFS~community+chemotherapyClass+anti_estrogen+anti_HER2,data=dataMatrix),
                                                  error = function(e) {
                                                    return(NA)
                                                  }
   )
   
  if(!is.na(  linearModels[["DFS_linearWithRx"]] )){
    linearModelsSumm[["DFS_linearWithRx"]] <- summary( linearModels[["DFS_linearWithRx"]])
   textOut <- capture.output( linearModelsSumm[["DFS_linearWithRx"]]$coefficients)
   cat(textOut,sep="\n",append=TRUE,file=survFile)
   
    linearModelsSumm[["DFS_WithRxMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
   
  }else{
    
    cat(paste("\nDFS with therapies analysis returned NA.\n",collapse="_"),append=TRUE,file=survFile)
    
  }
   
      }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))>1 &&
                length(unique(dataMatrix$anti_HER2))==1 ){
        
        
        cat(paste("\nFor DFS analysis, all anti_HER2 values the same: ",unique(unique(dataMatrix$anti_HER2)),"\n",collapse="_"),
            append=TRUE,file=survFile)
        linearModels[["DFS_linearWithRx"]] <- lm(DFS~community+chemotherapyClass+anti_estrogen,data=dataMatrix)
         linearModelsSumm[["DFS_linearWithRx"]] <- summary( linearModels[["DFS_linearWithRx"]])
        textOut <- capture.output( linearModelsSumm[["DFS_linearWithRx"]]$coefficients)
        cat(textOut,sep="\n",append=TRUE,file=survFile)
        
         linearModelsSumm[["DFS_WithRxMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
        
        
      }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))==1 
               && length(unique(dataMatrix$anti_HER2))>1 ){
        
        cat(paste("\nFor DFS analysis, all anti_estrogen values the same: ",unique(unique(dataMatrix$anti_estrogen)),"\n",collapse="_"),
            append=TRUE,file=survFile)
        
        linearModels[["DFS_linearWithRx"]] <- lm(DFS~community+chemotherapyClass+anti_HER2,data=dataMatrix)
         linearModelsSumm[["DFS_linearWithRx"]] <- summary( linearModels[["DFS_linearWithRx"]])
        textOut <- capture.output( linearModelsSumm[["DFS_linearWithRx"]]$coefficients)
        cat(textOut,sep="\n",append=TRUE,file=survFile)
        
         linearModelsSumm[["DFS_WithRxMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
        
        
        
      }else if(length(unique(dataMatrix$chemotherapyClass))==1 && length(unique(dataMatrix$anti_estrogen))>1 
               && length(unique(dataMatrix$anti_HER2))>1 ){
        
        
        cat(paste("\nFor DFS analysis, all chemo values the same: ",unique(unique(dataMatrix$chemotherapyClass)),"\n",collapse="_"),
            append=TRUE,file=survFile)
        
        linearModels[["DFS_linearWithRx"]] <- lm(DFS~community+anti_estrogen+anti_HER2,data=dataMatrix)
         linearModelsSumm[["DFS_linearWithRx"]] <- summary( linearModels[["DFS_linearWithRx"]])
        textOut <- capture.output( linearModelsSumm[["DFS_linearWithRx"]]$coefficients)
        cat(textOut,sep="\n",append=TRUE,file=survFile)
        
         linearModelsSumm[["DFS_WithRxMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
        
      
      }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))==1 
               && length(unique(dataMatrix$anti_HER2))==1 ){
        
        
        cat(paste("\nFor DFS analysis, the anti HER2 and estrogen variables have only one unique value each.\n",collapse="_"),append=TRUE,
            file=survFile)

        
        linearModels[["DFS_linearWithRx"]] <- lm(DFS~community+chemotherapyClass,data=dataMatrix)
         linearModelsSumm[["DFS_linearWithRx"]] <- summary( linearModels[["DFS_linearWithRx"]])
        textOut <- capture.output( linearModelsSumm[["DFS_linearWithRx"]]$coefficients)
        cat(textOut,sep="\n",append=TRUE,file=survFile)
        
         linearModelsSumm[["DFS_WithRxMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
        
        
      }else{
        
       # don't run models

      }
   

   #now try with stage, grade.
   cat(paste("\nLinear analysis with DFS and histological grade :\n",collapse="_"),append=TRUE,file=survFile)
   dataMatrix <- sampleClustCommPhenoData[which(!is.na((sampleClustCommPhenoData$DFS))), ]
   dataMatrix <- sampleClustCommPhenoData[which(!is.na((sampleClustCommPhenoData$hist_grade))),]
  dataMatrix <- data.frame(dataMatrix)
   cat(paste("\nMeta-clusters used: ",paste(unique(dataMatrix$community),collapse=","),"\n",collapse="_"),append=TRUE,file=survFile)
   
    linearModelsSumm[["DFS_withGrade_numPatients"]] <- nrow(dataMatrix)

  linearModels[["DFS_linearWithGrade"]] <- tryCatch( lm(DFS~community+hist_grade,data=dataMatrix),
                                                     error = function(e) {
                                                       return(NA)
                                                     }
  )
  if(!is.na(linearModels[["DFS_linearWithGrade"]])){
    linearModelsSumm[["DFS_linearWithGrade"]] <- summary( linearModels[["DFS_linearWithGrade"]])
    textOut <- capture.output( linearModelsSumm[["DFS_linearWithGrade"]]$coefficients)
    cat(textOut,sep="\n",append=TRUE,file=survFile) 
    linearModelsSumm[["DFS_WithGradeMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
    
  }else{
    
    cat(paste("\nHist grade with DFS analysis returned NA\n.",collapse="_"),append=TRUE,file=survFile)  
  }


#####pCR
cat(paste("\nLinear analysis with only pCR:\n",collapse="_"),append=TRUE,file=survFile)

dataMatrix <- sampleClustCommPhenoData[which(!is.na((sampleClustCommPhenoData$pCR))), ]
cat(paste("\nMeta-clusters used: ",paste(unique(dataMatrix$community),collapse=","),"\n",collapse="_"),append=TRUE,file=survFile)

#strongly predicts pCR.
linearModels[["pCR_linearAlone"]] <- lm(pCR~community,data=dataMatrix)
 linearModelsSumm[["pCR_linearAlone"]] <- summary(  linearModels[["pCR_linearAlone"]] ) 
textOut <- capture.output( linearModelsSumm[["pCR_linearAlone"]]$coefficients)
cat(textOut,sep="\n",append=TRUE,file=survFile)
 linearModelsSumm[["pCR_numPatients"]] <- nrow(dataMatrix)
 linearModelsSumm[["pCR_aloneMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
##with therapies
#interesting: community becomes insignificant if add in therapies now: this could be
#because community is most likely correlated with treatments; it's certainly correlated with pam50 status.
#low Rsquared.
cat(paste("\nLinear analysis with pCR and treatment:\n",collapse="_"),append=TRUE,file=survFile)
dataMatrix <- dataMatrix[which(!is.na(dataMatrix$chemotherapyClass)), ]
dataMatrix <- dataMatrix[which(!is.na(dataMatrix$anti_estrogen)), ]
dataMatrix <- dataMatrix[which(!is.na(dataMatrix$anti_HER2)), ]
dataMatrix <- data.frame(dataMatrix)
cat(paste("\nMeta-clusters used: ",paste(unique(dataMatrix$community),collapse=","),"\n",collapse="_"),append=TRUE,file=survFile)

 linearModelsSumm[["pCR_withRx_numPatients"]] <- nrow(dataMatrix)

if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))>1 &&
    length(unique(dataMatrix$anti_HER2))>1 ){
     
     linearModels[["pCR_linearWithRx"]] <- tryCatch(lm(pCR~community+chemotherapyClass+anti_estrogen+anti_HER2,data=dataMatrix),
                                                    error = function(e) {
                                                      return(NA)
                                                    }
     )
     
     if(!is.na(linearModels[["pCR_linearWithRx"]])){
      linearModelsSumm[["pCR_linearWithRx"]] <- summary( linearModels[["pCR_linearWithRx"]])
     textOut <- capture.output( linearModelsSumm[["pCR_linearWithRx"]]$coefficients)
     cat(textOut,sep="\n",append=TRUE,file=survFile)
     
      linearModelsSumm[["pCR_WithRxMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
     
     }else{
       
       
       linearModelsSumm[["pCR_WithRxMetaClusters_used"]] <- NA
       cat(paste("\nDFS analysis with therapies gave factor errors.\n",collapse="_"),append=TRUE,file=survFile)
       
     }
     
   }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))>1 
            && length(unique(dataMatrix$anti_HER2))==1 ){
     
     
     cat(paste("\nFor pCR analysis, all anti_HER2 values the same: ",unique(unique(dataMatrix$anti_HER2)),"\n",collapse="_"),
         append=TRUE,file=survFile)
     linearModels[["pCR_linearWithRx"]] <- lm(pCR~community+chemotherapyClass+anti_estrogen,data=dataMatrix)
      linearModelsSumm[["pCR_linearWithRx"]] <- summary( linearModels[["pCR_linearWithRx"]])
     textOut <- capture.output( linearModelsSumm[["pCR_linearWithRx"]]$coefficients)
     cat(textOut,sep="\n",append=TRUE,file=survFile)
     
      linearModelsSumm[["pCR_WithRxMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
     
     
   }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))==1 
            && length(unique(dataMatrix$anti_HER2))>1 ){
     
     cat(paste("\nFor pCR analysis, all anti_estrogen values the same: ",unique(unique(dataMatrix$anti_estrogen)),"\n",collapse="_"),
         append=TRUE,file=survFile)
     
     linearModels[["pCR_linearWithRx"]] <- lm(pCR~community+chemotherapyClass+anti_HER2,data=dataMatrix)
      linearModelsSumm[["pCR_linearWithRx"]] <- summary( linearModels[["pCR_linearWithRx"]])
     textOut <- capture.output( linearModelsSumm[["pCR_linearWithRx"]]$coefficients)
     cat(textOut,sep="\n",append=TRUE,file=survFile)
     
      linearModelsSumm[["pCR_WithRxMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
     
     
     
   }else if(length(unique(dataMatrix$chemotherapyClass))==1 && length(unique(dataMatrix$anti_estrogen))>1 
            && length(unique(dataMatrix$anti_HER2))>1 ){
     
     
     cat(paste("\nFor pCR analysis, all chemo values the same: ",unique(unique(dataMatrix$chemotherapyClass)),"\n",collapse="_"),
         append=TRUE,file=survFile)
     
     linearModels[["pCR_linearWithRx"]] <- lm(pCR~community+anti_estrogen+anti_HER2,data=dataMatrix)
      linearModelsSumm[["pCR_linearWithRx"]] <- summary( linearModels[["pCR_linearWithRx"]])
     textOut <- capture.output( linearModelsSumm[["pCR_linearWithRx"]]$coefficients)
     cat(textOut,sep="\n",append=TRUE,file=survFile)
     
      linearModelsSumm[["pCR_WithRxMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
     
     
   }else if(length(unique(dataMatrix$chemotherapyClass))>1 && length(unique(dataMatrix$anti_estrogen))==1 
            && length(unique(dataMatrix$anti_HER2))==1 ){
     
     
     cat(paste("\nFor pCR analysis, the anti HER2 and estrogen variables have only one unique value each.\n",collapse="_"),append=TRUE,
         file=survFile)
     
     
     linearModels[["pCR_linearWithRx"]] <- lm(pCR~community+chemotherapyClass,data=dataMatrix)
      linearModelsSumm[["pCR_linearWithRx"]] <- summary( linearModels[["pCR_linearWithRx"]])
     textOut <- capture.output( linearModelsSumm[["pCR_linearWithRx"]]$coefficients)
     cat(textOut,sep="\n",append=TRUE,file=survFile)
     
      linearModelsSumm[["pCR_WithRxMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
     
     
   }else{
     
     # don't run models
     
   }


#now try with stage, grade.
cat(paste("\nLinear analysis with pCR and histological grade :\n",collapse="_"),append=TRUE,file=survFile)
dataMatrix <- sampleClustCommPhenoData[which(!is.na((sampleClustCommPhenoData$pCR))), ]
dataMatrix <- sampleClustCommPhenoData[which(!is.na((sampleClustCommPhenoData$hist_grade))),]
dataMatrix <- data.frame(dataMatrix)
cat(paste("\nMeta-clusters used: ",paste(unique(dataMatrix$community),collapse=","),"\n",collapse="_"),append=TRUE,file=survFile)

 linearModelsSumm[["pCR_withGrade_numPatients"]] <- nrow(dataMatrix)
linearModels[["pCR_linearWithGrade"]] <- tryCatch( lm(pCR~community+hist_grade,data=dataMatrix),
                                                  error = function(e) {
                                                    return(NA)
                                                  }
)
 if(!is.na(linearModels[["pCR_linearWithGrade"]])){
 linearModelsSumm[["pCR_linearWithGrade"]] <- summary( linearModels[["pCR_linearWithGrade"]])
textOut <- capture.output( linearModelsSumm[["pCR_linearWithGrade"]]$coefficients)
cat(textOut,sep="\n",append=TRUE,file=survFile) 
 linearModelsSumm[["pCR_WithGradeMetaClusters_used"]] <- paste(unique(dataMatrix$community),collapse=",")
   
  }else{
    
  cat(paste("\nHist grade with pCR analysis returned NA.",collapse="_"),append=TRUE,file=survFile)  
  }
  
#end of if not ovarian.
}
#end of if survival analysis.
  
}
    output <- list(sampleClustCommPhenoData=sampleClustCommPhenoDataOrig,aggregateData=aggregateData,commInfo=commInfo,finalEdgeInfo=finalEdgeInfo,CoINcIDE_computeEdgesObject=CoINcIDE_output,
                            networkStats=networkStats,sigGenes=sigGenes,GSEA_out_unique=GSEA_out_unique,survivalResults=result,
                   linearModelSummaries=linearModelsSumm,linearSurvModels=linearModels,
                   clustSampleIndexList=clustSampleIndexList,clustFeatureIndexList=clustFeatureIndexList,
                   clinicalTables=clinicalTables)
  
  return(output)
}