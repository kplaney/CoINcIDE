community_metaGeneAnalysis <- function(community_membership,clusterMatrixList,genomeSize=10000,pvalueThresh=.05,
                                       dataMatrices,dataM_sampleCol=FALSE,geneCommIntersect,geneCommUnion,detailedBiclustNames,logScale=TRUE){
warning("Assumes that in each community, the intersecting or union genes are always in the same order");
source("/home/kplaney/gitRepos/IGP_network/igp_network/clust_robust.R");

if(missing(geneCommIntersect) || missing(geneCommUnion)){
  
  if(missing(detailedBiclustNames)){
    
    stop("Please provide the detailedBiclustNames variable.");
    
  }
  
  geneSumm <- communityGeneSummary(community_membership=community_membership,biclust=clusterMatrixList,detailedBiclustNames=detailedBiclustNames,
                                   genomeSize=genomeSize,
                                   sourceDir="/home/kplaney/gitRepos/RNAseq_pipeline/",
                                   default_refGeneListDir="/home/data/MSigDB/GSEA_base_MSigDB_lists_merged.RData.gzip");
  
  geneCommIntersect <-  geneSumm$geneCommIntersect;
  geneCommUnion <-   geneSumm$geneCommUnion;

}

comm_exprList <- communityHeatmaps(community_membership=community_membership,
                                   biclust=clusterMatrixList,dataMatrices=dataMatrices,dataM_sampleCol=dataM_sampleCol,geneCommIntersect=geneCommIntersect,
                          geneCommUnion=geneCommUnion,
                          fullHeatmaps=FALSE);



sourceDir="/home/kplaney/gitRepos/RNAseq_pipeline/";
load("/home/data/MSigDB/GSEA_base_MSigDB_lists_merged.RData.gzip");
  

  setwd(sourceDir);
  source("AO_GSEA.R");
  cat("\nUsing default MSigDB lists: MSigDB_onco_symbols, MSigDB_CanPath_symbols,MSigDB_TFT_symbols,MSigDB_immun_symbols,
            and MSigDB_cancerNeigh_symbols.\n");

  refGeneLists <- GSEA_base_MSigDB_lists_merged;

fishersMethod = function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower=FALSE);
sigGenesList <- list();
library("rmeta");

for(c in 1:length(comm_exprList$heatmaps_expr_union)){
  sigGenesList[[c]] <- list();
  #"up vs down"
  if(all(is.null(names(comm_exprList$heatmaps_expr_union[[c]])))){
    
    names(comm_exprList$heatmaps_expr_union[[c]]) <- c(1:length(comm_exprList$heatmaps_expr_union[[c]]));
    
  }

  index <- 0;
  for(i in 1:length(comm_exprList$heatmaps_expr_union[[c]])){
  #just need one index that is not null
    if(!is.null(comm_exprList$heatmaps_expr_union[[c]][[i]])){
  
      index <- i;
    break;
  
  }
  
  }
  if(index!=0){
    
  
  for(e in 1:length(comm_exprList$heatmaps_expr_union[[c]])){
    
    if(!is.null(comm_exprList$heatmaps_expr_union[[c]][[e]])){
      
      #run SAMR for each gene within each dataset?
      output <- samr()
      for(g in 1:nrow(comm_exprList$heatmaps_expr_union[[c]][[e]])){

        
        sigGeneMatrixUp[e,g] <-    
        sigGeneMatrixDown[e,g] <-    
                  
      }
      
    }
    #end of loop e.(studies in this community.) 
  }
  
  for each meta-cluster:
  for(g in 1:allGenesTestedAcrossAllStudies){
  ##now after do that for all studies: take s.d., calculate ES (relative difference), and take average q-value
  
    #The summary estimate is a weighted average. If weights are specified they are used, otherwise the
    #reciprocal of the estimated variance is used.
    ES <- meta.summaries(d=d, se=se, method=c("random"),
                         logscale=logscale,
                         conf.level=0.95);
  #summarize wilcoxon rank comm_exprList using fisher's p-values: combine across all studes
  qvalue_summ[g] <- fishersMethod(SAMR_qvalues[,g]);
  gene_up[g] <- fraction that were deemed to be overexpressed by SAMR for all studies.
    
  }
}


#remove any NA studies - ie where entire study was the cluster,
  #as these will have NA values then.
  if(any(is.na(sigGeneMatrix_hedgeG))){
    
  tmp_hedgeG <- sigGeneMatrix_hedgeG[-which(is.na(sigGeneMatrix_hedgeG[,1])),];
  tmp_hedgeG_se <- sigGeneMatrix_hedgeG_se[-which(is.na(sigGeneMatrix_hedgeG_se[,1])),];
  tmp_sigGeneMatrixUp <- sigGeneMatrixUp[-which(is.na(sigGeneMatrixUp[,1])), ];
  tmp_sigGeneMatrixDown <- sigGeneMatrixDown[-which(is.na(sigGeneMatrixDown[,1])), ];
  
  }else{
    
    tmp_hedgeG <-sigGeneMatrix_hedgeG;
    tmp_hedgeG_se <-  sigGeneMatrix_hedgeG_se;
    tmp_sigGeneMatrixUp <-  sigGeneMatrixUp;
    tmp_sigGeneMatrixDown <- sigGeneMatrixDown;
    
  }
  if(nrow(tmp_hedgeG) != nrow(tmp_hedgeG)){
    
    stop("\nAfter removing NAs, number of rows isn't matching up.")
  }
  summ_hedgeG_ES <- c();
  summ_hedgeG_ES_pvalue <- c();
  summ_hedgeG_ES_se <- c();
    sigGeneMatrixUp_summ <- c();
    sigGeneMatrixDown_summ <- c();
    for(g in 1:nrow(comm_exprList$heatmaps_expr_union[[c]][[e]])){
      
      #ours is random? A “group” effect is random if we can think of the levels we
      #observe in that group to be samples from a larger population.
      #■ Example: if collecting data from different medical centers,
      #“center” might be thought of as random.
      #■ Example: if surveying students on different campuses,
      #“campus” may be a random effect.
      #http://statweb.stanford.edu/~jtaylo/courses/stats203/notes/fixed+random.pdf
      d <- as.vector(tmp_hedgeG[,g]);
      se <- as.vector(tmp_hedgeG_se[,g]);
      #The summary estimate is a weighted average. If weights are specified they are used, otherwise the
      #reciprocal of the estimated variance is used.
      ES <- meta.summaries(d=d, se=se, method=c("random"),
                                  logscale=logscale,
                                  conf.level=0.95);
      
      summ_hedgeG_ES[g] <- ES$summary;
      summ_hedgeG_ES_se[g] <- ES$se.summary;
      #take from Purvesh's code: pool.inverseVar() function.
      #pval     <- 2*pnorm( abs(summ/se.summ), lower.tail=FALSE )
      summ_hedgeG_ES_pvalue[g] <- 2*pnorm( abs(summ_hedgeG_ES[g]/ summ_hedgeG_ES_se[g]), lower.tail=FALSE )
      #summarize wilcoxon rank comm_exprList using fisher's p-values: combine across all studes
      sigGeneMatrixUp_summ[g] <- fishersMethod(tmp_sigGeneMatrixUp[,g]);
      sigGeneMatrixDown_summ[g] <- fishersMethod(tmp_sigGeneMatrixDown[,g]);
                                                                                          
    }
    ##now: do fdr correction across all genes in this group.
  summ_hedgeG_ES_pvalue_fdr <- p.adjust( summ_hedgeG_ES_pvalue,method="fdr");
    sigGeneMatrixUp_summFDR  <- p.adjust(sigGeneMatrixUp_summ,method="fdr");
    sigGeneMatrixDown_summFDR  <- p.adjust(sigGeneMatrixDown_summ,method="fdr");
  
  names( summ_hedgeG_ES_pvalue_fdr ) <- rownames(comm_exprList$heatmaps_expr_union[[c]][[e]]);
  names(sigGeneMatrixUp_summFDR) <- rownames(comm_exprList$heatmaps_expr_union[[c]][[e]]);
  names(sigGeneMatrixDown_summFDR) <- rownames(comm_exprList$heatmaps_expr_union[[c]][[e]]);
  names(summ_hedgeG_ES_se) <- rownames(comm_exprList$heatmaps_expr_union[[c]][[e]]);
  names(summ_hedgeG_ES) <-  rownames(comm_exprList$heatmaps_expr_union[[c]][[e]]);
  
    sigGenesList[[c]] <- list(sigGeneMatrix_hedgeG=sigGeneMatrix_hedgeG,
                              sigGeneMatrix_hedgeG_se=sigGeneMatrix_hedgeG_se,
                              sigGeneMatrixDown=sigGeneMatrixDown,
                              sigGeneMatrixUp=sigGeneMatrixUp,
                              sigGeneMatrixUp_summFDR=sigGeneMatrixUp_summFDR,
                              sigGeneMatrixDown_summFDR=sigGeneMatrixDown_summFDR,summ_hedgeG_ES=summ_hedgeG_ES,
                              summ_hedgeG_ES_pvalue_fdr= summ_hedgeG_ES_pvalue_fdr,summ_hedgeG_ES_se=summ_hedgeG_ES_se
    );
    
   #if there is ANY data in this community. 
  }
    #end of loop c
  }
  

##NOW: summarize your results.
for(s in 1:length(sigGenesList)){
  
  if(s >1){
    
  totalGeneList <- union(totalGeneList,geneCommUnion[[c]]);
  
  }else{
    
    totalGeneList <- geneCommUnion[[c]];
    
  }
  
}

passSigThresh_matrix <- matrix(data=NA,nrow=length(sigGenesList),ncol=length( totalGeneList),
                               dimnames=list(c(1:length(sigGenesList)), totalGeneList));

wilcox_qvalue_matrix <- matrix(data=NA,nrow=length(sigGenesList),ncol=length( totalGeneList),
                               dimnames=list(c(1:length(sigGenesList)), totalGeneList));

ES_matrix <- matrix(data=NA,nrow=length(sigGenesList),ncol=length( totalGeneList),
                        dimnames=list(c(1:length(sigGenesList)), totalGeneList));
ES_se_matrix <- matrix(data=NA,nrow=length(sigGenesList),ncol=length( totalGeneList),
                    dimnames=list(c(1:length(sigGenesList)), totalGeneList));


ES_qvalue_matrix <- matrix(data=NA,nrow=length(sigGenesList),ncol=length( totalGeneList),
                    dimnames=list(c(1:length(sigGenesList)), totalGeneList));

GSEA_up <- list();
GSEA_down <- list();

genes_upList <- list();
genes_downList <- list();

for(c in 1:length(sigGenesList)){
  
  #intersect ES and wilcoxon test.
  up_genes <- intersect(names(sigGenesList[[c]]$sigGeneMatrixUp_summFDR[which(sigGenesList[[c]]$sigGeneMatrixUp_summFDR<=pvalueThresh)]), 
                        names(sigGenesList[[c]]$summ_hedgeG_ES_pvalue_fdr[which(sigGenesList[[c]]$summ_hedgeG_ES_pvalue_fdr<=pvalueThresh)]));
  down_genes <- intersect(names(sigGenesList[[c]]$sigGeneMatrixDown_summFDR[which(sigGenesList[[c]]$sigGeneMatrixDown_summFDR<=pvalueThresh)]),
                          names(sigGenesList[[c]]$summ_hedgeG_ES_pvalue_fdr[which(sigGenesList[[c]]$summ_hedgeG_ES_pvalue_fdr<=pvalueThresh)]));
  
  if(length(intersect(up_genes,down_genes))>0){
    
    matchGenes <- intersect(up_genes,down_genes);
    
    for(m in 1:length(matchGenes)){
      #break the tie.
      
      up_index <- grep(matchGenes[m],down_genes);
      down_index <- grep(matchGenes[m],up_genes);
      
      if(sigGenesList[[c]]$sigGeneMatrixUp_summFDR[up_index]<sigGenesList[[c]]$sigGeneMatrixDown_summFDR[ down_index]){
        
        #remove from down list.
        down_genes <- down_genes[-up_index];
        
      }else{
        #remove from up list.
        up_genes <- up_genes[- up_index];
        
      } 
    }
   #if any intersecting genes. 
  }
  
  if(length(up_genes)>0){
    
    passSigThresh_matrix[c,match(up_genes,colnames( passSigThresh_matrix))] <- "hi";
    #save mann-whitney/wilcoxon p-values
    wilcox_qvalue_matrix[c,match(up_genes,colnames(wilcox_qvalue_matrix))] <- up_genes;
    #run GSEA on up genes.
    #save up_genes list too.
    ES_qvalue_matrix[c,match(up_genes,colnames(ES_qvalue_matrix))] <-  
      sigGenesList[[c]]$summ_hedgeG_ES_pvalue_fdr[match(up_genes,names(sigGenesList[[c]]$summ_hedgeG_ES_pvalue_fdr))];
    ES_matrix[c,match(up_genes,colnames(ES_matrix))] <-
      sigGenesList[[c]]$summ_hedgeG_ES[match(up_genes,names(sigGenesList[[c]]$summ_hedgeG_ES))];
    ES_se_matrix[c,match(up_genes,colnames(ES_se_matrix))] <-
      sigGenesList[[c]]$summ_hedgeG_ES_se[match(up_genes,names(sigGenesList[[c]]$summ_hedgeG_ES_se))];
    
    genes_upList[[c]] <- up_genes;
    GSEA_up[[c]] <-  GSEA(testGeneVector=genes_upList[[c]],refGeneLists=refGeneLists,method=c("fisher"),
                          genomeSize=genomeSize)$qvalues;
    
    names(GSEA_up[[c]]) <- names(refGeneLists);
    GSEA_up[[c]] <- GSEA_up[[c]][which(GSEA_up[[c]]<=pvalueThresh)];

  
  }
  
  if(length(down_genes)>0){
    
  passSigThresh_matrix[c,match(down_genes,colnames( passSigThresh_matrix))] <- "lo";
  wilcox_qvalue_matrix[c,match(down_genes,colnames(wilcox_qvalue_matrix))] <- down_genes;
  ES_qvalue_matrix[c,match(down_genes,colnames(ES_qvalue_matrix))] <-  
    sigGenesList[[c]]$summ_hedgeG_ES_pvalue_fdr[match(down_genes,names(sigGenesList[[c]]$summ_hedgeG_ES_pvalue_fdr))]
  ES_matrix[c,match(down_genes,colnames(ES_matrix))] <-
    sigGenesList[[c]]$summ_hedgeG_ES[match(down_genes,names(sigGenesList[[c]]$summ_hedgeG_ES))];
  ES_se_matrix[c,match(down_genes,colnames(ES_se_matrix))] <-
    sigGenesList[[c]]$summ_hedgeG_ES_se[match(down_genes,names(sigGenesList[[c]]$summ_hedgeG_ES_se))];
  
  genes_downList[[c]] <- down_genes;
  GSEA_down[[c]] <-  GSEA(testGeneVector=genes_downList[[c]],refGeneLists=refGeneLists,method=c("fisher"),
                        genomeSize=genomeSize)$qvalues;
  
  names(GSEA_down[[c]]) <- names(refGeneLists);
  GSEA_down[[c]] <- GSEA_down[[c]][which(GSEA_down[[c]]<=pvalueThresh)];
 
  }
  
  #end of loop c (each community)
}

#make a binary matrix that's nice for heatmap plots.
  binarySigMatrix <- matrix(data=NA,nrow=length(sigGenesList),ncol=length( totalGeneList),
                            dimnames=list(c(1:length(sigGenesList)), totalGeneList));
binarySigMatrix[which(passSigThresh_matrix=="lo")] <- -1;
binarySigMatrix[which(passSigThresh_matrix=="hi")] <- 1;
binarySigMatrix[which(is.na(passSigThresh_matrix))] <- 0;


cat("\n indStudySigGeneData object contains the more granular statistics for all genes, from individual studies and summarized.
The ES,wilcox matrix objects are only for the significantly over or under-expressed genes found in each community.\n")
  sigGeneSummary <- list(passSigThresh_matrix=passSigThresh_matrix,
                         wilcox_qvalue_matrix=wilcox_qvalue_matrix,
                         ES_qvalue_matrix=ES_qvalue_matrix,binarySigMatrix=binarySigMatrix,ES_matrix=ES_matrix,
                         ES_se_matrix=ES_se_matrix, genes_downList=genes_downList,GSEA_down=GSEA_down,
                         GSEA_up=GSEA_up,genes_upList=genes_upList,indStudySigGeneData=sigGenesList,
                         comm_exprList=comm_exprList,geneCommIntersect=geneCommIntersect,geneCommUnion=geneCommUnion);



return(sigGeneSummary);

}
##################
#genome size: size of gene pool you started from
#can be a list?
#could be original size of your microarray...really just need dimensions
communityGeneSummary <- function(community_membership,biclust,detailedBiclustNames,genomeSize,
                                 refGeneLists,sourceDir="/home/kplaney/gitRepos/RNAseq_pipeline/",
                                 default_refGeneListDir="/home/data/MSigDB/GSEA_base_MSigDB_lists_merged.RData.gzip"){
  
  startTime <- proc.time();
  setwd(sourceDir);
  source("AO_GSEA.R");
  
  if(missing(refGeneLists)){
    
    warning("\nUsing default MSigDB lists: MSigDB_onco_symbols, MSigDB_CanPath_symbols,MSigDB_TFT_symbols,MSigDB_immun_symbols,
            and MSigDB_cancerNeigh_symbols.\n");
    load(default_refGeneListDir);
    refGeneLists <- GSEA_base_MSigDB_lists_merged;
    
  }
  #ALSO: track # genes, patients, orig dataset, for each node.
  biclustMasterList <- list();
  for(b in 1:length(biclust)){
    
    biclustMasterList <- append(biclustMasterList,biclust[[b]]);
    
  }
  names(biclustMasterList) <- c(1:length(biclustMasterList));
  
  #now create lists by communities. re-name to master biclust # (or underscore?)
  biclustCommunities <- list();
  communityNames <- unique(community_membership[,"community"]);
  for(c in 1:length(communityNames)){
    
    #community indices are same as that in the master list (same as those in IGP matrices)
    #CAREFUL: the levels make this weird...so just change to characters to match the names of the biclustMasterList.
    biclustCommunities[[c]] <- biclustMasterList[as.character(community_membership[which(community_membership[,"community"]==communityNames[c]), "clust"]) ];
    
  }
  
  
  qvalues_community <- list();
  geneCommIntersect <- list();
  geneCommUnion <- list();
  qvalues <- list();
  qvaluesDF <- list();
  for(c in 1:length(biclustCommunities)){
    
    qvalues[[c]] <- list();
    
    qvaluesDF[[c]] <- data.frame();
    for(g in 1:length(biclustCommunities[[c]])){
      #run hypergeometric tests for each gene
      #unlist biclusters?? hmm...or use detailed names??
      #assume genes are in columns
      testGeneVector  <- colnames(biclustCommunities[[c]][[g]]);
      
      #fisher is better for smaller genes lists (and can still handle large ones.)
      qvalues[[c]][[g]] <- GSEA(testGeneVector=testGeneVector,refGeneLists=refGeneLists,method=c("fisher"),genomeSize=genomeSize)$qvalues;  
      
      if(g >1){
        
        qvaluesDF[[c]] <- cbind(qvaluesDF[[c]], qvalues[[c]][[g]]);
        colnames(qvaluesDF[[c]])[g] <- names(biclustCommunities[[c]])[g];
        
      }else{
        
        qvaluesDF[[c]] <- data.frame(qvalues[[c]][[g]]);
        colnames(qvaluesDF[[c]])[g] <- names(biclustCommunities[[c]])[g];
      }
      
      if(g >1){
        
        geneCommIntersect[[c]] <- intersect(geneCommIntersect[[c]],testGeneVector);
        geneCommUnion[[c]] <- union(geneCommIntersect[[c]],testGeneVector);
        
        #initialize the lists g==1
      }else{
        
        geneCommIntersect[[c]] <- testGeneVector;
        geneCommUnion[[c]] <- testGeneVector;
        
      }
      
    }
    #run over-enrichmet again but just on the intersecting genes across all datasets
    if(length(geneCommIntersect[[c]])>0){
      
      #COME BACK: do I need to also correct for the # of clusters this came from?
      #ie what's the chance this gene list would pop up across N clusters to make this community
      qvalues_community[[c]] <- GSEA(testGeneVector=geneCommIntersect[[c]],refGeneLists=refGeneLists,method=c("fisher"),
                                     genomeSize=genomeSize)$qvalues;
      
    }else{
      
      warning("\nNo genes 100% overlapped in community ",c,"\n");
      qvalues_community[[c]] <- NA;
      
    }
    
    rownames(qvaluesDF[[c]]) <-  paste0(names(refGeneLists));
    
    if(c >1){
      
      qvalues_communityDF <- cbind(qvalues_communityDF,data.frame(unlist(qvalues_community[[c]])));
      
      
    }else{
      
      qvalues_communityDF  <- data.frame(unlist(qvalues_community[[c]]));
    }
    #end of loop c
  }
  
  rownames(qvalues_communityDF) <- paste0(names(refGeneLists));
  colnames(qvalues_communityDF) <- paste0("community_",c(1:length(biclustCommunities)));
  runTime <- startTime-proc.time();
  
  output <- list(qvalues=qvalues,qvalues_community=qvalues_community,geneCommIntersect=geneCommIntersect,geneCommUnion=geneCommUnion,biclustCommunities=biclustCommunities,runTime=runTime,
                 qvalues_communityDF=qvalues_communityDF,qvaluesDF=qvaluesDF,refGeneListNames=names(refGeneLists));
}


#######
#communityMembership=communityMembership$attrDF
communityHeatmaps <- function(community_membership,biclust,dataMatrices,dataM_sampleCol=FALSE,geneCommIntersect,
                              geneCommUnion,
                              fullHeatmaps=FALSE){
  
  #assumes studyNum is related to order of data matrices.
  names(dataMatrices) <- c(1:length(dataMatrices));

  warning("\nThis code assumes your *biclusters* have genes in the columns and patients in the rows.\n")
  library("RColorBrewer");
  biclustMasterList <- list();
  biclustMasterList_origStudies <- list();
  
  for(b in 1:length(biclust)){
    
    biclustMasterList <- append(biclustMasterList, biclust[[b]]);
    biclustMasterList_origStudies <- append(biclustMasterList_origStudies,rep(b,length(biclust[[b]])));
  }
  names(biclustMasterList) <- c(1:length(biclustMasterList));
  
  #now create lists by communities. re-name to master biclust # (or underscore?)
  biclustCommunities <- list();
  communityNames <- unique(community_membership[,"community"]);
  communityStudyNums <- list();
  
  heatmaps_expr <- list();
  heatmaps_corr <- list();
  heatmaps_expr_intersect <- list();
  heatmaps_corr_intersect <- list();
  heatmaps_expr_union <- list();
  heatmaps_corr_union <- list();
  GSMID_names <- list();
  expr_homogeneousClust <- list();

  for(c in 1:length(communityNames)){
    
    expr_homogeneousClust[[c]] <- list();
    #community indices are same as that in the master list (same as those in IGP matrices)
    #CAREFUL: the levels make this weird...so just change to characters to match the names of the biclustMasterList.
    biclustCommunities[[c]] <- biclustMasterList[as.character(community_membership[which(community_membership[,"community"]==communityNames[c]), "clust"]) ];
    
    communityStudyNums[[c]] <- community_membership[which(community_membership[,"community"]==communityNames[c]), "studyNum"];
    heatmaps_expr[[c]] <- list();
    heatmaps_corr[[c]] <- list();
    heatmaps_expr_intersect[[c]] <-list();
    heatmaps_corr_intersect[[c]] <- list();
    heatmaps_expr_union[[c]] <-list();
    heatmaps_corr_union[[c]] <- list();
    GSMID_names[[c]] <- list();
    
    for(e in 1:length(biclustCommunities[[c]])){
      
    
      #also only want genes in that biclust.
      if(dataM_sampleCol){
        
        if(ncol(dataMatrices[[as.character(communityStudyNums[[c]][e])]])>nrow(biclustCommunities[[c]][[e]])){
          
        if(length(setdiff(colnames(dataMatrices[[as.character(communityStudyNums[[c]][e])]]),
                          rownames(biclustCommunities[[c]][[e]])))==0){
          
          stop("\nNot intersecing bicluster and other patients correctly.")
        }
        
        #BUT some should intersect!
        if(length(intersect(colnames(dataMatrices[[as.character(communityStudyNums[[c]][e])]]),
                          rownames(biclustCommunities[[c]][[e]])))==0){
          
          stop("\nNot intersecing bicluster and other patients correctly.")
        }
        
        #DEBUG: is this OK?
        #take full gene list.
        comm_other <- dataMatrices[[as.character(communityStudyNums[[c]][e])]][colnames(biclustCommunities[[c]][[e]]),setdiff(colnames(dataMatrices[[as.character(communityStudyNums[[c]][e])]]),
                                                                                                               rownames(biclustCommunities[[c]][[e]]))];
        
        #keep full gene list in case need it later for union genes.
        comm_otherFull <- dataMatrices[[as.character(communityStudyNums[[c]][e])]][, setdiff(colnames(dataMatrices[[as.character(communityStudyNums[[c]][e])]]),
                                                                                                                     rownames(biclustCommunities[[c]][[e]]))];
        clustExprFull <- dataMatrices[[as.character(communityStudyNums[[c]][e])]][ ,rownames(biclustCommunities[[c]][[e]])];
        
        }else{
          #no "other" patients to test on!
          comm_other <- NULL;
          comm_otherFull <- NULL;
          clustExprFull <- NULL;
        }
        
      }else{
        
        #careful: need as.character to index these correctly! are we getting the EXACT same patients here? this should only be the case
        #if the bicluster is the ENTIRE set of patients in the matrix.
        if(nrow(dataMatrices[[as.character(communityStudyNums[[c]][e])]])>nrow(biclustCommunities[[c]][[e]])){
          
        if(length(setdiff(rownames(dataMatrices[[as.character(communityStudyNums[[c]][e])]]),rownames(biclustCommunities[[c]][[e]])))==0){
          
          stop("\nNot intersecting bicluster and other patients correctly.");
        }
        
        if(length(intersect(rownames(dataMatrices[[as.character(communityStudyNums[[c]][e])]]),rownames(biclustCommunities[[c]][[e]])))==0){
          
          stop("\nNot intersecting bicluster and other patients correctly.");
        }
        
        comm_other <- dataMatrices[[as.character(communityStudyNums[[c]][e])]][setdiff(rownames(dataMatrices[[as.character(communityStudyNums[[c]][e])]]),rownames(biclustCommunities[[c]][[e]])),
                                                                 colnames(biclustCommunities[[c]][[e]])];
        
        #keep full gene list in case need it later for union genes.
        comm_otherFull <- dataMatrices[[as.character(communityStudyNums[[c]][e])]][setdiff(rownames(dataMatrices[[as.character(communityStudyNums[[c]][e])]]),rownames(biclustCommunities[[c]][[e]])), ];
        
        clustExprFull <- dataMatrices[[as.character(communityStudyNums[[c]][e])]][rownames(biclustCommunities[[c]][[e]]), ];
        
        }else{
          #no "other" patients to test on!
          comm_other <- NULL;
          comm_otherFull <- NULL;
        }
      
      }

      
      if(!is.null(comm_other)){
        
      #keep track of true patient names.
      GSMID_names[[c]][[e]] <- append(rownames(clustExprFull),rownames(comm_otherFull));
      #make binary in cluster or not sample names.
      rownames(comm_other) <- rep(0,nrow(comm_other));
      rownames(comm_otherFull) <- rep(0,nrow(comm_otherFull));
      clustExpr <- biclustCommunities[[c]][[e]];
      rownames(clustExpr) <-  rep(1,nrow(clustExpr));
      rownames(clustExprFull) <-  rep(1,nrow(clustExprFull));
      names(GSMID_names[[c]][[e]]) <- append(rownames(clustExprFull),rownames(comm_otherFull));
      #want samples in column for visualization
      #rev(brewer.pal(11,"RdBu")) for red=overexpression. blue is under.
      
      
      if(fullHeatmaps){
        
      if(dataM_sampleCol){
        
       # heatmap(data.matrix(cbind(clustExpr,comm_other)), 
        #        col = rev(brewer.pal(11,"RdBu")),Rowv=NA,Colv="Rowv",scale="none",
        #        main=paste0("Expression heatmap for a cluster in community ",c ," \nvs rest of samples from study ",communityStudyNums[[c]][e]));
        heatmaps_expr[[c]][[e]] <- data.matrix(cbind(clustExpr,comm_other));
        
        names(heatmaps_expr[[c]])[e] <- paste0("comm_",c,"study_",communityStudyNums[[c]][e]);
                  
        #heatmap(cor(t(data.matrix(cbind(clustExpr,comm_other)))), 
        #        col = rev(brewer.pal(11,"RdBu")),Rowv=NA,Colv="Rowv",scale="none",
        #        main=paste0("Correlation heatmap for a cluster in community ",c ," \nvs rest of samples from study ",communityStudyNums[[c]][e]));
        heatmaps_corr[[c]][[e]] <- cor(t(data.matrix(cbind(clustExpr,comm_other)))); 
                                           
        names(heatmaps_corr[[c]])[e] <- paste0("comm_",c,"study_",communityStudyNums[[c]][e]);
        
      }else{
        
        heatmaps_expr[[c]][[e]] <- t(data.matrix(rbind(clustExpr,comm_other))); 
        names(heatmaps_expr[[c]])[e] <- paste0("comm_",c,"study_",communityStudyNums[[c]][e]);                                 
        heatmaps_corr[[c]][[e]] <- cor(t(data.matrix(rbind(clustExpr,comm_other))));
        names(heatmaps_corr[[c]])[e] <- paste0("comm_",c,"study_",communityStudyNums[[c]][e]);                                  
        
      }
 
      #end of if want full heatmaps.
    }
      if(!missing(geneCommIntersect)){
        #now only use the intersecting genes
        
        if(length(geneCommIntersect[[c]])>0){
          
        if(dataM_sampleCol){
          
          heatmaps_expr_intersect[[c]][[e]] <- data.matrix(cbind(clustExprFull,comm_otherFull)[geneCommIntersect[[c]], ]);
          names(heatmaps_expr_intersect[[c]])[e] <- paste0("comm_",c,"study_",communityStudyNums[[c]][e]);                                           
          heatmaps_corr_intersect[[c]][[e]] <-  data.matrix(cbind(clustExprFull,comm_otherFull)[geneCommIntersect[[c]], ]);
          names(heatmaps_corr_intersect[[c]][[e]] ) <- paste0("comm_",c,"study_",communityStudyNums[[c]][e]);                                             
        
          #run Sam R here? or in a separate function?
          }else{
          
          #cool...this finds better submatrices!!
            #transposing this will make more sense when plot later.
          heatmaps_expr_intersect[[c]][[e]] <- data.matrix(t(rbind(clustExprFull,comm_otherFull)[ ,geneCommIntersect[[c]]]));
          names(heatmaps_expr_intersect[[c]])[e] <- paste0("comm_",c,"study_",communityStudyNums[[c]][e]);                                          
          heatmaps_corr_intersect[[c]][[e]] <- cor(t(data.matrix(rbind(clustExprFull,comm_otherFull)[ ,geneCommIntersect[[c]]]))); 
          names(heatmaps_corr_intersect[[c]])[e] <- paste0("comm_",c,"study_",communityStudyNums[[c]][e]);
        
        }
        
      }
      #if there are no other patients in the matrix besides the one in this bicluster.
       names(heatmaps_expr_intersect[[c]])[e] <-  names(dataMatrices)[as.character(communityStudyNums[[c]][e])];
      names(heatmaps_corr_intersect[[c]])[e] <-  names(dataMatrices)[as.character(communityStudyNums[[c]][e])];
      
      #if not missing intersect genes
      }
    
      #NOW: do for UNION also.
      if(!missing(geneCommUnion)){
        #now only use the intersecting genes
        #COME BACK: IF CERTAIN GENES JUST DON'T EXIST IN FULL DATASET:  USE MATCH THEN INSTEAD OF JUST SUBSETTING BY ENETIRE GENE LIST.
        if(length(geneCommUnion[[c]])>0){
          
          if(dataM_sampleCol){
            
            heatmaps_expr_union[[c]][[e]] <- t(data.matrix(cbind(clustExprFull,comm_otherFull)[geneCommUnion[[c]], ]));
            names(heatmaps_expr_union[[c]])[e] <- paste0("comm_",c,"study_",communityStudyNums[[c]][e]);                                           
            heatmaps_corr_union[[c]][[e]] <-  data.matrix(cbind(clustExprFull,comm_otherFull)[geneCommUnion[[c]], ]);
            names(heatmaps_corr_union[[c]][[e]] ) <- paste0("comm_",c,"study_",communityStudyNums[[c]][e]);                                             
            
            
          }else{
          
            #cool...this finds better submatrices!!
            heatmaps_expr_union[[c]][[e]] <- t(data.matrix(rbind(clustExprFull,comm_otherFull)[ ,geneCommUnion[[c]]]));
            names(heatmaps_expr_union[[c]])[e] <- paste0("comm_",c,"study_",communityStudyNums[[c]][e]);                                          
            heatmaps_corr_union[[c]][[e]] <- cor(t(data.matrix(rbind(clustExprFull,comm_otherFull)[ ,geneCommUnion[[c]]]))); 
            names(heatmaps_corr_union[[c]])[e] <- paste0("comm_",c,"study_",communityStudyNums[[c]][e]);
            
          }
          
        }
     #if not missing union genes.
      }
 #if cluster is not entire dataset
    }else{
      
      #keep track of true patient names.
      GSMID_names[[c]][[e]] <- rownames(biclustCommunities[[c]][[e]])

      names(GSMID_names[[c]][[e]]) <- rep.int(1,times=nrow(biclustCommunities[[c]][[e]]));
      #want samples in column for visualization
      #rev(brewer.pal(11,"RdBu")) for red=overexpression. blue is under.
      #just use all genes here for now.
      #COME BACK: do intersecting genes too.
      expr_homogeneousClust[[c]][[e]] <- data.matrix(dataMatrices[[as.character(communityStudyNums[[c]][e])]][rownames(biclustCommunities[[c]][[e]]), ]);
      
      names(expr_homogeneousClust[[c]])[e] <- paste0("comm_",c,"study_",communityStudyNums[[c]][e]);
          
    
  
    }
    
 #end of loop e. 
  }
 #end of loop c
  }
  
  
  output <- list(studyNames=names(dataMatrices),GSMID_names=GSMID_names,heatmaps_expr_union=heatmaps_expr_union,
                 heatmaps_expr_intersect=heatmaps_expr_intersect,heatmaps_expr=heatmaps_expr,
                 heatmaps_corr_intersect=heatmaps_corr_intersect,heatmaps_corr=heatmaps_corr,
                 biclustCommunities=biclustCommunities,heatmaps_corr_union=heatmaps_corr_union,
                 expr_homogeneousClust=expr_homogeneousClust);
  
  return(output);
}

#######


assignCentroidSubtype <- function(origDataMatrix,minNumGenes=30,centroidRData="/home/data/breast_microarrayDB/pam50_centroids.RData"){
  
  warning("\nAssumes genes are in the columns of your data matrix, \nand for centroids, gene names are in the first columns.\n")
  load(centroidRData);
  #originally meant for use with pam50 centroids, but can be used on any matrix where first col is gene names, other columns are centroid groups.
  #RData object must be names centroidMatrix
  if(missing(centroidMatrix)){
    
    stop("\nRData object must have name 'centroidMatrix'.")
 
  }
  
  centroids <- centroidMatrix[, 2:ncol(centroidMatrix)];
  #only take the genes that are in this data matrix
  centroidMatrix_genes <- as.character(centroidMatrix[,1][centroidMatrix[,1] %in% colnames(origDataMatrix)]);
  
  if(length(centroidMatrix_genes)>minNumGenes){
    
    subtypes <- array(data=NA,dim=nrow(origDataMatrix));
    centroids <- centroids[centroidMatrix[,1] %in% colnames(origDataMatrix), ];
    dataMatrix <- origDataMatrix[ ,centroidMatrix_genes];
    corMatrix <- cor(t(dataMatrix),centroids);
    
    for(s in 1:length(subtypes)){
      #sometimes, if no variance for a list of genes: will return NA.
      subtypes[s] <- which(corMatrix[s,]==max(corMatrix[s,],na.rm=TRUE));
      
    }
    
    subtypeNames <- colnames(centroidMatrix)[2:ncol(centroidMatrix)];
    subtypeLabels <- array(data=NA,dim=nrow(dataMatrix));
    
    for(c in 1:length(unique(subtypeNames))){
      
      subtypeLabels[which(subtypes==c)] <- subtypeNames[c];
      
      
    }
    
    subtypes <- cbind(subtypes,subtypeLabels);
    output <- list(subtypes=subtypes,centroidMatrix_genes=centroidMatrix_genes);
    return(output);
    
  }else{
    
    return(NA);
    
  }

  return(output);
  
}
####
communitySubtypeBreakdown_plots <- function(community_membership,biclust,origDataMatrices,saveDir="./",clustMethodName="clust method",
                                             centroidSetName="pam50",centroidRData="/home/data/breast_microarrayDB/pam50_centroids.RData"){
  
  
  #specify colors:
  subtypeColorMatrix <- data.frame(c("#970eec","#ec0e58","#0e58ec","#7d7d7d","#0eec97"))
  #hmm can I order the colors in the same way each time?? may not really need to anyways.
  #must tranpose for ggplot to read it correctly and THEN add names for this to work.
  subtypeColorMatrix <- t(subtypeColorMatrix);
  names(subtypeColorMatrix) <- c("LumB","LumA","Her2","Basal","Normal");
  
  numTotalSamplesInStudy <- list();
  for(s in 1:length(origDataMatrices)){
    
    numTotalSamplesInStudy[[s]] <- nrow(origDataMatrices[[s]]);
    
  }
  library("limma");
  library("ggplot2");
  warning("\nThis code assumes your biclusters have genes in the columns and patients in the rows.\n")
  #library("RColorBrewer");
  biclustMasterList <- list();
  biclustMasterList_origStudies <- list();
  
  #IS THIS WORKING?? cluster 83 should have way more patients!!
  for(b in 1:length(biclust)){
    
    biclustMasterList <- append(biclustMasterList, biclust[[b]]);
    biclustMasterList_origStudies <- append(biclustMasterList_origStudies,rep(b,length(biclust[[b]])));
  }
  names(biclustMasterList) <- c(1:length(biclustMasterList));
  
  #now create lists by communities. re-name to master biclust # (or underscore?)
  biclustCommunities <- list();
  communityNames <- unique(community_membership[,"community"]);
  communityStudyNums <- list();
  
  
  subtype_plot <- list();
  subtype_dfList <- list();
  subtype_plot_stacked <- list();
  subtype_plot_fract <- list();
  
  #remove the levels.
  communityNames <- as.character(communityNames);
  for(c in 1:length(communityNames)){
    
    #community indices are same as that in the master list (same as those in IGP matrices)
    #CAREFUL: the levels make this weird...so just change to characters to match the names of the biclustMasterList.
    clust_nums <- as.numeric(as.character(community_membership[which(community_membership[,"community"]==communityNames[c]), "clust"]));
    biclustCommunities[[c]] <- biclustMasterList[clust_nums];
    
    #want to remove levels! otherwise indexing gets all messed up.
    communityStudyNums[[c]] <- as.numeric(as.character(community_membership[which(community_membership[,"community"]==communityNames[c]), "studyNum"]));
    
    subtype_dfList[[c]] <- data.frame();
    

    for(e in 1:length(biclustCommunities[[c]])){
      
      patIDs <- rownames(biclustCommunities[[c]][[e]]);
      
      studyNum <- rep.int(communityStudyNums[[c]][[e]],times=length(patIDs));
      subtype <- assignCentroidSubtype(origDataMatrix=origDataMatrices[[communityStudyNums[[c]][e]]][patIDs,],minNumGenes=30,centroidRData=centroidRData)$subtypes[,"subtypeLabels"];
      clust_numVar <- rep(clust_nums[e],length(subtype));
      temp <- data.frame(patIDs,subtype,studyNum,clust_numVar);
      subtype_dfList[[c]] <- rbind(subtype_dfList[[c]],temp); 

    }
    
    if(c >1){
      
      tmp <- data.frame(subtype_dfList[[c]],rep(c,times=nrow(subtype_dfList[[c]])));
      colnames(tmp)[ncol(tmp)] <- "community";
      
      subtype_dfMaster <- rbind(subtype_dfMaster,tmp);
      
    }else{
      
      tmp <- data.frame(subtype_dfList[[c]],rep(c,times=nrow(subtype_dfList[[c]])));
      colnames(tmp)[ncol(tmp)] <- "community";
      subtype_dfMaster <- tmp;
      #trying to create percentages.... a column for subtype, for the percentage, the community num,
      #and the studNum
      #subtype_dfSum <- data.frame(unique(subtype),table(subtype)/length(subtype),)
      
    }
 
    
    subtype_plot[[c]] <-  ggplot(data=subtype_dfList[[c]],aes(x=subtype))+geom_histogram(aes(fill=subtype))+
      labs(y="Number of samples" ,title=paste0(centroidSetName," subtypes for community ",c))+
      scale_fill_manual(values = subtypeColorMatrix)+ labs(fill="subtype")+
      theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
      theme(axis.text.x = element_text(colour = "black",size=18,angle=45,vjust=.9,hjust=1),axis.title.x=element_blank(),
            axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1))+
      theme(plot.title=element_text(colour="black",size=28,vjust=1,hjust=.2));

    #stacked bar for all studies.
    subtype_plot_stacked[[c]] <- ggplot(data=subtype_dfList[[c]],aes(x=factor(studyNum)))+geom_bar(aes(fill=factor(subtype)))+
      labs(y="Number of samples", x="Dataset number",title=paste0(centroidSetName," subtype by dataset for community ",c))+
      scale_fill_manual(values = subtypeColorMatrix)+labs(fill="subtype")+
      theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
      theme(axis.text.x = element_text(colour = "black",size=12,angle=45,vjust=1,hjust=1),
            axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1),
            axis.title.x = element_text(colour = "black",size=18,vjust=0))+
      theme(plot.title=element_text(colour="black",size=28,vjust=1,hjust=.2));
    
    #calculate percentage
    df <- split(subtype_dfList[[c]],f=subtype_dfList[[c]]$studyNum);
    subtype_fract <- lapply(df,FUN=function(e){
      
      output <- table(e$subtype)/nrow(e);
      
    });
    
    
    for(s in 1:length(subtype_fract)){
      
      for(t in 1:length(subtype_fract[[s]])){
        
        if(t >1){
          
          tmp <- rbind(tmp,data.frame(rep(names(subtype_fract)[s],round(subtype_fract[[s]][t]*100)),rep(names(subtype_fract[[s]])[t],round(subtype_fract[[s]][t]*100)),
                                      rep(c,round(subtype_fract[[s]][t]*100))));
        }else{
          
          tmp <- data.frame(rep(names(subtype_fract)[s],round(subtype_fract[[s]][t]*100)),rep(names(subtype_fract[[s]])[t],round(subtype_fract[[s]][t]*100)),
                            rep(c,round(subtype_fract[[s]][t]*100)));
        }
        
      }
      
      if(s >1){
        
        subtype_fractDF  <- rbind(subtype_fractDF,tmp);
        
      }else{
        
        subtype_fractDF  <- tmp;
        
      }
      
    }
    colnames( subtype_fractDF) <- c("studyNum","subtype","community");
    #now get percentages down to 100, not rounded to 101.
    tmp <- table(subtype_fractDF$studyNum);
    
    for(t in 1:length(tmp)){
      
      if(tmp[t]==101){
        #need to remove patients (should be just 1) to make 100.
        
        subtype_fractDF <-  subtype_fractDF[-which(subtype_fractDF$studyNum==names(tmp[t]))[1],];
      }
      
    }
    colnames( subtype_fractDF) <- c("studyNum","subtype","community");

    subtype_plot_fract[[c]]<-  ggplot(data=subtype_fractDF,aes(x=studyNum))+geom_bar(aes(fill=factor(subtype)))+
      labs(x="Dataset number",y="Composition of subtypes",title=paste0("Composition of ",centroidSetName," subtypes \nby dataset for community ",c))+
      labs(fill="subtype")+
    scale_fill_manual(values = subtypeColorMatrix)+  theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
      theme(axis.text.x = element_text(colour = "black",size=12,angle=45,vjust=1,hjust=1),
            axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1),
            axis.title.x = element_text(colour = "black",size=18,vjust=0))+
      theme(plot.title=element_text(colour="black",size=28,vjust=1,hjust=.5));
    
    if(c>1){
      
      
      subtype_fractMaster <- rbind(subtype_fractMaster,subtype_fractDF);
      
    }else{
      
      subtype_fractMaster <- subtype_fractDF;
      
    }
    
  }
  
  #do a facetted plot by all communities too
  subtype_plotALL <- ggplot(data=subtype_dfMaster,aes(x=subtype))+geom_histogram(aes(fill=subtype))+facet_grid(.~community,scales="free_x")+
    labs(y="Number of samples", x="Pam50 subtype breakdown",title=paste0(centroidSetName," subtype for each community"))+
    scale_fill_manual(values = subtypeColorMatrix)+
    labs(fill="subtype")+theme(strip.text.x = element_text(size = 17, colour = "black",face="bold"))+
    theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
    theme(axis.text.x = element_text(colour = "black",size=12,angle=45,vjust=1,hjust=1),axis.title.x= element_text(colour = "black",size=20,vjust=1),
          axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1))+
    theme(plot.title=element_text(colour="black",size=28,vjust=1,hjust=.2));
  
  png(filename=paste0(saveDir,"/",clustMethodName,"_subtype_plotALL_",Sys.Date(),".png"),
      width = 700, height = 1000);
  
  plot(subtype_plotALL);
  dev.off();
  
  library("grid");
  #stacked bar for all studies.
  #scales="free_x":  doesn't plot studies that are empty for that community.
  subtype_plot_stackedALL <- ggplot(data=subtype_dfMaster,aes(x=factor(studyNum)))+geom_bar(aes(fill=factor(subtype)))+
    facet_grid(facets=.~community,scales="free_x")+scale_fill_manual(values = subtypeColorMatrix)+
    labs(y="Number of samples",x="Dataset number",title=paste0(centroidSetName," subtype by community by dataset"))+
  labs(fill="subtype")+theme(strip.text.x = element_text(size = 17, colour = "black",face="bold"))+
    theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
    theme(axis.text.x = element_text(colour = "black",size=6,angle=45,vjust=1,hjust=1),axis.title.x= element_text(colour = "black",size=20,vjust=1),
          axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1))+
    theme(plot.title=element_text(colour="black",size=28,vjust=1,hjust=.2));
   
  
  
  png(filename=paste0(saveDir,"/",clustMethodName,"_subtype_plot_stackedALL_",Sys.Date(),".png"),
      width = 700, height = 1000);
  
  plot(subtype_plot_stackedALL);
  dev.off();

  #calculate percentage
  subtype_plot_fractALL <- ggplot(data=subtype_fractMaster,aes(x=studyNum))+geom_bar(aes(fill=factor(subtype)))+facet_grid(.~community,scales="free_x")+
    labs(xlim=100,x="Dataset number",y="Composition of subtypes",title=paste0("Composition of ",centroidSetName," subtype by community by dataset"))+
    scale_fill_manual(values = subtypeColorMatrix)+
    labs(fill="subtype")+theme(strip.text.x = element_text(size = 17, colour = "black",face="bold"))+
    theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
    theme(axis.text.x = element_text(colour = "black",size=6,angle=45,vjust=1,hjust=1),axis.title.x= element_text(colour = "black",size=20,vjust=1),
          axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1))+
    theme(plot.title=element_text(colour="black",size=28,vjust=1,hjust=.2));
    
  
  png(filename=paste0(saveDir,"/",clustMethodName,"_subtype_plot_fractALL_",Sys.Date(),".png"),
      width = 700, height = 1000);
  plot(subtype_plot_fractALL);
  dev.off();
  
  warning("Just becuase samples are from the same study, doesn't mean they weren't clustered into 2+ groups. 
          But these groups were deemed similar using my meta-cluster method.\n");
  
  output <- list(subtype_dfMaster=subtype_dfMaster,subtype_dfList=subtype_dfList,subtype_plot=subtype_plot,
                 subtype_plot_stacked=subtype_plot_stacked,subtype_plot_fract=subtype_plot_fract,
                 subtype_fractMaster =subtype_fractMaster,subtype_plot_fractALL=subtype_plot_fractALL,
                 subtype_plot_stackedALL=subtype_plot_stackedALL,subtype_plotALL=subtype_plotALL);
  return(output);
  
}

plot_communitiesWithSubtypes  <- function(subtype_dfMaster,community_membership,igraph_edgeDF,saveDir="./"){
    
    #add subtypes to communities
    community_membership <- cbind(community_membership,array(data=NA,dim=nrow(community_membership)));
    colnames(community_membership)[ncol(community_membership)] <- "subtype";
    #remove the levels. as.numeric messes it up further unless use as.character first...
    clustNum <- as.numeric(as.character(unique(community_membership[,"clust"])));
    
    if(length(clustNum)>0){
    for(s in 1:length(clustNum)){
      
      tmp <- subtype_dfMaster[which(subtype_dfMaster[,"clust_numVar"]==clustNum[s]), ];
      
      if(nrow(tmp)==0){
        
        stop("\nNot finding any samples that link to this cluster number.");
      }
      
     #assume each cluster can only be in 1 community.    
      subtypes <- table(tmp[,"subtype"]);
      #just take the max # as the subtype for this cluster.
      if(!all(is.null(names(subtypes)))){
        #if tie: just take first.
        community_membership[which(community_membership[,"clust"]==clustNum[s]),"subtype"] <- names(subtypes)[which(subtypes==max(subtypes))[1]]
          
      }else{
       
        community_membership[which(community_membership[,"clust"]==clustNum[s]),"subtype"] <- which(subtypes==max(subtypes))[1];
       
     }
     
  }
    }
  
  undirGraph <- graph.data.frame(igraph_edgeDF,directed=FALSE,vertices=community_membership);
  #color the nodes. V() gets the vertices. can do V()$ what added  with igraph_attrDF, like $clust.
  #vertex.size=4,vertex.label.dist=0.5,edge.weight=E(undirGraph)$sampleFract,
  #size is size of vertex
  finalNumCommunities <- length(unique(community_membership[,"community"]));
  png(filename=paste0(saveDir,clustMethodName,"_communitySubtypePlot_",Sys.Date(),".png"),
      width = 700, height = 1000);
  # vertex.color= V(undirGraph)$color,
  plot(undirGraph, layout=layout.fruchterman.reingold,vertex.label=V(undirGraph)$subtype,vertex.size=5,vertex.label.color="black",vertex.label.cex=.7,
                          vertex.color= V(undirGraph)$color, edge.arrow.size=3,main=paste0("Clusters deemed similar across ",length(unique(community_membership[,"studyNum"])), "studies.\n",
                                                                                           finalNumCommunities," pruned communities found using method ",clustMethodName, "."),xlab="color=community, label=subtype");
  
  
  dev.off();
  

  output <- list(undirGraph=undirGraph,community_membership=community_membership,finalNumCommunities=finalNumCommunities); 
  
  return(output);
  
  }

communityIHC_plots <- function(community_membership,biclust,receptorMatrix,origDataMatrices){
  
  numTotalSamplesInStudy <- list();
  for(s in 1:length(origDataMatrices)){
    
    numTotalSamplesInStudy[[s]] <- nrow(origDataMatrices[[s]]);
    
  }
  library("limma");
  library("ggplot2");
  warning("\nThis code assumes your biclusters have genes in the columns and patients in the rows.\n")
  #library("RColorBrewer");
  biclustMasterList <- list();
  biclustMasterList_origStudies <- list();
  
  for(b in 1:length(biclust)){
    
    biclustMasterList <- append(biclustMasterList, biclust[[b]]);
    biclustMasterList_origStudies <- append(biclustMasterList_origStudies,rep(b,length(biclust[[b]])));
  }
  names(biclustMasterList) <- c(1:length(biclustMasterList));
  
  #now create lists by communities. re-name to master biclust # (or underscore?)
  biclustCommunities <- list();
  communityNames <- unique(community_membership[,"community"]);
  communityStudyNums <- list();
  
  
  IHC_plot <- list();
  IHC_dfList <- list();
  
  for(c in 1:length(communityNames)){
    
    #community indices are same as that in the  list (same as those in IGP matrices)
    #CAREFUL: the levels make this weird...so just change to characters to match the names of the biclustMasterList.
    biclustCommunities[[c]] <- biclustMasterList[as.character(community_membership[which(community_membership[,"community"]==communityNames[c]), "clust"]) ];
    
    communityStudyNums[[c]] <- community_membership[which(community_membership[,"community"]==communityNames[c]), "studyNum"];
    
    IHC_dfList[[c]] <- data.frame();
    
    for(e in 1:length(biclustCommunities[[c]])){
      
      patIDs <- strsplit2(rownames(biclustCommunities[[c]][[e]]),"GSM")[,2];
      receptorStat <- receptorMatrix[na.omit(match(patIDs,receptorMatrix[,"patient_ID"])), "STATUS"];
      studyNum <- rep.int(communityStudyNums[[c]][[e]],times=length(patIDs));
      
      temp <- data.frame(patIDs,receptorStat,studyNum);
      IHC_dfList[[c]] <- rbind(IHC_dfList[[c]],temp); 
      
    }
    
    if(c >1){
      
      tmp <- data.frame(IHC_dfList[[c]],rep(c,times=nrow(IHC_dfList[[c]])));
      colnames(tmp)[ncol(tmp)] <- "community";
      
      IHC_dfMaster <- rbind(IHC_dfMaster,tmp);
      
    }else{
      
      tmp <- data.frame(IHC_dfList[[c]],rep(c,times=nrow(IHC_dfList[[c]])));
      colnames(tmp)[ncol(tmp)] <- "community";
      IHC_dfMaster <- tmp;
      
    }
    
    IHC_plot[[c]] <- ggplot(data=IHC_dfList[[c]],aes(x=receptorStat))+geom_histogram(aes(fill=receptorStat))+facet_grid(.~studyNum)+
      labs(x="Positive Receptor Status",title=paste0("Positive receptor status by dataset (batch) for community ",c));
    
    #COME BACK: we want the % of patients for this receptor status that are in this cluster, correct?
    #so we want geom_bar instead. basically a different data frame, right?
    #split() on each study. Then create a table() on each study to count # in each receptor status.
    #THEN: go back and do the same for each community/study cluster. this table / whole study table
    #will give you the percentages from each category. 
    #for each community: make a matrix of columns:  perc_ER, perc_,,,etc., along with the column: studyNum.
    #you'll fact on studyNum and make a bar plot for each one.
    #IHC_plot_perc[[c]] <- ggplot(data=IHC_dfList[[c]],aes(x=receptorStat))+geom_histogram(aes(fill=receptorStat))+facet_grid(.~studyNum)+
    # labs(x="Positive Receptor Status",title=paste0("Positive receptor status by dataset (batch) for community ",c));
    
  }
  
  
  output <- list(IHC_dfList=IHC_dfList,IHC_dfMaster=IHC_dfMaster,IHC_plot=IHC_plot,biclustCommunities=biclustCommunities,communityStudyNums=communityStudyNums);
  
  return(output);
}

#######wrapper analysis functions
#for now, requires you to have run the pam50 or intrinsic subtyping bc need "clust_numVar" - number of
#samples in each cluster.in the future, perhaps allow to just pass in the clustMatrixList and calculate it from that.
advanced_networkPlots <- function(analysisOutput=output$corr_analysis_restrictEdges,
                                  brewPal = c("Set3","Paired","Spectral","BrBG","PiYG","RdYlGn","RdYlBu","RdBu","PiYG","Set2"),
                                  saveDir="/home/kplaney/ISMB/",saveName="networks",colorCodes){
  
  #study summary.

  #number of edges
  #number of studies
  #number of communities.
  #number of clusters (nodes) and # started out with.
  #let's add size of each cluster
  
  network_stats <- c(length(unique(analysisOutput$communityMembership$attrDF$community)),
                     length(unique(analysisOutput$communityMembership$attrDF$clust)),
                     nrow(analysisOutput$adjMatricesList[[1]]),
                     nrow(analysisOutput$communityMembership$edgeDF),
                     length(unique(analysisOutput$communityMembership$attrDF$studyNum)));

  names(network_stats) <- c("numCommunities","numClusters","origNumClusters","numEdges",
                            "numStudies"); 
  
  if(!is.null(analysisOutput$subtypePlots$subtype_dfMaster)){
    
    clustSizes <- table(analysisOutput$subtypePlots$subtype_dfMaster[,"clust_numVar"]);
    #put in order of attributes
    clustSizes <- clustSizes[match(analysisOutput$communityMembership$attrDF$clust,names(clustSizes))];
    attrDF <- cbind(analysisOutput$communityMembership$attrDF,clustSizes);
    
  }else{
    
    attrDF <- analysisOutput$communityMembership$attrDF;
  }
  
  if(missing(colorCodes)){
    
  if(brewPal==FALSE){
    #make own color ramp.
    colorCodeF <- colorRampPalette(c("violet","blue","red"), bias = length(unique(attrDF[,"community"]))*2,
                                   space = "rgb", interpolate = "linear");
    #this will produce RGB representation in HEX, which cytoscape can take in.
    #if need plain old RGB: col2RGB(membership[,"colors"],alpha=FALSE);
    colorCodes <- colorCodeF(length(unique(attrDF[,"community"])));
    
  }else{
    #max is 11 colors
    colorCodes <- rev(brewer.pal(length(unique(attrDF[,"community"])),brewPal));
    
  }
  
  }
  attrDF$color <- as.character(attrDF$color);
  for(c in 1:length(unique(attrDF[,"community"]))){
    
    
    attrDF[which(attrDF$community==c),"color"] <- colorCodes[c];
    
  }
  
  undirGraph <- graph.data.frame(analysisOutput$communityMembership$edgeDF,directed=FALSE,vertices=attrDF)
  
  
  if(!is.null(analysisOutput$subtypePlots$subtype_dfMaster)){
    
  V(undirGraph)$size <- clustSizes;
  
  #save plots
  png(filename=paste0(saveDir,"/",saveName,"_communityPlot_scaledNodes_noLabels_",Sys.Date(),".png"),
      width = 700, height = 800);
  plot(undirGraph, layout=layout.fruchterman.reingold,vertex.label=rep.int("",times=nrow(attrDF)),vertex.size=V(undirGraph)$size/10,vertex.label.color="black",vertex.label.cex=.7,
       vertex.color= V(undirGraph)$color, edge.arrow.size=3);
  dev.off();
  

  #with study numbers
  png(filename=paste0(saveDir,"/",saveName,"_communityPlot_scaledNodes_studyNumlabels_",Sys.Date(),".png"),
      width = 700, height = 800);
  plot(undirGraph, layout=layout.fruchterman.reingold,vertex.label=V(undirGraph)$studyNum,vertex.size=V(undirGraph)$size/10,vertex.label.color="black",vertex.label.cex=.7,
       vertex.color= V(undirGraph)$color, edge.arrow.size=3);
  dev.off();
  
}
  #without relative sizes
  png(filename=paste0(saveDir,"/",saveName,"_communityPlot_unscaledNodes_nolabels_",Sys.Date(),".png"),
      width = 700, height = 800);
  plot(undirGraph, layout=layout.fruchterman.reingold,vertex.label=rep.int("",times=nrow(attrDF)),vertex.size=7,vertex.label.color="black",vertex.label.cex=.7,
       vertex.color= V(undirGraph)$color, edge.arrow.size=3);
  dev.off();
  #with study nums
  png(filename=paste0(saveDir,"/",saveName,"_communityPlot_unscaledNodes_studyNumlabels_",Sys.Date(),".png"),
      width = 700, height = 800);
  plot(undirGraph, layout=layout.fruchterman.reingold,vertex.label=V(undirGraph)$studyNum,vertex.size=7,vertex.label.color="black",vertex.label.cex=.7,
       vertex.color= V(undirGraph)$color, edge.arrow.size=3);
  dev.off();
  
  output <- list(undirGraph=undirGraph,attrDF=attrDF,network_stats=network_stats);
  return(output);
  
}


#####
computeEdgeMatrix <- function(clustRobust_output,restrictEdgesByStudy=TRUE,method=c("NN_distribution","diff_distribution"),
                             compareMetric=c("correlation","euclidean"),edge_trueCorr_thresh=.3,edge_fractFeatIntersect_thresh=.05,
                             edge_numFeatIntersect_thresh=20,edge_IGP_thresh=.7,sourceDir="/home/kplaney/gitRepos/IGP_network/igp_network/",
                             edge_simil_overNull_pvalue_thresh=.05,clustMethodName="clustRobust",saveDir="/home/kplaney/"){
  #setwd(sourceDir);
  #source("clust_robust.R");
  
  ###reusable code:
  clust_detailedNames <- clustRobust_output$clust_detailedNames;
  
  adjMatrix_fractFeat <- clustRobust_output$similMetrics$fractFeatIntersect;
  adjMatrix_numFeat <- clustRobust_output$similMetrics$numFeatIntersect;
  adjMatrix_simil <- clustRobust_output$similMetrics$simil;
  adjMatrix_pvalue <- clustRobust_output$simil_overNull_pvalue;
  
  if(restrictEdgesByStudy){
    library("limma");
    for(r in 1:nrow(adjMatrix_simil)){
      
      #we use > in the edge finding algorithm.
      compareIndices <- intersect(which(adjMatrix_simil[r,] > edge_trueCorr_thresh) , which(adjMatrix_pvalue[r,] < edge_simil_overNull_pvalue_thresh));
      
      if(length(compareIndices)>0){
        
        clust_names <- clust_detailedNames[compareIndices];
        clust_names_by_study <- split(clust_names,f=strsplit2(clust_names,"_")[,2],drop=TRUE);
        
        for(c in 1:length(clust_names_by_study)){
          
          if(length(clust_names_by_study[[c]])>1){
            #we have multiple matches
            orig_indices <- as.numeric(strsplit2(clust_names_by_study[[c]],"_")[,1]);
            #if a tie: then keep both - stay honest!
            bestMatch <- which(adjMatrix_simil[r,orig_indices]==max(adjMatrix_simil[r,orig_indices]));
            #remove this best index so that we don't set it to NA below.
            orig_indices <- orig_indices[-bestMatch];
            
            #if a tie, only 2 in orig_indices: length will be zero.
            if(length(orig_indices)>0){
              #null out these indices now so won't get compared
              adjMatrix_simil[r,orig_indices] <- NA;
              adjMatrix_pvalue[r,orig_indices] <- NA;
              
            }
            
          }
        }
        
        
      }
      # if(length(compareIndices)>0) 
    }
  }
  
  
  if(method=="NN_distribution"){
    
    adjMatrix_IGP <- clustRobust_output$IGP;
    
  }
  
  #if IGP: also have IGP.
  
  
  if(compareMetric=="correlation"){
    
    if(method=="diff_distribution"){
      
      adjMatricesList <- list(adjMatrix_simil,adjMatrix_fractFeat,adjMatrix_numFeat,
                              adjMatrix_pvalue);
      
      
      names(adjMatricesList) <- c("trueCorr","fractFeat","numFeat","corrOverNullpValue");
      
      
      thresholdVector <- c(edge_trueCorr_thresh,edge_fractFeatIntersect_thresh,edge_numFeatIntersect_thresh,
                           edge_simil_overNull_pvalue_thresh);
      
      #we want the edge pvalue to be insignificant so put a greater than...or just leave alone.
      threshDir <- c(">",">",">","<");
      
    }else if(method=="NN_distribution"){
      
      adjMatricesList <- list(adjMatrix_IGP,adjMatrix_simil,adjMatrix_fractFeat,adjMatrix_numFeat,
                              adjMatrix_pvalue);
      
      
      names(adjMatricesList) <- c("IGP","trueCorr","fractFeat","numFeat","corrOverNullpValue");
      
      
      thresholdVector <- c(edge_IGP_thresh,edge_trueCorr_thresh,edge_fractFeatIntersect_thresh,edge_numFeatIntersect_thresh,
                           edge_simil_overNull_pvalue_thresh);
      
      #we want the edge pvalue to be insignificant so put a greater than...or just leave alone.
      threshDir <- c(">",">",">",">","<",">");
      
    }else{
      
      stop("\nDid not pick an acceptable input for 'method' argument.")
    }
    
  }else if(compareMetric=="euclidean"){
    
    stop("haven't coded up yet!")
    #edge_distThresh <- 1000000;   
  }
  
  #check for erroneous/empty matrices fed in:
  for(a in 1:length(adjMatricesList)){
    
    if(is.null(dim(adjMatricesList[[a]]))){
      
      stop("\nYour adjacency matrix ",a, " is empty/is not 2D.\n");
    }
  }
  

 edgeResults <- computeEdges(adjMatricesList=adjMatricesList,thresholdVector=thresholdVector,threshDir=threshDir,saveDir=saveDir,fileTag=clustMethodName);

 cat("\n",nrow(edgeResults$edgeMatrix) ,"clust-clust linkages found after initial thresholding.\n");

 output <- list(adjMatricesList=adjMatricesList,edgeResults=edgeResults,thresholdVector=thresholdVector,threshDir=threshDir);

 return(output);

}

analyzeClustRobust_output <- function(clustRobust_output,fullMatrixList,clustMatrixList,clustMethodName,saveDir,method=c("diff_distribution","NN_dist"),compareMetric=c("correlation","euclidean"),
                                                  edge_trueCorr_thresh=.3,edge_fractFeatIntersect_thresh=.05,edge_numFeatIntersect_thresh=20,edge_IGP_thresh=.7,
                                                  edge_simil_overNull_pvalue_thresh=.05,genomeSize=20000,sourceDir="/home/kplaney/gitRepos/IGP_network/igp_network/",
                                                  commMethod="edgeBetween",minNumUniqueStudiesPerCommunity=3,
                                                  restrictEdgesByStudy=TRUE,centroidRData="/home/data/breast_microarrayDB/pam50_centroids.RData",centroidSetName="pam50",runGSEA=TRUE,
                                                  heatmaps=FALSE, makePlots=TRUE,saveGraphData=TRUE
){
  
  #setwd(sourceDir);
  #source("clust_robust.R");
  
  ###reusable code:
  clust_detailedNames <- clustRobust_output$clust_detailedNames;

  edgeOutput <- computeEdgeMatrix(clustRobust_output,restrictEdgesByStudy=restrictEdgesByStudy,method=method,
                                  compareMetric=compareMetric,edge_trueCorr_thresh=edge_trueCorr_thresh,
                                  edge_fractFeatIntersect_thresh=edge_fractFeatIntersect_thresh,clustMethodName=clustMethodName,
                                  edge_numFeatIntersect_thresh=edge_numFeatIntersect_thresh,edge_IGP_thresh=edge_IGP_thresh,
                                  edge_simil_overNull_pvalue_thresh=edge_simil_overNull_pvalue_thresh);
              
  adjMatricesList <-  edgeOutput$adjMatricesList;                                           
  edgeResults <- edgeOutput$edgeResults;
  thresholdVector <- edgeOutput$thresholdVector;
  threshDir <- edgeOutput$threshDir;
  #find communities 
  communityMembership <- findCommunities(edgeMatrix=edgeResults$edgeMatrix,edgeWeightMatrix=edgeResults$edgeWeightMatrix,
                                         clust_detailedNames=clust_detailedNames,fileTag=clustMethodName,saveDir=saveDir,
                                         clustMethodName=clustMethodName,minNumUniqueStudiesPerCommunity=minNumUniqueStudiesPerCommunity,commMethod=commMethod,
                                         makePlots=TRUE,saveGraphData=TRUE);
  
  cat("\n",communityMembership$numCommunities," final communities found.\n");
  
  if(runGSEA){
    #use default ref gene list.
    communityGeneSumm <- communityGeneSummary(community_membership=communityMembership$attrDF,
                                              biclust=clustMatrixList,detailedBiclustNames=clust_detailedNames,genomeSize=genomeSize);
    
    #pick the top scoring gene lists
    
    rowIndices <- c();
    
    for(c in 1:communityMembership$numCommunities){
      
      if(c==1){
        
        rowIndices <- which(communityGeneSumm$qvalues_communityDF[,c]<=.05);
        
      }else{
        
        rowIndices <- union(rowIndices,which(communityGeneSumm$qvalues_communityDF[,c]<=.05));
        
      }
      
      
    }
    
    bestRefGeneLists <- communityGeneSumm$refGeneListNames[rowIndices]; 
    finalGSEATable <- communityGeneSumm$qvalues_communityDF[rowIndices,];
    
  }else{
    
    communityGeneSumm <- NA;
    finalGSEATable <- NA;
    bestRefGeneLists <- NA;
    
  }
  
  if(heatmaps){
    
    community_hmaps <- communityHeatmaps(communityMembership$attrDF,biclust=clustMatrixList,
                                         dataMatrices=fullMatrixList,dataM_sampleCol=FALSE,
                                         geneCommIntersect=communityGeneSumm$geneCommIntersect);
    
  }else{
    
    community_hmaps=NA;
  }
  
  #if provide path to centroid RData object: also do centroid plots.
  if(!missing(centroidRData)){
    
    subtypePlots <- communitySubtypeBreakdown_plots(community_membership=communityMembership$attrDF,biclust=clustMatrixList,
                                                    origDataMatrices=fullMatrixList,saveDir=saveDir,clustMethodName=clustMethodName,
                                                    centroidRData="/home/data/breast_microarrayDB/pam50_centroids.RData");
    
    
    
    output <- list(adjMatricesList=adjMatricesList,edgeResults=edgeResults,thresholdVector=thresholdVector,threshDir=threshDir,
                   communityMembership=communityMembership,communityGeneSumm=communityGeneSumm,bestRefGeneLists=bestRefGeneLists,
                   finalGSEATable=finalGSEATable,subtypePlots=subtypePlots,community_hmaps=community_hmaps);
    
  }else{
    
    output <- list(adjMatricesList=adjMatricesList,edgeResults=edgeResults,thresholdVector=thresholdVector,threshDir=threshDir,
                   communityMembership=communityMembership,communityGeneSumm=communityGeneSumm,bestRefGeneLists=bestRefGeneLists,
                   finalGSEATable=finalGSEATable,community_hmaps=community_hmaps);
    
  }
  
  return(output);
  #save(clustRobust_output,file=paste0(saveDir,"/",clustMethodName,"_analysis.RData.gzip"),compress="gzip");
  
}

##############
#this code is also handy for some community numbers breakdowns.
#subtype data frame: outputted from clust robust analysis code.
#expects to find these EXACT column names: subtype studyNum clust_numVar community
pam50_pieCharts <- function(subtypeDataFrame,plotTitle,saveDir,saveName=plotTitle){
  
  df <- subtypeDataFrame;
  subtypeColorMatrix <- data.frame(c("#970eec","#ec0e58","#0e58ec","#7d7d7d","#0eec97"))
  #hmm can I order the colors in the same way each time?? may not really need to anyways.
  #must tranpose for ggplot to read it correctly and THEN add names for this to work.
  subtypeColorMatrix <- t(subtypeColorMatrix);
  names(subtypeColorMatrix) <- c("LumB","LumA","Her2","Basal","Normal");
  
  tmp <- split(df,f=df$community);
  
  
  
  for(t in 1:length(tmp)){
    
    if(t>1){
      
      subtypeBreakdowns <- rbind(subtypeBreakdowns,data.frame(table(tmp[[t]][,c("subtype")]),nrow(tmp[[t]]),t,length(unique(tmp[[t]][,c("studyNum")])),length(unique(tmp[[t]][,c("clust_numVar")]))));
      
    }else{
      
      subtypeBreakdowns <- data.frame(table(tmp[[t]][,c("subtype")]),nrow(tmp[[t]]),t,length(unique(tmp[[t]][,c("studyNum")])),length(unique(tmp[[t]][,c("clust_numVar")])));
    }
    
  }
  #COME BACK: add in # of studies?
  colnames(subtypeBreakdowns) <- c("subtype","number","communitySize","community","numStudyPerComm","numClustPerComm");
  subtypeBreakdowns <- cbind(subtypeBreakdowns,subtypeBreakdowns$number/subtypeBreakdowns$communitySize);
  colnames(subtypeBreakdowns) <- c("subtype","number","communitySize","community","numStudyPerComm","numClustPerComm","fract");
  library("ggplot2");
  
  p <- ggplot(subtypeBreakdowns, aes(x=1,y=fract, fill=subtype)) +facet_grid(.~community)+
    ggtitle(plotTitle) +
    geom_bar(stat="identity", color='black') +scale_fill_manual(values = subtypeColorMatrix)+
    # remove black diagonal line from legend
    guides(fill=guide_legend(override.aes=list(colour=NA))) + theme(strip.text.x = element_text(size = 17, colour = "black",face="bold"))+
    theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15));
  
  p <- p +coord_polar(theta='y') +theme(axis.ticks=element_blank(),  # the axis ticks
                                        axis.title=element_blank(),  # the axis labels
                                        axis.text.y=element_blank(), # the 0.75, 1.00, 1.25 labels.
                                        axis.text.x=element_blank())
  
  png(filename = paste0(saveDir,"/",saveName,".png"),
      width = 1000, height = 200);
  plot(p);
  dev.off();
  
  #just take first match for each community - values will be same across an entire community.
  communityStats <- subtypeBreakdowns[match(unique(subtypeBreakdowns$community),subtypeBreakdowns$community), ]
  communityStats <- communityStats[ ,c("communitySize","community","numStudyPerComm","numClustPerComm")];
  rownames(communityStats) <- communityStats$community;
  output <- list(pieChartPlot=p, communityStats=communityStats,subtypeBreakdowns=subtypeBreakdowns);
  
  return(output);
  
}