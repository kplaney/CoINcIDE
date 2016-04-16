

###no PAM50 analyses
#load("/home/data/genomeReferences/annotationHuman/druggableGenomeList_updatedGeneSymbols.RData.gzip");
load("~/druggableGenomeList_updatedGeneSymbols.RData.gzip");


#breast 50
saveDir <- "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/breast_50genesNoPAM50NoTop20_pear_meanCent_MM42015-12-18_2015-12-18/"
ES_topGenes <- read.table(paste0(saveDir,"breast_50genesNoPAM50NoTop20_pear_meanCent_MM42015-12-18_summaryGenes_ESpos_thresh_0.5_2015-12-18.txt"))

drugGenes <- druggableGenes[na.omit(match(rownames(ES_topGenes),druggableGenes))]

#look at ES above 1 only? no - might as well keep above .5 data.
thresh <- .5

ES_final <- ES_topGenes[drugGenes, ]
removeIndices <- c()
#first 1/2 of rows are effect size.
for(e in 1:nrow(ES_final)){
  
  if(all(ES_final[e,c(1:(ncol(ES_final)/2))] < thresh,na.rm=TRUE)){
    #remove this gene
    removeIndices <- c(removeIndices,e)
    
  }
  
}
if(length(removeIndices)>0){
  
  ES_final <- ES_final[-removeIndices, ]
  
}
#save this table.
write.table(ES_final, paste0(saveDir,"50F_NoPAM50_druggable.txt"),col.names=FALSE,quote=FALSE,row.names=TRUE)


#breast 200
saveDir <- "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/breast_200FNoPAM50_pear_meanCent_MM0.5_2015-12-18/"
ES_topGenes <- read.table(paste0(saveDir,"breast_200FNoPAM50_pear_meanCent_MM0.5_summaryGenes_ESpos_thresh_0.5_2015-12-18.txt"))

#Cool!! over 50 genes.
drugGenes <- druggableGenes[na.omit(match(rownames(ES_topGenes),druggableGenes))]

#look at ES above 1 only? no - might as well keep above .5 data.
thresh <- .5

ES_final <- ES_topGenes[drugGenes, ]
removeIndices <- c()
#first 1/2 of rows are effect size.
for(e in 1:nrow(ES_final)){
  
  if(all(ES_final[e,c(1:(ncol(ES_final)/2))] < thresh,na.rm=TRUE)){
    #remove this gene
    removeIndices <- c(removeIndices,e)
    
  }
  
}
if(length(removeIndices)>0){
  
  ES_final <- ES_final[-removeIndices, ]
  
}
#save this table.
write.table(ES_final, paste0(saveDir,"200F_NoPAM50_druggable.txt"),col.names=FALSE,quote=FALSE,row.names=TRUE)


##2000 genes:
saveDir <- "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/breast_2000FNoPAM50_pear_meanCent_MM0.5_2015-12-18/"
ES_topGenes <- read.table(paste0(saveDir,"breast_2000FNoPAM50_pear_meanCent_MM0.5_summaryGenes_ESpos_thresh_0.5_2015-12-18.txt"))

#Cool!! over 50 genes.
drugGenes <- druggableGenes[na.omit(match(rownames(ES_topGenes),druggableGenes))]

#look at ES above 1 only? no - might as well keep above .5 data.
thresh <- .5

ES_final <- ES_topGenes[drugGenes, ]
removeIndices <- c()
#first 1/2 of rows are effect size.
for(e in 1:nrow(ES_final)){
  
  if(all(ES_final[e,c(1:(ncol(ES_final)/2))] < thresh,na.rm=TRUE)){
    #remove this gene
    removeIndices <- c(removeIndices,e)
    
  }
  
}
if(length(removeIndices)>0){
  
  ES_final <- ES_final[-removeIndices, ]
  
}
#save this table.
write.table(ES_final, paste0(saveDir,"2000F_NoPAM50_druggable.txt"),col.names=FALSE,quote=FALSE,row.names=TRUE)


####old breast analyses:
#load("/home/data/genomeReferences/annotationHuman/druggableGenomeList_updatedGeneSymbols.RData.gzip");
load("~/druggableGenomeList_updatedGeneSymbols.RData.gzip");
#breast 200
ES_topGenes <- read.table("/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/breast_278genes_pear_meanCent_MM3_2015-07-21/breast_278genes_pear_meanCent_MM3_summaryGenes_ESpos_thresh_0.5_2015-07-21.txt")

#Cool!! over 50 genes.
drugGenes <- druggableGenes[na.omit(match(rownames(ES_topGenes),druggableGenes))]

#look at ES above 1 only? no - might as well keep above .5 data.
thresh <- .5

ES_final <- ES_topGenes[drugGenes, ]
removeIndices <- c()
#first 1/2 of rows are effect size.
for(e in 1:nrow(ES_final)){
  
  if(all(ES_final[e,c(1:(ncol(ES_final)/2))] < thresh,na.rm=TRUE)){
    #remove this gene
    removeIndices <- c(removeIndices,e)
    
  }
  
}
if(length(removeIndices)>0){
  
  ES_final <- ES_final[-removeIndices, ]
  
}
#save this table.
write.table(ES_final, "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/breast_278genes_pear_meanCent_MM3_2015-07-21/breast278_druggable.txt",col.names=FALSE,quote=FALSE,row.names=TRUE)

#breast 2000
ES_topGenes <- read.table("/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/breast_2052genes_pear_meanCent_MM4_2015-07-21/breast_2052genes_pear_meanCent_MM4_summaryGenes_ESpos_thresh_0.5_2015-07-21.txt")

#Cool!! 251 genes.
drugGenes <- druggableGenes[na.omit(match(rownames(ES_topGenes),druggableGenes))]
#look at ES above 1 only? no - might as well keep above .5 data.
thresh <- .5

ES_final <- ES_topGenes[drugGenes, ]
removeIndices <- c()
#first 1/2 of rows are effect size.
for(e in 1:nrow(ES_final)){
  
  if(all(ES_final[e,c(1:(ncol(ES_final)/2))] < thresh,na.rm=TRUE)){
    #remove this gene
    removeIndices <- c(removeIndices,e)
    
  }
  
}
if(length(removeIndices)>0){
  
  ES_final <- ES_final[-removeIndices, ]
  
}

write.table(ES_final, "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/breast_2052genes_pear_meanCent_MM4_2015-07-21/breast2052_druggable.txt",col.names=FALSE,quote=FALSE,row.names=TRUE)


##then look at ovarian:
#200
ES_topGenes <- read.table("/home/ywrfc09/ovarian_analysis/metaRankWithTop20Genes/240F_pear_meanCent_MM5_2015-07-21/240F_pear_meanCent_MM5_summaryGenes_ESpos_thresh_0.5_2015-07-21.txt")
#Cool!! 40 genes.
drugGenes <- druggableGenes[na.omit(match(rownames(ES_topGenes),druggableGenes))]


ES_final <- ES_topGenes[drugGenes, ]
removeIndices <- c()
#first 1/2 of rows are effect size.
for(e in 1:nrow(ES_final)){
  
  if(all(ES_final[e,c(1:(ncol(ES_final)/2))] < thresh,na.rm=TRUE)){
    #remove this gene
    removeIndices <- c(removeIndices,e)
    
  }
  
}
if(length(removeIndices)>0){
  
  ES_final <- ES_final[-removeIndices, ]
  
}

write.table(ES_final, "/home/ywrfc09/ovarian_analysis/metaRankWithTop20Genes/240F_pear_meanCent_MM5_2015-07-21/ovarian240_druggable.txt",col.names=FALSE,quote=FALSE,row.names=TRUE)

##2000
ES_topGenes <- read.table("/home/ywrfc09/ovarian_analysis/metaRankWithTop20Genes/2014F_pear_meanCent_MM5_2015-07-21/2014F_pear_meanCent_MM5_summaryGenes_ESpos_thresh_0.5_2015-07-21.txt")
#Cool!! 350 genes.
drugGenes <- druggableGenes[na.omit(match(rownames(ES_topGenes),druggableGenes))]


ES_final <- ES_topGenes[drugGenes, ]
removeIndices <- c()
#first 1/2 of rows are effect size.
for(e in 1:nrow(ES_final)){
  
  if(all(ES_final[e,c(1:(ncol(ES_final)/2))] < thresh,na.rm=TRUE)){
    #remove this gene
    removeIndices <- c(removeIndices,e)
    
  }
  
}
if(length(removeIndices)>0){
  
  ES_final <- ES_final[-removeIndices, ]
  
}
write.table(ES_final, "/home/ywrfc09/ovarian_analysis/metaRankWithTop20Genes/2014F_pear_meanCent_MM5_2015-07-21/ovarian2014_druggable.txt",col.names=FALSE,quote=FALSE,row.names=TRUE)
