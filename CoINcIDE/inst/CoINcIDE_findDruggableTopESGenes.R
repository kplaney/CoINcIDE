
load("/home/data/genomeReferences/annotationHuman/druggableGenomeList_updatedGeneSymbols.RData.gzip");

#breast 200
ES_topGenes <- read.table("/home/kplaney/breast_analysis_withTop20Genes/breast_200_features_2015-05-08/breast_200_features_summaryGenes_ESpos_thresh_0.5_2015-05-08.txt")

#Cool!! 50 genes.
drugGenes <- druggableGenes[na.omit(match(rownames(ES_topGenes),druggableGenes))]

#look at ES above 1 only?
thresh <- 1

ES_final <- ES_topGenes[drugGenes, ]
removeIndices <- c()
for(e in 1:nrow(ES_final)){
  
  if(all(ES_final[e,] < thresh,na.rm=TRUE)){
    #remove this gene
    removeIndices <- c(removeIndices,e)
    
  }
  
}
if(length(ES_final)>0){
  
  ES_final <- ES_final[-removeIndices, ]
  
}
#save this table.
write.table(ES_final, "/home/kplaney/breast200_druggable.txt",col.names=FALSE,quote=FALSE,row.names=TRUE)
#breast 2000
ES_topGenes <- read.table("/home/kplaney/breast_analysis_withTop20Genes/breast_2000_features_2015-05-18/breast_2000_features_summaryGenes_ESpos_thresh_0.5_2015-05-18.txt")

#Cool!! 251 genes.
drugGenes <- druggableGenes[na.omit(match(rownames(ES_topGenes),druggableGenes))]
ES_final <- ES_topGenes[drugGenes, ]
removeIndices <- c()
for(e in 1:nrow(ES_final)){
  
  if(all(ES_final[e,] < thresh,na.rm=TRUE)){
    #remove this gene
    removeIndices <- c(removeIndices,e)
    
  }
  
}
if(length(ES_final)>0){
  
  ES_final <- ES_final[-removeIndices, ]
  
}

write.table(ES_final, "/home/kplaney/breast2000_druggable.txt",col.names=FALSE,quote=FALSE,row.names=TRUE)
##then look at ovarian:
#200
ES_topGenes <- read.table("/home/kplaney/ovarian_analysis_withTop20Genes/ovarian_200_features_2015-05-05/ovarian_200_features_summaryGenes_ESpos_thresh_0.5_2015-05-05.txt")
#Cool!! 40 genes.
drugGenes <- druggableGenes[na.omit(match(rownames(ES_topGenes),druggableGenes))]

ES_final <- ES_topGenes[drugGenes, ]
removeIndices <- c()
for(e in 1:nrow(ES_final)){
  
  if(all(ES_final[e,] < thresh,na.rm=TRUE)){
    #remove this gene
    removeIndices <- c(removeIndices,e)
    
  }
  
}
if(length(ES_final)>0){
  
  ES_final <- ES_final[-removeIndices, ]
  
}

write.table(ES_final, "/home/kplaney/ovarian200_druggable.txt",col.names=FALSE,quote=FALSE,row.names=TRUE)
##2000
ES_topGenes <- read.table("/home/kplaney/ovarian_analysis_withTop20Genes/ovarian_2000_features_2015-05-08/ovarian_2000_features_summaryGenes_ESpos_thresh_0.5_2015-05-08.txt")
#Cool!! 350 genes.
drugGenes <- druggableGenes[na.omit(match(rownames(ES_topGenes),druggableGenes))]

ES_final <- ES_topGenes[drugGenes, ]
removeIndices <- c()
for(e in 1:nrow(ES_final)){
  
  if(all(ES_final[e,] < thresh,na.rm=TRUE)){
    #remove this gene
    removeIndices <- c(removeIndices,e)
    
  }
  
}
if(length(ES_final)>0){
  
  ES_final <- ES_final[-removeIndices, ]
  
}

write.table(ES_final, "/home/kplaney/ovarian2000_druggable.txt",col.names=FALSE,quote=FALSE,row.names=TRUE)