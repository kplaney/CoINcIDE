library("Biobase")
library("curatedOvarianData")
datasetNames <- data(package="curatedOvarianData")
datasetNames <- datasetNames[3]$results[,"Item"]



#odd way to get a dataset list...I just downloaded the source file
#and pull out their code.
#hmmm...is this the case for all genes?? WTF?
expandProbesets <- function (eset, sep = "///"){
  x <- lapply(featureNames(eset), function(x) strsplit(x, sep)[[1]])
  eset <- eset[order(sapply(x, length)), ]
  x <- lapply(featureNames(eset), function(x) strsplit(x, sep)[[1]])
  idx <- unlist(sapply(1:length(x), function(i) rep(i, length(x[[i]]))))
  xx <- !duplicated(unlist(x))
  idx <- idx[xx]
  x <- unlist(x)[xx]
  eset <- eset[idx, ]
  #Katie's notes: featureNames is from Biobase, for eset objects.
  featureNames(eset) <- x
  eset
}


add.surv.y <- function(X) Surv(X$days_to_death, X$vital_status=="deceased")
## -----------------------------------------------------------------------------
##load the esets
## -----------------------------------------------------------------------------

delim <- ":" 

data(list=data(package=package.name)[[3]][,3])

strEsets <- ls(pattern="^.*_eset$")
#ON SERVER: just download the tarbell to get package, code
#cd /home/kplaney/R/x86_64-redhat-linux-gnu-library
#sudo wget http://www.bioconductor.org/packages/release/data/experiment/src/contrib/curatedOvarianData_1.3.5.tar.gz
#tar xvzf curatedOvarianData_1.3.5.tar.gz
# strEsets <- list.files("/home/kplaney/R/x86_64-redhat-linux-gnu-library/curatedOvarianData/data/")
# #datalist is a text file...shouldn't be any other text files in this directory, only .rda files.
# strEsets <- strEsets[-which(strEsets=="datalist")]
# library("limma")
# names(strEsets) <- strsplit2(strEsets,split=".rda")[,1]
# #now load all of this data. equivalent of data() call above if can install package.
# #is having trouble reading from a connection?? what's up...
# #WTF..may have to manually step through this one..
# #weird...now suddently dataset 16 isn't working? why only THIS one?
# #makes me wonder if there is a space issue on the server.
# #I had to re-download the datasets...
# for (strEset in strEsets){
#   
#   load(paste0("/home/kplaney/R/x86_64-redhat-linux-gnu-library/curatedOvarianData/data/",strEset))
#   
# }
# #now call strEset again like above:
# strEsets <- ls(pattern="^.*_eset$")

##now back to common code for server or local computer:
esets <- list()
for (strEset in strEsets){
  
  eset <- get(strEset)
  ##Split out rows without unique gene name
  eset <- expandProbesets(eset)
  esets[[strEset]] <- eset
  
}



#pData(esets) already contains survival data.

#remove the extra TCGA datasets
#TCGA.RNASeqV2_eset                           Integrated genomic analyses of ovarian carcinoma.
#TCGA.mirna.8x15kv2_eset                      Integrated genomic analyses of ovarian carcinoma.
#TCGA_eset                                    Integrated genomic analyses of ovarian carcinoma.
indicesRemove <- na.omit(match(c("TCGA.RNASeqV2_eset","TCGA.mirna.8x15kv2_eset"),names(esets)))
esets <- esets[-indicesRemove]

#save for easier access next time:
#save(esets, file="/home/kplaney/curatedOvarianData_esetList.RData.gzip",compress="gzip")
#CHECK: remove any healthy samples. do this for ALL studies in case some have health tissue.
#for example: TCGA has some normal samples
table(pData(esets[["TCGA_eset"]])[,"batch"],pData(esets[["TCGA_eset"]])[,"sample_type"])
indicesRemove <- c()

for(e in 1:length(esets)){
  
  #are these all primary tumor? I wish they had kept the full TCGA labels.
  pDat <- pData(esets[[e]])[which(pData(esets[[e]])["sample_type"]=="tumor"),]
  expr <- exprs(esets[[e]])[ ,which(pData(esets[[e]])["sample_type"]=="tumor")]
  
  #if no samples left: remove
  if(ncol(expr)==0){
    
    indicesRemove <- append(indicesRemove,e)
    
  }else{
    
    exprs(esets[[e]]) <- expr
    pData(esets[[e]]) <- pDat
    protocolData(esets[[e]])@data <- data.frame(row.names=colnames(exprs(esets[[e]])))
    #check: make sure is still a valid object
    validObject(esets[[e]])
    
  }
  
}
message(paste("dataset indices to remove : ",indicesRemove))
#gut check: in TCGA, should be only tumor samples now!
if(any(unique(pData(esets[["TCGA_eset"]])[,"sample_type"]) != "tumor")){
  
  stop("\nNormal filtering appeared to not work.")
  
}

#check? any duplicated samples?
#it looks like unique_patient_ID ids duplicated samples, NOT the actual expression colnames?
#a little annoying - re-label the colnames to be this unique ID, if duplicated samples are found
#but esets can't have duplicated column names.
#in the end, looks like no samples were duplicated? 
for(e in 1:length(esets)){
  
  if(!all(is.na(pData(esets[[e]])[ ,"unique_patient_ID"]))){
    #for certain datasets: this is actually all NA values. so don't use then!
    if(any(duplicated(pData(esets[[e]])[,"unique_patient_ID"]))){
      
      #if(any(duplicated(colnames(exprs(esets[[e]]))))){
      
      
      message(paste0("Duplicated samples in study ",names(esets)[e]))
      
    }
    
  }
  
}
#}
#remove datasets that are too small? allow them for now.
#now process, then grab meta-features
#sounds like there SHOULD be some duplicated samples  in here.
#want very lowly varying genes to be removed for Combat.
#COME BACK: do variance filtering by batch for Combat.
esets <- processExpressionSetList(exprSetList=esets,outputFileDirectory="/home/kplaney/ovarian_analysis/",
                                  minVar=.001,featureDataFieldName="gene",uniquePDataID="unique_patient_ID")

#just strictly limit TCGA dataset now - bc of combat
#esets[["TCGA_eset"]] <- processExpressionSet(exprSet,esets[["TCGA_eset"]],
#                                 outputFileDirectory="/home/kplaney/ovarian_analysis/",
#                                  minVar=.1,featureDataFieldName="gene",uniquePDataID="unique_patient_ID")
#look for batches - break these apart?
#but it only looks like TCGA has actual site batches -
#the rest look like the date the samples were collected?
#in the R documentation it says 
#"Hybridization date when Affymetrix CEL files are available"
#perhaps TOO granular for our purposes.

for(e in 1:length(esets)){
  
  if(!all(is.na(pData(esets[[e]])[,"batch"]))){
    
    cat("\n",names(esets[[e]]))
    
  }
  
}


#only split up TCGA into batches
#na.omit: if any NA batch variables, then unique will return one "NA" 
TCGAbatches <- na.omit(unique(pData(esets[["TCGA_eset"]])[,"batch"]))
#any samples with NA batch? these will be thrown out.
which(is.na(pData(esets[["TCGA_eset"]])[,"batch"]))
#let's see if batches seem to have uneven distributions of known subtypes
#and/or grade...that will make me hesitant to use batch correction.
#hmm...so most are serious?
table(pData(esets[["TCGA_eset"]])[,"batch"],pData(esets[["TCGA_eset"]])[,"histological_type"])
unique(pData(esets[["TCGA_eset"]])[,"histological_type"])
#yep...only 10 a NA, rest are serous
length(which(is.na(pData(esets[["TCGA_eset"]])[,"histological_type"])))
#this matches the Nature paper saying the samples are only serous: http://www.nature.com/nature/journal/v474/n7353/full/nature10166.html

#and predominantly late-stage for all batches
table(pData(esets[["TCGA_eset"]])[,"batch"],pData(esets[["TCGA_eset"]])[,"summarystage"])
#same trend for substage
table(pData(esets[["TCGA_eset"]])[,"batch"],pData(esets[["TCGA_eset"]])[,"substage"])

outcomesAndCovariates <- pData(esets[["TCGA_eset"]])[,"batch",drop=FALSE]
rownames(outcomesAndCovariates) <- colnames(exprs(esets[["TCGA_eset"]]))
exprMatrix <- exprs(esets[["TCGA_eset"]])
rownames(exprMatrix) <- featureNames(esets[["TCGA_eset"]])
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_batchCorrection.R")
batchCorrect <- batchNormalization(countsMatrixNoNANoDup=exprMatrix,outcomesAndCovariates=outcomesAndCovariates,MinInBatch=4,combatModelFactorName=NULL,pvalueThresh=.05,batchColName="batch",outputFile="combatoutput.txt")

esets[["TCGA_eset"]] <- esets[["TCGA_eset"]][rownames(batchCorrect$GEN_Data_Corrected), colnames(batchCorrect$GEN_Data_Corrected)]
featureNames(esets[["TCGA_eset"]] ) <- rownames(batchCorrect$GEN_Data_Corrected)
validObject(esets[["TCGA_eset"]])




#now format just as a list of data matrices.
dataMatrixList <- exprSetListToMatrixList(esets,featureDataFieldName="gene")

names(dataMatrixList) <- names(esets)
#run meta-feature analysis.
metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),numFeatSelectByGlobalRank=1000,
                                                   numTopFeatFromEachDataset=20,fractNATopFeatAllowedPerDataset=.2,selectMethod=c("median"),
                                                   outputFile="/home/kplaney/ovarian_analysis//selectFeaturesMetaVarianceOut.txt")

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))
#remove datasets with too many missing top gene features
if(length(metaFeatures$datasetListIndicesToRemove)>0){
  
  dataMatrixList <- dataMatrixList[-metaFeatures$datasetListIndicesToRemove]
  
}


clustFeaturesList <- list()
for(d in 1:length(dataMatrixList)){
  
  clustFeaturesList[[d]] <- metaFeatures$finalFeatures
  
}

output <- list(dataMatrixList=dataMatrixList,clustFeaturesList=clustFeaturesList)
save(output,file="/home/kplaney/ovarian_analysis/curatedOvarianData_proc.RData.gzip",compress="gzip")

load("/home/kplaney/ovarian_analysis/curatedOvarianData_proc.RData.gzip")
dataMatrixList <- output$dataMatrixList
clustFeaturesList <- output$clustFeaturesList
source("/home/kplaney/gitRepos/CoINcIDE/coincide/CoINcIDE/R/CoINcIDE_cluster.R")

clusterOut_kmeansConsensus <- consensusClusterMatrixList(dataMatrixList,clustFeaturesList,maxNumClusters = 15, numSims=100, pItem=0.8, pFeature=1, clusterAlg=c("km"),
                                                         corUse=c("everything"),
                                                         distMethod=c("euclidean"),
                                                         minMeanClustConsensus=.7,outputFile=paste0("/home/kplaney/ovarian_analysis/kmeansConsensus_",Sys.Date(),".txt"))

#use average linkage for hclust: may work better with unevenly sized clusters?
#looks like ward's or average is best: http://www.stat.sc.edu/~hitchcock/compare_hier_fda.pdf, http://www.biomedcentral.com/1471-2105/13/S10/S7
#BUT ward.D and ward.D2 really overshoot k; so use average instead.
clusterOut_hlustEuclidConsensus <- consensusClusterMatrixList(dataMatrixList,clustFeaturesList,maxNumClusters = 15, numSims=100, pItem=0.8, pFeature=1, clusterAlg=c("hc"),
                                                         corUse=c("everything"),
                                                         distMethod=c("euclidean"),hclustAlgorithm="average",
                                                         minMeanClustConsensus=.7,outputFile=paste0("/home/kplaney/ovarian_analysis/hclustEuclidConsensus_",Sys.Date(),".txt"))


clusterOut_hlustSpearmanConsensus <- consensusClusterMatrixList(dataMatrixList,clustFeaturesList,maxNumClusters = 15, numSims=100, pItem=0.8, pFeature=1, clusterAlg=c("hc"),
                                                              corUse=c("everything"),
                                                              distMethod=c("spearman"),hclustAlgorithm="ward.D2",
                                                              minMeanClustConsensus=.7,outputFile=paste0("/home/kplaney/ovarian_analysis/hclustSpearmanConsensus_",Sys.Date(),".txt"))
