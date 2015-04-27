load("/home/kplaney/breast_analysis/kmeansConsensus_curatedBreastData.RData.gzip")
names(kmeansConsensus)

load("/home/kplaney/breast_analysis//curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
names(dataMatrixList)

load("/home/kplaney/breast_analysis//curatedBreastData_esets_proc.RData.gzip")


library("Biobase")
#this dataset has 3-5 clusters. NOTE: "normal" means normal as classified by pam50 status, NOT
#that it's actually normal tissue. pam50 centroids were trained on a bit of a skewed dataset.
table(pData(esets[["study_25055_GPL96_MDACC_M"]])$pam50)

#hmm...a lot are close to zero:

PAC <- unlist(kmeansConsensus$PAC[[which(names(dataMatrixList)=="study_GSE2226")]])

#perhaps round to hundreth decimal, take max?? but this would suggest seven clusters....
#but splitting hairs around the 100th decimal point makes no sense!
possibleKs <- which(round(PAC*100)==min(round(PAC*100),na.rm=TRUE))

maxPossibleK <- possibleKs[length(possibleKs)]

maxPossibleK

#now try I-SPY
table(pData(esets[["study_22226_GPL1708_all"]])$pam50)

PAC <- unlist(kmeansConsensus$PAC[[which(names(dataMatrixList)=="study_22226_GPL1708_all")]])

#perhaps round to hundreth decimal, take max?? but this would suggest seven clusters....
#but splitting hairs around the 100th decimal point makes no sense!
possibleKs <- which(round(PAC*100)==min(round(PAC*100),na.rm=TRUE))

maxPossibleK <- possibleKs[length(possibleKs)]

maxPossibleK
