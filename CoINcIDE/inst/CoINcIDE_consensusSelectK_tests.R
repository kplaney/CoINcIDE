#original run with gap test:
load(#original run:
  load("/home/data/breast_microarrayDB/output/pam50_subtypes/pam50Short_kmeans_allStudies.RData.gzip")

#pam50 run:
load("/home/kplaney/breast_analysis/kmeansConsensus_curatedBreastData.RData.gzip")
names(kmeansConsensus)


load("/home/kplaney/breast_analysis//curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
names(dataMatrixList)

#keep only matching studies - drops 1.
gapKmeans_pam50Short <- output$clustMatrixList_pam50Short[na.omit(match(names(dataMatrixList),names(output$clustMatrixList_pam50Short)))]
all(names(gapKmeans_pam50Short)==names(dataMatrixList))

load("/home/kplaney/breast_analysis//curatedBreastData_esets_proc.RData.gzip")


library("Biobase")

#look at MD-Anderson dataset. 
#how many K did gap test say? 7:
length(gapKmeans_pam50Short[["study_25055_GPL96_MDACC_M"]])
#this dataset has 3-5 clusters. NOTE: "normal" means normal as classified by pam50 status, NOT
#that it's actually normal tissue. pam50 centroids were trained on a bit of a skewed dataset.
table(pData(esets[["study_25055_GPL96_MDACC_M"]])$pam50)

#hmm...a lot are close to zero:
PAC <- unlist(kmeansConsensus$PAC[[which(names(dataMatrixList)=="study_25055_GPL96_MDACC_M")]])

PAC
#perhaps round to hundreth decimal, take max?? this would suggest seven clusters....
#splitting hairs around the 100th decimal point level makes no sense!
#look at all the ones that now match this rounded MIN
possibleKs <- which(round(PAC*100)==min(round(PAC*100),na.rm=TRUE))

maxPossibleK <- possibleKs[length(possibleKs)]

maxPossibleK

#here, we see that even k=9 is quite high!
meanConsensusByClust <- unlist(kmeansConsensus$meanConsensusClusterByK[[which(names(dataMatrixList)=="study_25055_GPL96_MDACC_M")]])

meanConsensusByClust

#look at all the ones that now match this rounded MAX
#looks like we need to go 2 decimal places for this one.
#need to add plus 1!
possibleKsMean <- which(round(meanConsensusByClust*10)==max(round(meanConsensusByClust*10),na.rm=TRUE)) +1
#need to add plus 1!
maxPossibleKMean <- possibleKsMean[length(possibleKsMean)]
#says 9!
maxPossibleKMean


#here, we see that even k=9 is quite high, but k=2 still wins out. perhaps too sensitive:
minConsensusByClust <- unlist(kmeansConsensus$minConsensusClusterByK[[which(names(dataMatrixList)=="study_25055_GPL96_MDACC_M")]])

minConsensusByClust

#look at all the ones that now match this rounded MAX
#looks like we need to go 4 decimal places for this one.
#need to add plus 1!
possibleKsMin <- which(round(minConsensusByClust*10)==max(round(minConsensusByClust*10),na.rm=TRUE)) +1
#need to add plus 1!
maxPossibleKMin <- possibleKsMin[length(possibleKsMin)]
#here will say 3: is fairly conservative.
maxPossibleKMin

##########
#now try I-SPY. this was not normalized by me.


#how many K did gap test say? 3:
length(gapKmeans_pam50Short[["study_22226_GPL1708_all"]])
#this dataset has 3-5 clusters. NOTE: "normal" means normal as classified by pam50 status, NOT
#that it's actually normal tissue. pam50 centroids were trained on a bit of a skewed dataset.
table(pData(esets[["study_22226_GPL1708_all"]])$pam50)

#hmm...a lot are close to zero:
PAC <- unlist(kmeansConsensus$PAC[[which(names(dataMatrixList)=="study_22226_GPL1708_all")]])

PAC
#perhaps round to hundreth decimal, take max?? this would suggest seven clusters....
#splitting hairs around the 100th decimal point level makes no sense!
#look at all the ones that now match this rounded MIN
possibleKs <- which(round(PAC*100)==min(round(PAC*100),na.rm=TRUE))

maxPossibleK <- possibleKs[length(possibleKs)]
#gives 4: about what we're looking for. in this scenario, 4 is pretty clearly the best cluster number.
maxPossibleK

#here, we see that even k=9 is quite high!
meanConsensusByClust <- unlist(kmeansConsensus$meanConsensusClusterByK[[which(names(dataMatrixList)=="study_22226_GPL1708_all")]])

meanConsensusByClust

#look at all the ones that now match this rounded MAX
#looks like we need to go 2 decimal places for this one.
#need to add plus 1!
possibleKsMean <- which(round(meanConsensusByClust*10)==max(round(meanConsensusByClust*10),na.rm=TRUE)) +1
#need to add plus 1!
maxPossibleKMean <- possibleKsMean[length(possibleKsMean)]
#says 2: overly conservative
maxPossibleKMean


#here, we see that even k=9 is quite high, but k=2 still wins out. perhaps too sensitive:
minConsensusByClust <- unlist(kmeansConsensus$minConsensusClusterByK[[which(names(dataMatrixList)=="study_22226_GPL1708_all")]])

minConsensusByClust

#look at all the ones that now match this rounded MAX
#looks like we need to go 4 decimal places for this one.
#need to add plus 1!
possibleKsMin <- which(round(minConsensusByClust*10)==max(round(minConsensusByClust*10),na.rm=TRUE)) +1
#need to add plus 1!
maxPossibleKMin <- possibleKsMin[length(possibleKsMin)]
#here says 4:
maxPossibleKMin

########
#now try one with little less heterogeneity: study_12093_GPL96_all
#how many K did gap test say? 1:
length(gapKmeans_pam50Short[["study_12093_GPL96_all"]])
#no pam50 info, but do have ER status:
#ALL are ER positive:
table(pData(esets[["study_12093_GPL96_all"]])$ER_preTrt)
#but we don't actually know about HER2 status: (try my pam50 centroids?)
table(pData(esets[["study_12093_GPL96_all"]])$HER2_preTrt)


#hmm...a lot are close to zero:
PAC <- unlist(kmeansConsensus$PAC[[which(names(dataMatrixList)=="study_12093_GPL96_all")]])

PAC
#perhaps round to hundreth decimal, take max?? this would suggest seven clusters....
#splitting hairs around the 100th decimal point level makes no sense!
#look at all the ones that now match this rounded MIN
possibleKs <- which(round(PAC*100)==min(round(PAC*100),na.rm=TRUE))

maxPossibleK <- possibleKs[length(possibleKs)]
#gives 3: (and it's pretty clear:)
maxPossibleK

#here, we see that even k=9 is quite high!
meanConsensusByClust <- unlist(kmeansConsensus$meanConsensusClusterByK[[which(names(dataMatrixList)=="study_12093_GPL96_all")]])

meanConsensusByClust

#look at all the ones that now match this rounded MAX
#looks like we need to go 2 decimal places for this one.
#need to add plus 1!
possibleKsMean <- which(round(meanConsensusByClust*10)==max(round(meanConsensusByClust*10),na.rm=TRUE)) +1
#need to add plus 1!
maxPossibleKMean <- possibleKsMean[length(possibleKsMean)]
#says 3:
maxPossibleKMean


#here, we see that even k=9 is quite high, but k=2 still wins out. perhaps too sensitive:
minConsensusByClust <- unlist(kmeansConsensus$minConsensusClusterByK[[which(names(dataMatrixList)=="study_12093_GPL96_all")]])

minConsensusByClust

#look at all the ones that now match this rounded MAX
#looks like we need to go 4 decimal places for this one.
#need to add plus 1!
possibleKsMin <- which(round(minConsensusByClust*10)==max(round(minConsensusByClust*10),na.rm=TRUE)) +1
#need to add plus 1!
maxPossibleKMin <- possibleKsMin[length(possibleKsMin)]
#here says 3:
maxPossibleKMin

######
#now try one with little heterogeneity: study_25065_GPL96_MDACC
#how many K did gap test say? 4:
length(gapKmeans_pam50Short[["study_25065_GPL96_MDACC"]])
#no pam50 info, but do have ER status:
#most of these are basal or luminal - looks like there is perhaps 3 main groups?
#but could be 2...
table(pData(esets[["study_25065_GPL96_MDACC"]])$pam50)


#hmm...a lot are close to zero:
PAC <- unlist(kmeansConsensus$PAC[[which(names(dataMatrixList)=="study_25065_GPL96_MDACC")]])

PAC
#perhaps round to hundreth decimal, take max?? this would suggest seven clusters....
#splitting hairs around the 100th decimal point level makes no sense!
#look at all the ones that now match this rounded MIN
possibleKs <- which(round(PAC*100)==min(round(PAC*100),na.rm=TRUE))

maxPossibleK <- possibleKs[length(possibleKs)]
#gives 2: (and it's pretty clear from PAC, although k=5 isn't bad, either)
#could be a bit by chance that got a "perfect" zero:
maxPossibleK

#here, we see that even k=9 is quite high!
meanConsensusByClust <- unlist(kmeansConsensus$meanConsensusClusterByK[[which(names(dataMatrixList)=="study_25065_GPL96_MDACC")]])

meanConsensusByClust

#look at all the ones that now match this rounded MAX
#looks like we need to go 2 decimal places for this one.
#need to add plus 1!
possibleKsMean <- which(round(meanConsensusByClust*10)==max(round(meanConsensusByClust*10),na.rm=TRUE)) +1
#need to add plus 1!
maxPossibleKMean <- possibleKsMean[length(possibleKsMean)]
#says 5:
maxPossibleKMean


#here, we see that even k=9 is quite high, but k=2 still wins out. perhaps too sensitive:
minConsensusByClust <- unlist(kmeansConsensus$minConsensusClusterByK[[which(names(dataMatrixList)=="study_25065_GPL96_MDACC")]])

minConsensusByClust

#look at all the ones that now match this rounded MAX
#looks like we need to go 4 decimal places for this one.
#need to add plus 1!
possibleKsMin <- which(round(minConsensusByClust*10)==max(round(minConsensusByClust*10),na.rm=TRUE)) +1
#need to add plus 1!
maxPossibleKMin <- possibleKsMin[length(possibleKsMin)]
#here says 2:
maxPossibleKMin

#########
#now try one with more heterogeneity: study_19615_GPL570_all
#how many K did gap test say? 4:
length(gapKmeans_pam50Short[["study_19615_GPL570_all"]])
#no pam50 info, but do have ER status:
#most of these are basal or luminal - looks like there IS actually 4 groups:
table(pData(esets[["study_19615_GPL570_all"]])$pam50)
table(pData(esets[["study_19615_GPL570_all"]])$ER_preTrt)
table(pData(esets[["study_19615_GPL570_all"]])$HER2_preTrt)

#hmm...a lot are close to zero:
PAC <- unlist(kmeansConsensus$PAC[[which(names(dataMatrixList)=="study_19615_GPL570_all")]])

PAC
#perhaps round to hundreth decimal, take max?? this would suggest seven clusters....
#splitting hairs around the 100th decimal point level makes no sense!
#look at all the ones that now match this rounded MIN
possibleKs <- which(round(PAC*100)==min(round(PAC*100),na.rm=TRUE))

maxPossibleK <- possibleKs[length(possibleKs)]
#gives 4
maxPossibleK

meanConsensusByClust <- unlist(kmeansConsensus$meanConsensusClusterByK[[which(names(dataMatrixList)=="study_19615_GPL570_all")]])

meanConsensusByClust

#look at all the ones that now match this rounded MAX
#looks like we need to go 2 decimal places for this one.
#need to add plus 1!
possibleKsMean <- which(round(meanConsensusByClust*10)==max(round(meanConsensusByClust*10),na.rm=TRUE)) +1
#need to add plus 1!
maxPossibleKMean <- possibleKsMean[length(possibleKsMean)]
#says 4:
maxPossibleKMean


#here, we see that even k=9 is quite high, but k=2 still wins out. perhaps too sensitive:
minConsensusByClust <- unlist(kmeansConsensus$minConsensusClusterByK[[which(names(dataMatrixList)=="study_19615_GPL570_all")]])

minConsensusByClust

#look at all the ones that now match this rounded MAX
#looks like we need to go 4 decimal places for this one.
#need to add plus 1!
possibleKsMin <- which(round(minConsensusByClust*10)==max(round(minConsensusByClust*10),na.rm=TRUE)) +1
#need to add plus 1!
maxPossibleKMin <- possibleKsMin[length(possibleKsMin)]
#here says 4:
maxPossibleKMinlength(gapKmeans_pam50Short[["study_20181_GPL96_all"]])

table(pData(esets[["study_20181_GPL96_all"]])$pam50)
#hmm: from publication, says are all ER + :http://www.ncbi.nlm.nih.gov/pubmed/20646288
#so expect lower # of clusters
table(pData(esets[["study_20181_GPL96_all"]])$ER_preTrt)
table(pData(esets[["study_20181_GPL96_all"]])$HER2_preTrt)

PAC <- unlist(kmeansConsensus$PAC[[which(names(dataMatrixList)=="study_20181_GPL96_all")]])

PAC
#perhaps round to hundreth decimal, take max?? this would suggest seven clusters....
#splitting hairs around the 100th decimal point level makes no sense!
#look at all the ones that now match this rounded MIN
possibleKs <- which(round(PAC*100)==min(round(PAC*100),na.rm=TRUE))

maxPossibleK <- possibleKs[length(possibleKs)]
#gives 2
maxPossibleK

meanConsensusByClust <- unlist(kmeansConsensus$meanConsensusClusterByK[[which(names(dataMatrixList)=="study_20181_GPL96_all")]])

meanConsensusByClust

#look at all the ones that now match this rounded MAX
#looks like we need to go 2 decimal places for this one.
#need to add plus 1!
possibleKsMean <- which(round(meanConsensusByClust*10)==max(round(meanConsensusByClust*10),na.rm=TRUE)) +1
#need to add plus 1!
maxPossibleKMean <- possibleKsMean[length(possibleKsMean)]
#says 3 - rounding bumped this one up.:
maxPossibleKMean


#here, we see that even k=9 is quite high, but k=2 still wins out. perhaps too sensitive:
minConsensusByClust <- unlist(kmeansConsensus$minConsensusClusterByK[[which(names(dataMatrixList)=="study_20181_GPL96_all")]])

minConsensusByClust

#look at all the ones that now match this rounded MAX
#looks like we need to go 4 decimal places for this one.
#need to add plus 1!
possibleKsMin <- which(round(minConsensusByClust*10)==max(round(minConsensusByClust*10),na.rm=TRUE)) +1
#need to add plus 1!
maxPossibleKMin <- possibleKsMin[length(possibleKsMin)]
#here says 2:
maxPossibleKMin

###conclusions: looks like if no PAC is below .1: let k=1.

