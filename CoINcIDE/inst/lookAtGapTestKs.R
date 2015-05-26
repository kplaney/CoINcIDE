load("/home/kplaney/breast_analysis_withTop20Genes/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.RData.gzip")
names(dataMatrixList)
studyNames <- c("study_2034_GPL96_all","study_25055_GPL96_MDACC_M","study_22226_GPL1708_all",
                "study_20181_GPL96_all","study_19615_GPL570_all" ,"study_16446_GPL570_all",
                "study_12093_GPL96_all","study_25065_GPL96_MDACC")

datasetIndices <- na.omit(match(studyNames,names(dataMatrixList)))
load("/home/kplaney/breast_analysis_withTop20Genes/gapTestKmeans_pam50Full_nstart12015-05-04.RData.gzip")

unlist(kmeansGapTest_pam50_full_Nstart1$bestK)[datasetIndices]

load("/home/kplaney/breast_analysis_withTop20Genes/gapTestKmeans_pam50Full_nstart252015-05-04.RData.gzip")
unlist(kmeansGapTest_pam50_full_Nstart25$bestK)[datasetIndices]

load("/home/kplaney/breast_analysis_withTop20Genes/gapTestKmeans_pam50Short_nstart12015-05-04.RData.gzip")

unlist(kmeansGapTest_pam50_short_Nstart1$bestK)[datasetIndices]

load("/home/kplaney/breast_analysis_withTop20Genes/gapTestKmeans_pam50Short_nstart252015-05-04.RData.gzip")

unlist(kmeansGapTest_pam50_short_Nstart25$bestK)[datasetIndices]


##now hclust

#(there's no actual "nstart" for hclust; it's all the same.)

unlist(kmeansGapTest_pam50_full_Nstart1$bestK)[datasetIndices]

load("/home/kplaney/breast_analysis_withTop20Genes/gapTesthclust_pam50Full_2015-05-15.RData.gzip")
unlist(hclustGapTest_pam50_full_Nstart25$bestK)[datasetIndices]

load("/home/kplaney/breast_analysis_withTop20Genes/hclustGapTest_pam50_short_2015-05-15.RData.gzip")

unlist(hclustGapTest_pam50_short$bestK)[datasetIndices]
