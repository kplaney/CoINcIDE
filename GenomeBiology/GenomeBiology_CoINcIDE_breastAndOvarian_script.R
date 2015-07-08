
saveDirVector <- c("/home/ywrfc09/breast_analysis/PAM50_analyses/","/home/ywrfc09/breast_analysis/PAM50_analyses/","/home/ywrfc09/breast_analysis/PAM50_analyses/",
             "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/","/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/","/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/",
             "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/","/home/ywrfc09/breast_analysis/metaRankNoTop20Genes/","/home/ywrfc09/breast_analysis/metaRankNoTop20Genes/",
             "/home/ywrfc09/breast_analysis/metaRankNoTop20Genes/","/home/ywrfc09/breast_analysis/metaRankNoTop20Genes/","/home/ywrfc09/breast_analysis/metaRankNoTop20Genes/",
             "/home/ywrfc09/ovarian_analysis/metaRankWithTop20Genes/","/home/ywrfc09/ovarian_analysis/metaRankWithTop20Genes/","/home/ywrfc09/ovarian_analysis/metaRankWithTop20Genes/",
             "/home/ywrfc09/ovarian_analysis/metaRankWithTop20Genes/")
             
#first object: save _PACR list indices just so can automatically grab the indices.
clusterDataVector <- c("/home/ywrfc09/breast_analysis/PAM50_analyses/pam50Full_centroidCluster.rds",
                       "/home/ywrfc09/breast_analysis/PAM50_analyses/kmeansConsensuspam50_short_Nstart1pItem9_pam50ShortFeatures_2015-04-28.rds",
                       "/home/ywrfc09/breast_analysis/PAM50_analyses//kmeansConsensuspam50_full_Nstart1pItem9_pam50FullFeatures_2015-04-28.rds",
                       "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/curatedbreastData_kmeansConsensus_nstart1pItem9200Features_2015-05-04.rds",
                       "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/curatedbreastData_kmeansConsensus_nstart1pItem9500Features_2015-05-04.rds",
                       "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/curatedbreastData_kmeansConsensus_nstart1pItem91000Features_2015-05-05.rds",
                       "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes/curatedbreastData_kmeansConsensus_nstart1pItem92000Features_2015-05-05.rds",
                       "/home/ywrfc09/breast_analysis/metaRankNoTop20Genes//curatedbreastData_kmeansConsensus_nstart1pItem9200Features_2015-05-19.rds",
                       "/home/ywrfc09/breast_analysis/metaRankNoTop20Genes//curatedbreastData_kmeansConsensus_nstart1pItem9300Features_2015-05-20.rds",
                       "/home/ywrfc09/breast_analysis/metaRankNoTop20Genes//curatedbreastData_kmeansConsensus_nstart1pItem9500Features_2015-05-19.rds",
                       "/home/ywrfc09/breast_analysis/metaRankNoTop20Genes//curatedbreastData_kmeansConsensus_nstart1pItem91000Features_2015-05-20.rds",
                       "/home/ywrfc09/breast_analysis/metaRankNoTop20Genes//curatedbreastData_kmeansConsensus_nstart1pItem92000Features_2015-05-20.rds",
                       "/home/ywrfc09/ovarian_analysis/metaRankWithTop20Genes/curatedOvarianData_kmeansConsensus_nstart1_200Features_2015-04-28.rds",
                       "/home/ywrfc09/ovarian_analysis/metaRankWithTop20Genes/curatedOvarianData_kmeansConsensus_nstart1_500Features_2015-04-29.rds",
                       "/home/ywrfc09/ovarian_analysis/metaRankWithTop20Genes/curatedOvarianData_kmeansConsensus_nstart1_1000Features_2015-04-29.rds",
                       "/home/ywrfc09/ovarian_analysis/metaRankWithTop20Genes/curatedOvarianData_kmeansConsensus_nstart1_2000Features_2015-04-30.RData.gzip")
                       
                       

dataMatrixListVector <- c("/home/ywrfc09/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.rds",
                          "/home/ywrfc09/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.rds",
                          "/home/ywrfc09/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.rds",
                          "/home/ywrfc09/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.rds",
                          "/home/ywrfc09/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.rds",
                          "/home/ywrfc09/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.rds",
                          "/home/ywrfc09/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.rds",
                          "/home/ywrfc09/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.rds",
                          "/home/ywrfc09/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.rds",
                          "/home/ywrfc09/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.rds",
                          "/home/ywrfc09/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.rds",
                          "/home/ywrfc09/breast_analysis/curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.rds",
                          "/home/ywrfc09/ovarian_analysis/curatedBreastData_dataMatrixList_proc_minVar001.rds",
                          "/home/ywrfc09/ovarian_analysis/curatedBreastData_dataMatrixList_proc_minVar001.rds",
                          "/home/ywrfc09/ovarian_analysis/curatedBreastData_dataMatrixList_proc_minVar001.rds",
                          "/home/ywrfc09/ovarian_analysis/curatedBreastData_dataMatrixList_proc_minVar001.rds")

experimentNameVector <- c("PAM50centroidCluster","PAM50kmeansShort","PAM50kmeansFull","breast278F","breast563F","breast1058F","breast2052F","breast200F","breast300F","breast500F","breast1000F",
                          "breast2000F","ovarian240F","ovarian527F","ovarian1019F","ovarian2014F")
centroidMethodVector <- c("mean","median")
edgeMethodVector <- c("pearson","spearman")

#run all pearson, and then all spearman, experiments
scriptDir <- "/home/ywrfc09/CoINcIDE/coincide/GenomeBiology/"
for(e in 1:length(edgeMethodVector)){
  
  for(c in 1:length(centroidMethodVector)){
    
    for(d in 1:length(clusterDataVector)){
      
      setwd(scriptDir)
      #kick off script!
      #you could also just paste this into the terminal if you're in the scriptDir already.
      system(paste0("Rscript GenomeBiology_baseAnalysisScript.R ",edgeMethodVector[e], " ",centroidMethodVector[c], " ", clusterDataVector[d], " ", 
                    dataMatrixListVector[d]," ", experimentNameVector[d], " ",saveDirVector[d]))
      
    }
    
  }
  
  
}