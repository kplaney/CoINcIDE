
#where do the baseAnalysis and baseAnalysisNullFDR scripts lie?
scriptDir <- "/home/ywrfc09/CoINcIDE/coincide/GenomeBiology/"

library("Coincide")

##CHANGE these paths to your user directory
##these variables will be used in scripts further down the pipeline
saveDirGlobal <- "/home/ywrfc09/breast_analysis/"
saveDirGlobal_ovarian <- "/home/ywrfc09/ovarian_analysis/"
saveDir_PAM50 <- "/home/ywrfc09/breast_analysis/PAM50_analyses/"
saveDir_20 <- "/home/ywrfc09/breast_analysis/metaRankWithTop20Genes"
saveDir_20_ovarian <-  "/home/ywrfc09/ovarian_analysis/metaRankNoTop20Genes"
outputFile <- "/home/kplaney/breast_analysis/clust_test_outMessages.txt"

saveDirVector <- c(saveDir_PAM50,saveDir_PAM50,saveDir_PAM50,
             saveDir_20,saveDir_20,saveDir_20,
             saveDir_20,saveDir_20_ovarian, saveDir_20_ovarian, 
             saveDir_20_ovarian,saveDir_20_ovarian)

#NOTE:  these are time-stamped: you will need to change the names of these file date tags
#to match the ones you outputted.
clusterDataVector <- c(paste0(saveDirVector, c("pam50Full_centroidCluster.rds",
                                                "kmeansConsensuspam50_short_Nstart1pItem9_pam50ShortFeatures_2015-04-28.rds",
                                               "kmeansConsensuspam50_full_Nstart1pItem9_pam50FullFeatures_2015-04-28.rds",
                                               "curatedbreastData_kmeansConsensus_nstart1pItem9200Features_2015-05-04.rds",
                                               "curatedbreastData_kmeansConsensus_nstart1pItem9500Features_2015-05-04.rds",
                                               "curatedbreastData_kmeansConsensus_nstart1pItem91000Features_2015-05-05.rds",
                                               "curatedbreastData_kmeansConsensus_nstart1pItem92000Features_2015-05-05.rds",
                                               "curatedOvarianData_kmeansConsensus_nstart1_200Features_2015-04-28.rds",
                                               "curatedOvarianData_kmeansConsensus_nstart1_500Features_2015-04-29.rds",
                                               "curatedOvarianData_kmeansConsensus_nstart1_1000Features_2015-04-29.rds",
                                               "curatedOvarianData_kmeansConsensus_nstart1_2000Features_2015-04-30.rds")))

       
             
dataMatrixListVector <- c(rep.int(paste0(saveDirGlobal,"curatedBreastData_dataMatrixList_proc_minVar001_min10kGenes_min40Samples.rds"),7),
                          rep.int(paste0(saveDirGlobal_ovarian,"curatedOvarianData_dataMatrixList_proc_minVar001.rds"),4))


experimentNameVector <- c("PAM50centroidCluster","PAM50kmeansShort","PAM50kmeansFull","breast278F","breast563F","breast1058F","breast2052F",
                          "ovarian240F","ovarian527F","ovarian1019F","ovarian2014F")
centroidMethodVector <- c("mean","median")
edgeMethodVector <- c("pearson")



for(e in 1:length(edgeMethodVector)){
  
  for(c in 1:length(centroidMethodVector)){
    
    for(d in 1:length(clusterDataVector)){
      
      setwd(scriptDir)
      #kick off script!
      #you could also just paste this into the terminal if you're in the scriptDir already.
      system(paste0("Rscript Coincide_baseAnalysisScript.R ",edgeMethodVector[e], " ",centroidMethodVector[c], " ", clusterDataVector[d], " ", 
                    dataMatrixListVector[d]," ", experimentNameVector[d], " ",saveDirVector[d]))
      
    }
    
  }
  
  
}

#separate loop for null FDR (as this takes longer, may want to run separately.)
for(e in 1:length(edgeMethodVector)){
  
  for(c in 1:length(centroidMethodVector)){
    
    for(d in 1:length(clusterDataVector)){
      
      setwd(scriptDir)
      #kick off script!
      #you could also just paste this into the terminal if you're in the scriptDir already.
      system(paste0("Rscript Coincide_baseAnalysisNullFDRScript.R ",edgeMethodVector[e], " ",centroidMethodVector[c], " ", clusterDataVector[d], " ", 
                    dataMatrixListVector[d]," ", experimentNameVector[d], " ",saveDirVector[d]))
      
    }
    
  }
  
  
}