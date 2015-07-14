 test_ov <- read.delim("~/Desktop/tmp_ov.txt",header=TRUE)
 test[na.omit(match(test_ov[,1],test[,1])),1]
# [1] SERPINF1 FAP      MMP2     COL6A3   C1S      TGFBR2   ODC1     ALDH1A1  ERBB3   
#[10] TFPI     FBP1  
 

 test[na.omit(match(test_ov[,1],test[,1])),]
#SERPINF1, FAP, MMP2, COL6A3, C1S, TGFBR2: breast meta-cluster 3, ovarian meta-cluster 2
 #ODC1: breast meta-cluster 2, ovarian meta-cluster 3 (only 4 genes heavily overexpressed)
 
 #ALDH1A1: breast meta-cluster 3, ovarian meta-cluster 5
 
 #ERBB3, FBP1: breast meta-cluster 1, ovarian meta-cluster 5
 
 #TFPI: breast meta-cluster 3, ovarian meta-cluster 5