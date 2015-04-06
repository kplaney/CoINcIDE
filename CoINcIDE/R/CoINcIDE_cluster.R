#use default Ward's method
library("cluster")
library("ConsensusClusterPlus")

clustF <- function(dataset,k){
  
  kmeans(dataset, centers=k,iter.max=30,nstart=25,algorithm="Hartigan-Wong");
  
}
library("cluster");
gapTest <- clusGap(dataset, FUNcluster=clustF,K.max=20, B = 1000, verbose = interactive());

#f(k) is the gap statistic.
# method = c("firstSEmax", "Tibs2001SEmax", "globalSEmax",
"firstmax", "globalmax")
bestK <- maxSE(f=gapTest$Tab[,"gap"], SE.f=gapTest$Tab[,"SE.sim"],
      method = c("Tibs2001SEmax"),
      SE.factor = 1)

select_metric <- gapTest$Tab[c(1:(K.max-1)),"gap"] > gapTest$Tab[c(2:K.max),"gap"] - gapTest$Tab[c(2:K.max),"SE.sim"];
#take first one that is true.
best_k <- which(select_metric)[1];
cat("\nYour chosen k is: ",output$kmeans_output[[1]]$bestK, "\ncluster library says: ",bestK, " ie ", best_k ,"\n");


#biobase:
consensusClusterPlus