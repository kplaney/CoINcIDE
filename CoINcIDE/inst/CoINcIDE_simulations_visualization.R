

###first: representative heatmap pictures for each of the 3 cases



##FIRST, SOME VALIDATION OF YOUR SIMULATED DATA
source("/home/kplaney/gitRepos//IGP_network/igp_network/clust_robust.R");
saveDir <- "/home/kplaney/ISMB/lung_sims/perfect/";

#,simType=c("highQualityClust","mixedClustQualityClust","unevenSizeClust")
lungData <- createLungMatrixList();
#this data already appears to be gene-centered:
head(rowMeans(lungData$lung200));

#good: samples are highly correlated to tissue type. These are our "clusters"
realLungCorData <- calcCorMeanMatrix(lungData$clustMatrixList,fileSaveName="/home/kplaney/ISMB/lung_sims/meanCorMatrix_realLungData.txt");

###simulate this covariance/correlation structure. zero noise at first.
numPerClust <- c(50,50,50,50);
eigenValueMin <- -.001;
lungSimData_perc <- createLungSimDatasets(sourceDir="/home/kplaney/gitRepos//IGP_network/igp_network/",
                                     saveDir = "/home/kplaney/ISMB/lung_sims/perfect/",numSimDatasets=11,
                                     eigenValueMin = -.001,simType=c("highQualityClust"),
                                     numPerClust = c(50,50,50,50),
                                     stddevNoise=0,numRows=200);

#inspect- looks similar to other results.
simLungCorData <- calcCorMeanMatrix(lungSimData$simClustList[[1]], fileSaveName="/home/kplaney/ISMB/lung_sims/meanCorMatrix_simLungData_200_evenSamples.txt");

#prep for ggplots.
library("reshape2");
lungOrigDF <- melt(lungData$lung200);
lungOrigDF <-  cbind(lungOrigDF,rep.int("orig",times=nrow(lungOrigDF)));
colnames(lungOrigDF )[4] <- "dataSource";
lungSimDF <- melt(lungSimData$simMatrixList[[1]]);
library("limma");
#make the tissue categories uniform so can facet.
lungSimDF$Var2 <- strsplit2(lungSimDF$Var2,"_")[,1];
lungSimDF <- cbind(lungSimDF,rep.int("sim",times=nrow(lungSimDF)));
colnames(lungSimDF )[4] <- "dataSource";

#will get a warning but we don't care - don't need the gene names.
masterDF <- rbind(lungOrigDF,lungSimDF);
masterDF$dataSource <- factor(masterDF$dataSource);

#remove gray background + theme(panel.background = element_blank())
library("ggplot2")
simData_density <- ggplot(data = masterDF, aes(x=value,linetype=dataSource))+
  #alpha in geom_density. controls transparency of color, if used fill= above.
  #if want transparent fill: geom_density(alpha=0). but different lines
  #are better for a black and white publication.
  geom_density() +
  labs(title = "Density curves for \noriginal and simulated data",
       y="Density",x="Logged gene expression")+
  theme(panel.background = element_rect(fill='white', colour='black')) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

png(filename=paste0(saveDir,"/simVsRealLungData_densityPlots_",Sys.Date(),".png"),
    width = 700, height = 1000)
origNetworkPlot <- plot(simData_density);
dev.off();

simData_density_facettedByTissue <- ggplot(data = masterDF, aes(x=value,linetype=dataSource))+
  #alpha in geom_density. controls transparency of color, if used fill= above.
  #if want transparent fill: geom_density(alpha=0). but different lines
  #are better for a black and white publication.
  geom_density() +
  #if want in a 2x2 grid: facet_wrap(~Var2,ncol=2)+
  facet_grid(Var2~.)+
  labs(title = "Density curves for \noriginal and simulated data",
       y="Density",x="Logged gene expression")+
  theme(panel.background = element_rect(fill='white', colour='black')) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

png(filename=paste0(saveDir,"/simVsRealLungData_densityPlots_byTissueType_",Sys.Date(),".png"),
    width = 700, height = 1000)
origNetworkPlot <- plot(simData_density_facettedByTissue);
dev.off();


########heatmaps
source("/home/kplaney/gitRepos/IGP_network/igp_network/clust_robust.R");
library("gplots");
library("RColorBrewer");
colorCodeF <-  colorRampPalette(c("black","blue","red"), bias = 4,
                                space = "rgb", interpolate = "linear");
colorCodes <- colorCodeF(4);


####perfect simes
###NOTE: these first parts use old code chunks, because the network plots are the same.
saveDir <- "/home/kplaney/ISMB/lung_sims/perfect/";

simData <- createLungSimDatasets(sourceDir="/home/kplaney/gitRepos//IGP_network/igp_network/",
                                   saveDir = saveDir,numSimDatasets=10,
                                   eigenValueMin = -.001,simType=c("highQualityClust"),
                                   numPerClust = c(50,50,50,50),
                                   stddevNoise=0,numRows=200);

simPerc_corr <- findSimilarClusters(clustMatrixList=simData$clustMatrixList,dataMatrixList=simData$dataMatrixList,minClustSampleSize=10,resampleWithClustSamples=TRUE,compareMetric=c("correlation"),
                                    normFeatureWise=c("none"),method=c("diff_distribution"),numSims=500,
                                    minFractFeatureIntersect=.05,minNumFeatureIntersect=5,fract_nullMax=.1,
                                    centroidType=c("mean"),outputFile="/home/kplaney/ISMB/lung_sims/tmp.txt",
                                    nullClustTag ="perf_temp" ,sourceDir="/home/kplaney/gitRepos/IGP_network/igp_network/",
                                    trueDiffMin=.3,intraSimilMax=.3,simil_pvalueMax=1,tempSaveDir=saveDir,
                                    runIntraStabilityTests=FALSE,compareIntraInterDistributions=FALSE
);
save(simPerc_corr,file="/home/kplaney/ISMB/lung_sims/perfect/simPerc_corr_clustRobustResults.RData.gzip");
load("/home/kplaney/ISMB/lung_sims/perfect/simPerc_corr_clustRobustResults.RData.gzip");
output_corr <- analyzeClustRobust_output(clustRobust_output=simPerc_corr,fullMatrixList=simPerc_corr$dataMatrixList
                                         ,clustMatrixList=simPerc_corr$clustMatrixList,
                                         clustMethodName="simulation 1 correlation similarity",saveDir=saveDir,method=c("diff_distribution"),compareMetric=c("correlation"),
                                      edge_trueCorr_thresh=.3,edge_fractFeatIntersect_thresh=.05,edge_numFeatIntersect_thresh=20,edge_IGP_thresh=.7,
                                      edge_simil_overNull_pvalue_thresh=.05,genomeSize=20000,sourceDir="/home/kplaney/gitRepos/IGP_network/igp_network/",
                                      commMethod="edgeBetween",minNumUniqueStudiesPerCommunity=3,
                                      restrictEdgesByStudy=TRUE,runGSEA=FALSE,
                                      heatmaps=FALSE
);

output <- list(simPerc_corr=simPerc_corr,output_corr=output_corr);
save(output,file="/home/kplaney/ISMB/lung_sims/perfect/simPerc_corr_clustRobustResults.RData.gzip",compress="gzip");

simPlots <- advanced_networkPlots(analysisOutput=output_corr,
                                  colorCodes=colorCodes,
                                  saveDir=saveDir,saveName="Simulation1_simCorr_noise0");

simNum <- 1;
simData <- createLungSimDatasets(sourceDir="/home/kplaney/gitRepos//IGP_network/igp_network/",
                                 saveDir = saveDir,numSimDatasets=1,
                                 eigenValueMin = -.001,simType=c("highQualityClust"),
                                 numPerClust = c(50,50,50,50),
                                 stddevNoise=0,numRows=200);
simSampleColors <- rownames(simData$dataMatrixList[[1]]);
simSampleColors[grep("Carcinoid",simSampleColors)] <- colorCodes[1];
simSampleColors[grep("Colon",simSampleColors)] <- colorCodes[2];
simSampleColors[grep("Normal",simSampleColors)] <- colorCodes[3];
simSampleColors[grep("SmallCell",simSampleColors)] <- colorCodes[4];


##do just heatmap, then corr data

#NOTE: if get error: Error in .External.graphics(C_layout, num.rows, num.cols, mat, as.integer(num.figures),  : 
#invalid graphics state: just save plot directly instead of trying to plot it in RStudio. Sometimes just a margin issue.
png(filename=paste0(saveDir,paste0("/sim",simNum,"_expr_heatmap_",Sys.Date(),".png")),
    width = 700, height = 1000)

heatmap.2(t(data.matrix(simData$dataMatrixList[[1]])),Rowv=F,Colv=F,scale = c("none"), 
          density.info="none", trace="none",symkey=TRUE,key.title="Color Key",key.xlab="Logged \nexpression",
          col=rev(brewer.pal(11, "RdBu")),key=TRUE,labRow=F,labCol=F,ColSideColors=simSampleColors,
          main=paste0("\n\nSimulation ",simNum," :\nexpression heatmap \nby tissue type"));

dev.off();

png(filename=paste0(saveDir,paste0("sim",simNum,"_corr_heatmap_",Sys.Date(),".png")),
    width = 700, height = 1000)

#if want no density plot: density.info="none",
heatmap.2(cor(t(simData$dataMatrixList[[1]])),Rowv=F,Colv=F,scale = c("none"), 
          density.info="none",trace="none",symkey=TRUE,key.title="Color Key",key.xlab="correlation",
          col=rev(brewer.pal(11, "RdBu")),key=TRUE,labRow=F,labCol=F,ColSideColors=simSampleColors,
          RowSideColors=simSampleColors,main=paste0("\n\nSimulation ",simNum," :\ncorrelation heatmap \nby tissue type"));

dev.off();

######mixed


saveDir <- "/home/kplaney/ISMB/lung_sims/mixed//";

simData <- createLungSimDatasets(sourceDir="/home/kplaney/gitRepos//IGP_network/igp_network/",
                                  saveDir = "/home/kplaney/ISMB/lung_sims/mixed/",numSimDatasets=10,
                                  eigenValueMin = -.001,simType=c("mixedClustQualityClust"),
                                  numPerClust = c(50,50,50,50),
                                  stddevNoise=0,numRows=200);


#only want red and black now.
colorCodes_mixed <- colorCodes[c(1,4)];

simMixed_corr <- findSimilarClusters(clustMatrixList=simData$clustMatrixList,dataMatrixList=simData$dataMatrixList,minClustSampleSize=10,resampleWithClustSamples=TRUE,compareMetric=c("correlation"),
                                    normFeatureWise=c("none"),method=c("diff_distribution"),numSims=500,
                                    minFractFeatureIntersect=.05,minNumFeatureIntersect=5,fract_nullMax=.1,
                                    centroidType=c("mean"),outputFile="/home/kplaney/ISMB/lung_sims/tmp.txt",
                                    nullClustTag ="perf_temp" ,sourceDir="/home/kplaney/gitRepos/IGP_network/igp_network/",
                                    trueDiffMin=.3,intraSimilMax=.3,simil_pvalueMax=1,tempSaveDir=saveDir,
                                    runIntraStabilityTests=FALSE,compareIntraInterDistributions=FALSE
);

save(simMixed_corr,file="/home/kplaney/ISMB/lung_sims/mixed/simMixed_corr_clustRobustResults.RData.gzip");
load("/home/kplaney/ISMB/lung_sims/mixed/simMixed_corr_clustRobustResults.RData.gzip");

output_corr <- analyzeClustRobust_output(clustRobust_output=simMixed_corr,fullMatrixList=simMixed_corr$dataMatrixList,
                                         clustMatrixList=simMixed_corr$clustMatrixList,
                                         clustMethodName="simulation 2 correlation similarity",saveDir=saveDir,method=c("diff_distribution"),compareMetric=c("correlation"),
                                         edge_trueCorr_thresh=.3,edge_fractFeatIntersect_thresh=.05,edge_numFeatIntersect_thresh=20,edge_IGP_thresh=.7,
                                         edge_simil_overNull_pvalue_thresh=.05,genomeSize=20000,sourceDir="/home/kplaney/gitRepos/IGP_network/igp_network/",
                                         commMethod="edgeBetween",minNumUniqueStudiesPerCommunity=3,
                                         restrictEdgesByStudy=TRUE,runGSEA=FALSE,
                                         heatmaps=FALSE
);

output <- list(simMixed_corr=simMixed_corr,output_corr=output_corr);
save(output,file="/home/kplaney/ISMB/lung_sims/mixed/simMixed_corr_clustRobustResults.RData.gzip",compress="gzip");

simPlots <- advanced_networkPlots(analysisOutput=output_corr,
                                  colorCodes=colorCodes_mixed,
                                  saveDir=saveDir,saveName="Simulation2_simCorr_noise0");

simNum <- 2;
simData <- createLungSimDatasets(sourceDir="/home/kplaney/gitRepos//IGP_network/igp_network/",
                                 saveDir = saveDir,numSimDatasets=1,
                                 eigenValueMin = -.001,simType=c("mixedClustQualityClust"),
                                 numPerClust = c(50,50,50,50),
                                 stddevNoise=0,numRows=200);

simSampleColors <- rownames(simData$dataMatrixList[[1]]);
simSampleColors[grep("Carcinoid",simSampleColors)] <- colorCodes[1];
#gray out the "bad" clusters.
simSampleColors[grep("Colon",simSampleColors)] <- "#D3D3D3"
simSampleColors[grep("Normal",simSampleColors)] <- "#7e7e7e"
simSampleColors[grep("SmallCell",simSampleColors)] <- colorCodes[4];


##do just heatmap, then corr data

#NOTE: if get error: Error in .External.graphics(C_layout, num.rows, num.cols, mat, as.integer(num.figures),  : 
#invalid graphics state: just save plot directly instead of trying to plot it in RStudio. Sometimes just a margin issue.
png(filename=paste0(saveDir,paste0("/sim",simNum,"_expr_heatmap_",Sys.Date(),".png")),
    width = 700, height = 1000)

heatmap.2(t(data.matrix(simData$dataMatrixList[[1]])),Rowv=F,Colv=F,scale = c("none"), 
          density.info="none", trace="none",symkey=TRUE,key.title="Color Key",key.xlab="Logged \nexpression",
          col=rev(brewer.pal(11, "RdBu")),key=TRUE,labRow=F,labCol=F,ColSideColors=simSampleColors,
          main=paste0("\n\nSimulation ",simNum," :\nexpression heatmap \nby tissue type"));

dev.off();

png(filename=paste0(saveDir,paste0("sim",simNum,"_corr_heatmap_",Sys.Date(),".png")),
    width = 700, height = 1000)

#if want no density plot: density.info="none",
heatmap.2(cor(t(simData$dataMatrixList[[1]])),Rowv=F,Colv=F,scale = c("none"), 
          density.info="none",trace="none",symkey=TRUE,key.title="Color Key",key.xlab="correlation",
          col=rev(brewer.pal(11, "RdBu")),key=TRUE,labRow=F,labCol=F,ColSideColors=simSampleColors,
          RowSideColors=simSampleColors,main=paste0("\n\nSimulation ",simNum," :\ncorrelation heatmap \nby tissue type"));

dev.off();

######uneven
saveDir <- "/home/kplaney/ISMB/lung_sims/uneven//";

simData <- createLungSimDatasets(sourceDir="/home/kplaney/gitRepos//IGP_network/igp_network/",
                                  saveDir = "/home/kplaney/ISMB/lung_sims/uneven/",numSimDatasets=10,
                                  eigenValueMin = -.001,simType=c("unevenSizeClust"),
                                  numPerClust = c(50,50,50,50),
                                  stddevNoise=0,numRows=200);


simUneven_corr <- findSimilarClusters(clustMatrixList=simData$clustMatrixList,dataMatrixList=simData$dataMatrixList,minClustSampleSize=10,resampleWithClustSamples=TRUE,compareMetric=c("correlation"),
                                     normFeatureWise=c("none"),method=c("diff_distribution"),numSims=500,
                                     minFractFeatureIntersect=.05,minNumFeatureIntersect=5,fract_nullMax=.1,
                                     centroidType=c("mean"),outputFile="/home/kplaney/ISMB/lung_sims/tmp.txt",
                                     nullClustTag ="uneven_temp" ,sourceDir="/home/kplaney/gitRepos/IGP_network/igp_network/",
                                     trueDiffMin=.3,intraSimilMax=.3,simil_pvalueMax=1,tempSaveDir=saveDir,
                                     runIntraStabilityTests=FALSE,compareIntraInterDistributions=FALSE
);

save(simUneven_corr,file="/home/kplaney/ISMB/lung_sims/uneven/simUneven_corr_clustRobustResults.RData.gzip",compress="gzip");
load("/home/kplaney/ISMB/lung_sims/uneven/simUneven_corr_clustRobustResults.RData.gzip");

colorCodes_uneven <- colorCodes[c(1,2)];

output_corr <- analyzeClustRobust_output(clustRobust_output=simUneven_corr,fullMatrixList=simUneven_corrr$dataMatrixList,
                                         clustMatrixList=simUneven_corr$clustMatrixList,
                                         clustMethodName="simulation 3 correlation similarity",saveDir=saveDir,method=c("diff_distribution"),compareMetric=c("correlation"),
                                         edge_trueCorr_thresh=.3,edge_fractFeatIntersect_thresh=.05,edge_numFeatIntersect_thresh=20,edge_IGP_thresh=.7,
                                         edge_simil_overNull_pvalue_thresh=.05,genomeSize=20000,sourceDir="/home/kplaney/gitRepos/IGP_network/igp_network/",
                                         commMethod="edgeBetween",minNumUniqueStudiesPerCommunity=3,
                                         restrictEdgesByStudy=TRUE,runGSEA=FALSE,
                                         heatmaps=FALSE
);

output <- list(simUneven_corr=simUneven_corr,output_corr=output_corr);
save(output,file="/home/kplaney/ISMB/lung_sims/uneven/simUneven_corr_clustRobustResults.RData.gzip",compress="gzip");

#interesting..creates 2 "blank" clusters...why? I guess because it missed these clusters?
simPlots <- advanced_networkPlots(analysisOutput=output_corr,
                                  colorCodes=colorCodes,
                                  saveDir=saveDir,saveName="Simulation3_simCorr_noise0");

simNum <- 3;
simData <- createLungSimDatasets(sourceDir="/home/kplaney/gitRepos//IGP_network/igp_network/",
                                 saveDir = saveDir,numSimDatasets=1,
                                 eigenValueMin = -.001,simType=c("unevenSizeClust"),
                                 numPerClust = c(50,50,50,50),
                                 stddevNoise=0,numRows=200);

simSampleColors <- rownames(simData$dataMatrixList[[1]]);
simSampleColors[grep("Carcinoid",simSampleColors)] <- colorCodes[1];

simSampleColors[grep("Colon",simSampleColors)] <- colorCodes[2];
#these clusters ARE in this dataset for zero noise:
simSampleColors[grep("Normal",simSampleColors)] <- colorCodes[3];
simSampleColors[grep("SmallCell",simSampleColors)] <- colorCodes[4];


##do just heatmap, then corr data

#NOTE: if get error: Error in .External.graphics(C_layout, num.rows, num.cols, mat, as.integer(num.figures),  : 
#invalid graphics state: just save plot directly instead of trying to plot it in RStudio. Sometimes just a margin issue.
png(filename=paste0(saveDir,paste0("/sim",simNum,"_expr_heatmap_",Sys.Date(),".png")),
    width = 700, height = 1000)

heatmap.2(t(data.matrix(simData$dataMatrixList[[1]])),Rowv=F,Colv=F,scale = c("none"), 
          density.info="none", trace="none",symkey=TRUE,key.title="Color Key",key.xlab="Logged \nexpression",
          col=rev(brewer.pal(11, "RdBu")),key=TRUE,labRow=F,labCol=F,ColSideColors=simSampleColors,
          main=paste0("\n\nSimulation ",simNum," :\nexpression heatmap \nby tissue type"));

dev.off();

png(filename=paste0(saveDir,paste0("sim",simNum,"_corr_heatmap_",Sys.Date(),".png")),
    width = 700, height = 1000)

#if want no density plot: density.info="none",
heatmap.2(cor(t(simData$dataMatrixList[[1]])),Rowv=F,Colv=F,scale = c("none"), 
          density.info="none",trace="none",symkey=TRUE,key.title="Color Key",key.xlab="correlation",
          col=rev(brewer.pal(11, "RdBu")),key=TRUE,labRow=F,labCol=F,ColSideColors=simSampleColors,
          RowSideColors=simSampleColors,main=paste0("\n\nSimulation ",simNum," :\ncorrelation heatmap \nby tissue type"));

dev.off();




#vs original data
lungData <- createLungMatrixList();
realSampleColors <- colnames(lungData$lung200);
realSampleColors[grep("Carcinoid",realSampleColors)] <- colorCodes[1];
realSampleColors[grep("Colon",realSampleColors)] <- colorCodes[2];
realSampleColors[grep("Normal",realSampleColors)] <- colorCodes[3];
realSampleColors[grep("SmallCell",realSampleColors)] <- colorCodes[4];

png(filename=paste0(saveDir,"/realLungData_expr_heatmap_",Sys.Date(),".png"),
    width = 700, height = 1000)

heatmap.2(data.matrix(lungData$lung200),Rowv=F,Colv=F,scale = c("none"), 
          density.info="none", trace="none",symkey=TRUE,key.title="Color Key",key.xlab="Logged expression",
          col=rev(brewer.pal(11, "RdBu")),key=TRUE,labRow=F,labCol=F,ColSideColors=realSampleColors,
          main="\n\n\n\n\n\n\n\nExpression heatmap \nof real lung data by tissue type");

dev.off();

png(filename=paste0(saveDir,"/realLungData_cor_heatmap_",Sys.Date(),".png"),
    width = 700, height = 1000)

heatmap.2(cor(lungData$lung200),Rowv=F,Colv=F,scale = c("none"), 
          density.info="none", trace="none",symkey=TRUE,key.title="Color Key",key.xlab="correlation",
          col=rev(brewer.pal(11, "RdBu")),key=TRUE,labRow=F,labCol=F,ColSideColors=realSampleColors,
          RowSideColors=realSampleColors,main="\n\n\n\n\n\n\n\nCorrelation heatmap \nof original lung data by tissue type");

dev.off();
###############TPR plots. 
saveDir <- "/home/kplaney/lungSims/"

####these were saved images() i.e. I saved the whole caboose to be on the safe side...
data1 <- "/home/kplaney/lungSims/highQuality_2015-05-12.RData.gzip";
data3 <- "/home/kplaney/lungSims/mixSims_2015-05-08.RData.gzip";
data5 <- "/home/kplaney/lungSims/unevenSims_2015-05-11.RData.gzip";

load(data1);
data1 <- highQuality
load(data3);
data3 <- mixQuality
load(data5);
data5 <- unevenQuality


#####REST DOES NOT CHANGE WITH DATASET NAMES.

#now create a ggplot data frame.
masterDF <- data.frame(data1$ROC_matrixFull,rep("highQuality",nrow(data1$ROC_matrixFull)));
colnames(masterDF)[ncol(masterDF)] <- "sim_type";

tmp <- data.frame(data3$ROC_matrixFull,rep("mixedQuality",nrow(data3$ROC_matrixFull)));
colnames(tmp)[ncol(tmp)] <- "sim_type";
masterDF <- rbind(masterDF,tmp);

tmp <- data.frame(data5$ROC_matrixFull,rep("unevenSize",nrow(data5$ROC_matrixFull)));
colnames(tmp)[ncol(tmp)] <- "sim_type";
masterDF <- rbind(masterDF,tmp);


library("ggplot2")
#masterDF$sim_type <- as.factor(masterDF$sim_type);
#TPR_plot <- ggplot(data =masterDF,aes(x=sd_noise,y=TPR,linetype=sim_type))+geom_line()+  labs(title = "TPR for simulations with increasing noise.",
#                                                                      y="TPR",x="sd noise")+
#  theme(panel.background = element_rect(fill='white', colour='black')) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                                                                               plot.title = element_text(size = rel(2)));

library("RColorBrewer")
colorCodes <- brewer.pal(3,"Dark2");

names(colorCodes) <- c("highQuality","mixedQuality","unevenSize");
TPR_plot_color <- ggplot(data =masterDF,aes(x=sd_noise,y=TPR,group=sim_type, colour=sim_type))  + geom_line(size=.7)+ 
  labs(title = "Network edge TPR for simulations with increasing noise",y="TPR",x="sd of additive random normal noise")+
  scale_color_manual(values=colorCodes)+
  theme(panel.background = element_rect(fill='white', colour='black')) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(colour = "black",size=18),axis.title.x = element_text(colour = "black",size=20,vjust=0),
        axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1),
        plot.title=element_text(colour="black",size=22,vjust=1,hjust=.4))+
  theme(legend.title=element_text(colour="black",size=12),legend.text=element_text(colour="black",size=15))+
  labs(colour="simulation\ntype");
                                                                 


TPR_plot_point <- ggplot(data =masterDF,aes(x=sd_noise,y=TPR,group=sim_type))  + geom_line(size=.4)+ geom_point(aes(shape=factor(sim_type)),size=5)+
  labs(title = "",y="True positive rate",x="Noise level")+
  scale_color_manual(values=colorCodes)+
  theme(panel.background = element_rect(fill='white', colour='black')) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(colour = "black",size=18),axis.title.x = element_text(colour = "black",size=20,vjust=0),
        axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1),
        plot.title=element_text(colour="black",size=22,vjust=1,hjust=.4))+
  theme(legend.title=element_text(colour="black",size=18),legend.text=element_text(colour="black",size=16))+
  labs(shape="simulation\ntype")



png(filename=paste0(saveDir,"/TPR_plot_",Sys.Date(),".png"),width=1400,height=1600,res=200);

plot(TPR_plot_point);

dev.off();