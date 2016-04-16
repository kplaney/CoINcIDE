
##code to analyze the output from Coincide run on simulated datasets from the 
#Coincide_tissueSimulation_script.R

####
simTypes <- c("highQuality", "highQualityUnevenSize", "highQualityUnevenNumClustMin2",
              "highQualityUnevenNumClustMin1","highQualityUnevenSizeUnevenNumClustMin2",
              "highQualityUnevenSizeUnevenNumClustMin1", "mixedQualityClust")

#run ROC curves for one simil level 
#usually zero above 0.7 
meanSimilVector <- c(0.0, 0.3, 0.5, 0.7)



#set theis path to where you saved the simulated data in the preceding 
#Coincide tissueSimulation script.

saveDirSims <-"/home/ywrfc09/simulations" 
library("RColorBrewer")
library("grid")
library("ggplot2")


#######mamking TPR, FPR pltos:
TPR_plot<- list()
FPR_plot <- list()

for(s in 1:length(simTypes)){

for(m in 1:length(meanSimilVector)){
  
  meanSimil <- meanSimilVector[m]
  
 message(paste0(saveDirSims,"/",simTypes[s],"_minSimil_",meanSimil,".rds"))
  data <- readRDS(paste0(saveDirSims,"/",simTypes[s],"_minSimil_",meanSimil,".rds"))
  #meanFractMaxTissueType[m,s] <- mean(data$commMeanMaxTissueType,na.rm=TRUE)

  
  if(m==1){
    
    masterDF <- data.frame(data$ROC_matrixFull,data$ROC_SDmatrixFull, rep( meanSimil,nrow(data$ROC_matrixFull)), rep(simTypes[s],nrow(data$ROC_matrixFull)));
    colnames(masterDF)[c(ncol(masterDF)-1, ncol(masterDF))] <- c("meanSimil","sim_type");

    
  }else{
    
    tmp <- data.frame(data$ROC_matrixFull,data$ROC_SDmatrixFull, rep( meanSimil,nrow(data$ROC_matrixFull)), rep(simTypes[s],nrow(data$ROC_matrixFull)));
    colnames(tmp)[c(ncol(masterDF)-1, ncol(tmp))] <- c("meanSimil","sim_type");
    masterDF <- rbind(masterDF,tmp);

    
  }
  
}

#dataFrames[[s]] <- masterDF

#make error bars:
TPR_limits <- aes(ymax = TPR + sd_TPR, ymin=TPR - sd_TPR)
FPR_limits <- aes(ymax = FPR + sd_FPR, ymin=FPR - sd_FPR)
#to make all same y limits
#+coord_cartesian(ylim=c(0,plotStackedYLimit)
#+scale_y_continuous(breaks=seq(0,1,0.25))


#masterDF$sim_type <- as.factor(masterDF$sim_type);
#TPR_plot <- ggplot(data =masterDF,aes(x=sd_noise,y=TPR,linetype=sim_type))+geom_line()+  labs(title = "TPR for simulations with increasing noise.",
#                                                                      y="TPR",x="sd noise")+
#  theme(panel.background = element_rect(fill='white', colour='black')) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                                                                               plot.title = element_text(size = rel(2)));


  dodge <- position_dodge(width=0.9)

#if(simTypes[s] !="mixedQualityClust"){
  
  yMax <- 1.1
  yMin <- -0.1
  
#}else{
  
#  yMax <- 0.6
#  yMin <- 0.0
  
#}


TPR_plot[[s]] <- ggplot(data =masterDF,aes(x=sd_noise,y=TPR))  + geom_line(size=.7)+ 
  labs(y="TPR",x="Noise level")+
  theme(panel.background = element_rect(fill='white', colour='black')) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(colour = "black",size=18),axis.title.x = element_text(colour = "black",size=22,vjust=0),
        axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=22,vjust=1))+
  coord_cartesian(ylim=c(yMin,yMax))+facet_grid(meanSimil~.)+
  geom_errorbar(TPR_limits, width=0.1)+ theme(panel.margin = unit(1.5, "lines"))+theme(strip.text.x = element_text(size = 14))

#title = "Network edge TPR for simulations with increasing noise",
#removed title font size from first theme:
#        plot.title=element_text(colour="black",size=22,vjust=1,hjust=.4))+
#removed title from theme:
#legend.text=element_text(colour="black",size=10))+
#  labs(colour="simulation\ntype")
options(bitmapType="cairo")
png(filename=paste0(saveDirSims,"/TPR_plot_",simTypes[s],"_",Sys.Date(),".png"),width=1400,height=1600,res=200);
#pdf(file=paste0(saveDirSims,"/TPR_plot_",simTypes[s],"_",Sys.Date(),".pdf"),width=6,height=9,paper="letter");
plot(TPR_plot[[s]]);

dev.off();

#on to FPR
  
  yMax <- 0.1
  yMin <- -0.1

FPR_plot[[s]] <- ggplot(data =masterDF,aes(x=sd_noise,y=FPR))  + geom_line(size=.7)+ 
  labs(y="FPR",x="Noise level")+
  theme(panel.background = element_rect(fill='white', colour='black')) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(colour = "black",size=18),axis.title.x = element_text(colour = "black",size=22,vjust=0),
        axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=22,vjust=1))+
  coord_cartesian(ylim=c(yMin,yMax))+facet_grid(meanSimil~.)+
  geom_errorbar(FPR_limits, width=0.1)+ theme(panel.margin = unit(1.5, "lines"))+theme(strip.text.x = element_text(size = 14))


options(bitmapType="cairo")
png(filename=paste0(saveDirSims,"/FPR_plot_",simTypes[s],"_",Sys.Date(),".png"),width=1400,height=1600,res=200);
#pdf(file=paste0(saveDirSims,"/FPR_plot_",simTypes[s],"_",Sys.Date(),".pdf"),width=6,height=9,paper="letter");
plot(FPR_plot[[s]]);

dev.off();

}


##now ROC plots - use all simil levels

meanFractMaxTissueType <- matrix(data=NA,ncol=length(simTypes),nrow=length(meanSimilVector),
                                 dimnames=list(meanSimilVector,simTypes))

#meanSimilVector <- seq(0, 1, length.out = 11)
meanSimilVector <- seq(1, 0, length.out = 11)


colorCodes <- brewer.pal(length(meanSimilVector),"Paired");
names(colorCodes) <- meanSimilVector
ROC_plot <- list()

dataFrames <- list()
for(s in 1:length(simTypes)){

for(m in 1:length(meanSimilVector)){
  
  meanSimil <- meanSimilVector[m]
  
 message(paste0(saveDirSims,"/",simTypes[s],"_minSimil_",meanSimil,".rds"))
  data <- readRDS(paste0(saveDirSims,"/",simTypes[s],"_minSimil_",meanSimil,".rds"))
 # meanFractMaxTissueType[m,s] <- mean(data$commMeanMaxTissueType,na.rm=TRUE)

  
  if(m==1){
    
    masterDF <- data.frame(data$ROC_matrixFull,data$ROC_SDmatrixFull, rep( meanSimil,nrow(data$ROC_matrixFull)), rep(simTypes[s],nrow(data$ROC_matrixFull)));
    colnames(masterDF)[c(ncol(masterDF)-1, ncol(masterDF))] <- c("meanSimil","sim_type");

    
  }else{
    
    tmp <- data.frame(data$ROC_matrixFull,data$ROC_SDmatrixFull, rep( meanSimil,nrow(data$ROC_matrixFull)), rep(simTypes[s],nrow(data$ROC_matrixFull)));
    colnames(tmp)[c(ncol(masterDF)-1, ncol(tmp))] <- c("meanSimil","sim_type");
    masterDF <- rbind(masterDF,tmp);

    
  }
  
}

dataFrames[[s]] <- masterDF

#hmm..remove certain noise levels to make it easier to plot?



#now ROC plots
#hmm...is data in the order you want, though?
#for ROC plots, we probably want it to be facetted by noise level
#and then we want min simil going from zero to 1.
#this is indeed how the data is already ordered.
#make the shape the mean simil level?

  yMax <- 1.1
  yMin <- -0.1

xMin <- -0.001
  xMax <- 0.02
 



#remove some noise levels for easier plotting
noiseLevelsKeep <- c(0.0, 0.8, 1.6, 2.4)
indicesKeep <- c()
for(i in 1:length(noiseLevelsKeep)){
  #some type casting issues...typical data frame!
  indicesKeep <- append(indicesKeep,which(as.numeric(as.character(masterDF$sd_noise))==noiseLevelsKeep[i]))
  
}
masterDF <- masterDF[indicesKeep,]

#removed line: geom_line(size=.7)+ 
#position = "jitter"
#do NOT jitter - on such a small scale, this really places values outside of their "zone"
ROC_plot[[s]] <- ggplot(data =masterDF,aes(x=FPR,y=TPR))  + geom_point(aes(color=factor(meanSimil),size=6))+ 
  labs(y="TPR",x="FPR")+ scale_color_manual(values=colorCodes)+
  theme(panel.background = element_rect(fill='white', colour='black'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(colour = "black",size=18),axis.title.x = element_text(colour = "black",size=22,vjust=0),
        axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=22,vjust=1))+
  theme(legend.text=element_text(colour="black",size=14))+
  coord_cartesian(ylim=c(yMin,yMax),xlim=c(xMin,xMax))+facet_wrap(~sd_noise)+
 theme(panel.margin = unit(3.5, "lines"))+  labs(color="Minimum\nsimilarity")+theme(strip.text.x = element_text(size = 14))


options(bitmapType="cairo")
png(filename=paste0(saveDirSims,"/ROC_plot_",simTypes[s],"_",Sys.Date(),".png"),width=1900,height=1700,res=200);
#pdf(file=paste0(saveDirSims,"/ROC_plot_",simTypes[s],"_",Sys.Date(),".pdf"),width=6,height=6,paper="letter");
plot(ROC_plot[[s]]);

dev.off();

#end of loops

}


names(dataFrames) <- simTypes
names(TPR_plot) <- simTypes
names(ROC_plot) <- simTypes
names(FPR_plot) <- simTypes
saveRDS(dataFrames[[s]],paste0(saveDirSims,"/dataFramesForPlotting_only4NoiseLevels.rds"),compress=TRUE)


