#library("limma");
#library("ggplot2");
#library("gplots")
#library("RColorBrewer")
#library("igraph");

#library("grid");
###correlation heatmaps
#png(filename=paste0(saveDir,paste0("sim",simNum,"_corr_heatmap_",Sys.Date(),".png")),
    #width = 700, height = 1000)

#if want no density plot: density.info="none",
#heatmap.2(cor(t(simData$dataMatrixList[[1]])),Rowv=F,Colv=F,scale = c("none"), 
 #         density.info="none",trace="none",symkey=TRUE,key.title="Color Key",key.xlab="correlation",
  #        col=rev(brewer.pal(11, "RdBu")),key=TRUE,labRow=F,labCol=F,ColSideColors=simSampleColors,
   #       RowSideColors=simSampleColors,main=paste0("\n\nSimulation ",simNum," :\ncorrelation heatmap \nby tissue type"));

#dev.off();

####



#if(!is.null(analysisOutput$subtypePlots$subtype_dfMaster)){
  
#  clustSizes <- table(analysisOutput$subtypePlots$subtype_dfMaster[,"clust_numVar"]);
  #put in order of attributes
#  clustSizes <- clustSizes[match(analysisOutput$communityMembership$attrDF$clust,names(clustSizes))];
#  attrDF <- cbind(analysisOutput$communityMembership$attrDF,clustSizes);
  
#}else{
  
#  attrDF <- analysisOutput$communityMembership$attrDF;
#}


##function to merge existing attributes in attrDF

plotMetaFeatureRank <- function(metaMatrix,saveFile=FALSE,plotToScreen=TRUE,
                             saveDir="./",fileTag="test",
                             plotTitle="Median rank\nacross all samples/studies",
                             key.xlab=""){

  #don't reverse: higher rank means higher value here. want higher rank values to be blue/cool,
  #as this means they had lower expression.
  colorCodeF <- brewer.pal(11,"RdBu")
  
  if(saveFile){
    
    jpeg(filename=paste0(saveDir,"/",fileTag,"_metaRankHeatmap.jpeg"),
         height=1000,width=700,quality=100,res=160);
    #    #ColSideColors=colColorCodes,
    heatmap.2(metaMatrix,
              Rowv=FALSE,Colv=FALSE,cexRow=1,srtRow=45, trace="none",
              key.title=plotTitle,key.xlab=key.xlab,scale="none",notecol="black",
              col=colorCodeF)
    
    dev.off()
  
  }
  
  if(plotToScreen){
    
    #ColSideColors=colColorCodes,
    heatmap.2(metaMatrix,
              Rowv=FALSE,Colv=FALSE,cexRow=1,trace="none",
              key.title=plotTitle,key.xlab=key.xlab,scale="none",notecol="black",
              col=colorCodeF,srtRow=45)
    
  }


}

plotMetaFeatureES <- function(ESMatrix,saveFile=FALSE,plotToScreen=TRUE,
                              saveDir="./",fileTag="test",
                              plotTitle="Significant Positive \nEffect size",
                              key.xlab=""){
  
  
  #set NA to zero values.
  ESMatrix[which(is.na(ESMatrix))] <- 0
  #let NA/zero to be light gray.
  colorCodeF <-  colorRampPalette(c("gray94","pink","red"), bias = 4,
                                  space = "rgb", interpolate = "linear");
  
  if(saveFile){
    
    jpeg(filename=paste0(saveDir,"/",fileTag,"_metaESHeatmap.jpeg"),
         height=1000,width=700,quality=100,res=160);
    #    #ColSideColors=colColorCodes,
    heatmap.2(ESMatrix,
              Rowv=FALSE,Colv=FALSE,cexRow=1,srtRow=45, trace="none",
              key.title=plotTitle,key.xlab=key.xlab,scale="none",notecol="black",
              col=colorCodeF)
    
    dev.off()
    
  }
  
  if(plotToScreen){
    
    #ColSideColors=colColorCodes,
    heatmap.2(ESMatrix,
              Rowv=FALSE,Colv=FALSE,cexRow=1,trace="none",
              key.title=plotTitle,key.xlab=key.xlab,scale="none",notecol="black",
              col=colorCodeF,srtRow=45)
    
  }
  
  
}

advancedNetworkPlots <- function(communityMembership,clustIndexMatrix,
                                  brewPal = c("Set3","Paired","Spectral","BrBG","PiYG","RdYlGn","RdYlBu","RdBu","PiYG","Set2"),
                                  saveDir="/home/kplaney/ISMB/",experimentName="networks",colorCodes,
                                 plotToScreen=FALSE,nodePlotSize=10,nodeFontSize=.7,plotEdgeWeight = TRUE,
                                 edgeWeightsColName="simil",nodeSizeScaleFactor=1){
  


   if(!plotToScreen){
     
     dir.create(saveDir)
     #dir.create(paste0(saveDir,"/",experimentName,"_",Sys.Date()),showWarnings=TRUE);
     #saveDir <- paste0(saveDir,"/",experimentName,"_",Sys.Date())
   }
  #study summary.
  #number of edges
  #number of studies
  #number of communities.
  #number of clusters (nodes) and # started out with.
  #let's add size of each cluster
  expName <- gsub("_"," ",experimentName)
  network_stats <- c(length(unique(communityMembership$attrDF$community)),
                     length(unique(communityMembership$attrDF$clust)),
                     nrow(clustIndexMatrix),
                     nrow(communityMembership$edgeDF),
                     length(unique(communityMembership$attrDF$studyNum)));
  
  names(network_stats) <- c("numCommunities","numClusters","origNumClusters","numEdges",
                            "numStudies"); 
  if(missing(colorCodes)){
    
    if(brewPal==FALSE){
      #make own color ramp.
      colorCodeF <- colorRampPalette(c("violet","blue","red"), bias = length(unique(communityMembership$attrDF[,"community"]))*2,
                                     space = "rgb", interpolate = "linear");
      #this will produce RGB representation in HEX, which cytoscape can take in.
      #if need plain old RGB: col2RGB(membership[,"colors"],alpha=FALSE);
      colorCodes <- colorCodeF(length(unique(communityMembership$attrDF[,"community"])));
      
    }else{
      #max is 11 colors
      colorCodes <- rev(brewer.pal(length(unique(communityMembership$attrDF[,"community"])),brewPal));
      
    }
    
  }
  #make sure is a character vector
  communityMembership$attrDF$color <- as.character(communityMembership$attrDF$color);
  
  #replace colors from user-selected color palette.
  for(c in 1:length(unique(communityMembership$attrDF[,"community"]))){
    
    
    communityMembership$attrDF[which(communityMembership$attrDF$community==c),"color"] <- colorCodes[c];
    
  }
  
  undirGraph <- graph.data.frame(communityMembership$edgeDF,directed=FALSE,vertices=communityMembership$attrDF)
  
  
  if(plotEdgeWeight){
    
    #advanced users can change this keyword if they have another data frame column they want to use.

 
      #multiple times ten- these original values will all be between -1 and 1
      
      if(max(get.edge.attribute(graph=undirGraph,name=edgeWeightsColName))<=1){
        
        #just scale up so can see easier.
        plotEdgeWeights <- 3*get.edge.attribute(graph=undirGraph,name=edgeWeightsColName)
        
      }else{
        
        stop("You selected to plot edge weights but your similarity (or distance) metric has large values above 1. This will not plot well.")
        
      }
    
    E(undirGraph)$weights <- plotEdgeWeights
    
  }
  
  if(!is.null(communityMembership$attrDF$size)){
    

    V(undirGraph)$size <- communityMembership$attrDF$size*nodeSizeScaleFactor;
    
    if(!plotToScreen){
      
    if(!plotEdgeWeight){
    #save plots
    png(filename=paste0(saveDir,"/",experimentName,"_communityPlot_scaledNodes_noLabels_",Sys.Date(),".png"),width=1000,height=1000,res=160)
        
        layout(matrix(c(1,2), 1,2), widths=c(3,1))
    plot(undirGraph, layout=layout.fruchterman.reingold,vertex.label=rep.int("",times=nrow(communityMembership$attrDF)),vertex.size=V(undirGraph)$size/10,vertex.label.color="black",vertex.label.cex=.7,
         vertex.color= V(undirGraph)$color,edge.arrow.size=3,main=expName,xlab="color=meta-cluster,node=cluster,size=#samples");
    
    
    plot.new()
    #want community numbers from smallest to largest.
    colorOrder <- sort.int(as.numeric(as.character(unique(V(undirGraph)$community))),index.return=TRUE,decreasing=FALSE)$ix
    colours = unique(V(undirGraph)$color)[colorOrder]
    #labels = paste(1:length(colours))
    labels = as.numeric(as.character(unique(V(undirGraph)$community)))[colorOrder]
    
    legend("center",legend=labels, col=colours, pch=19,pt.cex=2,
           title="Meta-clusters")
    dev.off();
    
    
    #with study numbers
    png(filename=paste0(saveDir,"/",experimentName,"_communityPlot_scaledNodes_studyNumlabels_",Sys.Date(),".png"),width=1000,height=1000,res=160)
    layout(matrix(c(1,2), 1,2), widths=c(3,1))
    plot(undirGraph, layout=layout.fruchterman.reingold,vertex.label=V(undirGraph)$studyNum,vertex.size=V(undirGraph)$size/10,vertex.label.color="black",vertex.label.cex=.7,
         vertex.color= V(undirGraph)$color, edge.arrow.size=3,main=expName,xlab="color=meta-cluster,node=cluster,#=dataset,node size=#samples");
    

    plot.new()
    #want community numbers from smallest to largest.
    colorOrder <- sort.int(as.numeric(as.character(unique(V(undirGraph)$community))),index.return=TRUE,decreasing=FALSE)$ix
    colours = unique(V(undirGraph)$color)[colorOrder]
    #labels = paste(1:length(colours))
    labels = as.numeric(as.character(unique(V(undirGraph)$community)))[colorOrder]
    
    legend("center",legend=labels, col=colours, pch=19,pt.cex=2,
           title="Meta-clusters")
    dev.off();
    
    #with edge weight
    }else{
      #edge.arrow.size=3,
      #save plots
      png(filename=paste0(saveDir,"/",experimentName,"_communityPlot_scaledNodes_noLabels_",Sys.Date(),".png"),width=1000,height=1000,res=160)
      
      layout(matrix(c(1,2), 1,2), widths=c(3,1))
      plot(undirGraph, layout=layout.fruchterman.reingold,vertex.label=rep.int("",times=nrow(communityMembership$attrDF)),vertex.size=V(undirGraph)$size/10,vertex.label.color="black",vertex.label.cex=.7,
           vertex.color= V(undirGraph)$color, edge.width = E(undirGraph)$weights,main=expName,xlab="color=meta-cluster,node=cluster,size=#samples");
      
      
      plot.new()
      #want community numbers from smallest to largest.
      colorOrder <- sort.int(as.numeric(as.character(unique(V(undirGraph)$community))),index.return=TRUE,decreasing=FALSE)$ix
      colours = unique(V(undirGraph)$color)[colorOrder]
      #labels = paste(1:length(colours))
      labels = as.numeric(as.character(unique(V(undirGraph)$community)))[colorOrder]
      
      legend("center",legend=labels, col=colours, pch=19,pt.cex=2,
             title="Meta-clusters")
      dev.off();
      
      
      #with study numbers
      png(filename=paste0(saveDir,"/",experimentName,"_communityPlot_scaledNodes_studyNumlabels_",Sys.Date(),".png"),width=1000,height=1000,res=160)
      layout(matrix(c(1,2), 1,2), widths=c(3,1))
      plot(undirGraph, layout=layout.fruchterman.reingold,vertex.label=V(undirGraph)$studyNum,vertex.size=V(undirGraph)$size/10,vertex.label.color="black",vertex.label.cex=.7,
           vertex.color= V(undirGraph)$color, edge.width = E(undirGraph)$weights,edge.arrow.size=3,main=expName,xlab="color=meta-cluster,node=cluster,#=dataset,node size=#samples");
      
      
      plot.new()
      #want community numbers from smallest to largest.
      colorOrder <- sort.int(as.numeric(as.character(unique(V(undirGraph)$community))),index.return=TRUE,decreasing=FALSE)$ix
      colours = unique(V(undirGraph)$color)[colorOrder]
      #labels = paste(1:length(colours))
      labels = as.numeric(as.character(unique(V(undirGraph)$community)))[colorOrder]
      
      legend("center",legend=labels, col=colours, pch=19,pt.cex=2,
             title="Meta-clusters")
      dev.off();
      
      
      #end if plot edge weights
    }

      #don't plot to screen
    }else{
      
      if(!plotEdgeWeight){
      
      layout(matrix(c(1,2), 1,2), widths=c(3,1))
          plot(undirGraph, layout=layout.fruchterman.reingold,vertex.label=rep.int("",times=nrow(communityMembership$attrDF)),vertex.size=V(undirGraph)$size/10,
               vertex.color= V(undirGraph)$color, edge.arrow.size=3,main=expName,xlab="color=meta-cluster,node=cluster,node size=#samples");
          
          
          plot.new()
          #want community numbers from smallest to largest.
          colorOrder <- sort.int(as.numeric(as.character(unique(V(undirGraph)$community))),index.return=TRUE,decreasing=FALSE)$ix
          colours = unique(V(undirGraph)$color)[colorOrder]
          #labels = paste(1:length(colours))
          labels = as.numeric(as.character(unique(V(undirGraph)$community)))[colorOrder]
          
          legend("center",legend=labels, col=colours, pch=19,pt.cex=2,
                 title="Meta-clusters")

          
          #with study numbers
      layout(matrix(c(1,2), 1,2), widths=c(3,1))
          plot(undirGraph, layout=layout.fruchterman.reingold,vertex.label=V(undirGraph)$studyNum,vertex.size=V(undirGraph)$size/10,vertex.label.color="black",vertex.label.cex=.7,
               vertex.color= V(undirGraph)$color, edge.arrow.size=3,main=expName,xlab="color=meta-cluster,node=cluster,#=dataset,node size=#samples");
          
          plot.new()
          #want community numbers from smallest to largest.
          colorOrder <- sort.int(as.numeric(as.character(unique(V(undirGraph)$community))),index.return=TRUE,decreasing=FALSE)$ix
          colours = unique(V(undirGraph)$color)[colorOrder]
          #labels = paste(1:length(colours))
          labels = as.numeric(as.character(unique(V(undirGraph)$community)))[colorOrder]
          
          legend("center",legend=labels, col=colours, pch=19,pt.cex=2,
                 title="Meta-clusters")
             
      }else{
        #plot edge weights
        layout(matrix(c(1,2), 1,2), widths=c(3,1))
        plot(undirGraph, layout=layout.fruchterman.reingold,vertex.label=rep.int("",times=nrow(communityMembership$attrDF)),vertex.size=V(undirGraph)$size/10,
             vertex.color= V(undirGraph)$color, edge.width = E(undirGraph)$weights,edge.arrow.size=3,main=expName,xlab="color=meta-cluster,node=cluster,node size=#samples");
        
        
        plot.new()
        #want community numbers from smallest to largest.
        colorOrder <- sort.int(as.numeric(as.character(unique(V(undirGraph)$community))),index.return=TRUE,decreasing=FALSE)$ix
        colours = unique(V(undirGraph)$color)[colorOrder]
        #labels = paste(1:length(colours))
        labels = as.numeric(as.character(unique(V(undirGraph)$community)))[colorOrder]
        
        legend("center",legend=labels, col=colours, pch=19,pt.cex=2,
               title="Meta-clusters")
        
        
        #with study numbers
        layout(matrix(c(1,2), 1,2), widths=c(3,1))
        plot(undirGraph, layout=layout.fruchterman.reingold,vertex.label=V(undirGraph)$studyNum,vertex.size=V(undirGraph)$size/10,vertex.label.color="black",vertex.label.cex=.7,
             vertex.color= V(undirGraph)$color, edge.width = E(undirGraph)$weights,edge.arrow.size=3,main=expName,xlab="color=meta-cluster,node=cluster,#=dataset,node size=#samples");
        
        plot.new()
        #want community numbers from smallest to largest.
        colorOrder <- sort.int(as.numeric(as.character(unique(V(undirGraph)$community))),index.return=TRUE,decreasing=FALSE)$ix
        colours = unique(V(undirGraph)$color)[colorOrder]
        #labels = paste(1:length(colours))
        labels = as.numeric(as.character(unique(V(undirGraph)$community)))[colorOrder]
        
        legend("center",legend=labels, col=colours, pch=19,pt.cex=2,
               title="Meta-clusters")
        
        
    }
    
  }
  
}
  #now plot without relative sizes.
  if(!plotToScreen){
    
    if(!plotEdgeWeight){
    
    png(filename=paste0(saveDir,"/",experimentName,"_communityPlot_unscaledNodes_nolabels_",Sys.Date(),".png"),width=1000,height=1000,res=160)
    
    layout(matrix(c(1,2), 1,2), widths=c(3,1))
    plot(undirGraph, layout=layout.fruchterman.reingold,vertex.label=rep.int("",times=nrow(communityMembership$attrDF)),vertex.size=nodePlotSize,vertex.label.color="black",vertex.label.cex=nodeFontSize,
         vertex.color= V(undirGraph)$color, edge.arrow.size=3,main=expName,xlab="color=meta-cluster,node=cluster");
    plot.new()
    #want community numbers from smallest to largest.
    colorOrder <- sort.int(as.numeric(as.character(unique(V(undirGraph)$community))),index.return=TRUE,decreasing=FALSE)$ix
    colours = unique(V(undirGraph)$color)[colorOrder]
    #labels = paste(1:length(colours))
    labels = as.numeric(as.character(unique(V(undirGraph)$community)))[colorOrder]
    
    legend("center",legend=labels, col=colours, pch=19,pt.cex=2,
           title="Meta-clusters")
    dev.off();
    #with study nums
    png(filename=paste0(saveDir,"/",experimentName,"_communityPlot_unscaledNodes_studyNumlabels_",Sys.Date(),".png"),width=1000,height=1000,res=160)
    layout(matrix(c(1,2), 1,2), widths=c(3,1))
    plot(undirGraph, layout=layout.fruchterman.reingold,vertex.label=V(undirGraph)$studyNum,vertex.size=nodePlotSize,vertex.label.color="black",vertex.label.cex=nodeFontSize,
         vertex.color= V(undirGraph)$color, edge.arrow.size=3,main=expName,xlab="color=meta-cluster,node=cluster,#=dataset");
    plot.new()
    #want community numbers from smallest to largest.
    colorOrder <- sort.int(as.numeric(as.character(unique(V(undirGraph)$community))),index.return=TRUE,decreasing=FALSE)$ix
    colours = unique(V(undirGraph)$color)[colorOrder]
    #labels = paste(1:length(colours))
    labels = as.numeric(as.character(unique(V(undirGraph)$community)))[colorOrder]
    
    legend("center",legend=labels, col=colours, pch=19,pt.cex=2,
           title="Meta-clusters")
    dev.off();
    
    #plot edgeWeights
    }else{
      
      png(filename=paste0(saveDir,"/",experimentName,"_communityPlot_unscaledNodes_nolabels_",Sys.Date(),".png"),width=1000,height=1000,res=160)
      
      layout(matrix(c(1,2), 1,2), widths=c(3,1))
      plot(undirGraph, layout=layout.fruchterman.reingold,vertex.label=rep.int("",times=nrow(communityMembership$attrDF)),vertex.size=nodePlotSize,vertex.label.color="black",vertex.label.cex=nodeFontSize,
           vertex.color= V(undirGraph)$color, edge.width = E(undirGraph)$weights, edge.arrow.size=3,main=expName,xlab="color=meta-cluster,node=cluster");
      plot.new()
      #want community numbers from smallest to largest.
      colorOrder <- sort.int(as.numeric(as.character(unique(V(undirGraph)$community))),index.return=TRUE,decreasing=FALSE)$ix
      colours = unique(V(undirGraph)$color)[colorOrder]
      #labels = paste(1:length(colours))
      labels = as.numeric(as.character(unique(V(undirGraph)$community)))[colorOrder]
      
      legend("center",legend=labels, col=colours, pch=19,pt.cex=2,
             title="Meta-clusters")
      dev.off();
      #with study nums
      png(filename=paste0(saveDir,"/",experimentName,"_communityPlot_unscaledNodes_studyNumlabels_",Sys.Date(),".png"),width=1000,height=1000,res=160)
      layout(matrix(c(1,2), 1,2), widths=c(3,1))
      plot(undirGraph, layout=layout.fruchterman.reingold,vertex.label=V(undirGraph)$studyNum,vertex.size=nodePlotSize,vertex.label.color="black",vertex.label.cex=nodeFontSize,
           vertex.color= V(undirGraph)$color,  edge.width = E(undirGraph)$weights,edge.arrow.size=3,main=expName,xlab="color=meta-cluster,node=cluster,#=dataset");
      plot.new()
      #want community numbers from smallest to largest.
      colorOrder <- sort.int(as.numeric(as.character(unique(V(undirGraph)$community))),index.return=TRUE,decreasing=FALSE)$ix
      colours = unique(V(undirGraph)$color)[colorOrder]
      #labels = paste(1:length(colours))
      labels = as.numeric(as.character(unique(V(undirGraph)$community)))[colorOrder]
      
      legend("center",legend=labels, col=colours, pch=19,pt.cex=2,
             title="Meta-clusters")
      dev.off();
      
      
      
      
      
    }
    
  }else{
    #plot to screen
    
    if(!plotEdgeWeight){
      
    layout(matrix(c(1,2), 1,2), widths=c(3,1))
    plot(undirGraph, layout=layout.fruchterman.reingold,vertex.label=rep.int("",times=nrow(communityMembership$attrDF)),vertex.size=nodePlotSize,vertex.label.color="black",vertex.label.cex=nodeFontSize,
         vertex.color= V(undirGraph)$color, edge.arrow.size=3,main=expName,xlab="color=meta-cluster,node=cluster");
    
    plot.new()
    #want community numbers from smallest to largest.
    colorOrder <- sort.int(as.numeric(as.character(unique(V(undirGraph)$community))),index.return=TRUE,decreasing=FALSE)$ix
    colours = unique(V(undirGraph)$color)[colorOrder]
    #labels = paste(1:length(colours))
    labels = as.numeric(as.character(unique(V(undirGraph)$community)))[colorOrder]
    
    legend("center",legend=labels, col=colours, pch=19,pt.cex=2,
           title="Meta-clusters")
    
    layout(matrix(c(1,2), 1,2), widths=c(3,1))
    plot(undirGraph, layout=layout.fruchterman.reingold,vertex.label=V(undirGraph)$studyNum,vertex.size=nodePlotSize,vertex.label.color="black",vertex.label.cex=nodeFontSize,
         vertex.color= V(undirGraph)$color, edge.arrow.size=3,main=expName,xlab="color=meta-cluster,node=cluster,#=dataset");
    
    plot.new()
    #want community numbers from smallest to largest.
    colorOrder <- sort.int(as.numeric(as.character(unique(V(undirGraph)$community))),index.return=TRUE,decreasing=FALSE)$ix
    colours = unique(V(undirGraph)$color)[colorOrder]
    #labels = paste(1:length(colours))
    labels = as.numeric(as.character(unique(V(undirGraph)$community)))[colorOrder]
    
    legend("center",legend=labels, col=colours, pch=19,pt.cex=2,
           title="Meta-clusters")
    
    #plot edge weights
    }else{
      
      layout(matrix(c(1,2), 1,2), widths=c(3,1))
      plot(undirGraph, layout=layout.fruchterman.reingold,vertex.label=rep.int("",times=nrow(communityMembership$attrDF)),vertex.size=nodePlotSize,vertex.label.color="black",vertex.label.cex=nodeFontSize,
           vertex.color= V(undirGraph)$color,  edge.width = E(undirGraph)$weights,edge.arrow.size=3,main=expName,xlab="color=meta-cluster,node=cluster");
      
      plot.new()
      #want community numbers from smallest to largest.
      colorOrder <- sort.int(as.numeric(as.character(unique(V(undirGraph)$community))),index.return=TRUE,decreasing=FALSE)$ix
      colours = unique(V(undirGraph)$color)[colorOrder]
      #labels = paste(1:length(colours))
      labels = as.numeric(as.character(unique(V(undirGraph)$community)))[colorOrder]
      
      legend("center",legend=labels, col=colours, pch=19,pt.cex=2,
             title="Meta-clusters")
      
      layout(matrix(c(1,2), 1,2), widths=c(3,1))
      plot(undirGraph, layout=layout.fruchterman.reingold,vertex.label=V(undirGraph)$studyNum,vertex.size=nodePlotSize,vertex.label.color="black",vertex.label.cex=nodeFontSize,
           vertex.color= V(undirGraph)$color,  edge.width = E(undirGraph)$weights,edge.arrow.size=3,main=expName,xlab="color=meta-cluster,node=cluster,#=dataset");
      
      plot.new()
      #want community numbers from smallest to largest.
      colorOrder <- sort.int(as.numeric(as.character(unique(V(undirGraph)$community))),index.return=TRUE,decreasing=FALSE)$ix
      colours = unique(V(undirGraph)$color)[colorOrder]
      #labels = paste(1:length(colours))
      labels = as.numeric(as.character(unique(V(undirGraph)$community)))[colorOrder]
      
      legend("center",legend=labels, col=colours, pch=19,pt.cex=2,
             title="Meta-clusters")
      
      
      
    }
    
  }
  output <- list(undirGraph=undirGraph,attrDF=communityMembership$attrDF,network_stats=network_stats);
  return(output);
  
}

#######


assignCentroidSubtype <- function(origDataMatrix,minNumGenes=30,centroidRData="/home/data/breast_microarrayDB/pam50_centroids.RData"){
  
  message("\nAssumes genes are in the columns of your data matrix, \nand for centroids, gene names are in the first columns.\n")
  load(centroidRData);
  #originally meant for use with pam50 centroids, but can be used on any matrix where first col is gene names, other columns are centroid groups.
  #RData object must be names centroidMatrix
  if(missing(centroidMatrix)){
    
    stop("\nRData object must have name 'centroidMatrix'.")
 
  }

  centroids <- centroidMatrix[, 2:ncol(centroidMatrix)];
  #only take the genes that are in this data matrix
  centroidMatrix_genes <- as.character(centroidMatrix[,1][centroidMatrix[,1] %in% colnames(origDataMatrix)]);
  
  if(length(centroidMatrix_genes)>minNumGenes){
    
    subtypes <- array(data=NA,dim=nrow(origDataMatrix));
    centroids <- centroids[centroidMatrix[,1] %in% colnames(origDataMatrix), ];
    dataMatrix <- origDataMatrix[ ,centroidMatrix_genes];
    corMatrix <- cor(t(dataMatrix),centroids);
    
    for(s in 1:length(subtypes)){
      #sometimes, if no variance for a list of genes: will return NA.
      subtypes[s] <- which(corMatrix[s,]==max(corMatrix[s,],na.rm=TRUE));
      
    }
    
    subtypeNames <- colnames(centroidMatrix)[2:ncol(centroidMatrix)];
    subtypeLabels <- array(data=NA,dim=nrow(dataMatrix));
    
    for(c in 1:length(unique(subtypeNames))){
      
      subtypeLabels[which(subtypes==c)] <- subtypeNames[c];
      
      
    }
    
    subtypes <- cbind(subtypes,subtypeLabels);
    output <- list(subtypes=subtypes,centroidMatrix_genes=centroidMatrix_genes);
    return(output);
    
  }else{
    
    return(NA);
    
  }

  return(output);
  
}
####

#pam50 colors:  t(data.frame(c("#970eec","#ec0e58","#0e58ec","#7d7d7d","#0eec97")))
#(and name for subtypes-look up order in old clust robust code.)
phenoVarCommunityBreakdownPlots <- function(sampleClustCommKey,saveDir="./",fileTag="clust method",
                                                variableColorMatrix=NULL,phenoVar){
  



  if(!is.data.frame(sampleClustCommKey)){
    
    stop("sampleClustCommKey must be a data frame for ggplot. \nMake sure that numeric and character variables are not altered if you force it to be a data frame.")
  }
  
  if(length(sampleClustCommKey[,phenoVar])==0){
    
    stop("phenoVar input does not appear to be a column name in the sampleClustCommKey data frame inputted.")
    
  }
  
  if(is.null(variableColorMatrix)){
    
    numColors <- length(unique(sampleClustCommKey[,phenoVar]))
    
    if(numColors >12){
      
      warning("You have more than 12 unique values in your phenoVar variable and thus graded colors and not an RColorBrewer palette will be used.")
      
      colorCodeF <- colorRampPalette(c("violet","blue","red"), bias = numColors*2,
                                     space = "rgb", interpolate = "linear");
      #this will produce RGB representation in HEX, which cytoscape can take in.
      #if need plain old RGB: col2RGB(membership[,"colors"],alpha=FALSE);
      colorCodes <- colorCodeF(numColors)
      
    }else{
      #this is a more distinct color palette
      variableColorMatrix <- brewer.pal(numColors,"Paired")
    }
 
    names(variableColorMatrix) <- unique(sampleClustCommKey[,phenoVar]);
  }
  #pam50 colors:  t(data.frame(c("#970eec","#ec0e58","#0e58ec","#7d7d7d","#0eec97")))

 
  subtype_plot <- list();
  subtype_dfList <- list();
  subtype_plot_stacked <- list();
  subtype_plot_fract <- list();
  
  #remove the levels.
  communityNames <- unique(sampleClustCommKey$community)
  
  for(c in 1:length(communityNames)){
    
    #community indices are same as that in the master list (same as those in IGP matrices)
    #CAREFUL: the levels make this weird...so just change to characters to match the names of the biclustMasterList.
    clust_nums <- as.numeric(as.character(community_membership[which(community_membership[,"community"]==communityNames[c]), "clust"]));
    biclustCommunities[[c]] <- biclustMasterList[clust_nums];
    
    #want to remove levels! otherwise indexing gets all messed up.
    communityStudyNums[[c]] <- as.numeric(as.character(community_membership[which(community_membership[,"community"]==communityNames[c]), "studyNum"]));
    
    subtype_dfList[[c]] <- data.frame();

    
    
    subtype_plot[[c]] <-  ggplot(data=sampleClustCommKey[which(sampleClustCommKey$community==c), ],aes(x=phenoVar))+geom_histogram(aes(fill=phenoVar))+
      labs(y="Number of samples" ,title=paste0(phenoVar," for community ",c))+
      scale_fill_manual(values = variableColorMatrix)+ labs(fill=phenoVar)+
      theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
      theme(axis.text.x = element_text(colour = "black",size=18,angle=45,vjust=.9,hjust=1),axis.title.x=element_blank(),
            axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1))+
      theme(plot.title=element_text(colour="black",size=28,vjust=1,hjust=.2));
    
    
    #stacked bar for all studies.
    subtype_plot_stacked[[c]] <- ggplot(data=sampleClustCommKey[which(sampleClustCommKey$community==c), ],aes(x=factor(studyNum)))+geom_bar(aes(fill=factor(phenoVar)))+
      labs(y="Number of samples", x="Dataset number",title=paste0(phenoVar," by dataset for community ",c))+
      scale_fill_manual(values = variableColorMatrix)+labs(fill="phenoVar")+
      theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
      theme(axis.text.x = element_text(colour = "black",size=12,angle=45,vjust=1,hjust=1),
            axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1),
            axis.title.x = element_text(colour = "black",size=18,vjust=0))+
      theme(plot.title=element_text(colour="black",size=28,vjust=1,hjust=.2));
    
    #calculate percentage
    df <- split(subtype_dfList[[c]],f=subtype_dfList[[c]]$studyNum);
    phenoVar_fract <- lapply(df,FUN=function(e){
      
      output <- table(e$subtype)/nrow(e);
      
    });
    
    
    for(s in 1:length(phenoVar_fract)){
      
      for(t in 1:length(phenoVar_fract[[s]])){
        
        if(t >1){
          
          tmp <- rbind(tmp,data.frame(rep(names(phenoVar_fract)[s],round(phenoVar_fract[[s]][t]*100)),rep(names(phenoVar_fract[[s]])[t],round(phenoVar_fract[[s]][t]*100)),
                                      rep(c,round(subtype_fract[[s]][t]*100))));
        }else{
          
          tmp <- data.frame(rep(names(phenoVar_fract)[s],round(phenoVar_fract[[s]][t]*100)),rep(names(phenoVar_fract[[s]])[t],round(phenoVar_fract[[s]][t]*100)),
                            rep(c,round(phenoVar_fract[[s]][t]*100)));
        }
        
      }
      
      if(s >1){
        
        phenoVar_fractDF  <- rbind(phenoVar_fractDF,tmp);
        
      }else{
        
        phenoVar_fractDF  <- tmp;
        
      }
      
    }
    colnames(phenoVar_fractDF) <- c("studyNum",phenoVar,"community");
    #now get percentages down to 100, not rounded to 101.
    tmp <- table(phenoVar_fractDF$studyNum);
    
    for(t in 1:length(tmp)){
      
      if(tmp[t]==101){
        #need to remove patients (should be just 1) to make 100.
        
        phenoVar_fractDF <-  phenoVar_fractDF[-which(phenoVar_fractDF$studyNum==names(tmp[t]))[1],];
      }
      
    }
    colnames(phenoVar_fractDF) <- c("studyNum",phenoVar,"community");
    
    phenoVar_plot_fract[[c]]<-  ggplot(data=phenoVar_fractDF,aes(x=studyNum))+geom_bar(aes(fill=factor(phenoVar)))+
      labs(x="Dataset number",y="Composition of subtypes",title=paste0("Composition of ",phenoVar," \nby dataset for community ",c))+
      labs(fill=phenoVar)+
      scale_fill_manual(values = variableColorMatrix)+  theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
      theme(axis.text.x = element_text(colour = "black",size=12,angle=45,vjust=1,hjust=1),
            axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1),
            axis.title.x = element_text(colour = "black",size=18,vjust=0))+
      theme(plot.title=element_text(colour="black",size=28,vjust=1,hjust=.5));
    
    if(c>1){
      
      
      phenoVar_fractMaster <- rbind(phenoVar_fractMaster,phenoVar_fractDF);
      
    }else{
      
      phenoVar_fractMaster <- phenoVar_fractDF;
      
    }
    
  }

  
 
  #stacked bar for all studies.
  #scales="free_x":  doesn't plot studies that are empty for that community.
  phenoVar_plot_stackedALL <- ggplot(data=subtype_dfMaster,aes(x=community))+geom_bar(aes(fill=factor(phenoVar)))+
    scale_fill_manual(values = variableColorMatrix)+
    labs(y="Number of samples",x=element_blank(),title=phenoVar)+theme_bw()+
    labs(fill="phenoVar")+theme(strip.text.x = element_text(size = 17, colour = "black",face="bold"))+
    theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
    theme(axis.text.x = element_text(colour = "black",size=18,angle=45,vjust=1,hjust=1),axis.title.x= element_text(colour = "black",size=20,vjust=1),
          axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1))+
    theme(plot.title=element_text(colour="black",size=28,vjust=1,hjust=.2))+theme(panel.background = element_rect(colour = "black"))+
    theme(panel.grid.major = element_line(colour = 0),
          panel.grid.minor = element_line(colour = 0));
  
  
  
  png(filename=paste0(saveDir,"/",fileTag,"_","plot_stacked_",Sys.Date(),".png"),
      width = 700, height = 1000,res=160);
  
  plot( phenoVar_plot_stackedALL);
  dev.off();
  
  
  ##make some subtype breakdown, community stats
 
  tmp <- split(sampleClustCommKey,f=df$community);
  
  for(t in 1:length(tmp)){
    
    if(t>1){
      
      phenoVarBreakdowns <- rbind(subtypeBreakdowns,data.frame(table(tmp[[t]][,phenoVar]),nrow(tmp[[t]]),t,length(unique(tmp[[t]][,c("studyNum")])),length(unique(tmp[[t]][,c("clustNum")]))));
      
    }else{
      
      phenoVarBreakdowns <- data.frame(table(tmp[[t]][,phenoVar]),nrow(tmp[[t]]),t,length(unique(tmp[[t]][,c("studyNum")])),length(unique(tmp[[t]][,c("clustNum")])));
    }
    
  }
  #COME BACK: add in # of studies?  this is already in network stats, though...
  colnames(phenoVarBreakdowns) <- c(phenoVar,"number","communitySize","community","numStudyPerComm","numClustPerComm");
  phenoVarBreakdowns <- cbind(phenoVarBreakdownsphenoVarBreakdowns$number/phenoVarBreakdowns$communitySize);
  colnames(phenoVarBreakdowns) <- c(phenoVar,"number","communitySize","community","numStudyPerComm","numClustPerComm","fract");

  #plot some piechars too.
  p <- ggplot(phenoVarBreakdowns, aes(x=1,y=fract, fill=subtype)) +facet_grid(.~community,scales="free_x")+

    geom_bar(stat="identity", color='black') +scale_fill_manual(values = t(phenoVar))+
    # remove black diagonal line from legend
    guides(fill=guide_legend(override.aes=list(colour=NA))) + theme(strip.text.x = element_text(size = 17, colour = "black",face="bold"))+
    theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15));
  
  p <- p +coord_polar(theta='y') +theme(axis.ticks=element_blank(),  # the axis ticks
                                        axis.title=element_blank(),  # the axis labels
                                        axis.text.y=element_blank(), # the 0.75, 1.00, 1.25 labels.
                                        axis.text.x=element_blank())
  
  png(filename = paste0(saveDir,"/",fileTag,"_pieCharts",Sys.Date(),".png"),
      width = 1000, height = 500,res=160);
  plot(p);
  dev.off();
  
  #just take first match for each community - values will be same across an entire community.
  communityStats <- phenoVarBreakdowns[match(unique(phenoVarBreakdowns$community),subtypeBreakdowns$community), ]
  communityStats <- communityStats[ ,c("communitySize","community","numStudyPerComm","numClustPerComm")];
  rownames(communityStats) <- communityStats$community;
  output <- list(phenoVar=phenoVar,phenoVar_plot_fract=phenoVar_plot_fract,
                 phenoVar_plot_stacked=phenoVar_plot_stacked,phenoVar_plot=phenoVar_plot,phenoVar_plot_stackedALL= phenoVar_plot_stackedALL,
                 pieChartPlot=p,communityStats=communityStats,phenoVarBreakdowns=phenoVarBreakdowns);
  return(output);
  
}
# communitySubtypeBreakdown_plots <- function(community_membership,biclust,origDataMatrices,saveDir="./",experimentName="clust method",
#                                              centroidSetName="pam50",centroidRData="/home/data/breast_microarrayDB/pam50_centroids.RData"){
#   
#   
#   #specify colors:
#   phenoVar <- data.frame(c("#970eec","#ec0e58","#0e58ec","#7d7d7d","#0eec97"))
#   #hmm can I order the colors in the same way each time?? may not really need to anyways.
#   #must tranpose for ggplot to read it correctly and THEN add names for this to work.
#   subtypeColorMatrix <- t(subtypeColorMatrix);
#   names(subtypeColorMatrix) <- c("LumB","LumA","Her2","Basal","Normal");
#   
#   numTotalSamplesInStudy <- list();
#   for(s in 1:length(origDataMatrices)){
#     
#     numTotalSamplesInStudy[[s]] <- nrow(origDataMatrices[[s]]);
#     
#   }
#   library("limma");
#   library("ggplot2");
#   warning("\nThis code assumes your biclusters have genes in the columns and patients in the rows.\n")
#   #library("RColorBrewer");
#   biclustMasterList <- list();
#   biclustMasterList_origStudies <- list();
#   
#   #IS THIS WORKING?? cluster 83 should have way more patients!!
#   for(b in 1:length(biclust)){
#     
#     biclustMasterList <- append(biclustMasterList, biclust[[b]]);
#     biclustMasterList_origStudies <- append(biclustMasterList_origStudies,rep(b,length(biclust[[b]])));
#   }
#   names(biclustMasterList) <- c(1:length(biclustMasterList));
#   
#   #now create lists by communities. re-name to master biclust # (or underscore?)
#   biclustCommunities <- list();
#   communityNames <- unique(community_membership[,"community"]);
#   communityStudyNums <- list();
#   
#   
#   subtype_plot <- list();
#   subtype_dfList <- list();
#   subtype_plot_stacked <- list();
#   subtype_plot_fract <- list();
#   
#   #remove the levels.
#   communityNames <- as.character(communityNames);
#   for(c in 1:length(communityNames)){
#     
#     #community indices are same as that in the master list (same as those in IGP matrices)
#     #CAREFUL: the levels make this weird...so just change to characters to match the names of the biclustMasterList.
#     clust_nums <- as.numeric(as.character(community_membership[which(community_membership[,"community"]==communityNames[c]), "clust"]));
#     biclustCommunities[[c]] <- biclustMasterList[clust_nums];
#     
#     #want to remove levels! otherwise indexing gets all messed up.
#     communityStudyNums[[c]] <- as.numeric(as.character(community_membership[which(community_membership[,"community"]==communityNames[c]), "studyNum"]));
#     
#     subtype_dfList[[c]] <- data.frame();
#     
# 
#     for(e in 1:length(biclustCommunities[[c]])){
#       
#       patIDs <- rownames(biclustCommunities[[c]][[e]]);
#       
#       studyNum <- rep.int(communityStudyNums[[c]][[e]],times=length(patIDs));
#       subtype <- assignCentroidSubtype(origDataMatrix=origDataMatrices[[communityStudyNums[[c]][e]]][patIDs,],minNumGenes=30,centroidRData=centroidRData)$subtypes[,"subtypeLabels"];
#       clust_numVar <- rep(clust_nums[e],length(subtype));
#       temp <- data.frame(patIDs,subtype,studyNum,clust_numVar);
#       subtype_dfList[[c]] <- rbind(subtype_dfList[[c]],temp); 
# 
#     }
#     
#     if(c >1){
#       
#       tmp <- data.frame(subtype_dfList[[c]],rep(c,times=nrow(subtype_dfList[[c]])));
#       colnames(tmp)[ncol(tmp)] <- "community";
#       
#       subtype_dfMaster <- rbind(subtype_dfMaster,tmp);
#       
#     }else{
#       
#       tmp <- data.frame(subtype_dfList[[c]],rep(c,times=nrow(subtype_dfList[[c]])));
#       colnames(tmp)[ncol(tmp)] <- "community";
#       subtype_dfMaster <- tmp;
#       #trying to create percentages.... a column for subtype, for the percentage, the community num,
#       #and the studNum
#       #subtype_dfSum <- data.frame(unique(subtype),table(subtype)/length(subtype),)
#       
#     }
#  
#     
#     subtype_plot[[c]] <-  ggplot(data=subtype_dfList[[c]],aes(x=subtype))+geom_histogram(aes(fill=subtype))+
#       labs(y="Number of samples" ,title=paste0(centroidSetName," subtypes for community ",c))+
#       scale_fill_manual(values = subtypeColorMatrix)+ labs(fill="subtype")+
#       theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
#       theme(axis.text.x = element_text(colour = "black",size=18,angle=45,vjust=.9,hjust=1),axis.title.x=element_blank(),
#             axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1))+
#       theme(plot.title=element_text(colour="black",size=28,vjust=1,hjust=.2));
# 
#     #stacked bar for all studies.
#     subtype_plot_stacked[[c]] <- ggplot(data=subtype_dfList[[c]],aes(x=factor(studyNum)))+geom_bar(aes(fill=factor(subtype)))+
#       labs(y="Number of samples", x="Dataset number",title=paste0(centroidSetName," subtype by dataset for community ",c))+
#       scale_fill_manual(values = subtypeColorMatrix)+labs(fill="subtype")+
#       theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
#       theme(axis.text.x = element_text(colour = "black",size=12,angle=45,vjust=1,hjust=1),
#             axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1),
#             axis.title.x = element_text(colour = "black",size=18,vjust=0))+
#       theme(plot.title=element_text(colour="black",size=28,vjust=1,hjust=.2));
#     
#     #calculate percentage
#     df <- split(subtype_dfList[[c]],f=subtype_dfList[[c]]$studyNum);
#     subtype_fract <- lapply(df,FUN=function(e){
#       
#       output <- table(e$subtype)/nrow(e);
#       
#     });
#     
#     
#     for(s in 1:length(subtype_fract)){
#       
#       for(t in 1:length(subtype_fract[[s]])){
#         
#         if(t >1){
#           
#           tmp <- rbind(tmp,data.frame(rep(names(subtype_fract)[s],round(subtype_fract[[s]][t]*100)),rep(names(subtype_fract[[s]])[t],round(subtype_fract[[s]][t]*100)),
#                                       rep(c,round(subtype_fract[[s]][t]*100))));
#         }else{
#           
#           tmp <- data.frame(rep(names(subtype_fract)[s],round(subtype_fract[[s]][t]*100)),rep(names(subtype_fract[[s]])[t],round(subtype_fract[[s]][t]*100)),
#                             rep(c,round(subtype_fract[[s]][t]*100)));
#         }
#         
#       }
#       
#       if(s >1){
#         
#         subtype_fractDF  <- rbind(subtype_fractDF,tmp);
#         
#       }else{
#         
#         subtype_fractDF  <- tmp;
#         
#       }
#       
#     }
#     colnames( subtype_fractDF) <- c("studyNum","subtype","community");
#     #now get percentages down to 100, not rounded to 101.
#     tmp <- table(subtype_fractDF$studyNum);
#     
#     for(t in 1:length(tmp)){
#       
#       if(tmp[t]==101){
#         #need to remove patients (should be just 1) to make 100.
#         
#         subtype_fractDF <-  subtype_fractDF[-which(subtype_fractDF$studyNum==names(tmp[t]))[1],];
#       }
#       
#     }
#     colnames( subtype_fractDF) <- c("studyNum","subtype","community");
# 
#     subtype_plot_fract[[c]]<-  ggplot(data=subtype_fractDF,aes(x=studyNum))+geom_bar(aes(fill=factor(subtype)))+
#       labs(x="Dataset number",y="Composition of subtypes",title=paste0("Composition of ",centroidSetName," subtypes \nby dataset for community ",c))+
#       labs(fill="subtype")+
#     scale_fill_manual(values = subtypeColorMatrix)+  theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
#       theme(axis.text.x = element_text(colour = "black",size=12,angle=45,vjust=1,hjust=1),
#             axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1),
#             axis.title.x = element_text(colour = "black",size=18,vjust=0))+
#       theme(plot.title=element_text(colour="black",size=28,vjust=1,hjust=.5));
#     
#     if(c>1){
#       
#       
#       subtype_fractMaster <- rbind(subtype_fractMaster,subtype_fractDF);
#       
#     }else{
#       
#       subtype_fractMaster <- subtype_fractDF;
#       
#     }
#     
#   }
#   
#   #do a facetted plot by all communities too
#   subtype_plotALL <- ggplot(data=subtype_dfMaster,aes(x=subtype))+geom_histogram(aes(fill=subtype))+facet_grid(.~community,scales="free_x")+
#     labs(y="Number of samples", x="Pam50 subtype breakdown",title=paste0(centroidSetName," subtype for each community"))+
#     scale_fill_manual(values = subtypeColorMatrix)+
#     labs(fill="subtype")+theme(strip.text.x = element_text(size = 17, colour = "black",face="bold"))+
#     theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
#     theme(axis.text.x = element_text(colour = "black",size=12,angle=45,vjust=1,hjust=1),axis.title.x= element_text(colour = "black",size=20,vjust=1),
#           axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1))+
#     theme(plot.title=element_text(colour="black",size=28,vjust=1,hjust=.2));
#   
#   png(filename=paste0(saveDir,"/",experimentName,"_subtype_plotALL_",Sys.Date(),".png"),
#       width = 700, height = 1000);
#   
#   plot(subtype_plotALL);
#   dev.off();
#   
#   library("grid");
#   #stacked bar for all studies.
#   #scales="free_x":  doesn't plot studies that are empty for that community.
#   subtype_plot_stackedALL <- ggplot(data=subtype_dfMaster,aes(x=factor(studyNum)))+geom_bar(aes(fill=factor(subtype)))+
#     facet_grid(facets=.~community,scales="free_x")+scale_fill_manual(values = subtypeColorMatrix)+
#     labs(y="Number of samples",x="Dataset number",title=paste0(centroidSetName," subtype by community by dataset"))+
#   labs(fill="subtype")+theme(strip.text.x = element_text(size = 17, colour = "black",face="bold"))+
#     theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
#     theme(axis.text.x = element_text(colour = "black",size=6,angle=45,vjust=1,hjust=1),axis.title.x= element_text(colour = "black",size=20,vjust=1),
#           axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1))+
#     theme(plot.title=element_text(colour="black",size=28,vjust=1,hjust=.2));
#    
#   
#   
#   png(filename=paste0(saveDir,"/",experimentName,"_subtype_plot_stackedALL_",Sys.Date(),".png"),
#       width = 700, height = 1000);
#   
#   plot(subtype_plot_stackedALL);
#   dev.off();
# 
#   #calculate percentage
#   subtype_plot_fractALL <- ggplot(data=subtype_fractMaster,aes(x=studyNum))+geom_bar(aes(fill=factor(subtype)))+facet_grid(.~community,scales="free_x")+
#     labs(xlim=100,x="Dataset number",y="Composition of subtypes",title=paste0("Composition of ",centroidSetName," subtype by community by dataset"))+
#     scale_fill_manual(values = subtypeColorMatrix)+
#     labs(fill="subtype")+theme(strip.text.x = element_text(size = 17, colour = "black",face="bold"))+
#     theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
#     theme(axis.text.x = element_text(colour = "black",size=6,angle=45,vjust=1,hjust=1),axis.title.x= element_text(colour = "black",size=20,vjust=1),
#           axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1))+
#     theme(plot.title=element_text(colour="black",size=28,vjust=1,hjust=.2));
#     
#   
#   png(filename=paste0(saveDir,"/",experimentName,"_subtype_plot_fractALL_",Sys.Date(),".png"),
#       width = 700, height = 1000);
#   plot(subtype_plot_fractALL);
#   dev.off();
#   
#   warning("Just becuase samples are from the same study, doesn't mean they weren't clustered into 2+ groups. 
#           But these groups were deemed similar using my meta-cluster method.\n");
#   
#   output <- list(subtype_dfMaster=subtype_dfMaster,subtype_dfList=subtype_dfList,subtype_plot=subtype_plot,
#                  subtype_plot_stacked=subtype_plot_stacked,subtype_plot_fract=subtype_plot_fract,
#                  subtype_fractMaster =subtype_fractMaster,subtype_plot_fractALL=subtype_plot_fractALL,
#                  subtype_plot_stackedALL=subtype_plot_stackedALL,subtype_plotALL=subtype_plotALL);
#   return(output);
#   
# }

#TO DO:  copy most of this function to allow to bucket by ANY dataset attribute you input
#(i.e. a column found in data frame you feed in - feed in a dataset attribute table,
#and the edge Attributes table, and then can ID what nodes are in what datasets.
#allow for more color schemes.
plot_communitiesWithSubtypes  <- function(subtype_dfMaster,community_membership,igraph_edgeDF,saveDir="./"){
    
    #add subtypes to communities
    community_membership <- cbind(community_membership,array(data=NA,dim=nrow(community_membership)));
    colnames(community_membership)[ncol(community_membership)] <- "subtype";
    #remove the levels. as.numeric messes it up further unless use as.character first...
    clustNum <- as.numeric(as.character(unique(community_membership[,"clust"])));
    
    if(length(clustNum)>0){
    for(s in 1:length(clustNum)){
      
      tmp <- subtype_dfMaster[which(subtype_dfMaster[,"clust_numVar"]==clustNum[s]), ];
      
      if(nrow(tmp)==0){
        
        stop("\nNot finding any samples that link to this cluster number.");
      }
      
     #assume each cluster can only be in 1 community.    
      subtypes <- table(tmp[,"subtype"]);
      #just take the max # as the subtype for this cluster.
      if(!all(is.null(names(subtypes)))){
        #if tie: just take first.
        community_membership[which(community_membership[,"clust"]==clustNum[s]),"subtype"] <- names(subtypes)[which(subtypes==max(subtypes))[1]]
          
      }else{
       
        community_membership[which(community_membership[,"clust"]==clustNum[s]),"subtype"] <- which(subtypes==max(subtypes))[1];
       
     }
     
  }
    }
  
  undirGraph <- graph.data.frame(igraph_edgeDF,directed=FALSE,vertices=community_membership);
  #color the nodes. V() gets the vertices. can do V()$ what added  with igraph_attrDF, like $clust.
  #vertex.size=4,vertex.label.dist=0.5,edge.weight=E(undirGraph)$sampleFract,
  #size is size of vertex
  finalNumCommunities <- length(unique(community_membership[,"community"]));
  png(filename=paste0(saveDir,experimentName,"_communitySubtypePlot_",Sys.Date(),".png"),
      width = 700, height = 1000);
  # vertex.color= V(undirGraph)$color,
  plot(undirGraph, layout=layout.fruchterman.reingold,vertex.label=V(undirGraph)$subtype,vertex.size=5,vertex.label.color="black",vertex.label.cex=.7,
                          vertex.color= V(undirGraph)$color, edge.arrow.size=3,main=paste0("Clusters deemed similar across ",length(unique(community_membership[,"studyNum"])), "studies.\n",
                                                                                           finalNumCommunities," pruned communities found using method ",experimentName, "."),xlab="color=community, label=subtype");
  
  
  dev.off();
  

  output <- list(undirGraph=undirGraph,community_membership=community_membership,finalNumCommunities=finalNumCommunities); 
  
  return(output);
  
  }

plot_subtypes_network <- function(clustRobust_output,analysisOutput,
                                  experimentName,saveDir="/home/kplaney/pre_proposal/",
                                  centroidRData="/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData"){
  
  #source("/home/kplaney/gitRepos/IGP_network/igp_network/pre_proposal_2015_new_communitySubtypeBreakdown_plots.R")
  
  
  commDF <- analysisOutput$communityMembership$attrDF;
  clustMatrixList <- clustRobust_output$clustMatrixList;
  origDataMatrices <- clustRobust_output$dataMatrixList;
  
  subtype_plots <- new_communitySubtypeBreakdown_plots(community_membership=commDF,biclust=clustMatrixList,
                                                       origDataMatrices=origDataMatrices,saveDir=saveDir,
                                                       experimentName=experimentName,
                                                       centroidSetName="",centroidRData=centroidRData);
  
  write.table(subtype_plots$communityStats,file=paste0(saveDir,"/",experimentName,"_communityStats_",Sys.Date(),".txt"),row.names=FALSE,
              quote=FALSE);
  network_plots <- advanced_networkPlots(analysisOutput=analysisOutput,
                                         brewPal <- c("Set1"),
                                         saveDir=saveDir,experimentName=experimentName);
  
  write.table(t(as.matrix(network_plots$network_stats)),file=paste0(saveDir,"/",experimentName,"_networkStats_",Sys.Date(),".txt"),row.names=FALSE,
              quote=FALSE);
  
  output <- list(subtype_plots=subtype_plots,network_plots=network_plots);
  
  return(output);
  
}

merged_communitySubtypeBreakdown_plots <- function(mergedMatrixData,experimentName,
                                                   saveDir,centroidRData="/home/data/breast_microarrayDB/pam50_centroids_updatedSymbols.RData",
                                                   clusterMatrixList){

   orig_fullSubtypeMatrix <- assignCentroidSubtype(origDataMatrix=t(mergedMatrixData$mergedExprMatrix),minNumGenes=30,centroidRData=centroidRData)
   
   minClustSize <- 10;
   subtypeDF <- data.frame();
   clusterCount <- 0;
  
   minNumGenes <- 15;
   #there are only 10 intrinsic genes in here! so just do intrinsic...
   for(c in 1:length(clusterMatrixList[[1]])){
     
     if(nrow(clusterMatrixList[[1]][[c]])>=minClustSize){
       clusterCount <- clusterCount+1;
       patIDs <- rownames(clusterMatrixList[[1]][[c]]);
       clustID <- rep(c,times=length(patIDs));
       studyID <- mergedMatrixData$study[match(patIDs,mergedMatrixData$GSMID)]
       subtypes <- assignCentroidSubtype(origDataMatrix=clusterMatrixList[[1]][[c]],minNumGenes=minNumGenes,centroidRData=centroidRData);
       
       if(clusterCount==1){
         
         subtypeDF <- data.frame(patIDs,clustID,studyID,subtypes$subtypes[,2] );
         
       }else{
         
         tmp <- data.frame(patIDs,clustID,studyID ,subtypes$subtypes[,2] );
         
         subtypeDF <- rbind(tmp,subtypeDF);
         
       }
       
     }
     
   }
   subtypeDF$subtypes <- factor(subtypeDF$subtypes)
   subtypeDF$studyFactor <- as.numeric(factor(subtypeDF$studyID));
   
   patientData <- data.frame(mergedMatrixData$study,factor(as.numeric(factor(mergedMatrixData$study))),mergedMatrixData$GSMID);
   colnames(patientData) <- c("fullStudyName","studyID","GSMID"); 

   colnames(subtypeDF)[5] <- "subtype";
   colnames(subtypeDF)[2] <- "community";
   saveDir <- "/home/kplaney/ISMB/FINAL_figs";
   experimentName <- "Kmeans BMC merged";
   centroidSetName <- "";
   #do a facetted plot by all communities too
   subtypeDF$community <- paste0("metacluster_",subtypeDF$community);
   
   #specify colors:
   subtypeColorMatrix <- data.frame(c("#970eec","#ec0e58","#0e58ec","#7d7d7d","#0eec97"))
   #hmm can I order the colors in the same way each time?? may not really need to anyways.
   #must tranpose for ggplot to read it correctly and THEN add names for this to work.
   #subtypeColorMatrix <- t(subtypeColorMatrix);
   rownames(subtypeColorMatrix) <- c("LumB","LumA","Her2","Basal","Normal");
   subtypeColorMatrix <- t(subtypeColorMatrix);
   #subtypeColorMatrix <- data.frame(c("#970eec","#ec0e58","#0e58ec","#7d7d7d","#0eec97"))
   #hmm can I order the colors in the same way each time?? may not really need to anyways.
   #must tranpose for ggplot to read it correctly and THEN add names for this to work.
   #subtypeColorMatrix <- t(subtypeColorMatrix);
   #names(subtypeColorMatrix) <- c("LumB","LumA","Her2","Basal","Normal");

   #stacked bar for all studies.
   #scales="free_x":  doesn't plot studies that are empty for that community.
   subtype_plot_stackedALL <- ggplot(data=subtypeDF,aes(x=community))+geom_bar(aes(fill=factor(subtype)))+
     scale_fill_manual(values = subtypeColorMatrix)+
     labs(y="Number of samples",x=element_blank(),title=paste0(centroidSetName,""))+theme_bw()+
     labs(fill="subtype")+theme(strip.text.x = element_text(size = 17, colour = "black",face="bold"))+
     theme(legend.title=element_text(colour="black",size=17),legend.text=element_text(colour="black",size=15))+
     theme(axis.text.x = element_text(colour = "black",size=8,angle=45,vjust=1,hjust=1),axis.title.x= element_text(colour = "black",size=20,vjust=1),
           axis.text.y = element_text(colour = "black",size=18),axis.title.y = element_text(colour = "black",size=20,vjust=1))+
     theme(plot.title=element_text(colour="black",size=28,vjust=1,hjust=.2))+theme(panel.background = element_rect(colour = "black"))+
     theme(panel.grid.major = element_line(colour = 0),
           panel.grid.minor = element_line(colour = 0));
   
   
   #make a bit wider..
   png(filename=paste0(saveDir,"/",experimentName,"_",Sys.Date(),".png"),
       width = 800, height = 1000,res=160);
   
   plot(subtype_plot_stackedALL);
   dev.off();
   
   mergedMatrixData <- list(subtypeDF=subtypeDF,subtype_plot_stackedALL=subtype_plot_stackedALL,patientData= patientData);
   
}

meanMetricDensityPlot <- function(meanMetricMatrix,saveDir="./",savePlot=TRUE,experimentName="expr",
                                  yLimit){

  if(missing(yLimit)){
    
    dataF <- as.vector(meanMetricMatrix)
    masterDF <- dataF[-which(is.na(dataF))]
    dens <- density(masterDF)
    masterDF <- as.data.frame(masterDF,stringsAsFactors=FALSE)
    densityPlot <- ggplot(data = masterDF, aes(x=masterDF))+
      #alpha in geom_density. controls transparency of color, if used fill= above.
      #if want transparent fill: geom_density(alpha=0). but different lines
      #are better for a black and white publication.
      geom_density() +labs(title = "Density curves ")+
      labs(title = "Density curve for \nCoINcIDE edge mean metric",
           y="Density",x="Mean metric")+
      
      theme(panel.background = element_rect(fill='white', colour='black')) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  }else{
    
    dataF <- as.vector(meanMetricMatrix)
    masterDF <- dataF[-which(is.na(dataF))]
    dens <- density(masterDF)
    masterDF <- as.data.frame(masterDF,stringsAsFactors=FALSE)
    densityPlot <- ggplot(data = masterDF, aes(x=masterDF))+
      #alpha in geom_density. controls transparency of color, if used fill= above.
      #if want transparent fill: geom_density(alpha=0). but different lines
      #are better for a black and white publication.
      geom_density() +labs(title = "Density curves ")+
      labs(title = "Density curve for \nCoINcIDE edge mean metric",
           y="Density",x="Mean metric")+
      theme(panel.background = element_rect(fill='white', colour='black')) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
       ylim(c(0, yLimit))
  }
  
  if(savePlot){
    
    options(bitmapType="cairo")
    png(filename=paste0(saveDir,"/",experimentName,"_densityPlot_",Sys.Date(),".png"),
        width = 700, height = 1000,res=160);
    plot(densityPlot)
    dev.off()
    
  }
  
  meanMetricAtMaxDensity <- dens$x[which(dens$y==max(dens$y))]
  message("The meanMetricAtMaxDensity is a global maximum; view the density plot to estimate a value at a local maximum that may have a higher mean metric.")
  output <- list(meanMetricAtMaxDensity=meanMetricAtMaxDensity,densityPlot=densityPlot)
  return(output)
  
}
