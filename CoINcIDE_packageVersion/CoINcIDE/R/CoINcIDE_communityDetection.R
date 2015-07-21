
#library("igraph");
#library("limma");
#library("RColorBrewer")
#######ANALYZE OUTPUT FUNCTIONS
assignFinalEdges <- function(computeTrueSimilOutput,pvalueMatrix,indEdgePvalueThresh=.1,
                             meanEdgePairPvalueThresh=.05,
                             minTrueSimilThresh=-Inf,maxTrueSimilThresh=Inf,
                             fractFeatIntersectThresh=0,numFeatIntersectThresh=0 ,
                             clustSizeThresh=0, clustSizeFractThresh=0,saveDir="./",fileTag="CoINcIDE_edges",
                             clustIndexMatrix,minFractNN =.7,minNumEdgesForCluster=1
                             ){
  
  clustSizeIndexRemove <- c()
  count <- 0
  meanEdgePvalueMatrix <- matrix(data=NA,ncol=ncol(pvalueMatrix),nrow=nrow(pvalueMatrix))
  meanFractNNmatrix <- matrix(data=NA,ncol=ncol(pvalueMatrix),nrow=nrow(pvalueMatrix))
  meanMeanMetricMatrix <- matrix(data=NA,ncol=ncol(pvalueMatrix),nrow=nrow(pvalueMatrix))
  
  for(r in 1:nrow(computeTrueSimilOutput$clustSizeMatrix)){
    
    if(computeTrueSimilOutput$clustSizeFractMatrix[r]==1){
      #ROW is when this cluster was the reference cluster. p-value will by default be insignificant.
      #just replace with p-value in the other direction (columnwise)
      #(if both clusters in a pair were entire size of matrix, will have NA p-value still.)
      #which(computeTrueSimilOutput$clustSizeFractMatrix!=1): do NOT want to compare these clusters if both have only 1 cluster in entire dataset.
      #too biased.
      pvalueMatrix[r, which(computeTrueSimilOutput$clustSizeFractMatrix != 1)] <- pvalueMatrix[which(computeTrueSimilOutput$clustSizeFractMatrix !=1), r]
      #set these cases to NA pvalue so that they are not compared
      pvalueMatrix[r, which(computeTrueSimilOutput$clustSizeFractMatrix == 1)] <- NA
      #must do same for NN fraction - otherwise too biased bc other way will be 1!
      computeTrueSimilOutput$trueFractNNmatrix[r,which(computeTrueSimilOutput$clustSizeFractMatrix!=1)] <- computeTrueSimilOutput$trueFractNNmatrix[which(computeTrueSimilOutput$clustSizeFractMatrix!=1) ,r]
      #set to NA if both clusters entire size of dataset.
      computeTrueSimilOutput$trueFractNNmatrix[r, which(computeTrueSimilOutput$clustSizeFractMatrix == 1)] <- NA
      
    }
    
    if(computeTrueSimilOutput$clustSizeMatrix[r] <= clustSizeThresh || computeTrueSimilOutput$clustSizeFractMatrix[r] <= clustSizeFractThresh){
      #remove these clusters - do not pass thresholds
      pvalueMatrix[r, ] <- NA
      count <- count + 1
      clustSizeIndexRemove[count] <- r
    }
    
  }
  
  message(length(clustSizeIndexRemove), " clusters dropped because they were below the clust size thresh clustSizeThresh
          threshold of ",clustSizeThresh)
  #assign mean p-value matrix

  for(r in 1:nrow(meanEdgePvalueMatrix)){
    
    for(c in 1:nrow(meanEdgePvalueMatrix)){
      #have we already computed the mean here?
      if(is.na(meanEdgePvalueMatrix[r,c])){
        
        #other choice: could combine via fisher's p-value, stouffer's (Z-transform) or brown's method.
        #for independent p-values: some claim Stouffer's is preferred over fisher's method: http://onlinelibrary.wiley.com/doi/10.1111/j.1420-9101.2005.00917.x/full
        #Brown's method, not stouffer or fisher, are correct for dependent p-values (which is the case here.)
        #But Brown's only defines the variance/standard devation and mean - i.e. not a readily interpretable p-value.
        #this sum method is simpler but much more conservative than fisher's p-value, and 
        #presumably more conservative than stouffer's or brown's method.
        #don't remove NA - set to NA if one value in pair is NA.
        meanEdgePvalueMatrix[r,c] <- mean(c(pvalueMatrix[r,c],pvalueMatrix[c,r]),na.rm = FALSE)
        meanEdgePvalueMatrix[c,r] <- meanEdgePvalueMatrix[r,c]
        
        meanFractNNmatrix[r,c] <- mean(c(computeTrueSimilOutput$trueFractNNmatrix[r,c], computeTrueSimilOutput$trueFractNNmatrix[c,r]),na.rm = FALSE)
        meanFractNNmatrix[c,r] <- meanFractNNmatrix[r,c]
        meanMeanMetricMatrix[r,c] <- mean(c(computeTrueSimilOutput$meanMetricMatrix[r,c], computeTrueSimilOutput$meanMetricMatrix[c,r]),na.rm = FALSE)
        meanMeanMetricMatrix[c,r] <- meanMeanMetricMatrix[r,c]
        
        
      }
      
    }
    
  }

  adjMatricesList <- list(meanMeanMetricMatrix,meanMeanMetricMatrix,
                          pvalueMatrix,meanEdgePvalueMatrix,computeTrueSimilOutput$fractFeatIntersectMatrix,
                          computeTrueSimilOutput$numFeatIntersectMatrix ,meanFractNNmatrix=meanFractNNmatrix)
  
  names(adjMatricesList) <- c("simil","similDup","indPvalue","meanPvalue","fractFeatIntersect","numFeatIntersect","minFractNN")
  
  thresholdVector <- c(minTrueSimilThresh,maxTrueSimilThresh,indEdgePvalueThresh,meanEdgePairPvalueThresh,
                       fractFeatIntersectThresh,numFeatIntersectThresh,minFractNN)
  
  threshDir <- c(">=","<=","<=","<=",">=",">=",">=")
  
  if(length(adjMatricesList) != length(thresholdVector) || length(adjMatricesList) != length(threshDir)){
    
    stop("Bug in setting up thresholding vectors.")
    
  }
  
  filterEdgeOutput <- filterEdges(adjMatricesList,thresholdVector,threshDir=threshDir,saveDir=saveDir,fileTag=fileTag,saveEdgeFile=TRUE)
  
  #now filter one step further: nodes that have below minimum number of edges
  clustFurtherRemoved <- setdiff(filterEdgeOutput$clustRemoved,clustSizeIndexRemove)
  
  allClustRemoved <- filterEdgeOutput$clustRemoved
  
  message("A total of ",length(allClustRemoved), " clusters removed because they have no edges that pass the thresholds.")

  allRemovedMinEdgeClusters <- c()
  
  if(minNumEdgesForCluster >1){
    
  while(nrow(filterEdgeOutput$edgeMatrix)>0){
    #keep looping: if remove a certain cluster, may leave another cluster with few edges.
    totalNumEdges <- c()
      
    clustNames <- unique(as.vector(filterEdgeOutput$edgeMatrix))
    
      for(m in 1:length(clustNames)){
  
        clustNum <-  clustNames[m]
        #take all edges
        edgeMatrixIndices <- union(which(filterEdgeOutput$edgeMatrix[,1]==clustNum), which(filterEdgeOutput$edgeMatrix[,2]==clustNum))
        totalNumEdges[m] <- length(edgeMatrixIndices)
        
      }
    
    names(totalNumEdges) <- clustNames
    
    #unweighted - care about how many other clusters this cluster was connected to.
    removeClusters <- names(which(totalNumEdges<=minNumEdgesForCluster))
    
    if(length(removeClusters)==0){
      #no more clusters to remove
      break;
      
    }
    
    allRemovedMinEdgeClusters <- append(allRemovedMinEdgeClusters,removeClusters)
    
    for(r in 1:length(removeClusters)){
      
      removeRowIndices <- union(which(filterEdgeOutput$edgeMatrix[,1]==removeClusters[r]), which(filterEdgeOutput$edgeMatrix[,2]==removeClusters[r]))
      
      filterEdgeOutput$edgeMatrix <- filterEdgeOutput$edgeMatrix[-removeRowIndices, ]
      filterEdgeOutput$edgeWeightMatrix <- filterEdgeOutput$edgeWeightMatrix[-removeRowIndices, ]
      
    }
    
  }
  
  }else{
  
    #keep looping: if remove a certain cluster, may leave another cluster with few edges.
    totalNumEdges <- c()
    
    clustNames <- unique(as.vector(filterEdgeOutput$edgeMatrix))
    
    for(m in 1:length(clustNames)){
      
      clustNum <-  clustNames[m]
      #take all edges
      edgeMatrixIndices <- union(which(filterEdgeOutput$edgeMatrix[,1]==clustNum), which(filterEdgeOutput$edgeMatrix[,2]==clustNum))
      totalNumEdges[m] <- length(edgeMatrixIndices)
      
    }
    
    names(totalNumEdges) <- clustNames
  }
  message("A total of ",length(allRemovedMinEdgeClusters), " clusters removed because they have too few edges.")
  
  allClustersRemoved <- append(allClustRemoved,allRemovedMinEdgeClusters)
  
  if(length(allClustersRemoved) != length(unique(allClustersRemoved))){
    
    stop("Error: double-counting edges (or recounting) clusters after they are removed.")
    
  }
  
  message("A total of ",length(union(unique(filterEdgeOutput$edgeMatrix[,2]),unique(filterEdgeOutput$edgeMatrix[,1]))), " clusters have edges that pass the thresholds.")
  
  output <- list(meanEdgePvalueMatrix=meanEdgePvalueMatrix,clustSizeIndexRemove=clustSizeIndexRemove,clustFurtherRemoved=clustFurtherRemoved,filterEdgeOutput=filterEdgeOutput,
                 clustRemovedBeforeMinNumEdges=allClustRemoved,adjMatricesList=adjMatricesList,meanFractNNmatrix=meanFractNNmatrix,meanMeanMetricMatrix=meanMeanMetricMatrix,
                 totalNumEdgesForCluster=totalNumEdges,minNumEdgesRemove=allRemovedMinEdgeClusters)
  return(output)
  
}

#create wrapper function for filterEdges (this is what user sees.)
filterEdges <- function(adjMatricesList,thresholdVector,threshDir=rep(">=",length(thresholdVector)),saveDir="./",fileTag="",saveEdgeFile=TRUE){
  
  #creates directory if there isn't one
  dir.create(saveDir,showWarnings = FALSE);

  numEdges <- 0;
  #create a "master" adjacency matrix.
  #if no edge: is a zero.
  adjMatrix <- matrix(data=NA,ncol=ncol(adjMatricesList[[1]]),nrow=nrow(adjMatricesList[[1]]));
  edgeMatrix <- matrix(data=NA,ncol=2);
  edgeWeightMatrix <- matrix(data=NA,ncol=length(adjMatricesList))
  
  if(all(!is.na(names(adjMatricesList)))){
    
    colnames(edgeWeightMatrix) <- names(adjMatricesList)
  
  }
  
  for(n in 1:nrow(adjMatricesList[[1]])){
    
    for(p in 1:ncol(adjMatricesList[[1]])){
      
      pass <- TRUE;
      #have we looked at this pair yet?
      if(is.na(adjMatrix[n,p]) && is.na(adjMatrix[p,n])){
        
        for(a in 1:length(adjMatricesList)){
          
          #must have non-NA adjacency valUes
          if(!is.na(adjMatricesList[[a]][p,n]) && !is.na(adjMatricesList[[a]][n,p])){
            #is it below threshold? must pass ALL thresholds.
            if(threshDir[a]==">="){
              
              if(adjMatricesList[[a]][p,n] >= thresholdVector[a] && adjMatricesList[[a]][n,p] >= thresholdVector[a]){
                
                pass <- TRUE;
                
              }else{
                pass <- FALSE;
                #break out of loop.
                break;
              }
              
            }else  if(threshDir[a]=="<="){
              
              if(adjMatricesList[[a]][p,n] <= thresholdVector[a] && adjMatricesList[[a]][n,p] <= thresholdVector[a]){
                
                pass <- TRUE;
                
              }else{
                
                #if one false pass: jump out of loop.
                pass <- FALSE;
                #break out of loop.
                break;
                
              }             
              
            }else{
              
              stop("\nPlease pick >= or <= for thresh directions.")
            }
            
          }else{
            #failed: had at least one NA adjacency value.
            pass <- FALSE;
            #break out of loop.
            break;
            
          }
          
          #end of loop a
        }
        
        if(pass){
          
          adjMatrix[n,p] <- 1;
          adjMatrix[p,n] <- 1;
          numEdges <- numEdges + 1;           
          edgeMatrix <- rbind(edgeMatrix, t(as.matrix(c(n,p))));
          
          temp <- c();
          #take average of each bidirectional weight to get final weight.
          for(a in 1:length(adjMatricesList)){
            
            temp[a] <- mean(adjMatricesList[[a]][p,n],adjMatricesList[[a]][n,p]);
            
            if(is.na(temp[a])){
              
              stop("Getting NA edge weights!");
              
            }
            
          }
          
          edgeWeightMatrix <- rbind(edgeWeightMatrix,t(temp));
          
        }else{
          
          adjMatrix[n,p] <- 0;
          adjMatrix[p,n] <- 0;
          
        }
        
        #end of if is not NA/already computed statement     
      }
      
    }
    #end of loops
  }
  
  #identify clusters who adjMatrix are ALL zeros (no edges.) return this info.
  clustRemoved <- which(rowSums(adjMatrix)==0)
  #remove first NA row
  edgeMatrix <- edgeMatrix[-1, ,drop=FALSE];
  edgeWeightMatrix <- edgeWeightMatrix[-1, ,drop=FALSE];
  
  if(saveEdgeFile){
  #can we attach intuitive column names to make this more readable?
  if(is.null(names(adjMatricesList)[1])){
    
    colnames(edgeWeightMatrix) <- names(adjMatricesList);
    write.table(edgeWeightMatrix,file=paste0(saveDir,"/IGPN_edgeWeights_",fileTag, "_",Sys.Date(),".txt"),sep="\t",quote=FALSE,row.names=FALSE,
                col.names=TRUE);
  }else{
    
    colnames(edgeWeightMatrix) <- names(adjMatricesList);
    write.table(edgeWeightMatrix,file=paste0(saveDir,"/edgeWeights_",fileTag, "_",Sys.Date(),".txt"),sep="\t",quote=FALSE,row.names=FALSE,
                col.names=FALSE);
    
  }

  #create edge and weight files for cytoscape or other graphing programs
  write.table(edgeMatrix,file=paste0(saveDir,"/edgeAssignments_",fileTag, "_",Sys.Date(),".txt"),sep="\t",quote=FALSE,row.names=FALSE,
              col.names=FALSE);
  
  }
  #create a weight column from each adjMatrix (use names() as column names.)
  output <- list(edgeMatrix=edgeMatrix,edgeWeightMatrix=edgeWeightMatrix,adjMatrix=adjMatrix,clustRemoved=clustRemoved);
  return(output);
  
}



###########################
#######
# IGP_communityDetection <- library("igraph")
# undirGraph <- graph.edgeMatrix(finalEdgeMatrix,directed=FALSE);
# plot(undirGraph)

findCommunities <- function(edgeMatrix,edgeWeightMatrix,clustIndexMatrix,fileTag="CoINcIDE_communityNodeAttributes_",
                            saveDir="./",minNumUniqueStudiesPerCommunity=3,experimentName="exp",
                            commMethod=c("edgeBetween","fastGreedy","walktrap","eigenvector","optimal","spinglass","multilevel"),
                            makePlots=TRUE,plotToScreen=FALSE,saveGraphData=TRUE,nodeFontSize=.7,nodePlotSize=10,
                            findCommWithWeights=FALSE, plotSimilEdgeWeight = TRUE,fractEdgesInVsOutComm=0,
                            fractEdgesInVsOutEdge=0){
  
  
  if(nrow(edgeMatrix)>0){
    
  if(findCommWithWeights || plotSimilEdgeWeight){
    #advanced users can change this keyword if they have another data frame column they want to use.
    edgeWeightsColName <- "simil"
    
    
  }else{
    
     edgeWeightsColName <- NULL
  }
  
  if(length(commMethod)>1){
    
    message("commMethod input variable longer than length 1; using default \'edgeBetween'\ method.")
    commMethod <- "edgeBetween"
  
  }
  dir.create(saveDir)
  #dir.create(paste0(saveDir,"/",experimentName,"_",Sys.Date()),showWarnings=TRUE);
  #saveDir <- paste0(saveDir,"/",experimentName,"_",Sys.Date())

  #  undirGraph <- graph.edgeMatrix(edgeMatrix,directed=FALSE);
  #igraph likes an edge list that is continous - i.e. if there's not an edge from node 10, it will still include node 10.
  # we don't wan't this! collapse down the node indices.
  graphKey <- cbind(unique(append(edgeMatrix[,1],edgeMatrix[,2])),
                    c(1:length(unique(append(edgeMatrix[,1],edgeMatrix[,2])))));
  
  colnames(graphKey) <- c("orig_id","new_id");
  igraph_edgeMatrix <- matrix(data=NA,nrow=nrow(edgeMatrix),ncol=ncol(edgeMatrix));
  igraph_weightMatrix <- matrix(data=NA,nrow=nrow(edgeWeightMatrix),ncol=ncol(edgeWeightMatrix));
  
  for(r in 1:nrow(graphKey)){
    #replace any id, whether it's in the first or second column
    igraph_edgeMatrix[which(edgeMatrix==graphKey[r,"orig_id"])] <- graphKey[r,"new_id"];
    
  }
  
  #add on weights: Additional columns in edgeMatrix are considered as edge attributes. 
  igraph_edgeMatrix <- data.frame(igraph_edgeMatrix, edgeWeightMatrix)

    
  clustNames <- union(unique(edgeMatrix[,1]),unique(edgeMatrix[,2]));
  studyNames <-  clustIndexMatrix[na.omit(match(clustNames,clustIndexMatrix[,1])) ,2];
  #order in same order as orig names.
  igraph_clustNames <- graphKey[match(clustNames,graphKey[,"orig_id"]),"new_id"];

  base_attr <- cbind(igraph_clustNames,studyNames);
  colnames(base_attr)[2] <- "studyNum";
  undirGraph_base <-  graph.data.frame(data.frame(igraph_edgeMatrix),directed=FALSE,vertices=data.frame(base_attr));

  #scaling factor
  if(!is.null(edgeWeightsColName)){
    
    #multiple times ten- these original values will all be between -1 and 1
    
    if(max(get.edge.attribute(graph=undirGraph_base,name=edgeWeightsColName))<=1){
      
      #don't scale up here just to make sure won't affect community detection.
      plotEdgeWeights <- get.edge.attribute(graph=undirGraph_base,name=edgeWeightsColName)
      
    }else{
      
      stop("You selected to plot edge weights but your similarity (or distance) metric has large values above 1. This will not plot well.")
      
    }
    
  
  if(length( plotEdgeWeights)==0){
    
    stop("\nEdge weights returned was a null lust.")
  
  }
  
  if(any(is.na(plotEdgeWeights))){
    
    stop("\nGetting NA edge weights.")
    
  }
  
  E(undirGraph_base)$weights <- plotEdgeWeights
  
  }
 
  if(makePlots){
    
    if(!plotToScreen){
      
    
    if(!is.null(edgeWeightsColName)){
      
    png(filename=paste0(saveDir,"/",experimentName,"_origNetworkPlot_edgesScaled",Sys.Date(),".png"),
      width = 700, height = 1000,res=160)
    
    origNetworkPlot <- plot(undirGraph_base, layout=layout.fruchterman.reingold,vertex.size=nodePlotSize,vertex.label=V(undirGraph_base)$studyNum,vertex.label.color="black",vertex.label.cex=nodeFontSize,
       edge.arrow.size=3,edge.width = E(undirGraph_base)$weights,edge.color="black", main=paste0(length(union(unique(edgeMatrix[,1]),unique(edgeMatrix[,2]))), " clusters deemed similar\n across ",length(unique(studyNames))," studies before community detection."
                                                                            ),xlab="node=cluster,#=dataset, edge weight = meanMetric");
   
    dev.off();
    
    }else{
      
      png(filename=paste0(saveDir,"/",experimentName,"_origNetworkPlot",Sys.Date(),".png"),
          width = 700, height = 1000)
      
      origNetworkPlot <- plot(undirGraph_base, layout=layout.fruchterman.reingold,vertex.size=nodePlotSize,vertex.label=V(undirGraph_base)$studyNum,vertex.label.color="black",vertex.label.cex=nodeFontSize,
                              edge.arrow.size=3,main=paste0(length(union(unique(edgeMatrix[,1]),unique(edgeMatrix[,2]))), " clusters deemed similar\n across ",length(unique(studyNames))," studies before community detection."
                              ),xlab="node=cluster,#=dataset");
      
      
      dev.off();
    }
    
    }else{
      
    if(!is.null(edgeWeightsColName)){
 
         origNetworkPlot <- plot(undirGraph_base, layout=layout.fruchterman.reingold,vertex.size=nodePlotSize,vertex.label=V(undirGraph_base)$studyNum,vertex.label.color="black",vertex.label.cex=nodeFontSize,
                                 edge.arrow.size=3,edge.width = plotEdgeWeights,edge.color="black", main=paste0(length(union(unique(edgeMatrix[,1]),unique(edgeMatrix[,2]))), " clusters deemed similar\n across ",length(unique(studyNames))," studies before community detection."
                                 ),xlab="#=Study Number, edge weight = meanMetric");
         
 
    }else{
         
      origNetworkPlot <- plot(undirGraph_base, layout=layout.fruchterman.reingold,vertex.size=nodePlotSize,vertex.label=V(undirGraph_base)$studyNum,vertex.label.color="black",vertex.label.cex=nodeFontSize,
                              edge.arrow.size=3,main=paste0(length(union(unique(edgeMatrix[,1]),unique(edgeMatrix[,2]))), " clusters deemed similar\n across ",length(unique(studyNames))," studies before community detection."
                              ),xlab="#=Study Number");
      
       }
       
       
       
    }
    
  }
  #commMethod=c("fastGreedy","edgeBetween","walktrap","eigenvector","optimal","spinglass","multilevel")

  if(!findCommWithWeights){
  
  if(commMethod=="fastgreedy"){
    
  comm <- fastgreedy.community(undirGraph_base,merges=TRUE,modularity=TRUE,
                                     membership=TRUE,weights=NULL);

  }else if(commMethod=="edgeBetween"){
    
    #note: can actually use weights here.
    #IF wanted weights: put it here, but only finding communities by weights does not always result in clear clusters.
    #weights=E(undirGraph_base)$weights
    comm <- edge.betweenness.community(undirGraph_base,directed=FALSE,merges=TRUE,modularity=TRUE,
                                 membership=TRUE,weights=NULL);
    
    
  }else if(commMethod=="walktrap"){
    #steps=4 is default
    comm <- walktrap.community(undirGraph_base,merges=TRUE,modularity=TRUE,
                               membership=TRUE,steps=4,weights=NULL);
    
  }else if(commMethod=="eigenvector"){
    
    
    comm <- leading.eigenvector.community(undirGraph_base,steps = -1,weights=NULL);
    
  }else if(commMethod=="optimal"){
    
    #maximizing the modularity measure over all possible partitions
    comm <- optimal.community(undirGraph_base,weights=NULL);
    
  }else if(commMethod=="spinglass"){
    
    comm <- spinglass.community(undirGraph_base,weights=NULL, vertex=NULL, spins=25,
                                parupdate=FALSE, start.temp=1, stop.temp=0.01,
                                cool.fact=0.99, update.rule=c("config"), gamma=1, implementation=c("orig"),
                                gamma.minus=1);
    
  }else if(commMethod=="multilevel"){
    
    comm <- multilevel.community(undirGraph_base,weights=NULL);
    
  }else{
    
    stop("\nPlease pick a viable community detection method for the 'commMethod' variable.")
  }
    
  #use weights, when can
  }else{
    
      
  if(commMethod=="fastgreedy"){
    
  comm <- fastgreedy.community(undirGraph_base,merges=TRUE,modularity=TRUE,
                                     membership=TRUE,weights=E(undirGraph_base)$weights);

  }else if(commMethod=="edgeBetween"){
    
    #note: can actually use weights here.
  #weights=E(undirGraph_base)$weights
    comm <- edge.betweenness.community(undirGraph_base,directed=FALSE,merges=TRUE,modularity=TRUE,
                                 membership=TRUE,weights=E(undirGraph_base)$weights);
    
    
  }else if(commMethod=="walktrap"){
    #steps=4 is default
    comm <- walktrap.community(undirGraph_base,merges=TRUE,modularity=TRUE,
                               membership=TRUE,steps=4,weights=E(undirGraph_base)$weights);
    
  }else if(commMethod=="eigenvector"){
    
    
    comm <- leading.eigenvector.community(undirGraph_base,steps = -1,weights=E(undirGraph_base)$weights);
    
  }else if(commMethod=="optimal"){
    
    #maximizing the modularity measure over all possible partitions
    comm <- optimal.community(undirGraph_base,weights=E(undirGraph_base)$weights);
    
  }else if(commMethod=="spinglass"){
    
    comm <- spinglass.community(undirGraph_base,weights=E(undirGraph_base)$weights, vertex=NULL, spins=25,
                                parupdate=FALSE, start.temp=1, stop.temp=0.01,
                                cool.fact=0.99, update.rule=c("config"), gamma=1, implementation=c("orig"),
                                gamma.minus=1);
    
  }else if(commMethod=="multilevel"){
    
    comm <- multilevel.community(undirGraph_base,weights=E(undirGraph_base)$weights);
    
  }else{
    
    stop("\nPlease pick a viable community detection method for the 'commMethod' variable.")
  }
    
    
    
  }
  #you can actually plot this later.
  unPrunedDendogram <- as.dendrogram(comm)
  modularity <- modularity(comm);
  #add back the old membership.
  membership <- cbind(graphKey[match(c(1:length(comm$membership)),graphKey[,"new_id"]),"orig_id"],comm$membership);
  colnames(membership) <- c("clust","community");
  numCommunitiesOrig <- length(unique(membership[,"community"]));
  membership[,"community"] <- as.factor(membership[,"community"]);
  membershipList <- split(membership[,"clust"],f=membership[,"community"]);
  
  if(makePlots){
    #make colors for the nodes
    #http://www.w3schools.com/html/html_colornames.asp
    # "blue", etc also works
    #turquoise/purple "#00FFFF","#FF1493"
    #this just creates a rainbow.
    #add x2 in the bias to make the colors a wee bit more distinct from each other.
    #bias: a positive number. Higher values give more widely spaced colors at the high end.
    colorCodeF <- colorRampPalette(c("pink","deepskyblue","yellow"), bias = length(unique(membership[,"community"]))*2,
                                   space = "rgb", interpolate = "linear");
    
    #this will produce RGB representation in HEX, which cytoscape can take in.
    #if need plain old RGB: col2RGB(membership[,"colors"],alpha=FALSE);
    colorCodes <- colorCodeF(length(unique(membership[,"community"])));
  
  }else{
    
    colorCodes <- rep.int(NA,times=length(unique(membership[,"community"])))
    
  }
  
  #let RcolorBrewer do the color mixing for you!
  #uggh..max is 11 colors!
  #colorCodes <- rev(brewer.pal(length(unique(membership[,"community"])),"Spectral"));
  
  membership <- cbind(membership,matrix(data=NA,nrow=nrow(membership),ncol=2));
  colnames(membership) <- c("clust","community","color","igraph_id");                   
  
  for(m in 1:length(unique(membership[,"community"]))){
    
    membership[which(membership[,"community"]==unique(membership[,"community"])[m]),"color"] <- colorCodes[m];
    
  }
  
  #add back the graphKey
  membership[ ,"igraph_id"] <-  graphKey[na.omit(match(membership[ ,"clust"],graphKey[,"orig_id"])), "new_id"];
  #put this one in front: it's the id for igraph
  membership <- membership[ ,c(4,1:3)];
  
  #add on study num
  #assumption: study number is in the second location. cluster number is in the first.:  x_x_..
  membership <- cbind(membership,clustIndexMatrix[na.omit(match(membership[,"clust"],clustIndexMatrix[,1])) ,2])
  colnames(membership)[ncol(membership)] <- "studyNum";

  #now see if this actually meets the thresholds
  commKeep <- c();
  k <- 0;
  numEdgesInComm <- c()
  totalNumEdges <- c()
  
  for(m in 1:nrow(membership)){
    
    commNum <- membership[m,"community"]
    clustNum <-  membership[m,"clust"]
    edgeMatrixIndices <- union(which(edgeMatrix[,1]==clustNum), which(edgeMatrix[,2]==clustNum))
    edgeMatrixShort <- edgeMatrix[edgeMatrixIndices, ]
    edgeMatrixVector <- as.vector(edgeMatrixShort)
    #remove nodes that are itself
    edgeMatrixVector <- edgeMatrixVector[-which(edgeMatrixVector==clustNum)]
    
    if(length(edgeMatrixVector) != length(unique(edgeMatrixVector))){
    
      stop("You're getting duplicated edges.")
      
    }
    
    totalNumEdges[m] <- length(unique(edgeMatrixVector))
    membershipShort <- membership[na.omit(match(edgeMatrixVector,membership[,"clust"])), ,drop=FALSE]
    numEdgesInComm[m] <- length(which(membershipShort[,"community"]==commNum))

  }
  
  fractEdgesInOut <- numEdgesInComm/totalNumEdges
  numEdgesNotInComm <- totalNumEdges - numEdgesInComm
  membership <- data.frame(membership, cbind(totalNumEdges,numEdgesInComm,fractEdgesInOut,numEdgesNotInComm))

  
  removeNodes <- which(fractEdgesInOut <= fractEdgesInVsOutEdge)
  
  if(length(removeNodes)>0){

    membership <- membership[-removeNodes, ]
    removeEdgesRowIndices <- c() 
    for(r in 1:length(removeNodes)){
     
    removeEdgesRowIndices <- append(removeEdgesRowIndices,union(which(igraph_edgeMatrix[,1]==removeNodes[r]), which(igraph_edgeMatrix[,2]==removeNodes[r])))
  
    }
    
    #these all share the same row indices.
    edgeMatrix <- edgeMatrix[-removeEdgesRowIndices, ]
    edgeWeightMatrix <- edgeWeightMatrix[-removeEdgesRowIndices, ]
    igraph_edgeMatrix <- igraph_edgeMatrix[-removeEdgesRowIndices, ]
    
  }
  #re-name nodes
  
  message(length(removeNodes), " nodes removed because they have a fraction of edges within a community over total number of edges below the threshold of ",fractEdgesInVsOutEdge)
  
  commEdgeInfo <- data.frame()
  k <- 0

  numMembersInComm <- table(membership[,"community"]);
  
  #want to remove the factors - otherwise can mess up final commKeep vector.
  commNum <- as.vector(unique(membership[,"community"]))
  if(!all(commNum == unique(membership[,"community"]))){
    
    stop("Potential bug in removing factor levels.")
    
  }
  for(c in 1:numCommunitiesOrig){
    
    #these will be double-counted as both nodes for each edge will be in the community
    numTotalEdgesInCommunity <- sum(membership[which(membership[,"community"]==commNum[c]),"numEdgesInComm"])/2
    #these will not be double-counted.
    numTotalEdgesNotInCommunity <- sum(membership[which(membership[,"community"]==commNum[c]),"numEdgesNotInComm"])
    
    fractEdgesInVsOut <- numTotalEdgesInCommunity/(numTotalEdgesNotInCommunity+numTotalEdgesInCommunity)
    
    numDatasetsPerCommunity <- length(unique(membership[which(membership[,"community"]==commNum[c]), "studyNum"]))
    #message(numDatasetsPerCommunity)
    tmp <- data.frame(numTotalEdgesInCommunity,numTotalEdgesNotInCommunity, fractEdgesInVsOut, numDatasetsPerCommunity,stringsAsFactors=FALSE)
    
    if( (numDatasetsPerCommunity>=minNumUniqueStudiesPerCommunity) && (fractEdgesInVsOut >= fractEdgesInVsOutComm) ){
        #message(unique(membership[,"community"])[c])
        #message(commNum[c])
        k <- 1 + k;
        commKeep[k] <- commNum[c];
        commEdgeInfo <- rbind(commEdgeInfo,tmp)
        rownames(commEdgeInfo)[k] <- k
      
    }
    
  }
  
  #then: bind together edge matrix and edge weight matrix (if certain nodes were not in a community) for igraph
  igraph_edgeDF_full <- data.frame(cbind(igraph_edgeMatrix,edgeWeightMatrix));
  #COME BACK: add column names
  #first column name must be vertex id
  igraph_attrDF_full <- data.frame(membership);
  #write over old graph object - can add edge, node attributes now that we've created.
  #df you feed in: A data frame containing a symbolic edge list in the first two columns. Additional columns are considered as edge attributes
  undirGraph_full <- graph.data.frame(igraph_edgeDF_full,directed=FALSE,vertices=igraph_attrDF_full);
  
  #scale up edge weights now that done with community detection.
  
  
  if(!is.null(edgeWeightsColName)){
  
    E(undirGraph_base)$weights <- 3*E(undirGraph_base)$weights
  
  }
  
  if(makePlots){
  
      if(!is.null(edgeWeightsColName)){
    
    if(!plotToScreen){
        
      png(filename=paste0(saveDir,"/",commMethod,"_communityPlotFull_",Sys.Date(),".png"),width=1000,height=1000,res=160)
      #,width = 700, height = 1000,res=160);
       
  
      layout(matrix(c(1,2), 1,2), widths=c(3,1))
      networkCommPlot_full <- plot(undirGraph_full, layout=layout.fruchterman.reingold,vertex.label=V(undirGraph_full)$studyNum,vertex.label.color="black",vertex.label.cex=nodeFontSize,vertex.size=nodePlotSize,
                                   vertex.color= V(undirGraph_full)$color, edge.arrow.size=3,main=paste0(nrow(igraph_attrDF_full)," clusters with edges across ",length(unique(membership[,ncol(membership)])), " datasets.\n",
                                                                                                         numCommunitiesOrig," meta-clusters of any size found ",                                                                                           "with\n community detection method ",commMethod,"."),xlab="color=meta-cluster,node=cluster,#=dataset",
                                   edge.width = E(undirGraph_base)$weights);
      
      plot.new()
      #want community numbers from smallest to largest.
      colorOrder <- sort.int(as.numeric(as.character(unique(V(undirGraph_full)$community))),index.return=TRUE,decreasing=FALSE)$ix
      colours = unique(V(undirGraph_full)$color)[colorOrder]
      #labels = paste(1:length(colours))
      labels = as.numeric(as.character(unique(V(undirGraph_full)$community)))[colorOrder]
      
      legend("center",legend=labels, col=colours, pch=19,pt.cex=2,
             title="Meta-clusters")
      dev.off();
    
    }else{
      
      
      layout(matrix(c(1,2), 1,2), widths=c(3,1))
      networkCommPlot_full <- plot(undirGraph_full, layout=layout.fruchterman.reingold,vertex.label=V(undirGraph_full)$studyNum,vertex.label.color="black",vertex.label.cex=nodeFontSize,vertex.size=nodePlotSize,
                                   vertex.color= V(undirGraph_full)$color, edge.arrow.size=3,main=paste0(nrow(igraph_attrDF_full)," clusters with edges across ",length(unique(membership[,ncol(membership)])), " datasets.\n",
                                                                                                         numCommunitiesOrig," meta-clusters of any size found ",
                                                                                                         "\nwith community detection method ",commMethod,"."),xlab="color=meta-cluster,node=cluster,#=dataset",
                                   edge.width = E(undirGraph_base)$weights);
      
      plot.new()
      #want community numbers from smallest to largest.
      colorOrder <- sort.int(as.numeric(as.character(unique(V(undirGraph_full)$community))),index.return=TRUE,decreasing=FALSE)$ix
      colours = unique(V(undirGraph_full)$color)[colorOrder]
      #labels = paste(1:length(colours))
      labels = as.numeric(as.character(unique(V(undirGraph_full)$community)))[colorOrder]
      
      legend("center",legend=labels, col=colours, pch=19,pt.cex=2,
             title="Meta-clusters")
      
    }
    
  #no edges
  }else{
    
        if(!plotToScreen){
        
      png(filename=paste0(saveDir,"/",commMethod,"_communityPlotFull_",Sys.Date(),".png"),width=1000,height=1000,res=160)
      #,width = 700, height = 1000,res=160);
       
  
      layout(matrix(c(1,2), 1,2), widths=c(3,1))
      networkCommPlot_full <- plot(undirGraph_full, layout=layout.fruchterman.reingold,vertex.label=V(undirGraph_full)$studyNum,vertex.label.color="black",vertex.label.cex=nodeFontSize,vertex.size=nodePlotSize,
                                   vertex.color= V(undirGraph_full)$color, edge.arrow.size=3,main=paste0(nrow(igraph_attrDF_full)," clusters with edges across ",length(unique(membership[,ncol(membership)])), " datasets.\n",
                                                                                                         numCommunitiesOrig," meta-clusters of any size found ",
                                                                                                         "with\n community detection method ",commMethod,"."),xlab="color=meta-cluster,node=cluster,#=dataset");
      
      plot.new()
      #want community numbers from smallest to largest.
      colorOrder <- sort.int(as.numeric(as.character(unique(V(undirGraph_full)$community))),index.return=TRUE,decreasing=FALSE)$ix
      colours = unique(V(undirGraph_full)$color)[colorOrder]
      #labels = paste(1:length(colours))
      labels = as.numeric(as.character(unique(V(undirGraph_full)$community)))[colorOrder]
      
      legend("center",legend=labels, col=colours, pch=19,pt.cex=2,
             title="Meta-clusters")
      dev.off();
    
    }else{
      
      
      layout(matrix(c(1,2), 1,2), widths=c(3,1))
      networkCommPlot_full <- plot(undirGraph_full, layout=layout.fruchterman.reingold,vertex.label=V(undirGraph_full)$studyNum,vertex.label.color="black",vertex.label.cex=nodeFontSize,vertex.size=nodePlotSize,
                                   vertex.color= V(undirGraph_full)$color, edge.arrow.size=3,main=paste0(nrow(igraph_attrDF_full)," clusters with edges across ",length(unique(membership[,ncol(membership)])), " datasets.\n",
                                                                                                         numCommunitiesOrig," meta-clusters of any size found ",
                                                                                                         "\nwith community detection method ",commMethod,"."),xlab="color=meta-cluster,node=cluster,#=dataset");
      
      plot.new()
      #want community numbers from smallest to largest.
      colorOrder <- sort.int(as.numeric(as.character(unique(V(undirGraph_full)$community))),index.return=TRUE,decreasing=FALSE)$ix
      colours = unique(V(undirGraph_full)$color)[colorOrder]
      #labels = paste(1:length(colours))
      labels = as.numeric(as.character(unique(V(undirGraph_full)$community)))[colorOrder]
      
      legend("center",legend=labels, col=colours, pch=19,pt.cex=2,
             title="Meta-clusters")
      
    }
    
  }
    
  }
  #NOW prune down...
  igraph_edgeDF <- igraph_edgeDF_full;
  igraph_attrDF <-  igraph_attrDF_full;
  
  
  commLose <- setdiff(unique(membership[ ,"community"]),commKeep);
  clustLose <- c()

  
  if(length(commLose)>0){
    
  for(l in 1:length(commLose)){
    
    nodesLose <- igraph_attrDF[which(igraph_attrDF[,"community"]==commLose[l]) ,"igraph_id"];
    clustLose <- append(clustLose,igraph_attrDF[which(igraph_attrDF[,"community"]==commLose[l]) ,"orig_id"])
    #remove these from your attribute data frame.
    igraph_attrDF <- igraph_attrDF[-which(igraph_attrDF[,"community"]==commLose[l]) ,];
    #and remove all edges that involve nodes in community commLose[l]...
    for(n in 1:length(nodesLose)){
      
      if(length(which(igraph_edgeDF[,1]==nodesLose[n]))>0){
        
        igraph_edgeDF <- igraph_edgeDF[-which(igraph_edgeDF[,1]==nodesLose[n]), ,drop=FALSE];
        
      }
      
      if(length(which(igraph_edgeDF[,2]==nodesLose[n]))>0){
        #and check end node of edge too.
        igraph_edgeDF <- igraph_edgeDF[-which(igraph_edgeDF[,2]==nodesLose[n]), ,drop=FALSE];
        
      }
      
    }
    
  }
  
  }



  if(nrow(igraph_edgeDF)==0){
    
    warning("all edges removed after pruning - no communities left.");
    finalNumCommunities <- NA;
    networkCommPlot <- NA;
    undirGraph <- NA;
    
  }else{
    
    if(length(commLose)>0){
      #NOW: update community numbers (don't want communit 1, 3, 6...but 1, 2, 3)
      commNums <- as.numeric(as.character(unique(igraph_attrDF$community)))
      newCommNums <- c(1:length(commNums))
      
      #remove factor status 
      igraph_attrDF$community <- as.numeric(as.character(igraph_attrDF$community))
      
      #renaming community numbers: index on original!
      commOrig <- igraph_attrDF$community
      for(c in 1:length(commNums)){
        
        igraph_attrDF$community[which(commOrig==commNums[c])] <- newCommNums[c]
        
      }
      
    }
  undirGraph <- graph.data.frame(igraph_edgeDF,directed=FALSE,vertices=igraph_attrDF);
  #color the nodes. V() gets the vertices. can do V()$ what added  with igraph_attrDF, like $clust.
  #vertex.size=4,vertex.label.dist=0.5,edge.weight=E(undirGraph)$sampleFract,
  #size is size of vertex
  finalNumCommunities <- length(unique(igraph_attrDF[,"community"]));
  
  if(finalNumCommunities >0){
    
  if(makePlots){
    
      
    if(!is.null(edgeWeightsColName)){
    
    if(!plotToScreen){
    
      
      png(filename=paste0(saveDir,"/",commMethod,"_communityPlotPruned_",Sys.Date(),".png"),width=1000,height=1000,res=160)
      #,width = 700, height = 1000,res=160);
      
      layout(matrix(c(1,2), 1,2), widths=c(3,1))
  
      networkCommPlot <- plot(undirGraph, layout=layout.fruchterman.reingold,vertex.label=V(undirGraph)$studyNum,vertex.label.color="black",vertex.label.cex=nodeFontSize,vertex.size=nodePlotSize,
                              vertex.color= V(undirGraph)$color, edge.arrow.size=3,main=paste0("Clusters deemed similar across ",length(unique(igraph_attrDF[,"studyNum"])), " datasets.\n",
                                                                                               finalNumCommunities," pruned meta-clusters found with\n community method ",commMethod,"."),xlab="color=meta-cluster,node=cluster,#=dataset",
                              edge.width = E(undirGraph_base)$weights);
      
      
      plot.new()
      #want community numbers from smallest to largest.
      colorOrder <- sort.int(as.numeric(as.character(unique(V(undirGraph)$community))),index.return=TRUE,decreasing=FALSE)$ix
      colours = unique(V(undirGraph)$color)[colorOrder]
      #labels = paste(1:length(colours))
      labels = as.numeric(as.character(unique(V(undirGraph)$community)))[colorOrder]
      
      legend("center",legend=labels, col=colours, pch=19,pt.cex=2,
             title="Meta-clusters")
      dev.off();
    
  
    }else{
      
      layout(matrix(c(1,2), 1,2), widths=c(3,1))
      networkCommPlot <- plot(undirGraph, layout=layout.fruchterman.reingold,vertex.label=V(undirGraph)$studyNum,vertex.label.color="black",vertex.label.cex=nodeFontSize,vertex.size=nodePlotSize,
                              vertex.color= V(undirGraph)$color, edge.arrow.size=3,main=paste0("Clusters deemed similar across ",length(unique(igraph_attrDF[,"studyNum"])), " datasets.\n",
                                                                                               finalNumCommunities," pruned meta-clusters found with\n community method ",commMethod,"."),xlab="color=meta-cluster,node=cluster,#=dataset",
                              edge.width = E(undirGraph_base)$weights);
      
      plot.new()
      #want community numbers from smallest to largest.
      colorOrder <- sort.int(as.numeric(as.character(unique(V(undirGraph)$community))),index.return=TRUE,decreasing=FALSE)$ix
      colours = unique(V(undirGraph)$color)[colorOrder]
      #labels = paste(1:length(colours))
      labels = as.numeric(as.character(unique(V(undirGraph)$community)))[colorOrder]
      
      legend("center",legend=labels, col=colours, pch=19,pt.cex=2,
             title="Meta-clusters")
      
    }
    
    #no edge weights
  }else{
    
      if(!plotToScreen){
    
      
      png(filename=paste0(saveDir,"/",commMethod,"_communityPlotPruned_",Sys.Date(),".png"),width=1000,height=1000,res=160)
      #,width = 700, height = 1000,res=160);
      
      layout(matrix(c(1,2), 1,2), widths=c(3,1))
  
      networkCommPlot <- plot(undirGraph, layout=layout.fruchterman.reingold,vertex.label=V(undirGraph)$studyNum,vertex.label.color="black",vertex.label.cex=nodeFontSize,vertex.size=nodePlotSize,
                              vertex.color= V(undirGraph)$color, edge.arrow.size=3,main=paste0("Clusters deemed similar across ",length(unique(igraph_attrDF[,"studyNum"])), " datasets.\n",
                                                                                               finalNumCommunities," pruned meta-clusters found with\n community method ",commMethod,"."),xlab="color=meta-cluster,node=cluster,#=dataset");
      
      
      plot.new()
      #want community numbers from smallest to largest.
      colorOrder <- sort.int(as.numeric(as.character(unique(V(undirGraph)$community))),index.return=TRUE,decreasing=FALSE)$ix
      colours = unique(V(undirGraph)$color)[colorOrder]
      #labels = paste(1:length(colours))
      labels = as.numeric(as.character(unique(V(undirGraph)$community)))[colorOrder]
      
      legend("center",legend=labels, col=colours, pch=19,pt.cex=2,
             title="Meta-clusters")
      dev.off();
    
  
    }else{
      
      layout(matrix(c(1,2), 1,2), widths=c(3,1))
      networkCommPlot <- plot(undirGraph, layout=layout.fruchterman.reingold,vertex.label=V(undirGraph)$studyNum,vertex.label.color="black",vertex.label.cex=nodeFontSize,vertex.size=nodePlotSize,
                              vertex.color= V(undirGraph)$color, edge.arrow.size=3,main=paste0("Clusters deemed similar across ",length(unique(igraph_attrDF[,"studyNum"])), " datasets.\n",
                                                                                               finalNumCommunities," pruned meta-clusters found with\n community method ",commMethod,"."),xlab="color=meta-cluster,node=cluster,#=dataset");
      
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
    
   #if make plots 
  }
    #zero communities after pruning.
  }else{
    
    message("no final communities")
    output <- list(numCommunities=0)
    return(output)
    
    
  }
  
  }
  
  }else{
    
    message("No final communities because zero edges in input edge matrix.")
    output <- list(numCommunities=0)
    return(output)
    
    
  }

  if(saveGraphData){
    
    write.table(igraph_attrDF,file=paste0(saveDir,"/nodeAttributes_pruned_",fileTag,"_",Sys.Date(),".txt"),sep="\t",quote=FALSE,row.names=FALSE,
                col.names=FALSE);
    
    
    write.table(igraph_attrDF_full,file=paste0(saveDir,"/nodeAttributes_full_",fileTag,"_",Sys.Date(),".txt"),sep="\t",quote=FALSE,row.names=FALSE,
                col.names=FALSE);
    
    write.table(igraph_edgeDF,file=paste0(saveDir,"/edgeMatrixWeights_pruned",fileTag,"_",Sys.Date(),".txt"),sep="\t",quote=FALSE,row.names=FALSE,
                col.names=FALSE);
    
    write.table(igraph_edgeDF_full,file=paste0(saveDir,"/edgeMatrixWeights_full",fileTag,"_",Sys.Date(),".txt"),sep="\t",quote=FALSE,row.names=FALSE,
                col.names=FALSE);
    
    
    #zero communities at initial input - no edges.
  }
  
  

  

  output <- list(numCommunities=finalNumCommunities,numCommunitiesOrig=numCommunitiesOrig,undirGraph=undirGraph,communityObjectUnpruned=comm,
                 ClustKey=graphKey,undirGraph_base=undirGraph_base,unPrunedDendogram=unPrunedDendogram,
                 edgeDF=igraph_edgeDF,attrDF=igraph_attrDF,undirGraph_full=undirGraph_full,edgeDF_full=igraph_edgeDF_full,attrDF_full=igraph_attrDF_full,
                 commMethod=commMethod,modularity=modularity,communityObject_full=comm,commLose=commLose,clustLose=clustLose,finalCommEdgeInfo=commEdgeInfo);
  return(output);
  
}

##return a membership matrix

returnSampleMemberMatrix <- function(clustSampleIndexList,dataMatrixList,communityInfo,
                                     noClustNA=FALSE){
  
sampleNames <- c()
for(d in 1:length(clustSampleIndexList)){
  
  for(c in 1:length(clustSampleIndexList[[d]])){
    #drop=FALSE: in case a cluster had only 1 patient.
    sampleNames <- append(sampleNames,colnames(dataMatrixList[[d]][ , clustSampleIndexList[[d]][[c]],drop=FALSE]))
    
  }
  
}

numTotalSamples <- length(sampleNames)

clustNum <- 1

globalClustNum <- c()
studyClustNum <- c()
studyNum <- c()
community <- c()
for(d in 1:length(clustSampleIndexList)){

  for(c in 1:length(clustSampleIndexList[[d]])){

    #just take first match - all 
    globalClustNum <- append(globalClustNum,rep.int(clustNum,times=length(clustSampleIndexList[[d]][[c]])))
    studyClustNum <- append(studyClustNum,rep.int(c,times=length(clustSampleIndexList[[d]][[c]])) )
    studyNum <- append(studyNum,rep.int(d,times=length(clustSampleIndexList[[d]][[c]])) )


    #was this in the final clustering?
    #==clustNum
    if(length(which(communityInfo$attrDF[,"clust"]==clustNum))>0){
      #message("commNum: ",commNum)
      commNum <- unique(communityInfo$attrDF[which(communityInfo$attrDF[,"clust"]==clustNum), "community"])
      
      if(length(commNum)!=1){
        
        stop("expected 1 unique community number per cluster, but returned greater than 1 or 0.")
        
      }
      
      community <- append(community,rep.int(commNum,times=length(clustSampleIndexList[[d]][[c]])) )
      
    }else{
      
      community <- append(community,rep.int(NA,times=length(clustSampleIndexList[[d]][[c]])) )
      
    }

   clustNum <- clustNum + 1

  }
  
}

sampleClustCommKey <- data.frame(sampleNames,globalClustNum,studyClustNum,studyNum,community,stringsAsFactors=FALSE)
colnames(sampleClustCommKey) <- c("sampleName","globalClustNum","studyClustNum","studyNum","community")

fullMemberMatrix <- matrix(data=0,ncol=numTotalSamples,nrow=numTotalSamples,dimnames=list(sampleNames,sampleNames))

#if for some reason sample not included in this run (i.e. we were resampling only a % of full datasets), set these samples to zero.
#also set to NA if the cluster was thresholded out or the community was too small based on user-defined thresholds and removed

if(noClustNA){
if(length(which(is.na(sampleClustCommKey[,"community"])))>0){
  #there can be no shared membership across any samples for these samples that were not included in the final communities.
  fullMemberMatrix[which(is.na(sampleClustCommKey[,"community"])), ] <- NA
  fullMemberMatrix[,which(is.na(sampleClustCommKey[,"community"])) ] <- NA

}

}#else: just set to zero if was never in a community. these are already set as zero in the original matrix

commNums <- unique(sampleClustCommKey[,"community"])
for(i in 1:length(commNums)){
   
  fullMemberMatrix[which(sampleClustCommKey[,"community"]==commNums[i]),which(sampleClustCommKey[,"community"]==commNums[i])] <- 1
  
}
  
output <- list(fullMemberMatrix=fullMemberMatrix,sampleClustCommKey=sampleClustCommKey)
 return(output)

}
   


