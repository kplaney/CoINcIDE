
#######ANALYZE OUTPUT FUNCTIONS
assignFinalEdges <- function(computeTrueSimilOutput,pvalueMatrix,indEdgePvalueThresh=.1,
                             meanEdgePairPvalueThresh=.05,
                             minTrueSimilThresh=-Inf,maxTrueSimilThresh=Inf,
                             minFractFeatureIntersect=0,fractFeatIntersectThresh=0,numFeatIntersectThresh=0 ,
                             clustSizeThresh=0, clustSizeFractThresh=0,saveDir="./",fileTag="CoINcIDE_edges"
){
  
  clustSizeIndexRemove <- c()
  count <- 0
  meanEdgePvalueMatrix <- matrix(data=NA,ncol=ncol(pvalueMatrix),nrow=nrow(pvalueMatrix))
  
  for(r in 1:nrow(computeTrueSimilOutput$clustSizeMatrix)){
    
    if(computeTrueSimilOutput$clustSizeFractMatrix[r]==1){
      #ROW is when this cluster was the reference cluster. p-value will by default be insignificant.
      #just replace with p-value in the other direction (columnwise)
      #(if both clusters in a pair were entire size of matrix, will have NA p-value still.)
      pvalueMatrix[, r] <-   pvalueMatrix[r, ]
      
    }
    if(computeTrueSimilOutput$clustSizeMatrix[r] <= clustSizeThresh || computeTrueSimilOutput$clustSizeFractMatrix[r] <= clustSizeFractThresh){
      #remove these clusters
      pvalueMatrix[r, ] <- NA
      count <- count + 1
      clustSizeIndexRemove[count] <- r
    }
    
  }
  
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
        meanEdgePvalueMatrix[r,c] <- mean(c(pvalueMatrix[r,c],pvalueMatrix[c,r]))
        meanEdgePvalueMatrix[c,r] <- meanEdgePvalueMatrix[r,c]
        
      }
      
    }
    
  }
  
  adjMatricesList <- list(computeTrueSimilOutput$similValueMatrix,computeTrueSimilOutput$similValueMatrix,
                          pvalueMatrix,meanEdgePvalueMatrix,computeTrueSimilOutput$fractFeatIntersectMatrix,
                          computeTrueSimilOutput$numFeatIntersectMatrix)
  
  names(adjMatricesList) <- c("simil","simil","meanPvalue","meanPvalue","fractFeatIntersect","numFeatIntersect")
  
  thresholdVector <- c(minTrueSimilThresh,maxTrueSimilThresh,indEdgePvalueThresh,meanEdgePairPvalueThresh,
                       fractFeatIntersectThresh,numFeatIntersectThresh)
  
  threshDir <- c(">=","<=","<=","<=",">=",">=")
  
  filterEdgeOutput <- filterEdges(adjMatricesList,thresholdVector,threshDir=threshDir,saveDir=saveDir,fileTag=fileTag,saveEdgeFile=TRUE)
  
  clustFurtherRemoved <- setdiff(filterEdgeOutput$clustRemoved,clustSizeIndexRemove)
  
  allClustRemoved <- filterEdgeOutput$clustRemoved
  
  output <- list(meanEdgePvalueMatrix=meanEdgePvalueMatrix,clustSizeIndexRemove=clustSizeIndexRemove,clustFurtherRemoved=clustFurtherRemoved,filterEdgeOutput=filterEdgeOutput,
                 allClustRemoved=allClustRemoved,adjMatricesList=adjMatricesList)
  return(output)
  
}

#create wrapper function for filterEdges (this is what user sees.)
filterEdges <- function(adjMatricesList,thresholdVector,threshDir=rep(">=",length(thresholdVector)),saveDir="/home/kplaney/",fileTag="",saveEdgeFile=TRUE){
  
  #creates directory if there isn't one
  dir.create(saveDir);

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
          
          #must have non-NA adjacency vales
          if(!is.na(adjMatricesList[[a]][p,n]) && !is.na(adjMatricesList[[a]][n,p])){
            #is it below threshold? must pass ALL thresholds.
            #NOTE: will need to figure something else if doing a Euc dist...want > then.
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

findCommunities <- function(edgeMatrix,edgeWeightMatrix,clustIndexMatrix,fileTag="IGPN_communityNodeAttributes_",
                            saveDir="./",minNumUniqueStudiesPerCommunity=3,clustMethodName="sparseBC",
                            commMethod=c("fastGreedy","edgeBetween","walktrap","eigenvector","optimal","spinglass","multilevel"),
                            makePlots=TRUE,saveGraphData=TRUE){
  
  dir.create(saveDir,showWarnings=TRUE);
  #add on community detection name.
  clustMethodName <- paste0(clustMethodName,"_commMethod_",commMethod);
  #library("network");
  library("igraph");
  library("limma");
  library("RColorBrewer")
  #  undirGraph <- graph.edgeMatrix(edgeMatrix,directed=FALSE);
  #igraph likes an edge list that is continous - i.e. if there's not an edge from node 10, it will still include node 10.
  # we don't wan't this! collapse down the node indices.
  graphKey <- cbind(unique(append(edgeMatrix[,1],edgeMatrix[,2])),
                    c(1:length(unique(append(edgeMatrix[,1],edgeMatrix[,2])))));
  
  colnames(graphKey) <- c("orig_id","new_id");
  igraph_edgeMatrix <- matrix(data=NA,nrow=nrow(edgeMatrix),ncol=ncol(edgeMatrix));
  
  for(r in 1:nrow(graphKey)){
    
    igraph_edgeMatrix[which(edgeMatrix==graphKey[r,"orig_id"])] <- graphKey[r,"new_id"];
    
  }
  
  
  clustNames <- union(unique(edgeMatrix[,1]),unique(edgeMatrix[,2]));
  studyNames <-  clustIndexMatrix[na.omit(match(clustNames,clustIndexMatrix[,1])) ,2];
  #order in same order as orig names.
  igraph_clustNames <- graphKey[match(clustNames,graphKey[,"orig_id"]),"new_id"];

  base_attr <- cbind(igraph_clustNames,studyNames);
  colnames(base_attr)[2] <- "studyNum";
  undirGraph_base <-  graph.data.frame(data.frame(igraph_edgeMatrix),directed=FALSE,vertices=data.frame(base_attr));
 
  if(makePlots){
    
    png(filename=paste0(saveDir,clustMethodName,"_origNetworkPlot_",Sys.Date(),".png"),
      width = 700, height = 1000)
    origNetworkPlot <- plot(undirGraph_base, layout=layout.fruchterman.reingold,vertex.size=5,vertex.label=V(undirGraph_base)$studyNum,vertex.label.color="black",vertex.label.cex=.5,
       edge.arrow.size=3,main=paste0(length(union(unique(edgeMatrix[,1]),unique(edgeMatrix[,2]))), " clusters deemed similar\n across ",length(unique(studyNames))," studies before community detection."
                                                                            ),xlab="#=Study Number");
   
    dev.off();
    
  }
  #commMethod=c("fastGreedy","edgeBetween","walktrap","eigenvector","optimal","spinglass","multilevel")
  #calculate the communities. fast greedy actually seems to work best for our cases. it can also worked on weighted graphs.
  if(commMethod=="fastgreedy"){
    
  comm <- fastgreedy.community(undirGraph_base,merges=TRUE,modularity=TRUE,
                                     membership=TRUE);

  }else if(commMethod=="edgeBetween"){
    
    #note: can actually use weights here.
    comm <- edge.betweenness.community(undirGraph_base,directed=FALSE,merges=TRUE,modularity=TRUE,
                                 membership=TRUE);
    
    
  }else if(commMethod=="walktrap"){
    #steps=4 is default
    comm <- walktrap.community(undirGraph_base,merges=TRUE,modularity=TRUE,
                               membership=TRUE,steps=4);
    
  }else if(commMethod=="eigenvector"){
    
    
    comm <- leading.eigenvector.community(undirGraph_base,steps = -1);
    
  }else if(commMethod=="optimal"){
    
    #maximizing the modularity measure over all possible partitions
    comm <- optimal.community(undirGraph_base);
    
  }else if(commMethod=="spinglass"){
    
    comm <- spinglass.community(undirGraph_base,weights=NULL, vertex=NULL, spins=25,
                                parupdate=FALSE, start.temp=1, stop.temp=0.01,
                                cool.fact=0.99, update.rule=c("config"), gamma=1, implementation=c("orig"),
                                gamma.minus=1);
    
  }else if(commMethod=="multilevel"){
    
    comm <- multilevel.community(undirGraph_base);
    
  }else{
    
    stop("\nPlease pick a viable community detection method for the 'commMethod' variable.")
  }
    
  modularity <- comm$modularity;
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
    colorCodeF <- colorRampPalette(c("violet","blue","red"), bias = length(unique(membership[,"community"]))*2,
                                   space = "rgb", interpolate = "linear");
    
    #this will produce RGB representation in HEX, which cytoscape can take in.
    #if need plain old RGB: col2RGB(membership[,"colors"],alpha=FALSE);
    colorCodes <- colorCodeF(length(unique(membership[,"community"])));
  
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
  
  numMembersInComm <- table(membership[,"community"]);
  
  #this is baseline...but those studies may not be all unique
  #commKeepTemp <- which(numMembersInComm>=minNumNodesPerCommunity);
  
  #now see if this actually meets the unique threshold
  commKeep <- c();
  k <- 0;
  for(c in 1:numCommunitiesOrig){
    
    if(length(unique(membership[which(membership[,"community"]==unique(membership[,"community"])[c]), "studyNum"]))>=minNumUniqueStudiesPerCommunity){
      
      k <- 1 + k;
      commKeep[k] <- unique(membership[,"community"])[c];
      
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
  
  if(makePlots){
  
    png(filename=paste0(saveDir,clustMethodName,"_communityPlotFull_",Sys.Date(),".png"),
      width = 700, height = 1000);
    networkCommPlot_full <- plot(undirGraph_full, layout=layout.fruchterman.reingold,vertex.label=V(undirGraph_full)$studyNum,vertex.size=5,vertex.label.color="black",vertex.label.cex=.5,
                               vertex.color= V(undirGraph_full)$color, edge.arrow.size=3,main=paste0(nrow(igraph_attrDF_full)," clusters deemed similar across ",length(unique(membership[,ncol(membership)])), " studies.\n",
                                                                                                     numCommunitiesOrig," communities of any size found \nusing method ",clustMethodName,
                                                                                                     "with community detection ",commMethod,"."),xlab="color=community, #=study");
  
    dev.off();
    
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
    
  undirGraph <- graph.data.frame(igraph_edgeDF,directed=FALSE,vertices=igraph_attrDF);
  #color the nodes. V() gets the vertices. can do V()$ what added  with igraph_attrDF, like $clust.
  #vertex.size=4,vertex.label.dist=0.5,edge.weight=E(undirGraph)$sampleFract,
  #size is size of vertex
  finalNumCommunities <- length(unique(igraph_attrDF[,"community"]));
  
  if(makePlots){
    #COME BACK: update resolutions, colors.
    png(filename=paste0(saveDir,clustMethodName,"_communityPlotPruned_",Sys.Date(),".png"),
      width = 700, height = 1000);
  
    networkCommPlot <- plot(undirGraph, layout=layout.fruchterman.reingold,vertex.label=V(undirGraph)$studyNum,vertex.size=5,vertex.label.color="black",vertex.label.cex=.5,
                          vertex.color= V(undirGraph)$color, edge.arrow.size=3,main=paste0("Clusters deemed similar across ",length(unique(membership[,ncol(membership)])), " studies.\n",
                                                                                           finalNumCommunities," pruned communities found using method ",clustMethodName, "with community detection ",
                                                                                           commMethod,"."),xlab="color=community, #=study");

  
    dev.off();
  
  }
  
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
    
  
  
  }
  output <- list(numCommunities=finalNumCommunities,numCommunitiesOrig=numCommunitiesOrig,undirGraph=undirGraph,communityObject=comm,networkCommPlot=networkCommPlot,
                 networkCommPlot_unpruned=networkCommPlot_full,ClustKey=graphKey,undirGraph_base=undirGraph_base,origNetworkPlot=origNetworkPlot,
                 edgeDF=igraph_edgeDF,attrDF=igraph_attrDF,undirGraph_full=undirGraph_full,edgeDF_full=igraph_edgeDF_full,attrDF_full=igraph_attrDF_full,
                 commMethod=commMethod,modularity=modularity,communityObject_full=comm,commLose=commLose,clustLose=clustLose);
  return(output);
  
}