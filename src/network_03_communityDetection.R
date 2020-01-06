#!/usr/bin/env Rscript

# Community Detection (Graph Clustering) Script
# Created by Kyla Collins
# Refactored and Parallelized by Isaac Robson

library('parallel')
library('purrr')
library("igraph")

# paste(..., sep = " ", collapse = NULL)

### Important Stuff ### 
dir_path = "/Users/isrobson/Github/kinotypeClusteringAndRegression/"
helper_file = "src/tools/communityHelpers.R"
graph_dir = "data/interactionNetworks/"

# clustering results filenames
fastgreedy_file = "fastgreedy_clusters.txt"
spinglass_file = "consensus_spinglass.txt"
lev_file = "consensus_eigenvector.txt"
lp_file = "consensus_label_propagation.txt"
wt_file = "consensus_walktrap.txt"
louv_file = "louvain_clusters.txt"
louv_small_file = "louvain_small_clusters.txt"
info_file = "consensus_infomap.txt"
eb_file = "edge_betweenness_community_clusters.txt"

# modularity results filename
mod_file = "clustering_modularity_results.txt"

# grab the helper file
source(paste(dir_path, helper_file, sep=""))

# This parameter controls the number of iterations for consensus stochastic algorithms
numiter <- 1000

### Actual Script   ### 

## algorithms present (description)
# fastgreedy.community (deterministic) -- run 1x
# spinglass.community (Potts model--statistical) -- run 1000x
# leading.eigenvector.community (top-down hierarchical) -- run 1000x
# label.propogation.community (relies on initial random seeds) -- run 1000x
# walktrap.community (random walks) -- run 1000x
# cluster_louvain (hierarchical modularity) -- run 1x
# cluster_infomap (random walks) -- run 1000x
# edge.betweenness.community  -- this will take a long time

for (network_config in c("anscombe_weighted", "unweighted", "both_weighted", "coals_weighted", "raw_weighted")){
  
  results_dir = "results/networkClusters/"
  
  # read in the graph and set results path
  if (network_config == "unweighted"){
    mainG <- read.graph(paste(dir_path, graph_dir, "kin_", network_config, ".csv", sep=""),format="ncol",names=TRUE,directed=FALSE)
    W = NULL
  }
  else {
    mainG <- read.graph(paste(dir_path, graph_dir, "kin_", network_config, ".csv", sep=""),format="ncol",names=TRUE,weights="yes",directed=FALSE)
    
    V(mainG)$comp <- components(mainG)$membership
    mainG <- induced_subgraph(mainG,V(mainG)$comp==1)
    W = E(mainG)$weight
  }

  # reconfigure results dir for subdir
  results_dir = paste(results_dir, network_config, "/", sep="")
  
  ### fastgreedy
  fc <- fastgreedy.community(mainG, weights=W)
  fg_clusts <- data.frame(names=fc$names, cluster=fc$membership)
  write.table(fg_clusts, paste(dir_path, results_dir, fastgreedy_file, sep=""),quote=FALSE,sep="\t",row.names=FALSE)

  # spinglass doesn't work well without asymmetric weights
  if (network_config != "unweighted"){
    ### spinglass.community
    sc <- spinglass.community(mainG, spins=100)
    
    # create a votes matrix to store individual cluster tallies
    numnodes <- length(sc$names)
    votes <- mat.or.vec(numnodes,numnodes)
    dim(votes)
    
    print(Sys.time())
    votes <- para_spinglass(g=mainG, numnodes = numnodes, numiter = numiter)
    print(Sys.time())
    
    thresh <- 0.9*numiter
    visited <- mat.or.vec(numnodes,1)
    groups <- mat.or.vec(numnodes,1)
    k <- 1
    for (i in 1:numnodes){
    	x <- which(votes[i,] > thresh)
    	if (visited[x[1]] == 0){
    		visited[x] <- 1
    		groups[x] <- k
    		k <- k + 1
    	}
    }
    
    sc_clusts <- data.frame(names=sc$names, cluster=groups)
    
    if (min(sc_clusts$cluster) == 0){
      sc_clusts$cluster = sc_clusts$cluster + 1
    }
    
    write.table(sc_clusts, paste(dir_path, results_dir, spinglass_file, sep=""),quote=FALSE,sep="\t",row.names=FALSE)
  }
  
  ### leading.eigenvector.community
  lev <- leading.eigenvector.community(mainG)
  
  numnodes <- length(lev$names)
  votes <- mat.or.vec(numnodes,numnodes)
  
  print(Sys.time())
  votes <- para_lev(g=mainG, numnodes = numnodes, numiter = numiter)
  print(Sys.time())
  
  thresh <- 0.9*numiter
  visited <- mat.or.vec(numnodes,1)
  groups <- mat.or.vec(numnodes,1)
  k <- 1
  for (i in 1:numnodes){
    x <- which(votes[i,] > thresh)
    if (visited[x[1]] == 0){
      visited[x] <- 1
      groups[x] <- k
      k <- k + 1
    }
  }
  
  lev_clusts <- data.frame(names=lev$names, cluster=groups)
  write.table(lev_clusts, paste(dir_path, results_dir, lev_file, sep=""),quote=FALSE,sep="\t",row.names=FALSE)
  
  ### label.propagation.community
  lp <- label.propagation.community(mainG)
  
  numnodes <- length(lp$names)
  votes <- mat.or.vec(numnodes,numnodes)
  
  print(Sys.time())
  votes <- para_lp(g=mainG, numnodes = numnodes, numiter = numiter)
  print(Sys.time())
  
  thresh <- 0.9*numiter
  visited <- mat.or.vec(numnodes,1)
  groups <- mat.or.vec(numnodes,1)
  k <- 1
  for (i in 1:numnodes){
    x <- which(votes[i,] > thresh)
    if (visited[x[1]] == 0){
      visited[x] <- 1
      groups[x] <- k
      k <- k + 1
    }
  }
  
  lp_clusts <- data.frame(names=lp$names, cluster=groups)
  write.table(lp_clusts, paste(dir_path, results_dir, lp_file, sep=""),quote=FALSE,sep="\t",row.names=FALSE)
  
  ### walktrap.community
  wt <- walktrap.community(mainG, modularity=TRUE, steps=gorder(mainG))
  
  numnodes <- length(wt$names)
  votes <- mat.or.vec(numnodes,numnodes)
  
  print(Sys.time())
  votes <- para_wt(g=mainG, numnodes = numnodes, numiter = numiter, steps=gorder(mainG))
  print(Sys.time())
  
  thresh <- 0.9*numiter
  visited <- mat.or.vec(numnodes,1)
  groups <- mat.or.vec(numnodes,1)
  k <- 1
  for (i in 1:numnodes){
    x <- which(votes[i,] > thresh)
    if (visited[x[1]] == 0){
      visited[x] <- 1
      groups[x] <- k
      k <- k + 1
    }
  }
  
  wt_clusts <- data.frame(names=wt$names, cluster=groups)
  write.table(wt_clusts, paste(dir_path, results_dir, wt_file, sep=""),quote=FALSE,sep="\t",row.names=FALSE)
  
  ### cluster_louvain
  louv <- cluster_louvain(mainG, weights = W)
  louv_clusts <- data.frame(names=louv$names, cluster=louv$memberships[2,])
  louv_small_clusts <- data.frame(names=louv$names, cluster=louv$memberships[1,])
  write.table(louv_clusts,  paste(dir_path, results_dir, louv_file, sep=""),quote=FALSE,sep="\t",row.names=FALSE)
  write.table(louv_small_clusts, paste(dir_path, results_dir, louv_small_file, sep=""),quote=FALSE,sep="\t",row.names=FALSE)
  
  ### cluster_infomap
  info <- cluster_infomap(mainG)
  
  numnodes <- length(info$names)
  votes <- mat.or.vec(numnodes,numnodes)
  
  print(Sys.time())
  votes <- para_info(g=mainG, numnodes = numnodes, numiter = numiter)
  print(Sys.time())
  
  thresh <- 0.9*numiter
  visited <- mat.or.vec(numnodes,1)
  groups <- mat.or.vec(numnodes,1)
  k <- 1
  for (i in 1:numnodes){
    x <- which(votes[i,] > thresh)
    if (visited[x[1]] == 0){
      visited[x] <- 1
      groups[x] <- k
      k <- k + 1
    }
  }
  
  info_clusts <- data.frame(names=info$names, cluster=groups)
  write.table(info_clusts, paste(dir_path, results_dir, info_file, sep=""),quote=FALSE,sep="\t",row.names=FALSE)
  
  ### edge.betweenness.community
  # edge betweenness treats weights as distance, but membership corrects for this
  # note - EB fails for zero weights -- skip if a zero present
  if (!is.null(W)){
    if (min(W) > 0){
      eb <- edge.betweenness.community(mainG, weights = W, membership=TRUE) 
      eb_clusts <- data.frame(names=eb$names, cluster=eb$membership)
      write.table(eb_clusts, paste(dir_path, results_dir, eb_file, sep=""),quote=FALSE,sep="\t",row.names=FALSE)
    }
  }
  ###
  
  #### collect modularity data
  mod <- data.frame(row.names = "modularity")
  mod$fast_greedy <- modularity(mainG,fg_clusts$cluster,weights=W)
  
  # modularity function doesn't accept the '0' cluster name, so we check if we
  # need to shift membership values up by one to get the modularity
  if (network_config != "unweighted"){
    mod$spinglass <- modularity(mainG,sc_clusts$cluster,weights=W)
  }
  mod$eigen <- modularity(mainG,lev_clusts$cluster,weights=W)
  mod$walktrap <- modularity(mainG,wt_clusts$cluster,weights=W)
  if (min(lp_clusts$cluster) == 0){
    lp_clusts$cluster = lp_clusts$cluster + 1
  }
  mod$label <- modularity(mainG,lp_clusts$cluster,weights=W)
  mod$louvain <- modularity(mainG,louv_clusts$cluster,weights=W)
  mod$small_louvain <- modularity(mainG,louv_small_clusts$cluster,weights=W)
  if (min(info_clusts$cluster) == 0){
    info_clusts$cluster = info_clusts$cluster + 1
  }
  mod$infomap <- modularity(mainG,info_clusts$cluster,weights=W)
  
  if (network_config != "unweighted"){
    if (min(W) > 0){
      mod$edge_between <- modularity(mainG,eb_clusts$cluster,weights=W)
    }
  }
  
  
  ### write modularity data out to a file
  
  write.table(mod,paste(dir_path, results_dir, mod_file, sep=""),quote=FALSE,sep="\t",row.names = FALSE)
}
