#!/usr/bin/env Rscript

# Cluster Rand Index Sensitivity Comparison Script
# Created by Isaac Robson

# install.packages('pdfCluster')

library("igraph")
library('pdfCluster')

source("~/GitHub/kinotypeClusteringAndRegression/src/tools/communityHelpers.R")

## do some parallelization
library('parallel')
library('purrr')

# read graph and weights
G <- read.graph("~/Github/KIN_ClusteringWithAnnotations/data/kin/kin_anscombe_weighted.csv",format="ncol",names=TRUE,weights="yes",directed=FALSE)
W <- read.table("~/Github/KIN_ClusteringWithAnnotations/data/kin/kin_anscombe_weighted.csv")$V3

## present algorithms
# cluster_louvain (hierarchical modularity) 

V(G)$comp <- components(G)$membership
mainG <- induced_subgraph(G,V(G)$comp==1)

hist(E(mainG)$weight)
hist(degree(mainG))

### cluster_louvain (deterministic, as in network_03)
louv <- cluster_louvain(mainG)
louv_clusts <- data.frame(names=louv$names, cluster=louv$memberships[2,])
louv_small_clusts <- data.frame(names=louv$names, cluster=louv$memberships[1,])

## reference info for comparisons
numnodes <- length(V(mainG))
mean_weight <- mean(E(mainG)$weight)

single_rand_compare_louv <- function(edge_to_add, g, new_weight, refclust){
  newG <- add_edges(g, c(edge_to_add[1], edge_to_add[2]), weight=new_weight)
  clust <- igraph::cluster_louvain(newG) # auto detects weights
  return(pdfCluster::adj.rand.index(refclust, clust$memberships[1,]))
}

para_rand_compare_louv <- function(g, numnodes, new_weight, refclust, cores=max(detectCores()-1,1)){
  # create a list of results and a list of parameters for the spinglasses
  new_edges <- combn(numnodes, 2, simplify = FALSE)
  
  # start up a cluster
  res <- parallel::mclapply(X = new_edges, FUN = single_rand_compare_louv, g=g, new_weight=new_weight, refclust=refclust, mc.cores = cores)
  
  # convert to upper triangular
  out <- mat.or.vec(numnodes, numnodes)
  res <- unlist(res)
  print(res)
  # unusual combn order does not align with upper triangular
  counter <- 1
  for (i in 1:(numnodes-1)){
    out[i, (i+1):numnodes] <- res[counter:(counter + numnodes - i - 1)]
    counter <- counter + numnodes - i - 1
  }
  return(out)
}

out <- para_rand_compare_louv(mainG, numnodes = length(V(mainG)), new_weight = mean_weight, refclust = louv_small_clusts$cluster)
which(V(mainG)$name == "TLK1")
## symmetrize
out[lower.tri(out)] <- t(out)[lower.tri(out)]

hist(colSums(out != 1), breaks = 50)
hist(out[1,], breaks = 100)
hist(out[133,], breaks = 100)

outfile="~/Github/kinotypeClusteringAndRegression/results/sensitivityNetworkClusters/sensitivity_randscores.tsv"
# write.table(out,outfile,quote=FALSE,sep="\t",row.names = FALSE)

# Rand Sensitivity Analysis
G_insrr <- G
G_insrr <- add_edges(G_insrr, c("INSRR", ""))

G_eif2ak2_amhr2 <- mainG
G_eif2ak2_amhr2 <- add_edges(G_eif2ak2_amhr2, c("MST1R", "FGR"), weight=mean_weight)
cluster_eif2ak2_amhr2 <- cluster_louvain(G_eif2ak2_amhr2)
eif2ak2_amhr2_clusts <- data.frame(names=cluster_eif2ak2_amhr2$names, cluster=cluster_eif2ak2_amhr2$memberships[2,])
eif2ak2_amhr2_small_clusts <- data.frame(names=cluster_eif2ak2_amhr2$names, cluster=cluster_eif2ak2_amhr2$memberships[1,])
E(G_eif2ak2_amhr2)[E(G_eif2ak2_amhr2)]
a <- E(G_eif2ak2_amhr2)[1]

V(mainG)
library(clues)
new_edges[[1]][[1]]
new_edges <- combn(numnodes, 2, simplify = FALSE)
for (i in 517:530){
  t<- single_rand_compare_louv(c(new_edges[[i]][[1]], new_edges[[i]][[2]]), mainG, mean_weight, louv_small_clusts$cluster)
  print(t)
}

new_edges[1:5]

adjustedRand(eif2ak2_amhr2_small_clusts$cluster, louv_small_clusts$cluster)

which(V(mainG)$name == "TYRO3")
E(mainG)[ from("TYRO3") ]
E(G_eif2ak2_amhr2)[ from("TYRO3") ]
V(mainG)[17]
louv_clusts

eif2ak2_amhr2_small_clusts$cluster == louv_small_clusts$cluster
