#!/usr/bin/env Rscript
source("~/GitHub/KIN_ClusteringWithAnnotations/src/tools/communityHelpers.R")
library("igraph")

# This parameter controls the number of iterations for consensus stochastic algorithms
numiter <- 1000

# read graph and weights
G <- read.graph("~/Github/KIN_ClusteringWithAnnotations/data/kin_anscombe_weighted.csv",format="ncol",names=TRUE,weights="yes",directed=FALSE)
W <- read.table("~/Github/KIN_ClusteringWithAnnotations/data/kin_anscombe_weighted.csv")$V3

## present algorithms
# fastgreedy.community (deterministic) -- run 1x
# spinglass.community (Potts model--statistical) -- run 1000x
# leading.eigenvector.community (top-down hierarchical) -- run 1000x
# label.propogation.community (relies on initial random seeds) -- run 1000x
# walktrap.community (random walks) -- run 1000x
# cluster_louvain (hierarchical modularity) -- run 1x
# cluster_infomap (random walks) -- run 1000x
# edge.betweenness.community  -- this will take a long time

V(G)$comp <- components(G)$membership
mainG <- induced_subgraph(G,V(G)$comp==1)

### fastgreedy
fc <- fastgreedy.community(mainG)
fg_clusts <- data.frame(names=fc$names, cluster=fc$membership)
write.table(fg_clusts, '~/GitHub/KIN_ClusteringWithAnnotations/results/weighted/fastgreedy_clusters.txt',quote=FALSE,sep="\t",row.names=FALSE)

### spinglass.community
sc <- spinglass.community(mainG, spins=100)

# create a votes matrix to store individual cluster tallies
numnodes <- length(sc$names)
votes <- mat.or.vec(numnodes,numnodes)
dim(votes)

## do some parallelization
library('parallel')
library('purrr')

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
write.table(sc_clusts, '~/Github/KIN_ClusteringWithAnnotations/results/weighted/consensus_spinglass.txt',quote=FALSE,sep="\t",row.names=FALSE)

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
write.table(lev_clusts, '~/Github/KIN_ClusteringWithAnnotations/results/weighted/consensus_eigenvector.txt',quote=FALSE,sep="\t",row.names=FALSE)

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
write.table(lp_clusts, '~/GitHub/KIN_ClusteringWithAnnotations/results/weighted/consensus_label_propagation.txt',quote=FALSE,sep="\t",row.names=FALSE)

### walktrap.community
wt <- walktrap.community(mainG, modularity=TRUE, steps=10)

numnodes <- length(wt$names)
votes <- mat.or.vec(numnodes,numnodes)

print(Sys.time())
votes <- para_wt(g=mainG, numnodes = numnodes, numiter = numiter)
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
write.table(wt_clusts, '~/Github/KIN_ClusteringWithAnnotations/results/weighted/consensus_walktrap.txt',quote=FALSE,sep="\t",row.names=FALSE)

### cluster_louvain
louv <- cluster_louvain(mainG)
louv_clusts <- data.frame(names=louv$names, cluster=louv$memberships[2,])
louv_small_clusts <- data.frame(names=louv$names, cluster=louv$memberships[1,])
write.table(louv_clusts, '~/GitHub/KIN_ClusteringWithAnnotations/results/weighted/louvain_clusters.txt',quote=FALSE,sep="\t",row.names=FALSE)
write.table(louv_small_clusts, '~/GitHub/KIN_ClusteringWithAnnotations/results/weighted/louvain_small_clusters.txt',quote=FALSE,sep="\t",row.names=FALSE)

### cluster_infomap
info <- walktrap.community(mainG)

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
write.table(info_clusts, '~/Github/KIN_ClusteringWithAnnotations/results/weighted/consensus_infomap.txt',quote=FALSE,sep="\t",row.names=FALSE)

### edge.betweenness.community
eb <- edge.betweenness.community(mainG)
eb_clusts <- data.frame(names=eb$names, cluster=eb$membership)
write.table(eb_clusts, '~/GitHub/KIN_ClusteringWithAnnotations/results/weighted/edge_betweenness_community_clusters.txt',quote=FALSE,sep="\t",row.names=FALSE)

###

#### collect modularity data
mod <- data.frame(row.names = "modularity")
mod$fast_greedy <- modularity(mainG,fg_clusts$cluster,weights=W)

# modularity function doesn't accept the '0' cluster name, so we shift all
# membership values up by one to get the modularity
mod$spinglass <- modularity(mainG,sc_clusts$cluster+1,weights=W)

mod$eigen <- modularity(mainG,lev_clusts$cluster,weights=W)
mod$walktrap <- modularity(mainG,wt_clusts$cluster,weights=W)
mod$label <- modularity(mainG,lp_clusts$cluster,weights=W)
mod$louvain <- modularity(mainG,louv_clusts$cluster,weights=W)
mod$small_louvain <- modularity(mainG,louv_small_clusts$cluster,weights=W)
mod$infomap <- modularity(mainG,info_clusts$cluster,weights=W)
mod$edge_between <- modularity(mainG,eb_clusts$cluster,weights=W)

### write modularity data out to a file

outfile="~/Github/KIN_ClusteringWithAnnotations/results/weighted/clustering_modularity_results.txt"
write.table(mod,outfile,quote=FALSE,sep="\t",row.names = FALSE)
