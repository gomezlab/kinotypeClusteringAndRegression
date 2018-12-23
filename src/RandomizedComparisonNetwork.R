library("igraph")
G <- read.graph("~/Github/KIN_ClusteringWithAnnotations/data/KIN_weighted_edges.txt",format="ncol",names=TRUE,weights="yes",directed=FALSE)
tmp <- read.table("~/Github/KIN_ClusteringWithAnnotations/data/KIN_weighted_edges.txt")
W <- tmp$V3
new_G <- read.graph('~/Github/KIN_ClusteringWithAnnotations/data/KIN_random_graph.txt',names=TRUE,format='ncol',directed=FALSE)
d <- degree(G)
new_G <- degree.sequence.game(out.deg = d, method=c("vl"))
V(new_G)$comp <- components(new_G)$membership
V(new_G)$name <- V(G)$name
E(new_G)$weight <- sample(W)

mainG <- induced_subgraph(new_G,V(new_G)$comp==1)

###1. fastgreedy.community
fc <- fastgreedy.community(mainG)
fg_clusts <- data.frame(names=fc$names, cluster=fc$membership)
write.table(fg_clusts, '~/GitHub/KIN_ClusteringWithAnnotations/randomizedResults/fastgreedy_clusters.txt',quote=FALSE,sep="\t",row.names=FALSE)

###2. spinglass.community - knows to use in weights
sc <- spinglass.community(mainG, spins=100)

# create a votes matrix to store individual cluster tallies
numnodes <- length(sc$names)
votes <- mat.or.vec(numnodes,numnodes)
numiter <- 1000
dim(votes)

## do some parallelization
library('parallel')
library('purrr')

single_spinglass <- function(numiter, g, numnodes, spins = 100){
  local_votes <- mat.or.vec(numnodes, numnodes)
  for (k in 1:numiter){
    sc <- igraph::spinglass.community(g, spins=100)
    for (i in 1:numnodes){
      local_votes[i,] = local_votes[i,] + (sc$membership == sc$membership[i])
    }
  }
  return(local_votes)
}

para_spinglass <- function(g, numnodes, numiter = 1000, spins = 100, cores=max(detectCores()-1,1)){
  # create a list of randomizedResults and a list of parameters for the spinglasses
  numiter_params = c(rep(numiter %/% cores, cores))
  
  # add one to the number of iterations per core if the total number of iterations
  # is not divisible by the number of cores
  if(numiter %% cores > 0){
    numiter_params[1:(numiter%%cores)] = numiter_params[1:(numiter%%cores)] + 1
  }
  print(numiter_params)
  
  # start up a cluster
  votes <- Reduce('+', parallel::mclapply(X = numiter_params, FUN = single_spinglass, g=mainG, numnodes=numnodes, mc.cores = cores))
  
  return(votes)
}

print(Sys.time())
votes <- para_spinglass(g=mainG, numnodes = numnodes, numiter = 1000)
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

sum()

sc_clusts <- data.frame(names=sc$names, cluster=groups)
write.table(sc_clusts, '~/Github/KIN_ClusteringWithAnnotations/randomizedResults/consensus_spinglass.txt',quote=FALSE,sep="\t",row.names=FALSE)

###3. leading.eigenvector.community -- knows to use weights
lev <- leading.eigenvector.community(mainG)

numnodes <- length(lev$names)
votes <- mat.or.vec(numnodes,numnodes)
numiter <- 1000

single_lev <- function(numiter, g, numnodes){
  local_votes <- mat.or.vec(numnodes, numnodes)
  for (k in 1:numiter){
    lev <- igraph::leading.eigenvector.community(g)
    for (i in 1:numnodes){
      local_votes[i,] = local_votes[i,] + (lev$membership == lev$membership[i])
    }
  }
  return(local_votes)
}

para_lev <- function(g, numnodes, numiter = 1000, cores=max(detectCores()-1,1)){
  # create a list of randomizedResults and a list of parameters for the spinglasses
  numiter_params = c(rep(numiter %/% cores, cores))
  
  # add one to the number of iterations per core if the total number of iterations
  # is not divisible by the number of cores
  if(numiter %% cores > 0){
    numiter_params[1:(numiter%%cores)] = numiter_params[1:(numiter%%cores)] + 1
  }
  print(numiter_params)
  
  # start up a cluster
  votes <- Reduce('+', parallel::mclapply(X = numiter_params, FUN = single_lev, g=mainG, numnodes=numnodes, mc.cores = cores))
  
  return(votes)
}

print(Sys.time())
votes <- para_lev(g=mainG, numnodes = numnodes, numiter = 1000)
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
write.table(lev_clusts, '~/Github/KIN_ClusteringWithAnnotations/randomizedResults/consensus_eigenvector.txt',quote=FALSE,sep="\t",row.names=FALSE)

###4. label.propagation.community
lp <- label.propagation.community(mainG)

numnodes <- length(lp$names)
votes <- mat.or.vec(numnodes,numnodes)
numiter <- 1000

single_lp <- function(numiter, g, numnodes){
  local_votes <- mat.or.vec(numnodes, numnodes)
  for (k in 1:numiter){
    lp <- igraph::label.propagation.community(g)
    for (i in 1:numnodes){
      local_votes[i,] = local_votes[i,] + (lp$membership == lp$membership[i])
    }
  }
  return(local_votes)
}

para_lp <- function(g, numnodes, numiter = 1000, cores=max(detectCores()-1,1)){
  # create a list of randomizedResults and a list of parameters for the spinglasses
  numiter_params = c(rep(numiter %/% cores, cores))
  
  # add one to the number of iterations per core if the total number of iterations
  # is not divisible by the number of cores
  if(numiter %% cores > 0){
    numiter_params[1:(numiter%%cores)] = numiter_params[1:(numiter%%cores)] + 1
  }
  print(numiter_params)
  
  # start up a cluster
  votes <- Reduce('+', parallel::mclapply(X = numiter_params, FUN = single_lp, g=mainG, numnodes=numnodes, mc.cores = cores))
  
  return(votes)
}

print(Sys.time())
votes <- para_lp(g=mainG, numnodes = numnodes, numiter = 1000)
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
write.table(lp_clusts, '~/GitHub/KIN_ClusteringWithAnnotations/randomizedResults/consensus_label_propagation.txt',quote=FALSE,sep="\t",row.names=FALSE)

###5. walktrap.community
wt <- walktrap.community(mainG, modularity=TRUE)

numnodes <- length(wt$names)
votes <- mat.or.vec(numnodes,numnodes)
numiter <- 1000
single_wt <- function(numiter, g, numnodes){
  local_votes <- mat.or.vec(numnodes, numnodes)
  for (k in 1:numiter){
    wt <- igraph::walktrap.community(g, modularity=TRUE)
    for (i in 1:numnodes){
      local_votes[i,] = local_votes[i,] + (wt$membership == wt$membership[i])
    }
  }
  return(local_votes)
}

para_wt <- function(g, numnodes, numiter = 1000, cores=max(detectCores()-1,1)){
  # create a list of randomizedResults and a list of parameters for the spinglasses
  numiter_params = c(rep(numiter %/% cores, cores))
  
  # add one to the number of iterations per core if the total number of iterations
  # is not divisible by the number of cores
  if(numiter %% cores > 0){
    numiter_params[1:(numiter%%cores)] = numiter_params[1:(numiter%%cores)] + 1
  }
  print(numiter_params)
  
  # start up a cluster
  votes <- Reduce('+', parallel::mclapply(X = numiter_params, FUN = single_wt, g=mainG, numnodes=numnodes, mc.cores = cores))
  
  return(votes)
}

print(Sys.time())
votes <- para_wt(g=mainG, numnodes = numnodes, numiter = 1000)
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
write.table(wt_clusts, '~/Github/KIN_ClusteringWithAnnotations/randomizedResults/consensus_walktrap.txt',quote=FALSE,sep="\t",row.names=FALSE)

###6. cluster_louvain
louv <- cluster_louvain(mainG)
louv_clusts <- data.frame(names=louv$names, cluster=louv$memberships[2,])
louv_small_clusts <- data.frame(names=louv$names, cluster=louv$memberships[1,])
write.table(louv_clusts, '~/GitHub/KIN_ClusteringWithAnnotations/randomizedResults/louvain_clusters.txt',quote=FALSE,sep="\t",row.names=FALSE)
write.table(louv_small_clusts, '~/GitHub/KIN_ClusteringWithAnnotations/randomizedResults/louvain_small_clusters.txt',quote=FALSE,sep="\t",row.names=FALSE)

###7. cluster_infomap (random walks) 
info <- walktrap.community(mainG)

numnodes <- length(info$names)
votes <- mat.or.vec(numnodes,numnodes)
numiter <- 1000

single_info <- function(numiter, g, numnodes){
  local_votes <- mat.or.vec(numnodes, numnodes)
  for (k in 1:numiter){
    info <- igraph::cluster_infomap(g)
    for (i in 1:numnodes){
      local_votes[i,] = local_votes[i,] + (info$membership == info$membership[i])
    }
  }
  return(local_votes)
}

para_info <- function(g, numnodes, numiter = 1000, cores=max(detectCores()-1,1)){
  # create a list of randomizedResults and a list of parameters for the spinglasses
  numiter_params = c(rep(numiter %/% cores, cores))
  
  # add one to the number of iterations per core if the total number of iterations
  # is not divisible by the number of cores
  if(numiter %% cores > 0){
    numiter_params[1:(numiter%%cores)] = numiter_params[1:(numiter%%cores)] + 1
  }
  print(numiter_params)
  
  # start up a cluster
  votes <- Reduce('+', parallel::mclapply(X = numiter_params, FUN = single_info, g=mainG, numnodes=numnodes, mc.cores = cores))
  
  return(votes)
}

print(Sys.time())
votes <- para_info(g=mainG, numnodes = numnodes, numiter = 1000)
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
write.table(info_clusts, '~/Github/KIN_ClusteringWithAnnotations/randomizedResults/consensus_infomap.txt',quote=FALSE,sep="\t",row.names=FALSE)

###8. edge.betweenness.community
eb <- edge.betweenness.community(mainG)
eb_clusts <- data.frame(names=eb$names, cluster=eb$membership)
write.table(eb_clusts, '~/GitHub/KIN_ClusteringWithAnnotations/randomizedResults/edge_betweenness_community_clusters.txt',quote=FALSE,sep="\t",row.names=FALSE)

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

outfile="~/Github/KIN_ClusteringWithAnnotations/randomizedResults/clustering_modularity_results.txt"
write.table(mod,outfile,quote=FALSE,sep="\t",row.names = FALSE)


