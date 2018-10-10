#!/usr/bin/env Rscript

#install.packages("ProNet")
#install.packages("igraph")
#install.packages('purrr')

#args<-commandArgs(TRUE) #args[1]==input network /path/to/file, args[2]==yes/no weighted, args[3]==output subnetworks path/to/file
library("igraph")
#library("ProNet")
#G <- read.graph("/Users/lutzka/Desktop/testgraph.txt",format="ncol",names=TRUE,weights="yes",directed=FALSE)
#G <- read.graph("/Users/lutzka/Desktop/Kinome+OneHopCompiled/All_Databases_Combined_Gene_Names_one_hop_to_kinases.txt",format="ncol",names=TRUE,weights="no",directed=FALSE)
G <- read.graph("~/Github/subnetclustering/data/Full_Kinome_Network_Compiled_weighted_pathways.txt",format="ncol",names=TRUE,weights="yes",directed=FALSE)
#G <- read.graph(args[1],format="ncol",names=TRUE,weights=args[2],directed=FALSE)
#dd <- read.table("/Users/lutzka/Dropbox/GomezLabShare/Kinome/NetworkData/All_Databases_Combined_Gene_Names_kinome_only.txt")
#G <- graph.data.frame(dd,directed=FALSE)
tmp <- read.table("~/Github/subnetclustering/data/Full_Kinome_Network_Compiled_weighted_pathways.txt")
W <- tmp$V3

##Algorithms to try:
###1. fastgreedy.community (deterministic) -- run 1x
###2. spinglass.community (Potts model--statistical) -- run 1000x
###3. leading.eigenvector.community (top-down hierarchical) -- run 1000x
###4. label.propogation.community (relies on initial random seeds) -- run 1000x
###5. walktrap.community (random walks) -- run 1000x

##Algorithms NOT to try:
###1. edge.betweenness.community -- size of the network makes this one not very feasible


V(G)$comp <- components(G)$membership
mainG <- induced_subgraph(G,V(G)$comp==1)

###1. fastgreedy.community
fc <- fastgreedy.community(mainG)


###2. spinglass.community
sc <- spinglass.community(mainG, spins=100)

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
      local_votes[i,] <- (sc$membership == sc$membership[i])
    }
  }
  return(local_votes)
}

para_spinglass <- function(g, numnodes, numiter = 1000, spins = 100, cores=max(detectCores()-1,1)){
  # create a list of results and a list of parameters for the spinglasses
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

compiled.votes <- data.frame(names=sc$names, cluster=groups)
fggroups <- data.frame(names=fc$names, cluster=fc$membership)
write.table(compiled.votes, '~/Github/subnetclustering/reproduced/consensusclusters_spinglass_greaterthan90percent.txt',quote=FALSE,sep="\t",row.names=FALSE)
write.table(fggroups, '~/Github/subnetclustering/reproduced/fastgreedy_clusters.txt',quote=FALSE,sep="\t",row.names=FALSE)


#mcode clustering
#mcodeclust <- mcode(mainG, vwp = 0.5, haircut = TRUE, fluff = FALSE, fdt = 0.9, loops = TRUE) #mcodeclust$COMPLEX[[i]] -- these are the lists of nodes (index numbers, not names) in each cluster
#mcodegroups <- mat.or.vec(numnodes,1)
#for (i in 1:length(mcodeclust$COMPLEX)){
#  mcodegroups[mcodeclust$COMPLEX[[i]]] <- i
#}

#mcgroups <- data.frame(names=sc$names, cluster=mcodegroups)
#write.table(mcgroups, '/Users/lutzka/Dropbox/Analyses/DrugComboMethod/Results/mcode_clusters.txt',quote=FALSE,sep="\t",row.names=FALSE)


###3. leading.eigenvector.community
lec <- leading.eigenvector.community(mainG)

numnodes <- length(lec$names)
votes <- mat.or.vec(numnodes,numnodes)
numiter <-20

for (k in 1:numiter){
  lec <- leading.eigenvector.community(mainG)
  for (i in 1:numnodes){
    for (j in 1:numnodes) {
      if (lec$membership[i] == lec$membership[j]) {
        votes[i,j] <- votes[i,j] + 1 #increase the count for this pair by 1
      }
    }
  }
}

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

lec.compiled.votes <- data.frame(names=lec$names, cluster=groups)
write.table(lec.compiled.votes, '~/Lab/subnetclustering/reproduced/consensusclusters_leadingeigenvector_greaterthan90percent.txt',quote=FALSE,sep="\t",row.names=FALSE)



###4. label.propagation.community
lp <- label.propagation.community(mainG)

numnodes <- length(lp$names)
votes <- mat.or.vec(numnodes,numnodes)
numiter <- 1000

for (k in 1:numiter){
  lp <- label.propagation.community(mainG)
  for (i in 1:numnodes){
    for (j in 1:numnodes) {
      if (lp$membership[i] == lp$membership[j]) {
        votes[i,j] <- votes[i,j] + 1 #increase the count for this pair by 1
      }
    }
  }
}

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

lp.compiled.votes <- data.frame(names=lp$names, cluster=groups)
write.table(lp.compiled.votes, '/Users/lutzka/Dropbox/Analyses/DrugComboMethod/Results/consensusclusters_labelpropogation_greaterthan90percent.txt',quote=FALSE,sep="\t",row.names=FALSE)


###5. walktrap.community
wt <- walktrap.community(mainG, modularity=TRUE)
#wmemb <- cutat(wt,steps=which.max(wt$modularity)-1)

numnodes <- length(wt$names)
votes <- mat.or.vec(numnodes,numnodes)
numiter <- 1000

for (k in 1:numiter){
  wt <- walktrap.community(mainG, modularity=TRUE)
  for (i in 1:numnodes){
    for (j in 1:numnodes) {
      if (wt$membership[i] == wt$membership[j]) {
        votes[i,j] <- votes[i,j] + 1 #increase the count for this pair by 1
      }
    }
  }
}

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

wt.compiled.votes <- data.frame(names=wt.$names, cluster=groups)
write.table(wt.compiled.votes, '/Users/lutzka/Dropbox/Analyses/DrugComboMethod/Results/consensusclusters_walktrap_greaterthan90percent.txt',quote=FALSE,sep="\t",row.names=FALSE)


### write output file with memberships for each community detection method
outfile="/Users/lutzka/Dropbox/Analyses/DrugComboMethod/Results/Kinase_all_clustering_membership_results.txt"
compiled.data = data.frame(KinaseNames=fc$names,
                           FastGreedyCommunity=fc$membership,
                           SpinGlassCommunity=sc$membership,
                           WalkTrapCommunity=wt$membership,
                           LeadingEigenvectorCommunity=lec$membership,
                           LabelPropogationCommunity=lp$membership,
                           MCODEcommunity=mcgroups$cluster)
write.table(compiled.data,outfile,quote=FALSE,sep="\t",row.names=FALSE)



mod <- list()
mod$fc <- modularity(mainG,fc$membership,weights=W)
mod$sc <- modularity(mainG,sc$membership,weights=W)
mod$mcode <- modularity(mainG,mcgroups,weights=W)
mod$lec <- modularity(mainG,lec$membership,weights=W)
mod$lp <- modularity(mainG,lp$membership,weights=W)
mod$wt <- modularity(mainG,wt$membership,weights=W)

outfile="~/Dropbox/Analyses/DrugComboMethod/Results/Kinase_all_clustering_modularity_results.txt"
write.table(mod,outfile,quote=FALSE,sep="\t",row.names = FALSE)

