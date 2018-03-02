#!/usr/bin/env Rscript

#install.packages("ProNet")
#install.packages("igraph")
library("igraph")
library("ProNet")

#G <- read.graph("/Users/lutzka/Desktop/testgraph.txt",format="ncol",names=TRUE,weights="yes",directed=FALSE)
#G <- read.graph("/Users/lutzka/Desktop/Kinome+OneHopCompiled/All_Databases_Combined_Gene_Names_one_hop_to_kinases.txt",format="ncol",names=TRUE,weights="no",directed=FALSE)
G <- read.graph("/Users/lutzka/Dropbox/Analyses/DrugComboMethod/Data/GenericKinome/Full_Kinome_Network_Compiled_weighted_pathways.txt",format="ncol",names=TRUE,weights="yes",directed=FALSE)
#G <- read.graph(args[1],format="ncol",names=TRUE,weights=args[2],directed=FALSE)
#dd <- read.table("/Users/lutzka/Dropbox/GomezLabShare/Kinome/NetworkData/All_Databases_Combined_Gene_Names_kinome_only.txt")
#G <- graph.data.frame(dd,directed=FALSE)
tmp <- read.table("/Users/lutzka/Dropbox/Analyses/DrugComboMethod/Data/GenericKinome/Full_Kinome_Network_Compiled_weighted_pathways.txt")
W <- tmp$V3

outfile="/Users/lutzka/Dropbox/Analyses/DrugComboMethod/Results/Kinase_all_clustering_membership_results.txt"

#Analyze largest CSpG cluster to split
tmp <- read.table(outfile,header=T)
freq <- table(tmp$SpinGlassCommunity)
x <- as.numeric(names(freq)[which.max(freq)])
lgclust <- which(tmp$SpinGlassCommunity == x)

subG <- induced_subgraph(G, lgclust, impl="auto")

###2. spinglass.community
sc <- spinglass.community(subG, spins=100)

numnodes <- length(sc$names)
votes <- mat.or.vec(numnodes,numnodes)
numiter <- 1000

for (k in 1:numiter){
  print(paste("Starting iteration",k))
  sc <- spinglass.community(subG, spins=100)
  for (i in 1:numnodes){
    for (j in 1:numnodes) {
      if (sc$membership[i] == sc$membership[j]) {
        votes[i,j] <- votes[i,j] + 1 #increase the count for this pair by 1
      }
    }
  }
  print(paste("Finished with iteration",k))
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

compiled.votes <- data.frame(names=sc$names, cluster=groups)
write.table(compiled.votes, '/Users/lutzka/Dropbox/Analyses/DrugComboMethod/Results/consensusclusters_spinglass_greaterthan90percent_largeclustersplit.txt',quote=FALSE,sep="\t",row.names=FALSE)

