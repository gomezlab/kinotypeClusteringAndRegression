install.packages("devtools")


#install and load ESSC
devtools::install_github("jdwilson4/ESSC")
library(ESSC, quietly = TRUE)


#load other required packages
library(Matrix, quietly = TRUE)
library(Rlab, quietly = TRUE)
library(devtools, quietly = TRUE)

# generate a random network
net <- stochastic.block(n = 1000, k = 3, P = cbind(c(0.1, 0.01, 0.01), c(0.01, 0.1, 0.01), c(0.01, 0.01, 0.1)), sizes = c(300, 300, 400))

#view the network
image(net$Adjacency)

# generate another random network
net2 <- stochastic.block(n = 1000, k = 3, P = cbind(c(0.1, 0.01, 0.01), c(0.01, 0.1, 0.01), c(0.01, 0.01, 0.1)), sizes = c(300, 300, 400), random.community.assignment = TRUE)

#view the network. notice that the community labels are assigned at random
image(net2$Adjacency)

#now let's rearrange the network according to the true community labels and view
image(net2$Adjacency[order(net2$Membership), order(net2$Membership)])

# try essc with a Poisson null
results.pois <- essc(net$Adjacency, alpha = 0.10, Null = "Poisson")
print(results.pois$Communities) # view communities

# try essc with a Binomial null
results.bin <- essc(net$Adjacency, alpha = 0.10, Null = "Binomial")
print(results.bin$Communities) # view communities

# read in the unweighted kinome network
unweighted_G <- read.graph("~/Github/subnetclustering/data/Full_Kinome_Network_Compiled_no_header.txt",format="ncol",names=TRUE,weights="no",directed=FALSE)
essc_adj <- as_adjacency_matrix(unweighted_G)

essc_results.pois <- essc(essc_adj, alpha = 0.05, Null = "Poisson")
print(essc_results.pois$Communities) # view communities

essc_results.bin <- essc(essc_adj, alpha = 0.05, Null = "Binomial")
print(essc_results.bin$Communities) # view communities

essc_results.bin$Background
essc_results.pois$Background

essc_results.pois$PValues
essc_results.bin$PValues
