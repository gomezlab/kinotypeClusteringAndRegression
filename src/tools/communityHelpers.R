### helper functions for community detection

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
  # create a list of results and a list of parameters for the spinglasses
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
  # create a list of results and a list of parameters for the spinglasses
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

single_wt <- function(numiter, g, numnodes, steps=10){
  local_votes <- mat.or.vec(numnodes, numnodes)
  for (k in 1:numiter){
    wt <- igraph::walktrap.community(g, modularity=TRUE, steps=steps)
    for (i in 1:numnodes){
      local_votes[i,] = local_votes[i,] + (wt$membership == wt$membership[i])
    }
  }
  return(local_votes)
}

para_wt <- function(g, numnodes, numiter = 1000, cores=max(detectCores()-1,1), steps=10){
  # create a list of results and a list of parameters for the spinglasses
  numiter_params = c(rep(numiter %/% cores, cores))

  # add one to the number of iterations per core if the total number of iterations
  # is not divisible by the number of cores
  if(numiter %% cores > 0){
    numiter_params[1:(numiter%%cores)] = numiter_params[1:(numiter%%cores)] + 1
  }
  print(numiter_params)

  # start up a cluster
  votes <- Reduce('+', parallel::mclapply(X = numiter_params, FUN = single_wt, g=mainG, numnodes=numnodes, mc.cores = cores, steps=steps))

  return(votes)
}

single_info <- function(numiter, g, numnodes, e.weights=NULL){
  local_votes <- mat.or.vec(numnodes, numnodes)
  for (k in 1:numiter){
    info <- igraph::cluster_infomap(g, e.weights=e.weights)
    for (i in 1:numnodes){
      local_votes[i,] = local_votes[i,] + (info$membership == info$membership[i])
    }
  }
  return(local_votes)
}

para_info <- function(g, numnodes, e.weights=NULL, numiter = 1000, cores=max(detectCores()-1,1)){
  # create a list of results and a list of parameters for the spinglasses
  numiter_params = c(rep(numiter %/% cores, cores))

  # add one to the number of iterations per core if the total number of iterations
  # is not divisible by the number of cores
  if(numiter %% cores > 0){
    numiter_params[1:(numiter%%cores)] = numiter_params[1:(numiter%%cores)] + 1
  }
  print(numiter_params)

  # start up a cluster
  votes <- Reduce('+', parallel::mclapply(X = numiter_params, FUN = single_info, g=mainG, e.weights=e.weights, numnodes=numnodes, mc.cores = cores))

  return(votes)
}
