
#' Evaluate individual nodes for their value in invasion detection (invasion starting from a single node) in one realization
#'
#' For a given introduction node (only one at this point), yields two outputs for one realization.  outmat has rows=time steps, columns = nodes, and entries = invasion status for each node at each time step (1 = invaded, 0 = not invaded).  sampnodes has rows = nodes, first column = time step at which invasion is detected at a node (Inf if node is never reached), second column = the number of nodes invaded at the time of detection (the total number that are ever invaded if the node is never reached), third column = the number of nodes not invaded at the time of detection.
#' @param adjmat adjacency matrix to be evaluated; if stoch=T, the entries are the probabilities that a link exists in any given realization
#' @param start.choice number (ID) of node where invasion starts
#' @param stoch logical var indicating whether adjacency matrix entries are fixed or probabilities (default is T)
#' @keywords prioritization sampling
#' @export

#' @examples
#' Amat <- matrix(c(1,0,0,0,1,1,0,0,0,1,1,0,1,1,1,1),nrow=4,ncol=4)
#' onestart(adjmat=Amat, start.choice=2, stoch=F) 
#' sAmat <- Amat * 0.7 # each potential link has probability 0.7 of existing in one realization
#' onestart(adjmat=sAmat, start.choice=2, stoch=T) 

onestart <- function(adjmat, start.choice, stoch){
   
  # the adjacency matrix is assumed to have link sources as rows, and link sinks as columns

  dimL <- dim(adjmat)[1] # number of rows of adjacency matrix
  t1 <- matrix(0 * 1:dimL, nrow=1)
  t1[,start.choice] <- 1 # starting node given status 1

  # if stochastic, generate the adjacency matrix for this realization - assumes the input adjacency matrix is made up of probabilities

    if(stoch){ 
    adjmat <- (matrix(runif(dimL^2), ncol=dimL) < adjmat)
  }

  # generate the infection status for each node at each time point until spread stops

  outmat <- t1
  infcount.pre <- 1
  infcount.post <- -99

  while(sum(t1) < dimL & infcount.pre != infcount.post){
    infcount.pre <- sum(t1) 
    t1 <- as.numeric(t1 %*% adjmat > 0)
    infcount.post <- sum(t1)
    outmat <- rbind(outmat,t1) # row i is the ith time point
  }

  # assign row names 
  rownames(outmat) <- 1:(dim(outmat)[1])

  # find the num infected for each individual sampling node, 
     # assuming invaded node stays invaded

  # find the number of nodes infected at each time
  inft <- rowSums(outmat)

  # find the first time each node is invaded, and number of nodes invaded at that time
  firstt <- -99 + 0*1:dimL # first time each node is invaded
  numinvt <- firstt # number of nodes invaded at that time

  for (i in 1:dimL){ # col i is the ith sample node

    if (sum(outmat[,i] > 0)){ # if node i ever gets invaded
      firstt[i] <- min(which(outmat[,i] > 0)) # first time
      numinvt[i] <- inft[firstt[i]] # how many nodes at that time
    }
    else {
      firstt[i] <- Inf
      numinvt[i] <- max(inft)
    }

  }

  # for each node, time first invaded and number of nodes invaded at that time, and number of nodes not invaded at that time
  numnotinvt <- dimL-numinvt
  sampnodes <- cbind(firstt, numinvt, numnotinvt)

  list(outmat=outmat, sampnodes=sampnodes)
}

