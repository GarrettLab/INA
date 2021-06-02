

#' Evaluate nodes for value in invasion detection (invasion starting from each node in turn)
#'
#' Evaluate nodes for value in invasion detection (invasion starting from each node in turn) using function onestart for one realization.  The output is a matrix with rows = starting nodes (introduction nodes), columns = samping nodes, and entries = number of nodes not invaded by time detected at the sampling node having started at the introduction node.  
#' @param adjmat adjacency matrix for evaluation
#' @param stoch logical var indicating whether adjacency matrix entries are fixed (1 or 0) or probabilities
#' @keywords prioritization sampling
#' @export

#' @examples
#' Amat <- matrix(c(1,0,0,0,1,1,0,0,0,1,1,0,1,1,1,1),nrow=4,ncol=4)
#' multistart(adjmat=Amat, stoch=F)
#' sAmat <- Amat * 0.7 # each potential link has probability 0.7 of existing in one realization
#' multistart(adjmat=sAmat, stoch=T)

multistart <- function(adjmat, stoch){

  dimL <- dim(adjmat)[1]

  alloutmat <- matrix(-99, ncol=dimL, nrow=dimL)

  # Changed to save 3rd col, number of nodes _not_ invaded

  for (i in 1:dimL) {
    alloutmat[i,] <- onestart(adjmat=adjmat, start.choice=i, stoch=stoch)$sampnodes[,3] 
  }

alloutmat
}


