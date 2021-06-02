  
#' Evaluate nodes for value in invasion detection (based on a specified number of realizations)
#'
#' Evaluate nodes for utility in invasion detection (based on a specified number of realizations) using the function multistart (which in turn uses the function onestart).  outarr saves the output from each realization from multistart, in a 3D array.  In meanarr, rows = starting nodes (introduction nodes), columns = sampling nodes, and entries = mean nodes free from invasion when invasion starting at the introduction node reaches the sampling node.  In vararr, for the same rows and columns, entries are the variance in the number of nodes free from invasion.  Note that nrealz=1 is all that is needed if there is no stochastic component.

#' @param adjmat adjacency matrix for evaluation
#' @param stoch logical var indicating whether adjacency matrix entries are fixed or probabilities
#' @param nrealz number of realizations to be analyzed (just 1 is needed if stoch=F)

#' @keywords prioritization sampling
#' @export
#' @examples
#' Amat <- matrix(c(1,0,0,0,1,1,0,0,0,1,1,0,1,1,1,1),nrow=4,ncol=4)
#' smartsurv(adjmat=Amat, stoch=F, nrealz=1)
#' sAmat <- Amat * 0.7
#' smartsurv(adjmat=sAmat, stoch=T, nrealz=10)



smartsurv <- function(adjmat, stoch, nrealz=1) {

  dimL <- dim(adjmat)[1]

  outarr <- array(-99, c(dimL,dimL,nrealz)) # all results
  meanarr <- matrix(-99, ncol=dimL, nrow=dimL) # mean of results
  vararr <- meanarr # variance of results

  for (i3 in 1:nrealz){
    outarr[,,i3]  <- multistart(adjmat=adjmat, stoch=stoch)
  }

  for (i3 in 1:dimL) {
    for (i4 in 1:dimL) {
      meanarr[i3,i4] <- mean(outarr[i3,i4,]) #mean across realizations
      vararr[i3,i4] <- var(outarr[i3,i4,]) #variance across realizations
  } }

  list(outarr=outarr, meanarr=meanarr, vararr=vararr)

}


