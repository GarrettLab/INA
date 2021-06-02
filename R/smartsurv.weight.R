

#' Evaluate nodes for value in invasion detection (based on a specified number of realizations, and including specified probabilities associated with each potential starting node)
#'
#' Uses weights that indicate how likely each node was to be the starting location (introduction point) for an invasion.  Uses output from smartsurv.  Finds the rate at which sampling success (other nodes free of invasion when invasion reaches sampling node) is expected to occur for each potential sampling node.

#' @param ss.out output object from function smartsurv
#' @param adjmat adjacency matrix for evaluation
#' @param wtvec vector of weights indicating the probability that each node would be the starting node for invasion
#' @param nodenam vector of node names

#' @keywords prioritization sampling weights
#' @export

#' @examples
#' Amat <- matrix(c(1,0,0,0,1,1,0,0,0,1,1,0,1,1,1,1),nrow=4,ncol=4)
#' smartsurv(adjmat=Amat, stoch=F, nrealz=1)
#' sAmat <- Amat * 0.7
#' ss.outex <- smartsurv(adjmat=sAmat, stoch=T, nrealz=10)
#' wtvec.ex <- c(1,10,100,1)
#' smartsurv.weight(ss.out=ss.outex, adjmat=sAmat, wtvec=wtvec.ex)
#' smartsurv.weight(ss.out=ss.outex, adjmat=sAmat, wtvec=wtvec.ex, nodenam=c("KS","NE","ND","SD"))
#' smartsurv.weight(ss.out=ss.outex, adjmat=sAmat, wtvec=c(1,100,10,1), nodenam=c("KS","NE","ND","SD"))



smartsurv.weight <- function(ss.out, adjmat, wtvec, nodenam=NA) {

  dimL <- dim(adjmat)[1]

  matop <- ss.out$meanarr
  wtarr <- wtvec * matop

  # find invasion free rate for each sampling node
  sampfree <- colSums(wtarr)
  tsampfree <- data.frame(1:dimL, sampfree, nodenam)
  
  list(wtarr=wtarr, tsampfree=tsampfree)

}

