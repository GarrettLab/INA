
#' Generate vector of status after one spread step
#'
#' This function generates a vector of status after one spread step.  It assumes that the input adjacency matrix is composed of 0s and 1s, and that the diagonal is all 1s.
#'
#' Updated 2020-09-05

#' @param amat the adjacency matrix, composed of 0s and 1s

#' @param vect1 status of nodes before spread, in terms of presence or absence of bioentity or information about management (a matrix with nrow=1 and ncol=number of nodes)

#' @keywords dispersal
#' @export 

#' @examples
#' x1 <- spreadstep(amat=matrix(c(1,1,0,0,1,1,0,0,1),ncol=3,byrow=T), vect1=matrix(c(1,0,0), ncol=3))

#' x2 <- spreadstep(amat=matrix(c(1,1,0,0,1,1,0,0,1),ncol=3,byrow=T), vect1=matrix(c(1,0,0), ncol=3))

#' x3 <- spreadstep(amat=matrix(c(1,1,0,0,1,1,0,0,1),ncol=3,byrow=T), vect1=matrix(c(1,1,0), ncol=3))

#' x9 <- spreadstep(amat=matrix(c(0,1,0,0,0,0, 0,0,1,0,0,0, 0,0,0,1,0,0, 0,0,0,0,1,0, 0,0,0,0,0,1, 1,0,0,0,0,0),byrow=T,ncol=6), vect1=c(1,1,1,0,0,0))


spreadstep <- function(amat, vect1) {

  diag(amat) <- 1

  vect2 <- vect1 %*% amat

  vect2 <- as.numeric(vect2 >= 1)
    
  vect2 # the resulting status of nodes 
}