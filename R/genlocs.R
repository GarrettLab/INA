
#' Generate geographic locations of nodes in network
#'
#' This function generates the geographic locations of nodes in network, for simulation studies where geographic locations are not observed.
#'
#' Updated 2020-09-05

#' @param xrange4 range of x coordinates
#' @param yrange4 range of y coordinates
#' @param numnodes4 the number of nodes
#' @param randgeo4 if TRUE then locations are randomly generated (other options for distributing locations will be added)
#' @keywords geography locations
#' @export

#' @examples
#' geolocs99 <- genlocs(xrange4=c(0,50), yrange4=c(0,50), numnodes4=100, randgeo4=TRUE)


genlocs <- function(xrange4, yrange4, numnodes4, randgeo4=TRUE){

  if (randgeo4) {

    xvec <- xrange4[1] + (xrange4[2] - xrange4[1])*runif(n=numnodes4)

    yvec <- yrange4[1] + (yrange4[2] - yrange4[1])*runif(n=numnodes4)

  }

  xyob <- cbind(xvec,yvec)
  xyob
}