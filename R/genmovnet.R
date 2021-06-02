
#' Generate network adjacency matrix for movement
#'
#' This function generates an adjacency matrix for movement, assumed symmetric (in this version).  It is used by functions including \code{INAscene}. The movement adjacency matrix is composed of 1s and 0s only if lktype="pa" option is used
#'
#' Updated 2020-09-05

#' @param amdist4 the function of distance used to estimate movement probability - 'random' (not related to distance) or 'powerlaw' (inverse power law) or 'exp' (negative exponential, to be added)
#' @param iplot if T, generates igraph plot of adjacency matrix
#' @param ampla4 inverse power law parameter a in ad^(-b)
#' @param amplb4 inverse power law parameter b in ad^(-b)
#' @param amrandp4 random matrix with entries binomial with probability p
#' @param geocoords4n the matrix of xy coordinates for node locations, used when the probability of a link is a function of distance (note that the distance between each pair of locations is assumed to be greater than 1)
#' @keywords dispersal
#' @export 
#' @import igraph

#' @examples

#' x1 <- genmovnet(j <- genlocs(xrange4=c(0,50), yrange4=c(0,50), numnodes4=50, randgeo4=TRUE), amdist4='random', amrandp4=0.01, iplot=T)

#' x2 <- genmovnet(j <- genlocs(xrange4=c(0,50), yrange4=c(0,50), numnodes4=100, randgeo4=TRUE), amdist4='random', amrandp4=0.02, iplot=T)

#' x7 <- genmovnet(geocoords4n=matrix(c(1,1, 1,2, 1,3, 2,1, 2,2, 2,3),ncol=2,byrow=T), amdist4='powerlaw', ampla4=2, amplb4=1, iplot=T)

#' x8 <- genmovnet(j <- genlocs(numnodes4=30, xrange4 = c(0, 10), yrange4 = c(0, 10)), amdist4='powerlaw', ampla4=2, amplb4=1, iplot=T)

#' x9 <- genmovnet(j <- genlocs(numnodes4=300, xrange4 = c(0, 10), yrange4 = c(0, 100)), amdist4='powerlaw', ampla4=2, amplb4=1, iplot=T)


genmovnet <- function(geocoords4n, amdist4, iplot=F, amrandp4, ampla4, amplb4){

  dimam <- dim(geocoords4n)[1]

  if (amdist4 == 'powerlaw') { # ad^(-b)

    tdist <- as.matrix(dist(geocoords4n, method = "euclidean", diag=T, upper=T))
    linkmat <- ampla4*tdist^(-amplb4)  
  }

  else if(amdist4 == 'random'){
  
    linkmat <- matrix(rbinom(n=dimam*dimam, size=1, prob=amrandp4), nrow=dimam)

    # Make the matrix symmetric

    linkmat[lower.tri(linkmat)] <- t(linkmat)[lower.tri(linkmat)]
    
  }

  
  # the diagonal is set to 1, so there are always self-loops

  # whether the spread to a location is maintained over time 
  #   is determined elsewhere

  diag(linkmat) <- 1  


  if(iplot){ 
    linkmat.diag0 <- linkmat  
    # make diagonal of temporary adjacency matrix 0 for
    #   clearer visualization when plotting  
    diag(linkmat.diag0) <- 0


    if (amdist4 == 'powerlaw') {
      linkmat.3 <- linkmat.diag0 > 0.3
      linkmati <- igraph::graph.adjacency(linkmat.3)
      plot(linkmati, edge.arrow.size=0, vertex.color='skyblue', main='Links with p > 0.3')
    }
    else if(amdist4 == 'random'){
      linkmati <- igraph::graph.adjacency(linkmat.diag0)
      plot(linkmati, edge.arrow.size=0, vertex.color='skyblue')
    }
    
  }

  linkmat
}