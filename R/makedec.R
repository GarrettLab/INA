
#' Decisions at each node about whether to adopt management
#'
#' This function generates decisions at each node about management
#'
#' Updated 2020-09-05

#' @param probadoptmean4 mean probability of adopting management if informed
#' @param probadoptsd4 sd for drawing probability of adoption, in truncated normal distribution
#' @param comvec vector of 1=info is present, 0=info is not present
#' @param probadoptvec4 vector of probabilities of adoption if informed for nodes in the network (\code{readprobadoptvec4} determines whether \code{probadoptvec4} is read in to \code{makedec} or generated within \code{makedec} based on \code{probadoptmean4} and \code{probadoptsd4 })
#' @param plotmp if T, then map of decision is plotted
#' @param readprobadoptvec4 if T, then a VECTOR of probadoptvec4 values is read in - if F, \code{probadoptmean4} and \code{probadoptsd4 } are used to generate the vector
#' @param geocoords4n matrix of x,y coordinates for the nodes, for mapping if \code{plotmp} = T (and for determining the number of nodes in current version)
#' @keywords decisions
#' @export
#' @import truncnorm

#' @examples

#' x5 <- makedec(comvec=c(1,1,1,0,0,0), geocoords4n=matrix(c(1,1, 1,2, 1,3, 2,1, 2,2, 2,3), byrow=T,ncol=2), readprobadoptvec4=F, probadoptmean4=0.1, probadoptsd4 =0.1, plotmp=T)

#' x6 <- makedec(comvec=c(1,1,1,0,0,0), geocoords4n=matrix(c(1,1, 1,2, 1,3, 2,1, 2,2, 2,3), byrow=T,ncol=2), readprobadoptvec4=F, probadoptmean4=0.5, probadoptsd4 =0.1, plotmp=T)

#' x7 <- makedec(comvec=c(1,1,1,0,0,0), geocoords4n=matrix(c(1,1, 1,2, 1,3, 2,1, 2,2, 2,3), byrow=T,ncol=2), readprobadoptvec4=T, probadoptvec4=c(0.1,0.1,0.9,0.1,0.2,0.3), plotmp=T)


makedec <- function(comvec, geocoords4n, probadoptmean4, probadoptsd4 , readprobadoptvec4, probadoptvec4, plotmp=F){

  nodlen <- dim(geocoords4n)[1] # number of nodes
 
############################## below repeated in ntsteps2???


  if(!readprobadoptvec4){
    probadoptvec4 <- truncnorm::rtruncnorm(n=nodlen, a=0, b=1, mean = probadoptmean4, sd=probadoptsd4 )
  }
  
  stochadopt <- runif(n=nodlen, min=0, max=1)
  hypdecvec <- stochadopt < probadoptvec4 # would adopt if informed
  decvec <- hypdecvec & comvec # does adopt

  if(plotmp){
    plot(geocoords4n, main='Decision to adopt, shaded = yes')
    if(!is.null(dim(geocoords4n[decvec==1,]))){
      points(geocoords4n[decvec==1,], pch=16)
    } 
    else if(is.null(dim(geocoords4n[decvec==1,]))){
      points(geocoords4n[decvec==1,][1], geocoords4n[decvec==1,][2], pch=16)
    }
  }

  decvec
}