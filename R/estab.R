
#' Bioentity establishment or not at each node
#'
#' This function determines whether a bioentity establishes or not at each node.  It is similar in some ways to the function \code{makedec}, which determines whether management is adopted.
#'
#' Updated 2020-09-05

#' @param decvec vector of 1=management, 0=no management
#' @param dispvec vector of 1=sp has dispersed to node, 0=sp has not dispersed to node
#' @param probestabmean4 mean probability of establishment (new or CONTINUED) in absence of management (used if \code{readprobestabvec4} =  F)
#' @param probestabsd4 sd in probability of establishment in absence of management (used if \code{readprobestabvec4} = F)

#' @param maneffdir4 if maneffdir4='decrease_estab', the management reduces the probability of establishment; if maneffdir4='increase_estab', the management reduces the probability that establishment does NOT occur

#' @param maneffmean4n mean effect of management (proportional change in estabp, where the combination of maneffmean4n = 1 and maneffsd4n = 0 makes establishment impossible)
#' @param maneffsd4n sd of management effect
#' @param probestabvec4 vector of probabilities of establishment (read in or generated when the function is run, depending on \code{readprobestabvec4} equal to T or F)
#' @param plotmp if T, then map of establishment is plotted
#' @param readprobestabvec4 if T, then a vector \code{probestabvec4} is read in; otherwise the vector is generated using \code{probestabmean4} and \code{probestabsd4}
#' @param geocoords4n matrix of x,y coordinates of nodes (used if \code{plotmp} = T)
#' @keywords establishment
#' @export
#' @import truncnorm

#' @examples
#' x6decrease <- estab(decvec=c(1,1,1,0,0,0), dispvec=c(1,1,1,1,1,0), geocoords4n=matrix(c(1,1, 1,2, 1,3, 2,1, 2,2, 2,3), byrow=T, ncol=2), readprobestabvec4=F, probestabmean4=0.5, probestabsd4=0.1, maneffdir4='decrease_estab', maneffmean4n=0.1, maneffsd4n=0.1, plotmp=T)

#' x7decrease <- estab(decvec=c(1,1,1,0,0,0), dispvec=c(1,1,1,1,1,0), geocoords4n=matrix(c(1,1, 1,2, 1,3, 2,1, 2,2, 2,3), byrow=T, ncol=2), readprobestabvec4=F, probestabmean4=0.5, probestabsd4=0.1, maneffdir4='decrease_estab', maneffmean4n=0.9, maneffsd4n=0.1, plotmp=T)

#' x8decrease <- estab(decvec=c(1,1,1,0,0,0), dispvec=c(1,1,1,1,1,0), geocoords4n=matrix(c(1,1, 1,2, 1,3, 2,1, 2,2, 2,3), byrow=T, ncol=2), readprobestabvec4=T, probestabvec4=c(0.9,0.1,0.1,0.1,0.1,0.9), maneffdir4='decrease_estab', maneffmean4n=0.9, maneffsd4n=0.1, plotmp=T)

#' x9decrease <- estab(decvec=c(1,1,1,0,0,0), dispvec=c(1,1,1,1,1,0), geocoords4n=matrix(c(1,1, 1,2, 1,3, 2,1, 2,2, 2,3), byrow=T, ncol=2), readprobestabvec4=T, probestabvec4=c(0.9,0.1,0.1,0.1,0.1,0.9), maneffdir4='decrease_estab', maneffmean4n=0.1, maneffsd4n=0.1, plotmp=T)

#' x6increase <- estab(decvec=c(1,1,1,0,0,0), dispvec=c(1,1,1,1,1,0), geocoords4n=matrix(c(1,1, 1,2, 1,3, 2,1, 2,2, 2,3), byrow=T, ncol=2), readprobestabvec4=F, probestabmean4=0.5, probestabsd4=0.1, maneffdir4='increase_estab', maneffmean4n=0.1, maneffsd4n=0.1, plotmp=T)

#' x7increase <- estab(decvec=c(1,1,1,0,0,0), dispvec=c(1,1,1,1,1,0), geocoords4n=matrix(c(1,1, 1,2, 1,3, 2,1, 2,2, 2,3), byrow=T, ncol=2), readprobestabvec4=F, probestabmean4=0.5, probestabsd4=0.1, maneffdir4='increase_estab', maneffmean4n=0.9, maneffsd4n=0.1, plotmp=T)

#' x8increase <- estab(decvec=c(1,1,1,0,0,0), dispvec=c(1,1,1,1,1,0), geocoords4n=matrix(c(1,1, 1,2, 1,3, 2,1, 2,2, 2,3), byrow=T, ncol=2), readprobestabvec4=T, probestabvec4=c(0.9,0.1,0.1,0.1,0.1,0.9), maneffdir4='increase_estab', maneffmean4n=0.9, maneffsd4n=0.1, plotmp=T)

#' x9increase <- estab(decvec=c(1,1,1,0,0,0), dispvec=c(1,1,1,1,1,0), geocoords4n=matrix(c(1,1, 1,2, 1,3, 2,1, 2,2, 2,3), byrow=T, ncol=2), readprobestabvec4=T, probestabvec4=c(0.9,0.1,0.1,0.1,0.1,0.9), maneffdir4='increase_estab', maneffmean4n=0.1, maneffsd4n=0.1, plotmp=T)




estab <- function(decvec, dispvec, probestabmean4, probestabsd4, maneffdir4, maneffmean4n, maneffsd4n, plotmp=F, geocoords4n, readprobestabvec4, probestabvec4){

  nodlen <- length(decvec)

############## ntsteps2 also does this?????

  # probability of establishment in absence of management
  if(!readprobestabvec4){
    probestabvec4 <- truncnorm::rtruncnorm(n=nodlen, a=0, b=1, mean=probestabmean4, sd=probestabsd4)
  }

  # how great would the management effect be where used
  obsman <- truncnorm::rtruncnorm(n=nodlen, a=0, b=1, mean=maneffmean4n, sd=maneffsd4n)

  ## probabiliy of estab adjusted for whether management used

  # if management reduces probability of establishing
  
  if (maneffdir4 == 'decrease_estab'){
    probestabvec4[decvec==1] <- (1-obsman)[decvec==1] * probestabvec4[decvec==1]

  # or if management reduces the probability of NOT establishing

  } else if (maneffdir4 == 'increase_estab') {
    probestabvec4[decvec==1] <- 1 - ((1-obsman)[decvec==1] * (1 - probestabvec4)[decvec==1])
  }

  # random draw used for comparison to probestabvec4
  stochestab <- runif(n=nodlen, min=0, max=1)

  # whether would establish if has dispersed
  estabvec <- stochestab < probestabvec4

  # has it also dispersed there?
  estabvec <- estabvec & dispvec

  if(plotmp){
    plot(geocoords4n, main='Establishment, shaded = yes')
    if(!is.null(dim(geocoords4n[estabvec==1,]))){
      points(geocoords4n[estabvec==1,], pch=16)
    } 
    else if(is.null(dim(geocoords4n[estabvec==1,]))){
      points(geocoords4n[estabvec==1,][1], geocoords4n[estabvec==1,][2], pch=16)
    }
  }

  estabvec
}