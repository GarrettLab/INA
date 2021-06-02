
#' Runs and summarizes multiple INA simulations with the same parameter values
#'
#' This function runs and summarizes multiple INA simulations with the same parameter values.  (It uses INA functions estinfo, genlocs, initvals, setup2, genmovnet, spreadstep, makedec, estab, and ntsteps2.)
#'
#' Updated 2020-09-05

#' @param nreals2 number of realizations to be evaluated
#' @param usethreshman2 if T, use the threshold for the management effect size estimate; thus, information is never present anywhere unless the management effect estimate exceeds the threshold
#' @param ntimesteps2 number of time steps to be evaluated

#' @param readseam2 if T, the communication adjacency matrix is read in rather than being generated from the outset
#' @param seam2 communication adjacency matrix, read in if readseam2=T
#' @param readbpam2 if T, the dispersal adjacency matrix is read in rather than being generated from the outset
#' @param bpam2 dispersal adjacency matrix, read in if readbpam2=T

#' @param maneffdir2 if maneffdir2='decrease_estab', the management reduces the probability of establishment; if maneffdir2='increase_estab', the management reduces the probability that establishment does NOT occur

#' @param maneffmean2 (estinfo, estab) the underlying mean change in establishment probability (as a proportion)
#' @param maneffsd2 (estinfo) the standard deviation of the effect
#' @param maneffthresh2 (estinfo) the threshold effect size for communicating about management (if maneffthresh2 = 0 there is no threshold so communication can always occur)
#' @param sampeffort2 (estinfo) sampling effort, where greater samping effort reduces the error in estimating the management effect

############ what about switch from sd to sampling effort - is maneffsd still needed somewhere???

#' @param maneffsd2 (estab) sd of management effect (same as maneffsd2 to function estinfo)


#' @param readgeocoords2 if T, read in geocoords2 - otherwise, generate it in each realization 
#' @param geocoords2 matrix of x,y coordinates of nodes, read in if readgeocoords2=T

#' @param xrange2 (genlocs) range of x coordinates
#' @param yrange2 (genlocs) range of y coordinates
#' @param numnodes2 (genlocs) the number of nodes
#' @param randgeo2 (genlocs) if TRUE then locations are randomly generated

#' @param readinitinfo2 info (setup2) if T, the initial values for the vector of starting locations for the presence of information are read in rather than generated
#' @param initinfo.s2 info (setup2) the vector of initial values read in if readinitinfo2 == T, a vector with lenth equal to the number of nodes and entries 1s or 0s with 1s indicating the initial presence of information

#' @param initinfo.dist2 info (initvals) the pattern of locations where initial presence occurs: 'random' indicates all equally likely, 'upedge' indicates that nodes closest to the upper edge have presence, 'rightedge' indicates that nodes closest to the right edge have presence
#' @param initinfo.n2 info (initvals) the number of initial locations for presence
#' @param initinfo.norp2 info (initvals) 'num' indicates initial number for presence, 'prop' indicates initial proportion
#' @param initinfo.p2 info (initvals) the proportion of initial locations for presence

#' @param readinitbio2 bio (setup2) if T, the initial values for the vector of starting locations for the presence of the bioentity are read in rather than generated
#' @param initbio2 bio (setup2) the vector of initial values read in if readinitbio2 == T, a vector with lenth equal to the number of nodes and entries 1s or 0s with 1s indicating the initial presence of the bioentity

#' @param initbio.dist2 estab (initvals) the type of locations where initial presence occurs: 'random' indicates all equally likely, 'upedge' indicates that nodes closest to the upper edge have presence, 'rightedge' indicates that nodes closest to the right edge have presence
#' @param initbio.n2 estab (initvals) the number of initial locations for presence
#' @param initbio.norp2 estab (initvals) 'num' indicates initial number for presence, 'prop' indicates initial proportion
#' @param initbio.p2 estab (initvals) the proportion of initial locations for presence

#### ending section of startup components

#' @param seamdist2 (genmovnet) the function of distance used to estimate movement probability - 'random' (not related to distance) or 'powerlaw' 
#' @param seampla2 (genmovnet) power law parameter a in ad^(-b)
#' @param seamplb2 (genmovnet) power law parameter b in ad^(-b)
#' @param seamrandp2 (genmovnet) random case, probability of link

#' @param bpamdist2 (genmovnet) the function of distance used to estimate movement probability - 'random' (not related to distance) or 'powerlaw' 
#' @param bpampla2 (genmovnet) power law parameter a in ad^(-b)
#' @param bpamplb2 (genmovnet) power law parameter b in ad^(-b)
#' @param bpamrandp2 (genmovnet) random case, probably of link

#' @param probadoptmean2 (makedec) mean probability of adopting management if informed
#' @param probadoptsd2 (makedec) sd in truncated normal distribution of probability of adoption
#' @param probadoptvec2 vector of probabilities of adoption if informed for nodes in the network (\code{readprobadoptvec2} determines whether \code{probadoptvec2} is read in to \code{makedec} or generated within \code{makedec} based on \code{probadoptmean2} and \code{probadoptsd2})
#' @param readprobadoptvec2 if T, then a VECTOR of probability of adoption values is read in; if F, \code{probadoptmean2} and \code{probadoptsd2} are used to generate the vector

#' @param probestabmean2 (estab) mean probability of establishment (new or CONTINUED) in absence of management
#' @param probestabsd2 (estab) sd in truncated normal distribution for probability of establishment in absence of management
#' @param probestabvec2 (estab) vector of probabilities of establishment (read in or generated when the function is run, depending on \code{readprobestabvec2} equal to T or F)
#' @param readprobestabvec2 (estab) if T, then a vector \code{probestabvec2} is read in; otherwise the vector is generated using \code{probestabmean2} and \code{probestabsd2}

#' @param doplot2 if true plots of resulting presence of information and species are generated

#' @keywords simulation experiments
#' @export
#' @examples


#' x5 <- multsame2(nreals2=10, ntimesteps2=3, usethreshman2=F, readgeocoords2=T, geocoords2=matrix(c(1,1, 1,2, 1,3, 2,1, 2,2, 2,3),byrow=T,ncol=2), maneffdir2='decrease_estab', maneffmean2=0.5, maneffsd2=0.1, maneffthresh2=0.5, sampeffort2=1, xrange2=NA, yrange2=NA, numnodes2=NA, randgeo2=NA, readinitinfo2=F, initinfo.dist2='random', initinfo.n2=5, initinfo.norp2='num', initinfo.p2=NA, readinitbio2=F, initbio.dist2='upedge', initbio.n2=5, initbio.norp2='num', initbio.p2=NA, readseam2=F, seam2=NA, seamdist2='powerlaw', seampla2=1, seamplb2=1, seamrandp2=NA, readbpam2=F, bpam2=NA, bpamdist2='powerlaw', bpampla2=1, bpamplb2=1, bpamrandp2=NA, probadoptmean2=0.5, probadoptsd2=0.1, probadoptvec2=NA, readprobadoptvec2=F, probestabmean2=0.5, probestabsd2=0.1, readprobestabvec2=F, probestabvec2=NA, doplot2=F)

#' x3 <- multsame2(nreals2=2,ntimesteps=3, usethreshman2=F, readgeocoords2=T, geocoords2=matrix(runif(n=100)*100,byrow=T,ncol=2), maneffdir2='decrease_estab', maneffmean2=0.5, maneffsd2=0.1, maneffthresh2=0.5, sampeffort2=1, xrange2=NA, yrange2=NA, numnodes2=NA, randgeo2=NA, readinitinfo2=F, initinfo.dist2='random', initinfo.n2=5, initinfo.norp2='num', initinfo.p2=0.05, readinitbio2=F, initbio.dist2='upedge', initbio.n2=5, initbio.norp2='num', initbio.p2=0.05, readseam2=F, seam2=NA, seamdist2='powerlaw', seampla2=1, seamplb2=1, seamrandp2=NA, readbpam2=F, bpam2=NA, bpamdist2='powerlaw', bpampla2=1, bpamplb2=1, bpamrandp2=NA, probadoptmean2=0.5, probadoptsd2=c(0,1), probadoptvec2=NA, readprobadoptvec2=F, probestabmean2=0.5, probestabsd2=0.1, readprobestabvec2=F, probestabvec2=NA, doplot2=F)



multsame2 <- function(nreals2, usethreshman2, ntimesteps2, geocoords2, readgeocoords2, maneffdir2, maneffmean2, maneffsd2, maneffthresh2, sampeffort2, xrange2, yrange2, numnodes2, randgeo2, readinitinfo2, initinfo2, initinfo.dist2, initinfo.n2, initinfo.norp2, initinfo.p2, readinitbio2, initbio2, initbio.dist2, initbio.n2, initbio.norp2, initbio.p2, readseam2, seam2, seamdist2, seampla2, seamplb2, seamrandp2, readbpam2, bpam2, bpamdist2, bpampla2, bpamplb2, bpamrandp2, probadoptmean2, probadoptsd2, probadoptvec2, readprobadoptvec2, probestabmean2, probestabsd2, readprobestabvec2, probestabvec2, doplot2=F){


  # determine the number of nodes, 
  #   if reading in the geographic coordinates

  if(readgeocoords2){numnodes2 <- nrow(geocoords2)}


  # prepare an output object with output from ntsteps2 
  #   for each realization

  multout <- as.list(1:nreals2)
  
  # prepare an output object with output from setup2 
  #   for each realization

  setupout <- multout

  # keep the proportion of nodes in each category:
  #   with com, with adoption, with disp, with establishment

  setocom <- 0*(1:nreals2)
  setodec <- setocom
  setodisp <- setocom
  setoestab <- setocom

  for(j in 1:nreals2) {


    tempo <- setup2(maneffmean3s=maneffmean2, maneffsd3s=maneffsd2, maneffthresh3=maneffthresh2, sampeffort3=sampeffort2, usethreshman3=usethreshman2, readgeocoords3s=readgeocoords2, geocoords3s=geocoords2, xrange3=xrange2, yrange3=yrange2, numnodes3=numnodes2, randgeo3=randgeo2, readinitinfo3=readinitinfo2, initinfo3=initinfo2, initinfo.dist3=initinfo.dist2, initinfo.n3=initinfo.n2, initinfo.norp3=initinfo.norp2, initinfo.p3=initinfo.p2, readinitbio3=readinitbio2, initbio3=initbio2, initbio.dist3=initbio.dist2, initbio.n3=initbio.n2, initbio.norp3=initbio.norp2, initbio.p3=initbio.p2, plotmp=doplot2)

# note that the following uses the output object tempo

    temp2 <- ntsteps2(nsteps=ntimesteps2, infon=tempo$com.yes, geocoords3n=tempo$geocoords3s, vect1cn=tempo$infovec, vect1dn=tempo$estabvec, readseam3=readseam2, seam3=seam2, seamdist3=seamdist2, seampla3=seampla2, seamplb3=seamplb2, seamrandp3=seamrandp2, readbpam3=readbpam2, bpam3=bpam2, bpamdist3=bpamdist2, bpampla3=bpampla2, bpamplb3=bpamplb2, bpamrandp3=bpamrandp2, readprobadoptvec3=readprobadoptvec2, probadoptvec3=probadoptvec2, probadoptmean3=probadoptmean2, probadoptsd3=probadoptsd2, probestabmean3=probestabmean2, probestabsd3=probestabsd2, maneffdir3=maneffdir2, maneffmean3n=maneffmean2, maneffsd3n=maneffsd2, readprobestabvec3=readprobestabvec2, probestabvec3=probestabvec2, plotmpn=doplot2)

    # save the output from setup2 and ntsteps2
    multout[[j]] <- temp2
    setupout[[j]] <- tempo

    # the proportion of nodes in each category
    setocom[j] <- sum(temp2$vect1cL[[ntimesteps2+1]])/numnodes2
    setodec[j] <- sum(temp2$decvecL[[ntimesteps2]])/numnodes2
    setodisp[j] <- sum(temp2$vect1dL[[ntimesteps2+1]])/numnodes2
    setoestab[j] <- sum(temp2$estabvecL[[ntimesteps2]])/numnodes2

  } ### end of nreals2 realizations

  meancom <- mean(setocom)
  meandec <- mean(setodec)
  meandisp <- mean(setodisp)
  meanestab <- mean(setoestab)

# add the 5th percentile and 95th percentile to the output

  com5 <- quantile(setocom, probs=0.05)
  dec5 <- quantile(setodec, probs=0.05)
  disp5 <- quantile(setodisp, probs=0.05)
  estab5 <- quantile(setoestab, probs=0.05)

  com95 <- quantile(setocom, probs=0.95)
  dec95 <- quantile(setodec, probs=0.95)
  disp95 <- quantile(setodisp, probs=0.95)
  estab95 <- quantile(setoestab, probs=0.95)


  list(multout=multout, setocom=setocom, setodec=setodec, setodisp=setodisp, setoestab=setoestab, meancom=meancom, meandec=meandec, meandisp=meandisp, meanestab=meanestab, com5=com5, dec5=dec5, disp5=disp5, estab5=estab5, com95=com95, dec95=dec95, disp95=disp95, estab95=estab95)
}