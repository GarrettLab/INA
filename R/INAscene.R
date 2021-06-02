
#' Evaluates scenarios in an impact network analysis (INA)
#'
#' This function implements and summarizes multiple simulations across a designated range of parameter values
#'
#' Updated 2020-10-13

#' @param nreals number of realizations to be evaluated
#' @param ntimesteps number of time steps to be evaluated
#' @param doplot if true, plots of resulting presence of information and bioentity are generated
#' @param outputvol output volume, where outputvol='less' excludes the list element $multdetails from output, and outputvol='more' includes it in output.  $multdetails becomes large quickly for many realizations for large matrices so outputvol='less' would be desirable for extensive analyses

#' @param readgeocoords if T, read in geocoords - otherwise, generate it in each realization
#' @param geocoords matrix of x,y coordinates of nodes

#' @param numnodes (genlocs) the number of nodes, can be a vector of different numbers of nodes for scenario comparisons (numnodes can be entered as a vector to evaluate multiple scenarios)
#' @param xrange (genlocs) range of x coordinates, e.g., c(0,50)
#' @param yrange (genlocs) range of y coordinates, e.g., c(0,20)
#' @param randgeo (genlocs) if TRUE then locations are randomly generated (the only location generation option for the moment)

#' @param readinitinfo info (setup2) if T, the initial values for the vector of starting locations for the presence of information are read in rather than generated
#' @param initinfo info (setup2) the vector of initial values read in if readinitinfo == T, a vector with lenth equal to the number of nodes and entries 1s or 0s with 1s indicating the initial presence of information

#' @param initinfo.norp info (initvals) 'num' indicates initial number for presence, 'prop' indicates initial proportion
#' @param initinfo.n info (initvals) the number of initial locations for presence (initinfo.n can be entered as a vector to evaluate multiple scenarios)
#' @param initinfo.p info (initvals) the proportion of initial locations for presence (initinfo.p can be entered as a vector to evaluate multiple scenarios)
#' @param initinfo.dist info (initvals) the type of locations where initial presence occurs: 'random' indicates all equally likely, 'upedge' indicates that nodes closest to the upper edge have presence, 'rightedge' indicates that nodes closest to the right edge have presence

#' @param readinitbio bio (setup2) if T, the initial values for the vector of starting locations for the presence of the bioentity are read in rather than generated
#' @param initbio bio (setup2) the vector of initial values read in if readinitbio.s == T, a vector with lenth equal to the number of nodes and entries 1s or 0s with 1s indicating the initial presence of the bioentity

#' @param initbio.norp estab (initvals) 'num' indicates initial number for presence, 'prop' indicates initial proportion
#' @param initbio.n estab (initvals) the number of initial locations for presence (initbio.n can be entered as a vector to evaluate multiple scenarios)
#' @param initbio.p estab (initvals) the proportion of initial locations for presence (initbio.p can be entered as a vector to evaluate multiple scenarios)
#' @param initbio.dist estab (initvals) the type of locations where initial presence occurs: 'random' indicates all equally likely, 'upedge' indicates that nodes closest to the upper edge have presence, 'rightedge' indicates that nodes closest to the right edge have presence

#' @param readseam if T, the socioeconomic network adjacency matrix is read in rather than being generated from the outset
#' @param seam socioeconnomic network adjacency matrix, read in if readseam=T (assumed to be structured so that rows are sources and columns are sinks)

#' @param seamdist com (genmovnet) the function of distance used to estimate movement probability - 'random' (not related to distance) or 'powerlaw'
#' @param seamrandp (genmovnet) probability of link existence in random network generated for the socioeconomic network (seamrandp can be entered as a vector to evaluate multiple scenarios)
#' @param seampla com (genmovnet) power law parameter a in ad^(-b) (seampla can be entered as a vector to evaluate multiple scenarios)
#' @param seamplb com (genmovnet) power law parameter b in ad^(-b) (seamplb can be entered as a vector to evaluate multiple scenarios)

#' @param readbpam if T, the biophysical network adjacency matrix, describing dispersal likelihoods, is read in rather than being generated from the outset
#' @param bpam biophysical adjacency matrix, read in if readbpam=T (assumed to be structured so that rows are sources and columns are sinks)

#' @param bpamdist disp (genmovnet) the function of distance used to estimate movement probability - 'random' (not related to distance) or 'powerlaw' 
#' @param bpamrandp (genmovnet) if bpamdist='random', the probability a link exists in the biophysical network adjacency matrix (bpamrandp can be entered as a vector to evaluate multiple scenarios)
#' @param bpampla disp (genmovnet) power law parameter a in ad^(-b) (bpampla can be entered as a vector to evaluate multiple scenarios)
#' @param bpamplb disp (genmovnet) power law parameter b in ad^(-b)  (bpamplb can be entered as a vector to evaluate multiple scenarios)

#' @param readprobadoptvec if T, read in a vector of probabilities of adoption for each node in socioeconomic network - if F, generate this vector
#' @param probadoptvec vector of probabilities of adoption for each node in socioeconomic network, read in if readprobadoptvec=T

#' @param probadoptmean (makedec) mean probability of adopting management if informed (probadoptmean can be entered as a vector to evaluate multiple scenarios)
#' @param probadoptsd (makedec) sd in truncated normal distribution of probability of adoption (probadoptsd can be entered as a vector to evaluate multiple scenarios)

#' @param readprobestabvec if T, read in a vector of probabilities of establishment for each node in biophysical network - if F, generate this vector
#' @param probestabvec vector of probabilities of establishment for each node in biophysical network, read in if readprobestabvec=T

#' @param probestabmean (estab) mean probability of establishment (new or CONTINUED) in absence of management (probestabmean can be entered as a vector to evaluate multiple scenarios)
#' @param probestabsd (estab) sd of probability of establishment in absence of management in truncated normal distribution (probestabsd can be entered as a vector to evaluate multiple scenarios)

#' @param maneffmean (estinfo, estab) the underlying mean change in establishment probability (as a proportion) resulting from management technology (maneffmean can be entered as a vector to evaluate multiple scenarios)
#' @param maneffsd (estinfo, estab) the standard deviation of the management technology effect (maneffsd can be entered as a vector to evaluate multiple scenarios)

#' @param maneffdir if maneffdir='decrease_estab', the management reduces the probability of establishment; if maneffdir='increase_estab', the management reduces the probability that establishment does NOT occur

#' @param usethreshman () if T, the threshold for management is used, so that information is never present anywhere unless the management effect estimate exceeds the threshold
#' @param maneffthresh (estinfo) the threshold effect size for communicating about management (maneffthresh can be entered as a vector to evaluate multiple scenarios)
#' @param sampeffort (estinfo) the sampling effort, where increasing sampling effort results in a more precise estimate of the management effect size (sampeffort can be entered as a vector to evaluate multiple scenarios)


#' @keywords simulation experiments
#' @export 

#' @examples

#' j25.readgeocoords <- INAscene(nreals=3, ntimesteps=3, doplot=F, readgeocoords=T, geocoords=matrix(c(1,1, 1,2, 1,3, 2,1, 2,2, 2,3),byrow=T,ncol=2), numnodes=NA, xrange=NA, yrange=NA, randgeo=F, readinitinfo=T, initinfo=c(1,1,1,0,0,0), initinfo.norp=NA, initinfo.n=NA, initinfo.p=NA, initinfo.dist=NA, readinitbio=T, initbio=c(0,0,0,1,1,1), initbio.norp=NA, initbio.n=NA, initbio.p=NA,  initbio.dist=NA, readseam=F, seam=NA, seamdist='random', seamrandp=c(0.01,0.05,0.1,0.5), seampla=NA, seamplb=NA, readbpam=F, bpam=NA, bpamdist='random', bpamrandp=0.1, bpampla=NA, bpamplb=NA, readprobadoptvec=F, probadoptvec=NA, probadoptmean=0.1, probadoptsd=0.1, readprobestabvec=F, probestabvec=NA, probestabmean=0.1, probestabsd=0.1, maneffdir='decrease_estab', maneffmean=0.5, maneffsd=0.1, usethreshman=F, maneffthresh=NA, sampeffort=NA) 

#' j25.baseline <- INAscene(nreals=3, ntimesteps=3, doplot=F, readgeocoords=F, geocoords=NA, numnodes=50, xrange=c(0,10), yrange=c(0,20), randgeo=T, readinitinfo=F, initinfo=NA, initinfo.norp='num', initinfo.n=5, initinfo.p=NA, initinfo.dist='random', readinitbio=F, initbio=NA, initbio.norp='num', initbio.n=5, initbio.p=NA,  initbio.dist='random', readseam=F, seam=NA, seamdist='powerlaw', seamrandp=NA, seampla=1, seamplb=0.5, readbpam=F, bpam=NA, bpamdist='powerlaw', bpamrandp=NA, bpampla=1, bpamplb=0.5, readprobadoptvec=F, probadoptvec=NA, probadoptmean=0.5, probadoptsd=0.2, readprobestabvec=F, probestabvec=NA, probestabmean=0.5, probestabsd=0.2, maneffdir='decrease_estab', maneffmean=0.5, maneffsd=0.2, usethreshman=T, maneffthresh=0.5, sampeffort=2) 

#' j25.baseline$multout

#' j25.readseam <- INAscene(nreals=3, ntimesteps=3, doplot=F, readgeocoords=T, geocoords=matrix(c(1,1, 1,2, 1,3, 2,1, 2,2, 2,3),byrow=T,ncol=2), numnodes=NA, xrange=NA, yrange=NA, randgeo=F, readinitinfo=T, initinfo=c(1,1,1,0,0,0), initinfo.norp=NA, initinfo.n=NA, initinfo.p=NA, initinfo.dist=NA, readinitbio=T, initbio=c(0,0,0,1,1,1), initbio.norp=NA, initbio.n=NA, initbio.p=NA,  initbio.dist=NA, readseam=T, seam=matrix(c(1,0,0,0,1,0, 0,1,0,0,1,0, 0,0,1,0,1,0, 0,1,0,1,0,0, 0,0,0,0,1,1, 0,1,0,0,0,1),byrow=T,nrow=6), seamdist=NA, seamrandp=NA, seampla=NA, seamplb=NA, readbpam=F, bpam=NA, bpamdist='random', bpamrandp=0.1, bpampla=NA, bpamplb=NA, readprobadoptvec=F, probadoptvec=NA, probadoptmean=0.1, probadoptsd=0.1, readprobestabvec=F, probestabvec=NA, probestabmean=0.1, probestabsd=0.1, maneffdir='decrease_estab', maneffmean=0.5, maneffsd=0.1, usethreshman=F, maneffthresh=NA, sampeffort=NA) 

#' sens.probadoptmean <- INAscene(nreals=15, ntimesteps=3, doplot=F, readgeocoords=F, geocoords=NA, numnodes=50, xrange=c(0,10), yrange=c(0,20), randgeo=T, readinitinfo=F, initinfo=NA, initinfo.norp='num', initinfo.n=5, initinfo.p=NA, initinfo.dist='random', readinitbio=F, initbio=NA, initbio.norp='num', initbio.n=5, initbio.p=NA,  initbio.dist='random', readseam=F, seam=NA, seamdist='powerlaw', seamrandp=NA, seampla=1, seamplb=0.5, readbpam=F, bpam=NA, bpamdist='powerlaw', bpamrandp=NA, bpampla=1, bpamplb=0.5, readprobadoptvec=F, probadoptvec=NA, probadoptmean=seq(0,1,0.1), probadoptsd=0.2, readprobestabvec=F, probestabvec=NA, probestabmean=0.5, probestabsd=0.2, maneffdir='decrease_estab', maneffmean=0.9, maneffsd=0.2, usethreshman=T, maneffthresh=0.3, sampeffort=2) 

#' jt <- sens.probadoptmean$multout
#' plot(jt$probadoptmean, jt$mestab, xlab='Mean probability of adopting technology if informed', ylab='Proportion nodes with bioentity', xlim=c(0,1), ylim=c(0,1))
#' plot(jt$probadoptmean, jt$mdec, xlab='Mean probability of adopting technology if informed', ylab='Proportion nodes with technology adoption', xlim=c(0,1), ylim=c(0,1))

#' sens.seamplb <- INAscene(nreals=15, ntimesteps=3, doplot=F, readgeocoords=F, geocoords=NA, numnodes=50, xrange=c(0,50), yrange=c(0,50), randgeo=T, readinitinfo=F, initinfo=NA, initinfo.norp='num', initinfo.n=5, initinfo.p=NA, initinfo.dist='random', readinitbio=F, initbio=NA, initbio.norp='num', initbio.n=5, initbio.p=NA,  initbio.dist='random', readseam=F, seam=NA, seamdist='powerlaw', seamrandp=NA, seampla=1, seamplb=seq(0,2,0.1), readbpam=F, bpam=NA, bpamdist='powerlaw', bpamrandp=NA, bpampla=1, bpamplb=0.5, readprobadoptvec=F, probadoptvec=NA, probadoptmean=0.7, probadoptsd=0.2, readprobestabvec=F, probestabvec=NA, probestabmean=0.5, probestabsd=0.2, maneffdir='decrease_estab', maneffmean=0.9, maneffsd=0.2, usethreshman=T, maneffthresh=0.3, sampeffort=2) 

#' jt2 <- sens.seamplb$multout
#' plot(jt2$seamplb, jt2$mestab, xlab='Power law parameter b in ad^(-b) for communication links', ylab='Proportion nodes with bioentity', xlim=c(0,2), ylim=c(0,1))
#' plot(jt2$seamplb, jt2$mdec, xlab='Power law parameter b in ad^(-b) for communication links', ylab='Proportion nodes with technology adoption', xlim=c(0,2), ylim=c(0,1))
#' plot(jt2$seamplb, jt2$mcom, xlab='Power law parameter b in ad^(-b) for communication links', ylab='Proportion nodes with information about technology', xlim=c(0,2), ylim=c(0,1))


#' j25.seamrandp <- INAscene(nreals=3, ntimesteps=3, doplot=F, readgeocoords=F, geocoords=NA, numnodes=50, xrange=c(0,10), yrange=c(0,20), randgeo=T, readinitinfo=F, initinfo=NA, initinfo.norp='num', initinfo.n=5, initinfo.p=NA, initinfo.dist='random', readinitbio=F, initbio=NA, initbio.norp='num', initbio.n=5, initbio.p=NA,  initbio.dist='random', readseam=F, seam=NA, seamdist='random', seamrandp=c(0.01,0.05,0.1,0.5), seampla=NA, seamplb=NA, readbpam=F, bpam=NA, bpamdist='random', bpamrandp=0.1, bpampla=NA, bpamplb=NA, readprobadoptvec=F, probadoptvec=NA, probadoptmean=0.1, probadoptsd=0.1, readprobestabvec=F, probestabvec=NA, probestabmean=0.1, probestabsd=0.1, maneffdir='decrease_estab', maneffmean=0.5, maneffsd=0.1, usethreshman=F, maneffthresh=NA, sampeffort=NA) 

#' sens.exp2a <- INAscene(nreals=10, ntimesteps=10, doplot=F, readgeocoords=F, geocoords=NA, numnodes=50, xrange=c(0,50), yrange=c(0,50), randgeo=T, readinitinfo=F, initinfo=NA, initinfo.norp='prop', initinfo.n=NA, initinfo.p=0.05, initinfo.dist='random', readinitbio=F, initbio=NA, initbio.norp='prop', initbio.n=NA, initbio.p=0.05,  initbio.dist='upedge', readseam=F, seam=NA, seamdist='powerlaw', seamrandp=NA, seampla=1, seamplb=0.5, readbpam=F, bpam=NA, bpamdist='powerlaw', bpamrandp=NA, bpampla=1, bpamplb=0.5, readprobadoptvec=F, probadoptvec=NA, probadoptmean=0.5, probadoptsd=0.2, readprobestabvec=F, probestabvec=NA, probestabmean=0.5, probestabsd=0.2, maneffdir='decrease_estab', maneffmean=0.5, maneffsd=0.2, usethreshman=F, maneffthresh=NA, sampeffort=2) 

#' sens.exp2a <- INAscene(nreals=10, ntimesteps=10, doplot=F, outputvol='less', readgeocoords=F, geocoords=NA, numnodes=50, xrange=c(0,50), yrange=c(0,50), randgeo=T, readinitinfo=F, initinfo=NA, initinfo.norp='prop', initinfo.n=NA, initinfo.p=0.05, initinfo.dist='random', readinitbio=F, initbio=NA, initbio.norp='prop', initbio.n=NA, initbio.p=0.05,  initbio.dist='upedge', readseam=F, seam=NA, seamdist='powerlaw', seamrandp=NA, seampla=1, seamplb=0.5, readbpam=F, bpam=NA, bpamdist='powerlaw', bpamrandp=NA, bpampla=1, bpamplb=0.5, readprobadoptvec=F, probadoptvec=NA, probadoptmean=0.5, probadoptsd=0.2, readprobestabvec=F, probestabvec=NA, probestabmean=0.5, probestabsd=0.2, maneffdir='decrease_estab', maneffmean=0.5, maneffsd=0.2, usethreshman=F, maneffthresh=NA, sampeffort=2) 



INAscene <- function(nreals, ntimesteps, doplot=F, outputvol='more', readgeocoords, geocoords=NA, numnodes=NA, xrange=NA, yrange=NA, randgeo=NA, readinitinfo, initinfo=NA, initinfo.norp=NA, initinfo.n=NA, initinfo.p=NA, initinfo.dist=NA, readinitbio, initbio=NA, initbio.norp=NA, initbio.n=NA, initbio.p=NA,  initbio.dist=NA, readseam, seam=NA, seamdist=NA, seamrandp=NA, seampla=NA, seamplb=NA, readbpam, bpam=NA, bpamdist=NA, bpamrandp=NA, bpampla=NA, bpamplb=NA, readprobadoptvec, probadoptvec=NA, probadoptmean=NA, probadoptsd=NA, readprobestabvec, probestabvec=NA, probestabmean=NA, probestabsd=NA, maneffdir=NA, maneffmean=NA, maneffsd=NA, usethreshman, maneffthresh=NA, sampeffort=NA) {

### error message(s)

if(maneffdir != 'decrease_estab' & maneffdir != 'increase_estab'){stop('Analysis stopped because of the parameter maneffdir which must be either increase_estab or decrease_estab')}

### prepare the output matrix, taking into account which
###   variables are used, as a function of which options
###   were selected

# create initial output matrix, starting with only the 12
#   output variables

vn <- c('mcom', 'mdec', 'mestab', 'vcom', 'vdec', 'vestab', 'com5', 'dec5', 'estab5', 'com95', 'dec95', 'estab95')

multout <- data.frame(matrix(-99,ncol=length(vn),nrow=1))

names(multout) <- vn


# function expandoutmat is used to expand the output matrix 
#   (add rows) as more combinations of parameter levels are
#   considered

expandoutmat <- function(startoutmat, newvar){
  endmat <- startoutmat
  if (length(newvar) > 1){
    for (i2 in 2:length(newvar)) {
      endmat <- rbind(endmat,startoutmat)
    }
  } 
  endmat
}

# 1. Add to output matrix, taking into account whether
#   reading in geographic locations -or- generating them
#
# numnodes can be a vector of values

if (readgeocoords==F) {
  nrowsbefore <- nrow(multout)
  multout <- expandoutmat(startoutmat=multout,newvar=numnodes)
  multout$numnodes <- rep(numnodes, each=nrowsbefore)
}

# 2.1. Add to output matrix, taking into account whether
#   reading in which nodes have info initially -or- generating
#   the node locations
#
# initinfo.n or initinfo.p can be a vector of values

if (readinitinfo==F) {
  if (initinfo.norp=='num') {
    nrowsbefore <- nrow(multout)
    multout <- expandoutmat(multout,newvar=initinfo.n)
    multout$initinfo.n <- rep(initinfo.n, each=nrowsbefore)
  } else if (initinfo.norp=='prop'){
    nrowsbefore <- nrow(multout)
    multout <- expandoutmat(multout,newvar=initinfo.p)
    multout$initinfo.p <- rep(initinfo.p, each=nrowsbefore)
  }
}

# 2.2. Add to output matrix, taking into account whether
#   reading in which nodes have the bioentity initially -or-
#   generating the node locations
#
# initbio.n or initbio.p can be a vector of values

if (readinitbio==F) {
  if (initbio.norp=='num') {
    nrowsbefore <- nrow(multout)
    multout <- expandoutmat(startoutmat=multout,newvar=initbio.n)
    multout$initbio.n <- rep(initbio.n, each=nrowsbefore)
  } else if (initbio.norp=='prop'){
    nrowsbefore <- nrow(multout)
    multout <- expandoutmat(multout,newvar=initbio.p)
    multout$initbio.p <- rep(initbio.p, each=nrowsbefore)
  }
}

# 3.1. Add to output matrix, taking into account whether
#   reading in the adjacency matrix for the socioeconomic
#     network -or- generating this adjacency matrix
#
# seamrandp or seampla and seamplb
#   can be vectors of values

if (readseam==F) {
  if (seamdist=='random') {
    nrowsbefore <- nrow(multout)
    multout <- expandoutmat(startoutmat=multout,newvar=seamrandp)
    multout$seamrandp <- rep(seamrandp, each=nrowsbefore)
  } else if (seamdist=='powerlaw'){
    nrowsbefore <- nrow(multout)
    multout <- expandoutmat(multout,newvar=seampla)
    multout$seampla <- rep(seampla, each=nrowsbefore)
    nrowsbefore <- nrow(multout)
    multout <- expandoutmat(multout,newvar=seamplb)
    multout$seamplb <- rep(seamplb, each=nrowsbefore)
  }
}

# 3.2. Add to output matrix, taking into account whether
#   reading in the adjacency matrix for the biophysical
#     network -or- generating this adjacency matrix
#
# bpamrandp or bpampla and bpamplb
#   can be vectors of values

if (readbpam==F) {
  if (bpamdist=='random') {
    nrowsbefore <- nrow(multout)
    multout <- expandoutmat(startoutmat=multout,newvar=bpamrandp)
    multout$bpamrandp <- rep(bpamrandp, each=nrowsbefore)
  } else if (bpamdist=='powerlaw'){
    nrowsbefore <- nrow(multout)
    multout <- expandoutmat(multout,newvar=bpampla)
    multout$bpampla <- rep(bpampla, each=nrowsbefore)
    nrowsbefore <- nrow(multout)
    multout <- expandoutmat(multout,newvar=bpamplb)
    multout$bpamplb <- rep(bpamplb, each=nrowsbefore)
  }
}

# 4.1. Add to output matrix, taking into account whether
#   reading in a vector of the probability of adoption
#   -or- generating this vector
#
# probadoptmean and probadoptsd can be vectors of values

if (readprobadoptvec==F) {
  nrowsbefore <- nrow(multout)
  multout <- expandoutmat(multout,newvar=probadoptmean)
  multout$probadoptmean <- rep(probadoptmean, each=nrowsbefore)
  nrowsbefore <- nrow(multout)
  multout <- expandoutmat(startoutmat=multout,newvar=probadoptsd)
  multout$probadoptsd <- rep(probadoptsd, each=nrowsbefore)
}

# 4.2. Add to output matrix, taking into account whether
#   reading in a vector of the probability of establishment in
#   the absence of management -or- generating this vector
#
# probestabmean, probestabsd, maneffmean, and maneffsd can be
#   vector of values

if (readprobestabvec==F) {
  nrowsbefore <- nrow(multout)
  multout <- expandoutmat(multout,newvar=probestabmean)
  multout$probestabmean <- rep(probestabmean, each=nrowsbefore)
  nrowsbefore <- nrow(multout)
  multout <- expandoutmat(startoutmat=multout,newvar=probestabsd)
  multout$probestabsd <- rep(probestabsd, each=nrowsbefore)
}
nrowsbefore <- nrow(multout)
multout <- expandoutmat(startoutmat=multout,newvar=maneffmean)
multout$maneffmean <- rep(maneffmean, each=nrowsbefore)
nrowsbefore <- nrow(multout)
multout <- expandoutmat(startoutmat=multout,newvar=maneffsd)
multout$maneffsd <- rep(maneffsd, each=nrowsbefore)


# 5. Add to output matrix, taking into account whether
#   the estimated management effect must exceed a 
#   threshold for communication to take place
#
# maneffthresh and sampeffort can be a vector of values

if (usethreshman==TRUE) {
  nrowsbefore <- nrow(multout)
  multout <- expandoutmat(multout,newvar=maneffthresh)
  multout$maneffthresh <- rep(maneffthresh, each=nrowsbefore)
  nrowsbefore <- nrow(multout)
  multout <- expandoutmat(multout,newvar=sampeffort)
  multout$sampeffort <- rep(sampeffort, each=nrowsbefore)
}


# number of rows in multout = number of parameter combinations
ncombs <- nrow(multout) 

# number of columns in multout = number of parameters varying 
#   plus 12 response variables
ncmultout <- ncol(multout) 

# set up an object to keep all the output
multdetails <- as.list(1:(ncombs))

# keep track of which row is being evaluated
rowcount <- 1

### Evaluate the scenario for each parameter combination

for (j22 in 1:ncombs) {

  # prepare input values for multsame2 taking the values
  #   from the current row of multout (params start at 13)

  for (i22 in 13:ncmultout) {
    assign(names(multout)[i22], value=multout[j22,i22])
  }


  temp <- multsame2(nreals2=nreals, ntimesteps2=ntimesteps, doplot2=doplot, readgeocoords2=readgeocoords, geocoords2=geocoords, numnodes2=numnodes, xrange2=xrange, yrange2=yrange, randgeo2=randgeo, readinitinfo2=readinitinfo, initinfo2=initinfo, initinfo.norp2=initinfo.norp, initinfo.n2=initinfo.n, initinfo.p2=initinfo.p, initinfo.dist2=initinfo.dist, readinitbio2=readinitbio, initbio2=initbio, initbio.norp2=initbio.norp, initbio.n2=initbio.n, initbio.p2=initbio.p, initbio.dist2=initbio.dist, readseam2=readseam, seam2=seam,  seamdist2=seamdist, seamrandp2=seamrandp, seampla2=seampla, seamplb2=seamplb, readbpam2=readbpam, bpam2=bpam, bpamdist2=bpamdist, bpamrandp2=bpamrandp, bpampla2=bpampla, bpamplb2=bpamplb, readprobadoptvec2=readprobadoptvec, probadoptvec2=probadoptvec, probadoptmean2=probadoptmean, probadoptsd2=probadoptsd, readprobestabvec2=readprobestabvec, probestabvec2=probestabvec, probestabmean2=probestabmean, probestabsd2=probestabsd, maneffdir2=maneffdir, maneffmean2=maneffmean, maneffsd2=maneffsd, usethreshman2=usethreshman, maneffthresh2=maneffthresh, sampeffort2=sampeffort)


  multdetails[[rowcount]] <- temp

  multout$mcom[rowcount] <- temp$meancom
  multout$mdec[rowcount] <- temp$meandec
  multout$mestab[rowcount] <- temp$meanestab
  multout$vcom[rowcount] <- var(temp$setocom)
  multout$vdec[rowcount] <- var(temp$setodec)
  multout$vestab[rowcount] <- var(temp$setoestab)

  multout$com5[rowcount] <- temp$com5
  multout$dec5[rowcount] <- temp$dec5
  multout$estab5[rowcount] <- temp$estab5

  multout$com95[rowcount] <- temp$com95
  multout$dec95[rowcount] <- temp$dec95
  multout$estab95[rowcount] <- temp$estab95

  print(paste("ending parameter combination",rowcount))
  rowcount <- rowcount + 1

}

if (outputvol=='more') {
    list(multdetails=multdetails, multout=multout)
} else if (outputvol=='less'){
    list(multout=multout)
}

}