
#' Set up starting conditions for INA scenario analysis
#'
#' This function sets up all starting conditions, including the estimated effect size (but not direction), the x,y coordinates of the nodes, the initial information locations, and the initial species locations. It runs functions \code{estinfo}, \code{genlocs}, and \code{initvals}. The next step after this is to address the adjacency matrices for communication and for dispersal.
#'
#' Updated 2020-09-05

#' @param maneffmean3s (estinfo) the underlying mean change in establishment probability (as a proportion)
#' @param maneffsd3s (estinfo) the standard deviation of the effect
#' @param maneffthresh3 (estinfo) the threshold effect size for communicating about management (if maneffthresh3 = 0 there is no threshold so communication can always occur)
#' @param sampeffort3 (estinfo) sampling effort, where greater samping effort reduces the error in estimating the management effect

#' @param usethreshman3 () if T, the threshold for management maneffthresh3 is used to determine whether communication occurs; communication only occurs if the observed management effect is greater than maneffthresh3

#' @param readgeocoords3s read in the xy coordinates for locations (readgeocoords3s = T) or generate the xy coordinates by calling genlocs (readgeocoords3s = F)
#' @param geocoords3s coordinates for node geographic locations, ncol= 2 (x,y) and nrow=number of nodes (to be read in if readorgenxy = 'read')

#' @param xrange3 (genlocs) range of x coordinates
#' @param yrange3 (genlocs) range of y coordinates
#' @param numnodes3 (genlocs) the number of nodes (nn in genlocs)
#' @param randgeo3 (genlocs) if TRUE then locations are randomly generated

#' @param readinitinfo3 if T, the initial values for the vector of starting locations for the presence of information are read in rather than generated
#' @param initinfo3 the vector of initial values read in if readinitinfo3 == T, a vector with lenth equal to the number of nodes and entries 1s or 0s with 1s indicating the initial presence of information

#' @param initinfo.dist3 (initvals) the type of locations where initial presence occurs: 'random' indicates all equally likely, 'upedge' indicates that nodes closest to the upper edge have presence, 'rightedge' indicates that nodes closest to the right edge have presence
#' @param initinfo.n3 (initvals) the number of initial locations for presence
#' @param initinfo.norp3 (initvals) 'num' indicates initial number for presence, 'prop' indicates initial proportion
#' @param initinfo.p3 (initvals) the proportion of initial locations for presence (ip in initvals)

#' @param readinitbio3 if T, the initial values for the vector of starting locations for the presence of the bioentity are read in rather than generated
#' @param initbio3 the vector of initial values read in if readinitbio3 == T, a vector with lenth equal to the number of nodes and entries 1s or 0s with 1s indicating the initial presence of the bioentity

#' @param initbio.dist3 bioentity (initvals) the type of locations where initial presence occurs: 'random' indicates all equally likely, 'upedge' indicates that nodes closest to the upper edge have presence, 'rightedge' indicates that nodes closest to the right edge have presence
#' @param initbio.n3 bioentity (initvals) the number of initial locations for presence
#' @param initbio.norp3 bioentity (initvals) 'num' indicates initial number for presence, 'prop' indicates initial proportion
#' @param initbio.p3 bioentity (initvals) the proportion of initial locations for presence (ip in initvals)

#' @param plotmp if T plots of initial presence of information and species are generated

#' @keywords simulation experiments
#' @export 

#' @examples
#' x2 <- setup2(maneffmean3s=0.5, maneffsd3s=0.5, maneffthresh3=0.5, sampeffort3=1, usethreshman3=T, readgeocoords3s=F, xrange3=c(0,50), yrange3=c(0,50), numnodes3=100, randgeo3=TRUE, readinitinfo3=F, initinfo.dist3='random', initinfo.n3=5, initinfo.norp3='num', readinitbio3=F, initbio.dist3='upedge', initbio.n3=5, initbio.norp3='num', plotmp=T)

#' x3 <- setup2(maneffmean3s=0.5, maneffsd3s=0.5, maneffthresh3=0.5, sampeffort3=1, usethreshman3=T, readgeocoords3s=F, xrange3=c(0,50), yrange3=c(0,50), numnodes3=100, randgeo3=TRUE, readinitinfo3=F, initinfo.dist3='random', initinfo.n3=15, initinfo.norp3='num', readinitbio3=F, initbio.dist3='upedge', initbio.n3=15, initbio.norp3='num', plotmp=T)

#' x4.readgeocoords3s <- setup2(maneffmean3s=0.5, maneffsd3s=0.5, maneffthresh3=0.5, sampeffort3=1, usethreshman3=T, readgeocoords3s=T, geocoords3s=matrix(c(1,1, 1,2, 1,3, 2,1, 2,2, 2,3),byrow=T,ncol=2), readinitinfo3=F, initinfo.dist3='random', initinfo.n3=2, initinfo.norp3='num', readinitbio3=F, initbio.dist3='upedge', initbio.n3=3, initbio.norp3='num', plotmp=T)


setup2 <- function(maneffmean3s, maneffsd3s, maneffthresh3, sampeffort3, usethreshman3, readgeocoords3s, geocoords3s, xrange3, yrange3, numnodes3, randgeo3, readinitinfo3, initinfo3, initinfo.dist3, initinfo.n3, initinfo.norp3, initinfo.p3, readinitbio3, initbio3, initbio.dist3, initbio.n3, initbio.norp3, initbio.p3, plotmp=F){


if (usethreshman3) {
  infout <- estinfo(maneffmean4s=maneffmean3s, maneffsd4s=maneffsd3s, maneffthresh4=maneffthresh3, sampeffort4=sampeffort3)
} else {
  infout <- list()
  infout$com.yes <- T
  infout$obschange <- NA
}


if (readgeocoords3s == F){
  geocoords3s <- genlocs(xrange4=xrange3, yrange4=yrange3, numnodes4=numnodes3, randgeo4=randgeo3)
}

# generate initial locations for info if not supplied
if (readinitinfo3 == FALSE) {
  infovec <- initvals(geocoords4s=geocoords3s, dist4=initinfo.dist3, init.n4=initinfo.n3, norp4=initinfo.norp3, init.p4=initinfo.p3) 
} else if (readinitinfo3 == TRUE) {
  infovec <- initinfo3
}

# generate initial locations for bioentity establishment
#   if not supplied
if (readinitbio3== FALSE) {
  estabvec <- initvals(geocoords4s=geocoords3s, dist4=initbio.dist3, init.n4=initbio.n3, norp4=initbio.norp3, init.p4=initbio.p3) 
} else if (readinitbio3 == TRUE) {
  estabvec <- initbio3
}

if(plotmp){
    plot(geocoords3s, main='Information, shaded = yes')
    if(!is.null(dim(geocoords3s[infovec==1,]))){
      points(geocoords3s[infovec==1,], pch=16)
    } 
    else if(is.null(dim(geocoords3s[infovec==1,]))){
      points(geocoords3s[infovec==1,][1], geocoords3s[infovec==1,][2], pch=16)
    }
  }

if(plotmp){
    plot(geocoords3s, main='Bioentity establishment, shaded = yes')
    if(!is.null(dim(geocoords3s[estabvec==1,]))){
      points(geocoords3s[estabvec==1,], pch=16)
    } 
    else if(is.null(dim(geocoords3s[estabvec==1,]))){
      points(geocoords3s[estabvec==1,][1], geocoords3s[estabvec==1,][2], pch=16)
    }
  }

  list(com.yes=infout$com.yes, obschange=infout$obschange, geocoords3s=geocoords3s, infovec=infovec, estabvec=estabvec)
}