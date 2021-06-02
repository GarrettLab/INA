
#' Function for plotting multilayer networks
#'
#' This function plots a multilayer network.  It assumes that the x and y coordinates are in [0,1]

#' @param upmat adjacency matrix for upper network
#' @param lowmat adjacency matrix for lower network
#' @param upcoords x-y coordinates for upper network
#' @param lowcoords x-y coordinates for lower network
#' @param yflat2 proportional flattening of map
#' @param xover2 proportional right movement of map
#' @param htdif2 distance between upper and lower network
#' @keywords plotting networks multilayer
#' @export 
#' @examples
#' y1 <- matrix(c(0,0,0,1, 0,0,1,0, 0,0,0,1, 0,0,0,0),byrow=T, ncol=4)
#' y2 <- matrix(c(0,1,0,0,0, 0,0,1,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,1,0),byrow=T, ncol=5)
#' x1 <- matrix(c(0.5,0.8, 0.9,0.9, 0.1,0.1, 0.5,0.5), byrow=T, ncol=2)
#' x2 <- matrix(c(0.5,0.8, 0.9,0.9, 0.1,0.1, 0.5,0.5, 0.3,0.8), byrow=T, ncol=2)
#' layerplot(upmat=y1,lowmat=y2,upcoords=x1,lowcoords=x2,yflat=0.5,xover=0.5, htdif=1)

# to do - GENERAL TESTING


### layerplot function

layerplot <- function(upmat, lowmat, upcoords, lowcoords, yflat2, xover2, htdif, lwdhere=1){
 
# position the corners of the two plots

  lowcorners <- matrix(c(0,0, 1,0, 0,1, 1,1), byrow=T, ncol=2)
  upcorners <- matrix(c(0,0, 1,0, 0,1, 1,1), byrow=T, ncol=2)

  lowcorners2 <- squish(lowcorners, yflat=yflat2, xover=xover2)
  
  upcorners2 <- squish(upcorners, yflat=yflat2, xover=xover2)
  upcorners2[,2] <- upcorners2[,2] + htdif

# find the limits and produce an empty plot for filling in

  maxx <- max(c(lowcorners2[,1],upcorners2[,1]))
  minx <- min(c(lowcorners2[,1],upcorners2[,1]))
  maxy <- max(c(lowcorners2[,2],upcorners2[,2]))
  miny <- min(c(lowcorners2[,2],upcorners2[,2]))

  plot(1, type="n", axes=F, xlab="", ylab="", xlim=c(minx,maxx), ylim=c(miny,maxy))

# add the lines outlining the plot areas
  lines(c(lowcorners2[1,1],lowcorners2[2,1],lowcorners2[4,1],lowcorners2[3,1],lowcorners2[1,1]), c(lowcorners2[1,2],lowcorners2[2,2],lowcorners2[4,2],lowcorners2[3,2],lowcorners2[1,2]), lwd=lwdhere)

lines(c(upcorners2[1,1],upcorners2[2,1],upcorners2[4,1],upcorners2[3,1],upcorners2[1,1]), c(upcorners2[1,2],upcorners2[2,2],upcorners2[4,2],upcorners2[3,2],upcorners2[1,2]), lwd=lwdhere)  

# add the points in each plot

  lowcoords2 <- squish(lowcoords, yflat=yflat2, xover=xover2)
  
  upcoords2 <- squish(upcoords, yflat=yflat2, xover=xover2)
  upcoords2[,2] <- upcoords2[,2] + htdif
  
  points(lowcoords2)
  points(upcoords2)

# add lines within plots

for (i in 1:dim(upmat)[1]) {
  for (j in 1:dim(upmat)[1]) {
    if(upmat[i,j] == 1) {
      lines(c(upcoords2[i,1],upcoords2[j,1]), c(upcoords2[i,2],upcoords2[j,2]), lwd=lwdhere)
    }
  }
}

for (i in 1:dim(lowmat)[1]) {
  for (j in 1:dim(lowmat)[1]) {
    if(lowmat[i,j] == 1) {
      lines(c(lowcoords2[i,1],lowcoords2[j,1]), c(lowcoords2[i,2],lowcoords2[j,2]), lwd=lwdhere)
    }
  }
}

# add lines between layers - assuming link exists if same coordinates

for (i in 1:dim(lowmat)[1]) {
  for (j in 1:dim(upmat)[1]) {
    if(lowcoords[i,1] == upcoords[j,1] & lowcoords[i,2] == upcoords[j,2]) {
      lines(c(lowcoords2[i,1],upcoords2[j,1]), c(lowcoords2[i,2],upcoords2[j,2]),lty=2, lwd=lwdhere)
    }
  }
}

  list(lowcoords2=lowcoords2, upcoords2=upcoords2)
}