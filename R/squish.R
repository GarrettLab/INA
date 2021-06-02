
#' Function for plotting network on an angle
#'
#' This function plots a network on the angle specified, for layering plots

#' @param coords matrix of coordinates, with each row a node, the first column the x coordinate, and the second column the y 
#' @param yflat proportional flattening of map
#' @param xover proportional right movement of map
#' @keywords plotting networks multilayer
#' @export 
#' @examples
#' Z1 <- matrix(c(0.5,0.8, 0.9,0.9, 0.1,0.1, 0.5,0.5), byrow=T, ncol=2)
#' Z1
#' plot(Z1, xlim=c(0,2), ylim=c(0,1))
#' Z2 <- squish(Z1)
#' plot(Z2, xlim=c(0,2), ylim=c(0,1))

# to do - GENERAL TESTING

squish <- function(coords, yflat=0.5, xover=0.5){

  coords2 <- coords
  coords2[,1] <- coords[,1] + xover*coords[,2] # x
  coords2[,2] <- yflat * coords[,2] # y

  coords2
}
