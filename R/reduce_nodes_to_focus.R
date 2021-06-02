
#' Reduces the number of nodes in a network to those which are linked to a set of key nodes
#'
#' There is the option to select all nodes linked to the set of key nodes within a certain number of steps, or only in one direction
#'
#' Created 2020-05-29

#' @param startmat the starting adjacency matrix (NOT an igraph object), to be reduced
#' @param maxdist the greatest distance (in terms of the number of links traversed) at which a node is considered connected to the key nodes
#' @param keynodes the set of key nodes (e.g., c(1,5,9))
#' @param tofrom choice between keeping nodes that are (1) linked towards the key nodes, 'to', (2) linked away from the key nodes, 'from', or (3) either to or from, 'either'
#' @keywords visualization
#' @export 
#' @examples
#' j <- matrix(rbinom(n=100,p=0.1,size=1),nrow=10)
#' diag(j) <- 0
#' jnames <- letters[1:10]
#' library(igraph)
#' ji <- graph_from_adjacency_matrix(j)
#' V(ji)$names <- jnames
#' plot.igraph(ji, vertex.label=V(ji)$names, vertex.color='light green', main='Starting matrix')
#' 
#' jout <- reduce_nodes_to_focus(startmat=j, maxdist=2, keynodes=c(1,5), tofrom='to')
#' jouti <- graph_from_adjacency_matrix(jout$outmat)
#' V(jouti)$names <- jnames[jout$keepers]
#' plot(jouti, vertex.label=V(jouti)$names, vertex.color='light green', main='2 steps TO key nodes 1 and 5')
#'
#' jout2 <- reduce_nodes_to_focus(startmat=j, maxdist=2, keynodes=c(1,5), tofrom='from')
#' jout2i <- graph_from_adjacency_matrix(jout2$outmat)
#' V(jout2i)$names <- jnames[jout2$keepers]
#' plot(jout2i, vertex.label=V(jout2i)$names, vertex.color='light green', main='2 steps FROM key nodes 1 and 5')
#'
#' jout3 <- reduce_nodes_to_focus(startmat=j, maxdist=2, keynodes=c(1,5), tofrom='either')
#' jout3i <- graph_from_adjacency_matrix(jout3$outmat)
#' V(jout3i)$names <- jnames[jout3$keepers]
#' plot(jout3i, vertex.label=V(jout3i)$names, vertex.color='light green', main='2 steps TO OR FROM key nodes 1 and 5')

# Note - need to confirm how igraph function distances deals with 'infinite' distances - this will influence how the parameter maxdist needs to be set up, when nodes with any possible path to or from key nodes are desired

reduce_nodes_to_focus <- function(startmat, maxdist, keynodes, tofrom){
  
  library(igraph)
  # make a version of the starting adjacency matrix as
  #   an igraph adjacency matrix
  startmati <- graph_from_adjacency_matrix(j)

  # construct matrix of distances between nodes
  #  (in number of steps)
  distmat <- distances(startmati, mode="out")

  # construct matrix of indicators for whether the distance
  #	between node pairs is <= the maximum distance being
  #	considered
  distmat2 <- (distmat <= maxdist) * 1

  # obtain the rows of the matrix corresponding to movement
  #   FROM the key nodes and TO the key nodes
  distmat2FROM <- distmat2[keynodes,,drop=F]
  distmat2TO <- distmat2[,keynodes,drop=F]

  # determine which nodes to keep based on movement 
  #   FROM key nodes and TO key nodes
  keepersFROM <- apply(distmat2FROM,2,sum)
  keepersTO <- apply(distmat2TO,1,sum)

  # if only considering the nodes on paths FROM the key nodes,
  #   these are the nodes to keep
  if (tofrom == 'from'){
    keepers <- keepersFROM > 0
  
  # if only considering the nodes on paths TO the key nodes,
  #   these are the nodes to keep
  } else if (tofrom == 'to'){
    keepers <- keepersTO > 0
  
  # if consider either from or to the key nodes,
  #   these are the nodes to keep
  } else if (tofrom == 'either'){
    keeperseither <- keepersTO + keepersFROM
    keepers <- keeperseither > 0
  }
  
  # make output matrix with only those nodes close enough to 
  #   key nodes
  outmat <- startmat[keepers,keepers,drop=F]

  list(outmat=outmat,keepers=keepers)

}