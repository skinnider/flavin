#' Convert a data frame to an adjacency matrix
#' 
#' Convert a data frame of pairwise interactions into an adjacency matrix
#' between all objects. 
#' 
#' @param dat a data frame with the objects that are interacting in the 
#' first two columns
#' 
#' @return a square adjacency matrix with all unique values of the first 
#' two columns of \code{dat} as row and column names 
#' 
#' @export
to_adjacency_matrix <- function(dat) {
  nodes <- unique(c(dat[[1]], dat[[2]]))
  adj <- matrix(0, nrow = length(nodes), ncol = length(nodes),
                   dimnames = list(nodes, nodes))
  adj[cbind(dat[[1]], dat[[2]])] <- 1
  adj[cbind(dat[[2]], dat[[1]])] <- 1
  return(adj)
}