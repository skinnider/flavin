#' Calculate the number of shared annotations between each pair of genes or 
#' proteins in a network. 
#' 
#' @param network a data.frame with nodes (genes) in the first two columns
#' @param shared a matrix containing the number of shared annotations, e.g. 
#' GO terms, between each pair of proteins (as output by `shared_annotations`)
#' @param keep.all if false, only results for network edges where both nodes
#' were found in the shared annotations matrix will be reported 
#' @return the input network, with an extra column containing the number of 
#' shared annotations (e.g., GO terms)
#' @export
#' @examples 
#' # read network
#' net <- read.delim("mouse-network.tsv.gz")
#' # get proteome-wide shared annotations
#' shared <- shared_annotations(ann)
#' # calculate number of shared GO terms for each edge the network
#' net.shared <- network_shared_go(net, shared)
network_shared_annotations <- function(network, shared, keep.all = F) {
  inSharedIndices <- network[,1] %in% rownames(shared) & 
    network[,2] %in% rownames(shared)
  if (sum(inSharedIndices) == 0)
    stop("no nodes in the network were found in the shared GO terms matrix")
  in.shared <- network[inSharedIndices,]
  rowIndices <- match(in.shared[,1], rownames(shared))
  colIndices <- match(in.shared[,2], colnames(shared))
  in.shared$shared <- shared[cbind(rowIndices, colIndices)]
  if (keep.all) {
    not.in.shared <- network[!inSharedIndices,]
    not.in.shared$shared <- 0
    return(rbind(in.shared, not.in.shared))
  } else {
    return(in.shared)  
  }
}
