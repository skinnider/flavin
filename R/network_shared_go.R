#' Calculate the number of shared GO terms between each pair of genes or 
#' proteins in a network. 
#' 
#' @param network a data.frame with nodes (genes) in the first two columns
#' @param shared.go a matrix containing the number of shared GO terms between 
#' each pair of proteins (as output by `proteome_shared_go`)
#' @param keep.all if false, only results for network edges where both nodes
#' were found in the shared GO matrix will be reported 
#' @return the input network, with an extra column containing the number of 
#' shared GO terms
#' @export
#' @examples 
#' # read network
#' net <- read.delim("mouse-network.tsv.gz")
#' # get proteome-wide shared GO terms
#' shared <- proteome_shared_go("goa_mouse.gpa.gz", database = org.Mm.eg.db,
#'    accession = "ENSEMBL")
#' # calculate number of shared GO terms for each edge the network
#' net.shared <- network_shared_go(net, shared)
network_shared_go <- function(network, shared.go, keep.all = F) {
  inSharedIndices <- network[,1] %in% rownames(shared.go) & 
    network[,2] %in% rownames(shared.go)
  if (sum(inSharedIndices) == 0)
    stop("no nodes in the network were found in the shared GO terms matrix")
  in.shared <- network[inSharedIndices,]
  rowIndices <- match(in.shared[,1], rownames(shared.go))
  colIndices <- match(in.shared[,2], colnames(shared.go))
  in.shared$shared <- shared.go[cbind(rowIndices, colIndices)]
  if (keep.all) {
    not.in.shared <- network[!inSharedIndices,]
    not.in.shared$shared <- 0
    return(rbind(in.shared, not.in.shared))
  } else {
    return(in.shared)  
  }
}
