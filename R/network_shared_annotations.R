#' Calculate shared annotations between network nodes
#' 
#' Calculate the number of shared annotations between each pair of genes or 
#' proteins in a network, and append this value to the input data frame
#' as an extra column. 
#' 
#' @param network a data.frame with nodes (genes) in the first two columns
#' @param shared a matrix containing the number of shared annotations, e.g. 
#' GO terms, between each pair of proteins (as output by `shared_annotations`)
#' @param col_name the name of the new column to create; defaults to "coexpr"
#' 
#' @return the input network, with an extra column containing the number of 
#' shared annotations (e.g., GO terms)
#' 
#' @export
#'
#' @examples 
#' # read network
#' net <- read.delim("mouse-network.tsv.gz")
#' # get proteome-wide shared annotations
#' shared <- shared_annotations(ann)
#' # calculate number of shared GO terms for each edge the network
#' net.shared <- network_shared_go(net, shared)
network_shared_annotations <- function(network, shared, col_name = "shared") {
  # make sure there is some overlap between network and expression data
  inSharedIndices <- network[,1] %in% rownames(shared) & 
    network[,2] %in% rownames(shared)
  if (sum(inSharedIndices) == 0)
    stop("no nodes in the network were found in the shared annotations matrix")
  # create empty column
  network[[col_name]] <- NA
  # add values for nodes in shared matrix
  in.shared <- network[inSharedIndices,]
  rowIndices <- match(in.shared[,1], rownames(shared))
  colIndices <- match(in.shared[,2], colnames(shared))
  network[[col_name]][inSharedIndices] <- shared[cbind(rowIndices, colIndices)]
  return(network)
}
