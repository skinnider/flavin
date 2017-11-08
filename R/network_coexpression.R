#' Calculate coexpression of network nodes
#' 
#' Calculate the mRNA or protein coexpression of each pair of genes or 
#' proteins in a network, and append this value to the input data frame
#' as an extra column.
#' 
#' @param network a matrix or data frame, with nodes (genes) in the first 
#' two columns
#' @param expr a matrix containing mRNA or protein expression data, with 
#' genes as columns
#' @param col_name the name of the new column to create; defaults to "coexpr"
#' @param ... further arguments passed directly to the cor function, such as
#' coefficient (e.g., Pearson or Spearman)
#' 
#' @return the input network, with an extra column containing the 
#' coexpression of each interacting gene or protein pair
#' 
#' @export
network_coexpression <- function(network, expr, col_name = "coexpr", ...) {
  # make sure there is some overlap between network and expression data  
  exprGenes <- colnames(expr)
  overlap <- network[,1] %in% exprGenes & network[,2] %in% exprGenes
  if (sum(overlap) == 0)
    stop("no nodes in the network were found in the expression matrix")
  # get the network genes in the expression matrix, for faster calculations
  subset <- network[overlap,]
  networkGenes <- unique(c(subset[, 1], subset[, 2]))
  # calculate full coexpression matrix
  expr <- expr[, genes %in% networkGenes]
  coexpr <- cor(expr, ...)
  # find network edges in coexpression matrix
  network[[col_name]] <- NA
  network[[col_name]][overlap] <- coexpr[cbind(subset[,1], subset[,2])]
  return(network)
}
