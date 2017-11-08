#' Calculate enrichment for coexpression among interacting proteins
#' 
#' This function calculates the difference between the median coexpression 
#' of interacting protein pairs within a network, and random pairs from the same 
#' expression matrix. 
#' 
#' @param network a matrix or data frame, with nodes (genes) in the first 
#' two columns
#' @param expr mRNA or protein expression data for 
#' @param ... further arguments passed directly to the cor function, such as
#' coefficient (e.g., Pearson or Spearman)
#' @return the difference in medians between interacting and random protein 
#' pairs
#' 
#' @export
network_coexpression <- function(network, expr, ...) {
  # make sure there is some overlap between network and genes 
  genes <- colnames(expr)
  overlap <- network[,1] %in% genes & network[,2] %in% genes
  if (sum(overlap) == 0)
    stop("no nodes in the network were found in the expression matrix")
  # calculate coexpression matrix
  coexpr <- cor(expr, ...)
  # calculate median coexpression
  median <- median(coexpr, na.rm = T)
  # calculate median coexpression of network interactions
  subset <- network[network[,1] %in% genes & network[,2] %in% genes,]
  networkCoexpr <- coexpr[cbind(subset[,1], subset[,2])]
  networkMedian <- median(networkCoexpr, na.rm = T)
  # return the difference
  return(networkMedian - median)
}
