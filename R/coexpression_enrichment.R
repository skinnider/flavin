#' Calculate a Z score for enrichment for coexpression
#' 
#' Calculate the enrichment among a set of interacting proteins for mRNA or 
#' protein co-expression, relative to a rewired version of the same network. 
#' 
#' @param network a data frame containing interactions between protein pairs
#' in the first two columns
#' @param expr a matrix containing mRNA or protein expression data, with 
#' genes as columns
#' @param bootstraps the number of times to rewire the network in order
#' to calculate the Z score; defaults to 100
#' @param seed the number to seed the rewiring algorithm with; defaults to zero
#' @param colname the column name to assign to the shared annotations 
#' column created in the network data frame; defaults to "coexpr"
#' 
#' @return a list with three objects: the median coexpression in the 
#' observed network; the median coexpression in each randomized network;
#' and the Z score for the enrichment
#' 
#' @export
coexpression_enrichment <- function(network, expr, bootstraps = 100,
                                         seed = 0, colname = "coexpr") {
  # calculate observed median coexpression
  network <- network_coexpression(network, expr, col_name = colname)
  obs <- median(network[[colname]], na.rm = T)
  
  # convert network to graph
  g <- igraph::graph_from_data_frame(network, directed = F)
  m <- length(igraph::E(g))
  n_iterations <- 6.9077 * m ## from https://arxiv.org/pdf/1202.3473.pdf
  set.seed(seed) ## set seed for reproducibility
  rewires <- map(seq_len(bootstraps), ~ igraph::rewire(
    g, with = igraph::keeping_degseq(niter = n_iterations)))
  rnd_shared <- map(rewires, ~ network_coexpression(
    igraph::as_data_frame(.), expr))
  rnds <- purrr:: map(rnd_shared, colname)
  rnd <- purrr::map_dbl(rnds, ~ median(., na.rm = T))
  
  # calculate Z score
  Z <- (obs - mean(rnd)) / sd(rnd)
  
  # return results
  results <- list(
    observed = obs,
    rewired = rnd,
    Z = Z)
  return(results)
}