#' Calculate a Z score for enrichment for shared annotations
#' 
#' Calculate the enrichment among a set of interacting proteins for shared
#' annotations, relative to a rewired version of the same network. 
#' 
#' @param network a data frame containing interactions between protein pairs
#' in the first two columns
#' @param shared a matrix denoting the number of shared annotations between
#' each protein pair
#' @param bootstraps the number of times to rewire the network in order
#' to calculate the Z score; defaults to 100
#' @param seed the number to seed the rewiring algorithm with; defaults to zero
#' @param colname the column name to assign to the shared annotations 
#' column created in the network data frame; defaults to "shared"
#' 
#' @return a list with three objects: the proportion of interacting 
#' protein pairs that share GO terms in the observed network; the proportions
#' in each randomized network; and the Z score for the enrichment
#' 
#' @export
shared_annotation_enrichment <- function(network, shared, bootstraps = 100,
                                         seed = 0, colname = "shared") {
  # calculate observed shared annotations
  network <- network_shared_annotations(network, shared, col_name = colname)
  obs <- mean(network[[colname]] > 0, na.rm = T)
  
  # convert network to graph
  g <- igraph::graph_from_data_frame(network, directed = F)
  m <- length(igraph::E(g))
  n_iterations <- 6.9077 * m ## from https://arxiv.org/pdf/1202.3473.pdf
  set.seed(seed) ## set seed for reproducibility
  rewires <- map(seq_len(bootstraps), ~ igraph::rewire(
    g, with = igraph::keeping_degseq(niter = n_iterations)))
  rnd_shared <- map(rewires, ~ tryCatch(
    network_shared_annotations(igraph::as_data_frame(.), shared),
    error = function(e) return(0)))
  rnds <- purrr::map(rnd_shared, colname)
  rnd <- purrr::map_dbl(rnds, ~ mean(. > 0, na.rm = T))
  
  # calculate Z score
  Z <- (obs - mean(rnd)) / sd(rnd)
  
  # return results
  results <- list(
    observed = obs,
    rewired = rnd,
    Z = Z)
  return(results)
} 