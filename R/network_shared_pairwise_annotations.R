#' Calculate shared pairwise annotations between network nodes
#' 
#' Calculate the number of shared pairwise annotations between each pair of 
#' genes or proteins in a network, and append this value to the input data frame
#' as an extra column. By "paired," it is meant that the annotations apply
#' to pairs of proteins: for instance, interactions between Pfam domains 
#' derived from analysis of crystal structures (as curated by 3did and other
#' databases). 
#' 
#' @param network a data.frame with nodes (genes) in the first two columns
#' @param annotations an annotation list, created by 
#' \code{\link{as_annotation_list}}, with nodes as names and annotations as 
#' entries in the list; for example, a list associating proteins with their
#' Pfam domains 
#' @param annotation_pairs a logical or binary matrix containing pairs of 
#' annotations that are considered to be shared; for example, a matrix of
#' Pfam domains, where known domain-domain interactions are marked with
#' ones
#' @param col_name the name of the new column to create; defaults to 
#' "shared_pairs"
#' 
#' @return the input network, with an extra column containing the number of 
#' shared pairwise annotations (e.g., domain-domain interactions)
#' 
#' @export
network_shared_pairwise_annotations <- function(network, annotations, 
                                                annotation_pairs,
                                                col_name = "shared_pairs") {
  # make sure there is some overlap between network and annotations data
  annotatedNodes <- names(annotations)
  sharedIdxs <- network[[1]] %in% annotatedNodes & 
    network[[2]] %in% annotatedNodes
  if (sum(sharedIdxs) == 0)
    stop("no nodes in the network were found in the annotation list")
  # create empty column
  network[[col_name]] <- NA
  # calculate values for nodes in pairwise annotation matrix
  subset <- network[sharedIdxs,]
  grid <- data.frame(col1 = I(annotations[subset[[1]]]),
                     col2 = I(annotations[subset[[2]]]))
  shared_pairs <- purrr::pmap_dbl(grid, ~ sum(annotation_pairs[.x, .y])) 
  network[[col_name]][sharedIdxs] <- shared_pairs
  return(network)
}

