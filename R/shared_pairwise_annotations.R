#' Construct a matrix of shared pairwise annotation counts
#' 
#' This function is related to \code{\link{shared_annotations}}, but whereas 
#' that function calculates the total number of annotations shared by each
#' pair of genes or proteins, this function calculates the total number of 
#' pairwise annotations shared by each gene or protein pair. These are 
#' annotations that are only interpretable in a pairwise manner: for instance,
#' interactions between Pfam domains. 
#' 
#' @param annotations a list of annotations, with the names of the list 
#' corresponding to nodes in the network and entries corresponding to 
#' annotations
#' @param annotation_pairs a binary or logical matrix 
#' 
#' @return a matrix containing the total number of shared pairwise annotations
#' between each unique gene or protein pair
#' 
#' @export
shared_pairwise_annotations <- function(annotations, 
                                        annotation_pairs) {
  # remove any annotations not present in the pairwise matrix
  annotations <- annotations %>%
    stack() %>%
    dplyr::filter(values %in% colnames(annotation_pairs)) %>%
    unstack()
  if (length(unlist(annotations)) == 0)
    stop("no annotations were found in the pairwise annotation matrix")
  # calculate shared pairwise annotations
  n <- length(annotations)
  shared <- matrix(nrow = n, ncol = n)
  grid <- expand.grid(annotations, annotations)
  shared[] <- purrr::pmap_dbl(grid, ~ sum(annotation_pairs[.x, .y])) 
  colnames(shared) <- rownames(shared) <- names(annotations)
  return(shared)
}
