#' Filter GO terms from a particular ontological category (BP, CC, or MF) from 
#' a GOA file. 
#'
#' @param go the GO annotation data frame
#' @param ontology the GO ontology, as read with the \code{ontologyIndex} 
#' package 
#' @param root the ontological category to restrict the GOA file to; one of
#' "BP", "CC", or "MF"
#' @return the filtered GO annotation file 
#' @export 
filter_roots <- function(goa, ontology, root = c("BP", "CC", "MF")) {
  message("filtering GOA data for terms with root ", root, "...")
  root <- match.arg(root)
  # define roots
  rootNames <- c(BP = "GO:0008150", CC = "GO:0005575", MF = "GO:0003674")
  # get roots 
  col_name <- ifelse("GO.ID" %in% colnames(goa), "GO.ID", "GO_ID")
  ids <- unique(goa[[col_name]])
  ancestors <- sapply(ids, function(x) ontology$ancestors[x])
  roots <- lapply(ancestors, function(x) x[x %in% rootNames])
  names(roots) <- ids
  # remove terms not in the target root
  terms <- names(roots)[roots %in% rootNames[root]]
  # filter GOA
  goa <- goa[goa[[col_name]] %in% terms,]
  return(goa)
}
