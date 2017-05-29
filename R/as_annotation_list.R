# as_annotation_list(goa, "GO.ID", "UNIPROT")

#' Convert an annotation data.frame to a named list, where the names of the 
#' items are the annotated targets and the entries are unique annotations. 
#'
#' @param annotations a data frame containing some annotations
#' @param termCol the name of the column that contains the annotation terms 
#' (e.g. Gene Ontology terms)
#' @param keyCol the name of the column that contains the annotated targets 
#' (e.g. UniProt accessions)
#' @return a named list where the names are found in keyCol and the entries are 
#' unique items from termCol
#' @export
#' @examples 
#' # read a GOA file 
#' go <- read_gpa("goa_human.gpa.gz")
#' # convert to a named annotation list 
#' ann <- as_annotation_list(go, "GO.ID", "UNIPROT")
as_annotation_list <- function(annotations, termCol, keyCol) {
  message("aggregating annotations to the same source node ...")
  annList <- unstack(annotations[, c(termCol, keyCol)])
  ann <- sapply(annList$x, function(x) unique(x))
  names(ann) <- annList[,1]
  return(ann)
}
