#' Filter GO terms from a particular ontological category (BP, CC, or MF) from 
#' a GOA file. 
#'
#' @param go the GPA file, read with `read_gpa`
#' @param ontology.file location of the GO ontology file 
#' @param root the ontological category to restrict the GOA file to; one of
#' "BP", "CC", or "MF"
#' @return the filtered GPA file 
#' @export 
#' @examples 
#' # read a GOA file 
#' go <- read_gpa("goa_human.gpa.gz")
#' # get BP terms only
#' bp <- filter_roots(go, ontology.file = "go-basic.obo.gz", root = "BP")
filter_roots <- function(goa, ontology.file, root = c("BP", "CC", "MF")) {
  message("filtering GOA data for terms with root ", root, "...")
  if (!root %in% names(root))
    stop("invalid GO root ", root)
  # define roots
  rootNames <- c(BP = "GO:0008150", CC = "GO:0005575", MF = "GO:0003674")
  # read the ontology
  ontology <- ontologyIndex::get_ontology(ontology.file, 
                                          extract_tags = "minimal")
  # get roots 
  ids <- unique(goa$GO.ID)
  ancestors <- sapply(ids, function(x) ontology$ancestors[x])
  roots <- lapply(ancestors, function(x) x[x %in% rootNames])
  names(roots) <- ids
  # remove terms not in the target root
  terms <- names(roots)[roots %in% rootNames[root]]
  # filter GOA
  goa <- goa[goa$GO.ID %in% terms,]
  return(goa)
}
