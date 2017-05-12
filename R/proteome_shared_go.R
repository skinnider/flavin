#' Rapidly identify shared GO terms for all pairs of proteins in a proteome. 
#'
#' @param filepath the location of the GPA file. Files can be downloaded from
#' ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/
#' @param database optionally, provide a Bioconductor `AnnotationDbi` database
#' to map identifiers from UniProt to some other accession 
#' @param accession optionally, the type of identifier to map UniProt accessions
#' to using a Bioconductor `AnnotationDbi` database. Must be a keytype in the
#' database. 
#' @param filter_NOT if true, filter annotations with the qualifier NOT
#' @param filter_evidence optionally, specify evidence codes to filter. By 
#' default, evidence codes ND, IPI, IEA and NAS are filtered. 
#' @param filter.breadth if true, filter annotations that are very specific or 
#' very broad 
#' @param min.breadth minimum size of a GO group to consider
#' @param max.breadth maximum size of a GO group to consider
#' @param ontology.file optionally, provide a GO ontology file in order to 
#' restrict analysis to a specific ontological category (BP, CC, or MF)
#' @param ontology if an ontology file is specified, the GO root to consider 
#' (BP, CC, or MF)
#' @return a `dsCMatrix` containing all 
#' @export
#' @examples 
#' # get all shared GO terms between mouse Ensembl genes, considering only 
#' # BP groups of size 20-100 proteins
#' shared <- proteome_shared_go("goa_mouse.gpa.gz", database = org.Mm.eg.db,
#'    accession = "ENSEMBL", filter.breadth = T, max.breadth = 100),
#'    ontology.file = "go-basic.obo.gz", root = "BP")
proteome_shared_go <- function(gpa.file,
                      database = NULL, accession = "UNIPROT", 
                      filter.NOT = T, 
                      filter.evidence = c("ND", "IPI", "IEA", "NAS"),
                      filter.breadth = T,
                      min.breadth = 20, max.breadth = 1000, 
                      ontology.file = NULL, ontology = c("BP", "CC", "MF")) {
  # read GPA file 
  goa <- read_gpa(gpa.file, accession, database, filter.NOT, filter.evidence)
  
  # filter based on term breadth
  if (filter.breadth) {
    goa <- filter_breadth(goa, min.breadth, max.breadth)
  }
  # filter based on ontology
  if (!is.null(ontology.file)) {
    goa <- filter_roots(goa, ontology.file, ontology)
  }
  
  # aggregate by term
  go <- as_annotation_list(goa, "GO.ID", accession)
  
  # get shared terms 
  message("calculating shared terms for ", length(go), " x ", length(go),
          " protein pairs...")
  shared <- Matrix::crossprod(xtabs( ~ values + ind, stack(go), 
                                     sparse = T)) 
  return(shared)
}
