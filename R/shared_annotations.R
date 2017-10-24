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
#'    


#' Construct a matrix of shared annotation counts
#' 
#' Given an annotation list, where the name of the gene or protein corresponds
#' to the names of the list, and the annotations (e.g., GO term) are items
#' within the list, this function constructs a square matrix wherein the gene or 
#' protein identifiers are the row and column names of the matrix, and the 
#' entries within the matrix are the total number of annotations (e.g., 
#' GO terms) shared by each pair of genes or proteins. 
#' 
#' @param ann a list of annotations, as described above
#' @return a sparse matrix containing the total number of shared annotations
#' between each unique gene or protein pair
shared_annotations <- function(ann) {
  Matrix::crossprod(xtabs( ~ values + ind, stack(ann), sparse = T))
}
