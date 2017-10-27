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
  Matrix::crossprod(xtabs( ~ values + ind, utils::stack(ann), sparse = T))
}
