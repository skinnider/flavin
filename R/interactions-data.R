#' Sample protein-protein interaction data
#' 
#' Interactions predicted by PRInCE, with a precision greater than or equal to 
#' 90%, from Scott et al.'s apoptosis dataset. Columns 1 and 2 contain the 
#' UniProt accessions for each interacting protein pair. Column 3 contains the 
#' calculated precision of the interaction (see Stacey et al. for details on 
#' calculation of precision in PRInCE).
#' 
#' @docType data 
#' 
#' @usage data(interactions)
#' 
#' @format a data frame with three columns 
#' 
#' @references Scott et al.
#' Stacey et al. 
#' 
"interactions"