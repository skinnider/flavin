#' Filter very broad or very specific GO terms from a GOA file. 
#'  
#' @param go the GOA file, read with `read_gpa`
#' @param min.breadth minimum size of a GO group to consider
#' @param max.breadth maximum size of a GO group to consider
#' @return the filtered GPA file
#' @export
#' @examples 
#' # read a GOA file 
#' go <- read_gpa("goa_human.gpa.gz")
#' # get GO terms annotated to between 30 and 100 proteins 
#' bp <- filter_breadth(go, min.breadth = 30, max.breadth = 100)
filter_breadth <- function(goa, min.breadth = 20, max.breadth = 1000) {
  message("filtering GOA data for terms with minimum breadth ", min.breadth, 
          " and maximum breadth ", max.breadth, "...")
  # calculate GO term frequency
  if (!"GO.ID" %in% colnames(goa) & !"GO_ID" %in% colnames(goa))
    stop("Couldn't find GO term column (expecting one of: GO.ID, GO_ID)")
  colname <- ifelse("GO.ID" %in% colnames(goa), "GO.ID", "GO_ID")
  freq <- table(goa[[colname]])
  # exclude narrow or broad terms
  freq <- freq[freq >= min.breadth & freq <= max.breadth]
  # now filter the GOA
  goa <- goa[goa[[colname]] %in% names(freq),]
  return(goa)  
}
