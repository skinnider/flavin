#' Identify shared subcellular localization terms between interacting proteins.
#' 
#' Identify interacting protein pairs that share one or more subcellular 
#' localization annotations. Optionally, the maximum breadth of subcellular
#' localization annotations within a set of interactions can be restricted, to
#' avoid considering very universal terms. 
#' 
#' @param interactions A list of interactions, where the first two columns
#' correspond to the accessions of the interacting proteins. 
#' @param localization A dataset of protein subcellular localization, where the 
#' first column contains the accession of the protein and the second column 
#' contains one or more localization annotations. 
#' @param breadth Optional: specify the maximum fraction of proteins in the 
#' network that can be annotated with a subcellular localization term in order 
#' to consider it.
#' @param delim Optional: specify the delimiter separating subcellular 
#' localization annotations in the second column. Defaults to space. 
#' 
#' @return The input table of interactions, with the shared annotations and 
#' number of shared annotations in new columns
#' 
#' @examples 
#' interactions <- read.csv("interactions.csv")
#' localization <- read.csv("subcellular_localization.csv")
#' shared.localization <- GetSharedSubcellularLocalization(interactions, 
#' expression)
#' # Exclude all terms that are annotated to more than 25% of input proteins
#' shared.localization.25 <- GetSharedSubcellularLocalization(interactions, 
#' expression, breadth = 0.25)
#' 
#' @export 
GetSharedSubcellularLocalization <- function(interactions, localization, 
                                                   breadth = 0, delim = " ") {
  proteins <- unique(c(interactions[,1], interactions[,2]))
  # Ignore proteins that can't be mapped to subcellular localization
  n <- length(proteins)
  proteins <- proteins[proteins %in% localization$UniProt]
  message("Retrieved subcellular localization for ", length(proteins), " of ", 
          n, " proteins")
  # Ignore annotations that don't correspond to network proteins 
  localization <- subset(localization, UniProt %in% proteins)
  if (breadth > 0) {
    # Get term breadth
    terms <- unique(unlist(strsplit(localization[,2], delim)))
    b <- data.frame(term = terms)
    b$count <- unlist(sapply(b$term, function(term) 
      sum(grepl(term, localization[,2]))))
    # Define breadth cutoff 
    cutoff <- breadth * length(proteins)
    # Get terms that are above breadth cutoff 
    too.broad <- b$term[b$count >= cutoff]
    # Remove terms from 
    localization[,2] <- unlist(sapply(localization[,2], function(l) {
      terms <- unlist(strsplit(l, delim))
      paste(terms[!(terms %in% too.broad)], collapse = " ")))
  }
  
  
  # Apply breadth cutoff
  if (BREADTH > 0) {
    cutoff <- BREADTH * nrow(df)
    too.broad <- breadth$term[breadth$count >= cutoff]
    breadth.subset <- function(go) {
      terms <- strsplit(go, " ")[[1]]
      terms <- terms[!(terms %in% too.broad)]
      return(paste(terms, collapse=" "))
    }
    df$localization <- unlist(sapply(df$localization, breadth.subset))
    df <- df[df$localization != "",]
    message(nrow(df), " had a subcellular localization term present in < ",
            100*BREADTH, "% of the network")
  }
}