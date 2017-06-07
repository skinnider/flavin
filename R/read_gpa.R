#' Read a GO annotation file in GPA format. 
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
#' @param propatate if true, and an ontology file is provided, all ancestors of 
#' a given term are associated with each protein.
#' @return a data.frame containing the filtered GPA file  
#' @export
#' @examples 
#' # read human GOA and map to Ensembl, filtering IEA annotations
#' go <- read_gpa("goa_human.gpa.gz", database = org.Hs.eg.db,
#'    accession = "ENSEMBL", filter.evidence = "IEA")
#' # read mouse GOA and map to gene symbol, filtering IEA and IPI annotations
#' go <- read_gpa("goa_mouse.gpa.gz", database = org.Mm.eg.db,
#'    accession = "SYMBOL", filter.evidence = c("IEA", "IPI"))
read_gpa <- function(filepath,
                     database = NULL, accession = "UNIPROT", 
                     filter.NOT = T, 
                     filter.evidence = c("ND", "IPI", "IEA", "NAS"),
                     ontology.file = NULL, propagate = T) {
  message("reading GOA file ", filepath, "...")
  # read GOA file 
  gpa.colnames <- c("DB", "DB_Object_ID", "Qualifier", "GO.ID", "DB.Reference",
                    "ECO.evidence.code", "With.From", "Interacting.taxon.ID",
                    "Date", "Assigned_by", "Annotation.Extension", 
                    "Annotation.Properties")
  goa <- read.table(filepath, header = F, sep = "\t", comment.char = "!",
                    quote = "", col.names = gpa.colnames)
  # optionally, filter annotations with the qualifier NOT
  if (filter.NOT) 
    goa <- goa[!grepl("NOT", goa$Qualifier),]
  # optionally, filter out evidence 
  evidence <- gsub("go_evidence=", "", goa$Annotation.Properties)
  goa <- goa[!grepl(paste(filter.evidence, collapse = "|"), evidence),]
  # map to accession 
  if (accession != "UNIPROT") {
    map <- suppressMessages(AnnotationDbi::select(
      database, keys = goa$DB_Object_ID, keytype = "UNIPROT", 
      columns = accession))
    goa[[accession]] <- map[[accession]][match(goa$DB_Object_ID, map$UNIPROT)]
    # remove NAs
    goa <- goa[!is.na(goa[[accession]]),]
  } else {
    goa$UNIPROT <- goa$DB_Object_ID
  }
  # read the ontology
  if (!is.null(ontology.file) & propagate) {
    ontology <- ontologyIndex::get_ontology(ontology.file, 
                                            extract_tags = "minimal")
    goa$ancestors <- ontology$ancestors[goa$GO.ID]
    goa <- unnest(goa, ancestors)
    # replace column
    goa[["GO.ID"]] <- goa[["ancestors"]]
    goa <- goa[, -ncol(goa)]
  }
  return(goa)
}
