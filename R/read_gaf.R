#' Read a GAF file 
#' 
#' Read a GO annotation file in GAF format. Further information about the GAF
#' format is available from the 
#' \href{http://geneontology.org/page/go-annotation-file-gaf-format-21}{
#' Gene Ontology Consortium}. 
#' 
#' @param filepath the location of the GAF file. Files can be downloaded from
#' ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/
#' @param database optionally, provide a Bioconductor `AnnotationDbi` database
#' to map identifiers from UniProt to some other accession 
#' @param accession optionally, the type of identifier to map UniProt accessions
#' to using a Bioconductor `AnnotationDbi` database. Must be a keytype in the
#' database. 
#' @param filter_NOT if true, filter annotations with the qualifier NOT
#' @param filter_evidence optionally, specify evidence codes to filter. By 
#' default, evidence codes ND, IPI, IEA and NAS are filtered. 
#' @param propagate if true, and an ontology file is provided, all ancestors of 
#' a given term are associated with each protein.
#' @return a data.frame containing the filtered GAF file  
#' @export
#' @examples 
#' # read human GOA and map to Ensembl, filtering IEA annotations
#' go <- read_gaf("goa_human.gpa.gz", database = org.Hs.eg.db,
#'    accession = "ENSEMBL", filter.evidence = "IEA")
#' # read mouse GOA and map to gene symbol, filtering IEA and IPI annotations
#' go <- read_gaf("goa_mouse.gpa.gz", database = org.Mm.eg.db,
#'    accession = "SYMBOL", filter.evidence = c("IEA", "IPI"))
read_gaf <- function(filepath,
                     database = NULL, accession = "UNIPROT", 
                     filter.NOT = T, 
                     filter.evidence = c("ND", "IPI", "IEA", "NAS"),
                     ontology = NULL, propagate = T) {
  message("reading GOA file ", filepath, "...")
  # read GOA file 
  gaf.colnames <- c("DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", 
                    "GO_ID", "DB_Reference", "Evidence_Code", "With_From",
                    "Aspect", "DB_Object_Name", "DB_Object_Synonym",
                    "DB_Object_Type", "Taxon", "Date", "Assigned_By", 
                    "Annotation_Extension", "Gene_Product_Form_ID")
  goa <- suppressMessages(
    readr::read_tsv(filepath, comment = "!", col_names = gaf.colnames))
  # optionally, filter annotations with the qualifier NOT
  if (filter.NOT) 
    goa <- goa[!grepl("NOT", goa$Qualifier),]
  # optionally, filter out evidence 
  if (!is.null(filter.evidence) & length(filter.evidence) > 0) 
    goa <- goa[!goa$Evidence_Code %in% filter.evidence,]
  # map to accession 
  if (!is.null(accession) & accession != "UNIPROT") {
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
  if (!is.null(ontology) & propagate) {
    if (!"ontology_index" %in% class(ontology))
      stop("Ontology must be of class ontology_index")
    goa$ancestors <- ontology$ancestors[goa$GO_ID]
    # filter out terms missing ancestors (deprecated)
    goa <- goa[lengths(goa$ancestors) > 0,]
    goa <- tidyr::unnest(goa, ancestors)
    # replace column
    goa[["GO_ID"]] <- goa[["ancestors"]]
    goa <- goa[, -ncol(goa)]
  }
  return(goa)
}
