#' Get a collection of gene sets to perform enrichment on
#'
#' This function allows users to select libraries and specific 
#' gene.sets to form a GeneSetCollection that is a list of gene sets.
#
#' @param species The scientific name of the species of interest in 
#' order to get correcent gene nomenclature
#' @param library Individual collection(s) of gene sets, e.g. c("H", "C5").
#' See \href{https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp}{msigdbr}for
#' all MSigDB collections.
#' @param subcategory MSigDB sub-collection abbreviation, such as CGP or BP.
#' @param gene.sets Select gene sets or pathways, using specific names, 
#' example: pathways = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB"). Will only be
#' honored if library is set, too.
#'
#' @examples 
#' GS <- getGeneSets(library = "H")
#' 
#' @export
#' 
#' @importFrom GSEABase GeneSet GeneSetCollection
#' @importFrom msigdbr msigdbr msigdbr_species
#' 
#' @author Nick Borcherding, Jared Andrews
#' @return A \code{GeneSetCollection} object containing the requested \code{GeneSet} objects.
getGeneSets <- function(species = "Homo sapiens", 
                        library = NULL, 
                        subcategory = NULL,
                        gene.sets = NULL) {
  spec <- msigdbr_species()
  spec_check <- unlist(spec[spec$species_name %in% species,][,1])
  if (length(spec_check) == 0) {
    message(paste0("Please select a compatible species: ", 
                   paste(spec, collapse = ", ")))
  }
  if(!is.null(library)) {
    if (length(library) == 1) {
      if (is.null(subcategory)) {
        m_df = msigdbr(species = spec_check, category = library)
      } else {
        m_df = msigdbr(species = spec_check, category = library, subcategory = subcategory)
      }
    }
    m_df <- NULL
    for (x in seq_along(library)) {
      if (is.null(subcategory)) {
        tmp2 = msigdbr(species = spec_check, category = library[x])
      } else {
        tmp2 = msigdbr(species = spec_check, category = library, subcategory = subcategory)
      }
      m_df <- rbind(m_df, tmp2)
    }
    if(!is.null(gene.sets)) {
      m_df <- m_df[m_df$gs_name %in% gene.sets,]
    }    
  }
  gs <- unique(m_df$gs_name)
  ls <- list()
  for (i in seq_along(gs)) {
    tmp <- m_df[m_df$gs_name == gs[i],]
    tmp <- tmp$gene_symbol
    tmp <- unique(tmp)
    tmp <- GeneSet(tmp, setName=paste(gs[i]))
    ls[[i]] <- tmp
  }
  gsc <- GeneSetCollection(ls)
  return(gsc)
}