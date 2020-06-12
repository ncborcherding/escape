#' Calculate Principal Components for the Enrichment Scores
#'
#' Using all or selected enrichment scores of individual 
#' single-cells, this function will calculate 
#' principal components using scaled values and attach 
#' to the output columns to use to graph later.
#'
#' @param enriched The output of \code{\link{enrichIt}}.
#' @param groups The column headers to use in future graphing functions.
#'
#' @importFrom dplyr select_if
#' @importFrom stats prcomp
#' 
#' @examples 
#' ES2 <- readRDS(url(
#' "https://ncborcherding.github.io/vignettes/escape_enrichment_results.rds"))
#' PCA <- performPCA(enriched = ES2, groups = c("Type", "Cluster"))
#'
#' @export
#' @return Data frame of principal compoenents
#'
#' @author Nick Borcherding
#'
performPCA <- function(enriched, groups) {
    groups <- enriched[,colnames(enriched) %in% groups]
    input <- select_if(enriched, is.numeric)
    PCA <- prcomp(input, scale. = TRUE)
    merged <- merge(PCA$x, groups, by = "row.names")
    rownames(merged) <- merged[,1]
    merged <- merged[,-1]

    return(merged)
}

#' Get a collection of gene sets to perform enrichment on
#'
#' This function allows users to select libraries and specific 
#' gene.sets to form a GeneSetCollection that is a list of gene sets.
#
#' @param species The scientific name of the species of interest in 
#' order to get correcent gene nomenclature
#' @param library Individual libraries or multiple libraries to select, 
#' example: library = c("H", "C5").
#' @param gene.sets Select gene sets or pathways, using specific names, 
#' example: pathways = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB").
#'
#' @examples 
#' GS <- getGeneSets(library = "H")
#' 
#' @export
#' 
#' @importFrom GSEABase GeneSet GeneSetCollection
#' @importFrom msigdbr msigdbr msigdbr_show_species
#' 
#' @author Nick Borcherding, Jared Andrews
#' @return List of GeneSets in collection format
getGeneSets <- function(species = "Homo sapiens", 
                    library = NULL, gene.sets = NULL) {
    spec <- msigdbr_show_species()
    spec_check <- spec[spec %in% species]
    if (length(spec_check) == 0) {
        stop(paste0("Please select a compatible species: ", 
                    paste(spec, collapse = ", ")))
    }
    m_df = msigdbr(species = spec_check)
    if(!is.null(library)) {
        m_df <- m_df[m_df$gs_cat %in% library,]
    }
    if(!is.null(gene.sets)) {
        m_df <- m_df[m_df$gs_name %in% gene.sets,]
    }
    gs <- unique(m_df$gs_name)
    ls <- list()
    for (i in seq_along(gs)) {
        tmp <- m_df[m_df$gs_name == gs[i],]
        tmp <- tmp$human_gene_symbol
        tmp <- unique(tmp)
        tmp <- GeneSet(tmp, setName=paste(gs[i]))
        ls[[i]] <- tmp
    }
    gsc <- GeneSetCollection(ls)
    return(gsc)
}



