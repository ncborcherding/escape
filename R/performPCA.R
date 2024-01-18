#' Calculate Principal Components for the Enrichment Scores
#'
#' Using all or selected enrichment scores of individual 
#' single-cells, this function will calculate 
#' principal components using scaled values and attach 
#' to the output columns to use to graph later.
#'
#' @param enriched The output of \code{\link{enrichIt}}.
#' @param gene.sets Names of gene sets to include in the PCA
#' @param groups The column headers to use in future graphing functions.
#'
#' @importFrom dplyr select_if
#' @importFrom stats prcomp
#' 
#' @examples 
#' ES2 <- readRDS(url(
#' "https://ncborcherding.github.io/vignettes/escape_enrichment_results.rds"))
#' 
#' PCA <- performPCA(enriched = ES2, groups = c("Type", "Cluster"), 
#' gene.sets = colnames(ES2))
#'
#' @export
#' @return Data frame of principal components
#'
#' @author Nick Borcherding
#'
performPCA <- function(enriched, gene.sets = NULL, groups) {
  groups <- data.frame(enriched[,colnames(enriched) %in% c(groups)])
  input <- select_if(enriched, is.numeric)
  if (!is.null(gene.sets)) {
    input <- input[,colnames(input) %in% gene.sets]
  }
  PCA <- prcomp(input, scale. = TRUE)
  merged <- cbind(PCA$x, groups)
  return(merged)
}