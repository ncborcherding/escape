#' Visualize the components of the PCA analysis of the enrichment results
#'
#' Graph the major gene set contributors to the \code{\link{pcaEnrichment}}.
#'
#' @param enriched The output of \code{\link{enrichIt}}.
#' @param gene.sets Names of gene sets to include in the PCA
#' @param PCx The principal component graphed on the x-axis.
#' @param PCy The principal component graphed on the y-axis.
#' @param top.contribution The number of gene sets to graph, organized 
#' by PCA contribution.
#'
#' @importFrom ggplot2 ggplot
#' @importFrom stats prcomp
#' @import dplyr
#' 
#' @examples 
#' ES2 <- readRDS(url(
#' "https://ncborcherding.github.io/vignettes/escape_enrichment_results.rds"))
#' 
#' masterPCAPlot(ES2, PCx = "PC1", PCy = "PC2", gene.sets = colnames(ES2), 
#' top.contribution = 10)
#'
#' @export
#'
#' @seealso \code{\link{enrichIt}} for generating enrichment scores.
#' @return ggplot2 object sumamrizing the PCA for the enrichment scores
masterPCAPlot <- function(enriched, 
                          gene.sets, 
                          PCx, 
                          PCy, top.contribution = 10) {
  input <- select_if(enriched, is.numeric)
  if (!is.null(gene.sets)) {
    input <- input[,colnames(input) %in% gene.sets]
  }
  PCA <- prcomp(input, scale. = TRUE)
  var_explained <- PCA$sdev^2/sum(PCA$sdev^2)
  
  tbl <- data.frame(names = rownames(PCA$rotation), 
                    factors.y = PCA$rotation[,PCy]^2/sum(PCA$rotation[,PCy]^2),
                    factors.x = PCA$rotation[,PCx]^2/sum(PCA$rotation[,PCx]^2)) 
  names <- tbl %>% 
    slice_max(n = top.contribution, 
              order_by = (factors.x + factors.y)/2)
  names <- names$names
  df <- as.data.frame(PCA$rotation)
  df <- df[rownames(df) %in% names,]
  df$names <- rownames(df)
  
  plot <- df %>%
    ggplot(aes(x=df[,PCx],y=df[,PCy])) + 
    geom_point() + 
    geom_text(aes_string(label = "names"), 
              size=2, hjust = 0.5, nudge_y = -0.01) + 
    geom_hline(yintercept = 0, lty=2) + 
    geom_vline(xintercept = 0, lty=2) +
    labs(x=paste0(PCx,": ",round(var_explained[1]*100,1),"%"),
         y=paste0(PCy, ": ",round(var_explained[2]*100,1),"%")) +
    theme_classic()
  return(plot)
}