#' Density plot of the principal components
#'
#' @param PCAout The output of \code{\link{performPCA}}
#' @param PCx The principal component graphed on the x-axis
#' @param PCy The principal component graphed on the y-axis
#' @param palette Colors to use in visualization - input any 
#' \link[grDevices]{hcl.pals}.
#' @param contours Binary classifier to add contours to the density plot
#' @param facet A parameter to separate the graph
#'
#' @import ggplot2
#' @importFrom grDevices colorRampPalette
#' 
#' @examples 
#' ES2 <- readRDS(url(
#' "https://ncborcherding.github.io/vignettes/escape_enrichment_results.rds"))
#' PCA <- performPCA(enriched = ES2, groups = c("Type", "Cluster"))
#' pcaEnrichment(PCA, PCx = "PC1", PCy = "PC2", contours = TRUE)
#'
#' @export
#'
#' @seealso \code{\link{performPCA}} for generating PCA results.
#' @return ggplot2 object of the results of PCA for the enrichment scores
pcaEnrichment <- function(PCAout, 
                          PCx, 
                          PCy, 
                          palette = "inferno", 
                          contours = TRUE, facet = NULL) 
{
  plot <- ggplot(PCAout, aes(x=PCAout[,PCx], y=PCAout[,PCy])) +
    stat_binhex() +
    scale_fill_gradientn(name = "Percentage", colors = .colorizer(palette,21)) +
    theme_classic() +
    ylab(PCy) +
    xlab(PCx)
  
  if (contours == TRUE) {
    plot <- plot + stat_density_2d(color = "black")
  }
  else {
    plot <- plot
  }
  
  if (!is.null(facet)) {
    plot <- plot + facet_wrap(as.formula(paste('~', facet)))
  }
  plot <- plot +
    geom_hline(yintercept = 0, lty=2) + 
    geom_vline(xintercept = 0, lty=2)
  
  return(plot)
}