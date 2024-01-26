#' Generate a ridge plot to examine enrichment distributions
#' 
#' This function allows to the user to examine the distribution of 
#' enrichment across groups by generating a ridge plot.
#'
#' @param enriched The output of \code{\link{enrichIt}}
#' @param group The parameter to group, displayed on the y-axis.
#' @param gene.set The gene set to graph on the x-axis. 
#' @param scale.bracket This will filter the enrichment scores to remove 
#' extreme outliers. Values entered (1 or 2 numbers) will be the filtering 
#' parameter using z-scores of the selected gene.set. If only 1 value is given, 
#' a seocndary bracket is autommatically selected as the inverse of the number.
#' @param palette Colors to use in visualization - input any 
#' \link[grDevices]{hcl.pals}.
#' @param facet A parameter to separate the graph.
#' @param add.rug Binary classifier to add a rug plot to the x-axis.
#'
#' @import ggplot2
#' @importFrom ggridges geom_density_ridges geom_density_ridges2 position_points_jitter
#' 
#' @examples
#' ES2 <- readRDS(url(
#' "https://ncborcherding.github.io/vignettes/escape_enrichment_results.rds"))
#' ridgeEnrichment(ES2, gene.set = "HALLMARK_DNA_REPAIR", group = "cluster", 
#' facet = "Type", add.rug = TRUE)
#'
#' @export
#'
#' @seealso \code{\link{enrichIt}} for generating enrichment scores.
#' @return ggplot2 object with ridge-based distributions of selected gene.set
ridgeEnrichment <- function(enriched, 
                            group.by = "cluster", 
                            gene.set = NULL, 
                            #color.by = "group.by"
                            scale = NULL, 
                            facet.by = NULL, 
                            add.rug = FALSE,
                            palette = "inferno") {

  enriched <- .prepData(input.data, assay, group.by, split.by, facet.by) 
  
  if(scale) {
    enriched[,gene.set] <- scale(enriched[,gene.set])
  }
  
  if(!is.null(order.by) && !is.null(group.by)) {
    enriched <- .orderFunction(enriched, order.by, group.by)
  }
  
  col <- length(unique(enriched[,split.by]))
  plot <- ggplot(enriched, aes(x = enriched[,gene.set], 
                               y = enriched[,group.by], 
                               fill = enriched[,group]))
  
  if (add.rug == TRUE) {
    plot <- plot + geom_density_ridges(
      jittered_points = TRUE,
      position = position_points_jitter(width = 0.05, height = 0),
      point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7) 
    
  } else {
    plot <- plot + 
      geom_density_ridges2(alpha = 0.8) 
  }
  
  plot <- plot + ylab(group) +
    xlab(paste0(gene.set, " Enrichment Score")) +
    labs(fill = ) + #############
    scale_fill_manual(values = .colorizer(palette, col))
    theme_classic() +
    guides(fill = "none")
  
  if (!is.null(facet.by)) {
    plot <- plot + facet_grid(as.formula(paste('. ~', facet.by))) 
    }
  
  return(plot)
}