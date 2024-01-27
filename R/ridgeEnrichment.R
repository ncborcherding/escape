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
                            assay = NULL,
                            group.by =NULL, 
                            gene.set = NULL, 
                            color.by = "group",
                            order.by = NULL,
                            scale = FALSE, 
                            facet.by = NULL, 
                            add.rug = FALSE,
                            palette = "inferno") {

  if(is.null(group.by)) {
    group.by <- "ident"
  }
  
  if(color.by == "group") {
    color.by <- group.by
  }
  
  enriched <- .prepData(input.data, assay, gene.set, group.by, NULL, facet.by) 
  
  if(inherits(enriched[,color.by], "numeric") && gene.set == color.by) {
    gradient.format <- TRUE
  } else {
    gradient.format <- FALSE
  }
  
  if(scale) {
    enriched[,gene.set] <- as.numeric(scale(enriched[,gene.set]))
  }
  
  if(!is.null(order.by) && !is.null(group.by)) {
    enriched <- .orderFunction(enriched, order.by, group.by)
  } 
  

  
  if(gradient.format) {
    plot <- ggplot(enriched, aes(x = enriched[,gene.set], 
                               y = enriched[,group.by], 
                               fill = after_stat(x)))
  } else {
    plot <- ggplot(enriched, aes(x = enriched[,gene.set], 
                                 y = enriched[,group.by], 
                                 fill = enriched[,group.by]))
  }
  
  if (add.rug) {
    if(gradient.format) {
      plot <- plot + geom_density_ridges_gradient(jittered_points = TRUE,
                                   position = position_points_jitter(width = 0.05, height = 0),
                                   point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7,
                                   quantile_lines = TRUE, quantile_fun = median,
                                   vline_width = 1)
    } else {
      plot <- plot + geom_density_ridges(
        jittered_points = TRUE,
        position = position_points_jitter(width = 0.05, height = 0),
        point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7,
        quantile_lines = TRUE, quantile_fun = median,
        vline_width = 1) 
    }
    
  } else {
    if(gradient.format) {
      plot <- plot + 
        geom_density_ridges_gradient(alpha = 0.8, 
                                     quantile_lines = TRUE, 
                                     quantile_fun = median,
                                     vline_width = 1) 
    } else {
      plot <- plot + 
        geom_density_ridges2(alpha = 0.8,
                             quantile_lines = TRUE, 
                             quantile_fun = median,
                             vline_width = 1) 
    }
  }
  
  plot <- plot + 
            ylab(group.by) +
            xlab(paste0(gene.set, " Enrichment Score")) +
            labs(fill = color.by) + #############
            theme_classic() +
            guides(fill = "none")
  
  plot <- .colorby(enriched,
                   plot, 
                   color.by)
  
  if (!is.null(facet.by)) {
    plot <- plot + facet_grid(as.formula(paste('. ~', facet.by))) 
  }
  
  return(plot)
}
