#' Visualize enrichment results with a ridge plot
#' 
#' This function allows to the user to examine the distribution of 
#' enrichment across groups by generating a ridge plot.
#'
#' @param input.data Enrichment output from \code{\link{escape.matrix}} or
#' \code{\link{runEscape}}.
#' @param assay Name of the assay to plot if data is a single-cell object.
#' @param group.by Categorical parameter to plot along the x.axis. If input is
#' a single-cell object the default will be cluster.
#' @param gene.set Gene set to plot (on y-axis).
#' @param color.by How the color palette applies to the graph - can 
#' be \strong{"group"} for a categorical color palette based on the 
#' \strong{group.by} parameter or use the \strong{gene.set} name if wanting to 
#' apply a gradient palette.
#' @param order.by Method to organize the x-axis: \strong{"mean"} will arrange
#' the x-axis by the mean of the gene.set, while \strong{"group"} will arrange
#' the x-axis by in alphanumerical order. Using \strong{NULL} will not reorder
#' the x-axis.
#' @param facet.by Variable to facet the plot into n distinct graphs.
#' @param scale Visualize raw values \strong{FALSE} or Z-transform 
#' enrichment values \strong{TRUE}.
#' @param add.rug Add visualization of the discrete cells along
#' the ridge plot (\strong{TRUE}).
#' @param palette Colors to use in visualization - input any 
#' \link[grDevices]{hcl.pals}.
#'
#' @import ggplot2
#' @importFrom ggridges geom_density_ridges geom_density_ridges2 position_points_jitter geom_density_ridges_gradient
#' 
#' @examples 
#' GS <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
#'            Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))
#' pbmc_small <- SeuratObject::pbmc_small
#' pbmc_small <- runEscape(pbmc_small, 
#'                         gene.sets = GS, 
#'                         min.size = NULL)
#'                         
#' ridgeEnrichment(pbmc_small, 
#'                 assay = "escape",
#'                 gene.set = "Tcells")
#'                 
#' ridgeEnrichment(pbmc_small, 
#'                 assay = "escape",
#'                 gene.set = "Tcells", 
#'                 color.by = "Tcells")
#'
#' @export
#'
#' @return ggplot2 object with ridge-based distributions of selected gene.set
ridgeEnrichment <- function(input.data, 
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
                                                  point_shape = '|', 
                                                  point_size = 3, 
                                                  point_alpha = 1, 
                                                  alpha = 0.7,
                                                  quantile_lines = TRUE, 
                                                  quantile_fun = median,
                                                  vline_width = 1)
    } else {
      plot <- plot + geom_density_ridges(jittered_points = TRUE,
                                         position = position_points_jitter(width = 0.05, height = 0),
                                         point_shape = '|', 
                                         point_size = 3, 
                                         point_alpha = 1, 
                                         alpha = 0.7,
                                         quantile_lines = TRUE, 
                                         quantile_fun = median,
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
                   color.by,
                   palette)
  
  if (!is.null(facet.by)) {
    plot <- plot + 
            facet_grid(as.formula(paste('. ~', facet.by))) 
  }
  
  return(plot)
}
