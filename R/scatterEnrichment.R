#' Generate a density-based scatter plot
#' 
#' This function allows to the user to examine the distribution of 
#' 2 gene sets along the x.axis and y.axis. The color gradient
#' is generated using the a density estimate. See 
#' \href{https://github.com/LKremer/ggpointdensity}{ggpointdensity})
#' for more information.
#'
#' @param input.data Enrichment output from \code{\link{escape.matrix}} or
#' \code{\link{runEscape}}.
#' @param assay Name of the assay to plot if data is a single-cell object.
#' @param x.axis Gene set to plot on the x.axis.
#' @param y.axis Gene set to plot on the y.axis.
#' \strong{group.by} parameter or use the \strong{gene.set} name if wanting to 
#' apply a gradient palette.
#' @param order.by Method to organize the x-axis: \strong{"mean"} will arrange
#' the x-axis by the mean of the gene.set, while \strong{"group"} will arrange
#' the x-axis by in alphanumerical order. Using \strong{NULL} will not reorder
#' the x-axis.
#' @param facet.by Variable to facet the plot into n distinct graphs.
#' @param scale Visualize raw values \strong{FALSE} or Z-transform 
#' enrichment values \strong{TRUE}.
#' @param palette Colors to use in visualization - input any 
#' \link[grDevices]{hcl.pals}.
#'
#' @import ggplot2
#' @importFrom ggpointdensity geom_pointdensity
#' 
#' @examples
#' GS <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
#'            Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))
#' pbmc_small <- SeuratObject::pbmc_small
#' pbmc_small <- runEscape(pbmc_small, 
#'                         gene.sets = GS, 
#'                         min.size = NULL)
#'                         
#' scatterEnrichment(pbmc_small, 
#'                   assay = "escape",
#'                   x.axis = "Tcells",
#'                   y.axis = "Bcells")
#'
#' @export
#' 
#' @return ggplot2 object with a scatter plot of selected gene.sets
scatterEnrichment <- function(enriched, 
                              assay = NULL,
                              x.axis = NULL, 
                              y.axis = NULL,
                              scale = FALSE, 
                              facet.by = NULL, 
                              palette = "inferno") {
  
  gene.set <- c(x.axis, y.axis)
  
  enriched <- .prepData(input.data, assay, gene.set, group.by, NULL, facet.by) 
  
  if(scale) {
    enriched[,gene.set] <- apply(enriched[,gene.set], 2, scale)
  }
  
  plot <- ggplot(data = enriched, aes(x = enriched[,x.axis], 
                              y = enriched[,y.axis])) +
                geom_pointdensity() +
                scale_color_gradientn(colors = .colorizer(palette, 11)) + 
                ylab(paste0(y.axis, " Enrichment Score")) +
                xlab(paste0(x.axis, " Enrichment Score")) +
                labs(color = "Relative Density") + 
                theme_classic()
  
  if (!is.null(facet.by)) {
    plot <- plot + 
      facet_grid(as.formula(paste('. ~', facet.by))) 
  }
  return(plot)
}
