#' Generate a ridge plot to examine enrichment distributions
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
#' @param palette Colors to use in visualization - input any 
#' \link[grDevices]{hcl.pals}.
#'
#' @import ggplot2
#' @importFrom ggdist stat_pointinterval
#' 
#' @examples
#' GS <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
#'            Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))
#' pbmc_small <- SeuratObject::pbmc_small
#' pbmc_small <- runEscape(pbmc_small, 
#'                         gene.sets = GS, 
#'                         min.size = NULL)
#'                         
#' geyserEnrichment(pbmc_small, 
#'                  assay = "escape",
#'                  gene.set = "Tcells")
#'
#' @export
#'
#' @return ggplot2 object with geyser-based distributions of selected gene.set
geyserEnrichment <- function(enriched, 
                             assay = NULL,
                             group.by =NULL, 
                             gene.set = NULL, 
                             color.by = "group",
                             order.by = NULL,
                             scale = FALSE, 
                             facet.by = NULL, 
                             palette = "inferno") {
  
  if(is.null(group.by)) {
    group.by <- "ident"
  }
  
  if(color.by == "group") {
    color.by <- group.by
  }
  
  enriched <- .prepData(input.data, assay, gene.set, group.by, NULL, facet.by) 
  
  if(!is.null(order.by) && !is.null(group.by)) {
    enriched <- .orderFunction(enriched, order.by, group.by)
  } 
  
  if(scale) {
    enriched[,gene.set] <- as.numeric(scale(enriched[,gene.set]))
  }
  
  if(inherits(enriched[,color.by], "numeric") && gene.set == color.by) {
    gradient.format <- TRUE
  } else {
    gradient.format <- FALSE
  }
  
  plot <- ggplot(data = enriched,
                 mapping = aes(x = enriched[,group.by],
                               y = enriched[,gene.set],
                               color = enriched[,color.by])) 
  
  plot <- plot + 
          geom_jitter(size = 2,
                      na.rm = TRUE) + 
          stat_pointinterval(interval_size_range = c(2, 3),
                             fatten_point = 1.5,
                             interval_color = "white",
                             point_color = "white",
                             position = position_dodge(width = 1),
                             na.rm = TRUE,
                             show.legend = FALSE) +
          stat_pointinterval(interval_size_range = c(1, 2),
                             interval_color = "black",
                             point_color = "black",
                             position = position_dodge(width = 1),
                             na.rm = TRUE,
                             show.legend = FALSE)
  
  plot <- plot + 
          xlab(group.by) +
          ylab(paste0(gene.set, " Enrichment Score")) +
          theme_classic() +
          guides(fill = "none")
  
  plot <- .colorby(enriched,
                   plot, 
                   color.by,
                   palette,
                   type = "color")
  
  if (!is.null(facet.by)) {
    plot <- plot + 
            facet_grid(as.formula(paste('. ~', facet.by))) 
  }
  return(plot)
}