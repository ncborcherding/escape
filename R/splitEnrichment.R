#Developing split violin plot
#Code from: https://stackoverflow.com/a/45614547
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), 
                                               xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = 
                                                                  if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], 
                                              newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- 
                               round(newdata[1, "x"])
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), 
                                         all(draw_quantiles <= 1))
                               quantiles <- 
                                 ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), 
                                                  setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", 
                                                grid::grobTree(GeomPolygon$draw_panel(newdata, ...), 
                                                               quantile_grob))
                             } else {
                               ggplot2:::ggname("geom_split_violin", 
                                                GeomPolygon$draw_panel(newdata, ...))}
                           })

#Defining new geometry
#Code from: https://stackoverflow.com/a/45614547
geom_split_violin <- 
  function(mapping = NULL, data = NULL, 
           stat = "ydensity", position = "identity", ..., draw_quantiles = NULL, 
           trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, 
           inherit.aes = TRUE) {
    layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
          position = position, show.legend = show.legend, 
          inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, 
                                                   draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
  }

#' Visualize enrichment results with a split violin plot
#' 
#' This function allows to the user to examine the distribution of 
#' enrichment across groups by generating a split violin plot.
#'
#' @param input.data Enrichment output from \code{\link{escape.matrix}} or
#' \code{\link{runEscape}}.
#' @param assay Name of the assay to plot if data is a single-cell object.
#' @param split.by Variable to form the split violin, must have 2 levels.
#' @param group.by Categorical parameter to plot along the x.axis. If input is
#' a single-cell object the default will be cluster.
#' @param gene.set Gene set to plot (on y-axis).
#' @param order.by Method to organize the x-axis - \strong{"mean"} will arrange
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
#' 
#' @examples
#' GS <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
#'            Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))
#' pbmc_small <- SeuratObject::pbmc_small
#' pbmc_small <- runEscape(pbmc_small, 
#'                         gene.sets = GS, 
#'                         min.size = NULL)
#'                         
#' splitEnrichment(pbmc_small, 
#'                 assay = "escape",
#'                 split.by = "groups",
#'                 gene.set = "Tcells")
#'
#' @export
#'
#' @return ggplot2 object violin-based distributions of selected gene.set
splitEnrichment <- function(input.data,
                            assay = NULL,
                            split.by = NULL,
                            group.by = NULL, 
                            gene.set = NULL,
                            order.by = NULL,
                            facet.by = NULL,
                            scale = TRUE,
                            palette = "inferno") {
  if(is.null(split.by)){
    stop("Please select a variable with 'split.by' to generate the splitEnrichment() plots")
  } 
  
  if(is.null(group.by)) {
    group.by <- "ident"
  }
  
  enriched <- .prepData(input.data, assay, gene.set, group.by, split.by, facet.by) 
  
  if (length(unique(enriched[,split.by])) != 2) {
    message("SplitEnrichment() can only work for binary variables - reselect 'split.by'")
  }
  
  if(scale) {
    enriched[,gene.set] <- scale(enriched[,gene.set])
  }
  
  if(!is.null(order.by) && !is.null(group.by)) {
    enriched <- .orderFunction(enriched, order.by, group.by)
  }
  
  col <- length(unique(enriched[,split.by]))
  if (is.null(group.by)) {
    plot <- ggplot(enriched, aes(x = ".", 
                                 y = enriched[,gene.set], 
                                 fill = enriched[,split.by])) 
    check = 1
  } else {
    plot <- ggplot(enriched, aes(x = enriched[,group.by], 
                                 y = enriched[,gene.set], 
                                 fill = enriched[,split.by])) + 
                  xlab(group.by) 
    check = NULL
    }
  
  plot <- plot + 
          geom_split_violin(alpha=0.8, lwd= 0.25) +
          geom_boxplot(width=0.1, 
                       fill = "grey", 
                       alpha=0.5, 
                       outlier.alpha = 0)  + 
          ylab(paste0(gene.set, " Enrichment Score")) +
          labs(fill = split.by) + 
          scale_fill_manual(values = .colorizer(palette, col))+
          theme_classic() 

    if (!is.null(check)) {
      plot <- plot + 
              theme(axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                    axis.ticks.x = element_blank())}
  if (!is.null(facet.by)) {
    plot <- plot + 
            facet_grid(as.formula(paste('. ~', facet.by))) 
  }
  return(plot)
}
