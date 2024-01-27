#' Generate a ridge plot to examine enrichment distributions
#' 
#' This function allows to the user to examine the distribution of 
#' enrichment across groups by generating a ridge plot.
#'
#' @param enriched The output of \code{\link{enrichIt}}
#' @param group The parameter to group, displayed on the y-axis.
#' @param gene.set The gene set to graph on the x-axis. 
#' @param scale This will filter the enrichment scores to remove 
#' extreme outliers. Values entered (1 or 2 numbers) will be the filtering 
#' parameter using z-scores of the selected gene.set. If only 1 value is given, 
#' a seocndary bracket is autommatically selected as the inverse of the number.
#' @param palette Colors to use in visualization - input any 
#' \link[grDevices]{hcl.pals}.
#' @param facet A parameter to separate the graph.
#' @param add.rug Binary classifier to add a rug plot to the x-axis.
#'
#' @import ggplot2
#' @importFrom ggdist stat_pointinterval
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
#' 
#' 
geyserEnrichment <- function(enriched, 
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
  
  plot <- plot + geom_jitter(size = 2) + 
    stat_pointinterval(interval_size_range = c(2, 3),
                               fatten_point = 1.5,
                               interval_color = "white",
                               point_color = "white",
                               position = ggplot2::position_dodge(width = 1),
                               na.rm = TRUE,
                               show.legend = FALSE) +
    stat_pointinterval(interval_size_range = c(1, 2),
                               interval_color = "black",
                               point_color = "black",
                               position = ggplot2::position_dodge(width = 1),
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
    plot <- plot + facet_grid(as.formula(paste('. ~', facet.by))) 
  }
  return(plot)
}