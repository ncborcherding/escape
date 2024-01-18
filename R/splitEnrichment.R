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

#' Generate a split violin plot examine enrichment distributions
#' 
#' This function allows to the user to examine the distribution of 
#' enrichment across groups by generating a split violin plot.
#'
#' @param enriched The output of \code{\link{enrichIt}}
#' @param x.axis Optional parameter for seperation.
#' @param gene.set The gene set to graph on the y-axis. 
#' @param scale.bracket This will filter the enrichment scores to remove 
#' extreme outliers. Values entered (1 or 2 numbers) will be the filtering 
#' parameter using z-scores of the selected gene.set. If only 1 value is given, 
#' a secondary bracket is automatically selected as the inverse of the number.
#' @param split The parameter to split, must be binary.
#' @param palette Colors to use in visualization - input any 
#' \link[grDevices]{hcl.pals}.
#'
#' @import ggplot2
#' 
#' @examples
#' ES2 <- readRDS(url(
#' "https://ncborcherding.github.io/vignettes/escape_enrichment_results.rds"))
#' splitEnrichment(ES2, x.axis = "cluster", split = "Type", 
#' gene.set = "HALLMARK_DNA_REPAIR")
#'
#' @export
#'
#' @seealso \code{\link{enrichIt}} for generating enrichment scores.
#' @return ggplot2 object violin-based distributions of selected gene.set
splitEnrichment <- function(enriched, 
                            x.axis = NULL, 
                            scale.bracket = NULL,
                            split = NULL, 
                            gene.set = NULL, 
                            palette = "inferno") {
  
  if (length(unique(enriched[,split])) != 2) {
    message("SplitEnrichment() can only work for binary classification")}
  
  if (!is.null(scale.bracket)) {
    if (length(scale.bracket) != 1 | length(scale.bracket) != 1) {
      message("Please indicate one or two values for the scale.bracket 
                parameter, such as scale.bracket = c(-2,2)")
    }
    scale.bracket <- order(scale.bracket)
    if(length(scale.bracket) == 1) {
      scale.bracket <- c(scale.bracket, -scale.bracket)
      scale.bracket <- order(scale.bracket)
    } 
    tmp <- enriched
    tmp[,gene.set]<- scale(tmp[,gene.set])
    rows_selected <- rownames(tmp[tmp[,gene.set] >= scale.bracket[1] & 
                                    tmp[,gene.set] <= scale.bracket[2],])
    enriched <- enriched[rownames(enriched) %in% rows_selected,]
  }
  cols <- length(unique(enriched[,split]))
  if (is.null(x.axis)) {
    plot <- ggplot(enriched, aes(x = ".", 
                                 y = enriched[,gene.set], 
                                 fill = enriched[,split])) 
    check = 1
  } else {
    plot <- ggplot(enriched, aes(x = enriched[,x.axis], 
                                 y = enriched[,gene.set], 
                                 fill = enriched[,split])) + 
      xlab(x.axis) 
    check = NULL}
  plot <- plot + 
    geom_split_violin(alpha=0.8) +
    geom_boxplot(width=0.1, fill = "grey", alpha=0.5, 
                 outlier.alpha = 0)  + 
    ylab(paste0(gene.set, " (NES)")) +
    labs(fill = split) + 
    scale_fill_manual(values = .colorizer(palette, col))+
    theme_classic() +
    #guides(fill = FALSE)
    if (!is.null(check)) {
      plot <- plot + theme(axis.title.x = element_blank(),
                           axis.text.x = element_blank(),
                           axis.ticks.x = element_blank())}
  return(plot)
}