# Assigning color pallette
#' @importFrom grDevices colorRampPalette
assignColor <- function(x, enriched, group) {
    if (length(x) != length(unique(enriched[,group]))) {
        x <- colorRampPalette(x)(length(unique(enriched[,group])))
    } else { x <- x }
    return(x)
}


#' Density plot of the principal components
#'
#' @param PCAout The output of \code{\link{performPCA}}
#' @param PCx The principal component graphed on the x-axis
#' @param PCy The principal component graphed on the y-axis
#' @param colors The color palette for the density plot
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
pcaEnrichment <- function(PCAout, PCx, PCy, 
    colors = c("#0348A6", "#7AC5FF", "#C6FDEC", "#FFB433", "#FF4B20"), 
    contours = TRUE, facet = NULL) 
{
    plot <- ggplot(PCAout, aes(x=PCAout[,PCx], y=PCAout[,PCy])) +
        stat_binhex() +
        scale_fill_gradientn(colours = colorRampPalette(colors)(50)) +
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

#' Visualize the components of the PCA analysis of the enrichment results
#'
#' Graph the major gene set contributors to the \code{\link{pcaEnrichment}}.
#'
#' @param enriched The output of \code{\link{enrichIt}}.
#' @param PCx The principal component graphed on the x-axis.
#' @param PCy The principal component graphed on the y-axis.
#' @param top.contribution The number of gene sets to graph, organized 
#' by PCA contribution.
#'
#' @importFrom factoextra fviz_pca_var
#' @importFrom ggplot2 ggplot
#' @importFrom ggrepel geom_text_repel
#' @importFrom stats prcomp
#' @importFrom rlang .data
#' @import dplyr
#' 
#' @examples 
#' ES2 <- readRDS(url(
#' "https://ncborcherding.github.io/vignettes/escape_enrichment_results.rds"))
#' masterPCAPlot(ES2, PCx = "PC1", PCy = "PC2", top.contribution = 10)
#'
#' @export
#'
#' @seealso \code{\link{enrichIt}} for generating enrichment scores.
#' @return ggplot2 object sumamrizing the PCA for the enrichment scores
masterPCAPlot <- function(enriched, PCx, PCy, top.contribution = 10) {
    input <- select_if(enriched, is.numeric)
    PCA <- prcomp(input, scale. = TRUE)
    PCx1 <- which(colnames(PCA$x) == PCx)
    PCy1 <- which(colnames(PCA$x) == PCy)
    p <- fviz_pca_var(PCA, axes = c(PCx1, PCy1))
    output <- data.frame(p$data$name, p$data$contrib, p$data$x, p$data$y)
    output <- output %>% top_n(n = top.contribution, wt = .data$p.data.contrib)

    ggplot(output, aes_string(x = "p.data.x", y = "p.data.y")) +
        geom_point() +
        geom_text_repel(data = output, aes_string(label = "p.data.name"), 
            size=2) +
        theme_classic() +
        xlab(PCx) +
        ylab(PCy) + 
        geom_hline(yintercept = 0, lty=2) + 
        geom_vline(xintercept = 0, lty=2)
}

#' Generate a ridge plot to examine enrichment distributions
#' 
#' This function allows to the user to examine the distribution of 
#' enrichment across groups by generating a ridge plot.
#'
#' @param enriched The output of \code{\link{enrichIt}}
#' @param group The parameter to group, displayed on the y-axis.
#' @param gene.set The gene set to graph on the x-axis. 
#' @param scale.bracket The convert enrichment to z-scores filter 
#' ithin the selected range. Must include lower and upper limit.
#' @param colors The color palette for the ridge plot.
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
ridgeEnrichment <- function(enriched, group = "cluster", gene.set = NULL, 
            scale.bracket = NULL, facet = NULL, add.rug = FALSE,
            colors = c("#0348A6", "#7AC5FF", "#C6FDEC", "#FFB433", "#FF4B20")) 
             {
    if (!is.null & length(scale.bracket) == 2) {
        enriched[,gene.set] <- scale(enriched[,gene.set])
        enirched <- enriched[enriched[,gene.set] >= brackets[1] & enriched[,gene.set] <= brackets[2], ]
    }
    colors <- assignColor(colors, enriched, group) 
    plot <- ggplot(enriched, aes(x = enriched[,gene.set], 
                    y = enriched[,group], fill = enriched[,group]))
    
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
        xlab(paste0(gene.set, " (NES)")) +
        labs(fill = group) + 
        scale_fill_manual(values = colors) + 
        theme_classic() +
        guides(fill = FALSE)
    
    if (!is.null(facet)) {
        plot <- plot + facet_grid(as.formula(paste('. ~', facet))) }
    
    return(plot)
}

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
                    grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
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
#' @param scale.bracket The convert enrichment to z-scores filter 
#' ithin the selected range. Must include lower and upper limit.
#' @param split The parameter to split, must be binary.
#' @param colors The color palette for the ridge plot.
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
splitEnrichment <- function(enriched, x.axis = NULL, 
                            split = NULL, gene.set = NULL, 
                            colors = c("#0348A6", "#7AC5FF", "#C6FDEC", 
                                "#FFB433", "#FF4B20")) {
    
    if (length(unique(enriched[,split])) != 2) {
        stop("SplitEnrichment() can only work for binary classification")}
    
    if (!is.null & length(scale.bracket) == 2) {
        enriched[,gene.set] <- scale(enriched[,gene.set])
        enirched <- enriched[enriched[,gene.set] >= brackets[1] & enriched[,gene.set] <= brackets[2], ]
    }
    colors <- assignColor(colors, enriched, split) 
    if (is.null(x.axis)) {
        plot <- ggplot(enriched, aes(x = ".", y = enriched[,gene.set], 
                    fill = enriched[,split])) 
        check = 1
    } else {
        plot <- ggplot(enriched, aes(x = enriched[,x.axis], 
                    y = enriched[,gene.set], 
                    fill = enriched[,split])) + 
            xlab(x.axis) 
        check = NULL
    }
    plot <- plot + 
        geom_split_violin(alpha=0.8) +
        geom_boxplot(width=0.1, fill = "grey", alpha=0.5, 
            outlier.alpha = 0)  + 
        ylab(paste0(gene.set, " (NES)")) +
        labs(fill = split) + 
        scale_fill_manual(values = colors) + 
        theme_classic() +
        guides(fill = FALSE)
    if (!is.null(check)) {
        plot <- plot + theme(axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                    axis.ticks.x = element_blank())
    }
    return(plot)
}


