#' Visualize the PCA of enrichment values
#' 
#' This function allows to the user to examine the distribution
#' of principal components run on the enrichment values.
#'
#' @param input.data PCA from \code{\link{performPCA}}.
#' @param dimRed Name of the dimensional reduction to plot if data is a single-cell object.
#' @param x.axis Component to plot on the x.axis.
#' @param y.axis Component set to plot on the y.axis.
#' @param facet.by Variable to facet the plot into n distinct graphs.
#' @param style Return a \strong{"hex"} bin plot or a \strong{"point"}-based plot.
#' @param add.percent.contribution Add the relative percent of contribution of the 
#' selected components to the axis labels.
#' @param display.factors Add an arrow overlay to show the direction and magnitude of individual 
#' gene sets on the PCA dimensions.
#' @param number.of.factors The number of gene.sets to display on the overlay.
#' @param palette Colors to use in visualization - input any 
#' \link[grDevices]{hcl.pals}.
#'
#' @import ggplot2
#' @importFrom dplyr slice_max %>%
#' 
#' @examples 
#' GS <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
#'            Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))
#' pbmc_small <- SeuratObject::pbmc_small
#' pbmc_small <- runEscape(pbmc_small, 
#'                         gene.sets = GS, 
#'                         min.size = NULL)
#'                         
#' pbmc_small <- performPCA(pbmc_small, 
#'                          assay = "escape")
#'                          
#' pcaEnrichment(pbmc_small,
#'               x.axis = "PC1",
#'               y.axis = "PC2",
#'               dimRed = "escape.PCA")
#'
#' @export
#'
#' @return ggplot2 object with PCA distribution
pcaEnrichment <- function(input.data, 
                          dimRed = NULL,
                          x.axis = "PC1", 
                          y.axis = "PC2",
                          facet.by = NULL,
                          style = "point",
                          add.percent.contribution = TRUE,
                          display.factors = FALSE,
                          number.of.factors = 10,
                          palette = "inferno") {
  
  
  if (is_seurat_or_se_object(input.data)) {
    pca.values <- .grabDimRed(input.data, dimRed) 
  } else if (inherits(input.data, "list") & length(input.data) == 4) {
    pca.values <- input.data
    if(!is.null(facet.by)) {
      stop("group.by parameter requires input.data to be a single-cell object.")
    }
  } else {
    stop("input.data does not seem to be a single-cell object or a product of performPCA().")
  }
  
  x.axis.dim <- as.numeric(substring(x.axis, 3, nchar(x.axis)))
  y.axis.dim <- as.numeric(substring(y.axis, 3, nchar(y.axis)))
  
  if(add.percent.contribution & length(pca.values) == 4) {
    x.axis.title <- paste0(x.axis, "\n (", pca.values[[3]][x.axis.dim],"%)")
    y.axis.title <- paste0(y.axis, "\n (", pca.values[[3]][y.axis.dim],"%)")
  } else {
    x.axis.title <- x.axis
    y.axis.title <- y.axis
  }
  
  plot.df <- as.data.frame(pca.values[[1]])
  
  if(!is.null(facet.by)) {
    meta <- .grabMeta(input.data)
    if(facet.by %!in% colnames(meta)) {
      stop("Please select a variable in your meta data to use for facet.by.")
    }
    col.pos <- ncol(plot.df)
    plot.df <- cbind.data.frame(plot.df, meta[,facet.by])
    colnames(plot.df)[col.pos+1] <- facet.by
  }
  
  plot <- ggplot(data = plot.df,
                 mapping = aes(x = plot.df[,x.axis.dim], 
                               y = plot.df[,y.axis.dim])) 
  
  if(style == "point") {
    plot <- plot + 
      geom_pointdensity() + 
      scale_color_gradientn(colors = .colorizer(palette, 11)) + 
      labs(color = "Relative Density")
  } else if (style == "hex") {
    plot <- plot + 
      stat_binhex() +
      scale_fill_gradientn(colors = .colorizer(palette, 11))
      labs(fill = "Relative Density")
  }
  
  plot <- plot + 
    ylab(y.axis.title) +
    xlab(x.axis.title) +
    theme_classic()
  
  if (!is.null(facet.by)) {
    plot <- plot + 
      facet_grid(as.formula(paste('. ~', facet.by))) 
  }
  
  if(display.factors) {
    x.range <- range(plot.df[,x.axis.dim])
    
    y.range <- range(plot.df[,y.axis.dim])
    
    tbl <- data.frame(names = row.names(pca.values[[4]]), 
                      factors.y = pca.values[[4]][,y.axis.dim]^2/sum(pca.values[[4]][,y.axis.dim]^2),
                      factors.x = pca.values[[4]][,x.axis.dim]^2/sum(pca.values[[4]][,x.axis.dim]^2)) %>%
      slice_max(n = number.of.factors, 
                order_by = (factors.x + factors.y)/2)
    names <- tbl$names
    
    df <- as.data.frame(pca.values[[4]])
    df <- df[rownames(df) %in% names,]
    df$names <- rownames(df)
    if(!is.null(facet.by)) {
      facets <- sort(unique(plot.df[,facet.by]))
      df[,facet.by] <- facets[1]
    }
    
    plot <- plot + 
            geom_hline(yintercept = 0, lty=2) + 
            geom_vline(xintercept = 0, lty=2) +
            geom_segment(data = df, 
                         aes(x = 0,
                             y = 0,
                             xend = .scale.variable(df[,x.axis.dim], x.range),
                             yend = .scale.variable(df[,y.axis.dim], y.range)),
                         arrow = arrow(length = unit(0.25, "cm"))) + 
            geom_label(data = df, aes(label = names,
                                      x = .scale.variable(df[,x.axis.dim], x.range),
                                      y = .scale.variable(df[,y.axis.dim], y.range)), 
                    size=2, 
                    hjust = 0.5, 
                    nudge_y = -0.01,
                    label.padding = unit(0.1, "lines")) 
  }
  
  
  return(plot)
}

# Function to scale the new variable
.scale.variable <- function(new_var, existing_range) {
  new_range <- range(new_var)
  existing_range <- existing_range* 0.8
  normalized <- (new_var - min(new_range)) / (max(new_range) - min(new_range))
  scaled <- normalized * (max(existing_range) - min(existing_range)) + min(existing_range)
  return(scaled)
}
