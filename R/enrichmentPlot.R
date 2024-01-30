#' Gene Rank Enrichment Plot
#' Display the rank order density for individual gene sets 
#' by group identify. This function will use the group variable 
#' to take the mean rank order across all individual cells.
# 
#' @param obj The Seurat or SingleCellExperiment object.
#' @param gene.set The name of the specific gene set to visualize
#' @param gene.sets Gene sets from \code{\link{getGeneSets}} to use 
# for the enrichment analysis. 
#' @param group The header in the meta data that will be used for the comparison
#' @param colors The color palette for the enrichment plot
#' @examples 
#'  \dontrun{
#'  GS <- list(Housekeeping = c("ACTA1", "ACTN1", "GAPDH"),
#'  Cancer = c("TP53","BRCA2","ERBB2","MYC"))
#'  pbmc_small <- suppressWarnings(SeuratObject::pbmc_small)
#'  
#'  enrichmentPlot(pbmc_small gene.set = "Cancer",
#'                 gene.sets = GS, group = "groups")
#'  }
#' @import patchwork
#' @importFrom utils getFromNamespace
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
#'
#' @return ggplot2 object mean rank gene density
densityEnrichment <- function(input.data, 
                              gene.set.use = NULL, 
                              gene.set.reference = NULL, 
                              group.by = NULL, 
                              palette = "inferno") {
  if (!inherits(x=input.data, what ="Seurat") || 
      !inherits(x=input.data, what ="SummarizedExperiment")) {
    stop("Currently this function only support single-cell objects")
  }
  
  if(is.null(group.by)) {
    group.by <- "ident"
  }
  
  compute.gene.density<-utils::getFromNamespace("compute.gene.density", "GSVA")
  
  
  
  cnts <- .cntEval(input.data, 
                   assay = "RNA", 
                   type = "counts")
  cnts.filter <- .filterFeatures(cnts)
  grouping <- as.vector(.grabMeta(input.data)[,group.by])
  
  groups <- unique(grouping)
  lapply(seq_len(length(groups)), function(x) {
    tmp <- cnts.filter[,which(grouping == groups[x])]
    density <- compute.gene.density(tmp, seq_len(ncol(tmp)), TRUE, TRUE)
    rank.scores <- rep(0, nrow(tmp))
    sort.sgn.idxs <- apply(density, 2, order, decreasing=TRUE)
    gsva_rnk2 <- apply(sort.sgn.idxs, 2, compute_rank_score.mod, nrow(cnts))
    means <- rowMeans(gsva_rnk2)
    rank <- round(order(means)/2)
    rank
  }) -> ranks
  
  output <- do.call(cbind, ranks)
  output <- as.data.frame(output)
  colnames(output) <- paste0(group.by, ".", groups)
  rownames(output) <- rownames(cnts.filter)
  
  mapped.gset.idx.list <- na.omit(match(gene.set, rownames(cnts.filter)))
  
  output$gene.set.query <- NA
  output$gene.set.query[mapped.gset.idx.list] <- "yes"
  melted.data.frame <- suppressMessages(melt(output))
  col <- .colorizer(palette, length(groups))
  
  
  plot1 <- ggplot(melted.data.frame, aes(x = value)) + 
    geom_density(data = subset(melted.data.frame, gene.set.query == "yes"), 
                 #linetype="dashed",
                 aes(fill = variable), 
                 alpha = 0.4,
                 color = "black") + 
    theme_classic() + 
    scale_fill_manual(values = col) + 
    labs(fill = "Group") + 
    ylab("Rank Density") + 
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(), 
          axis.text.x = element_blank())
  melted.data.frame$segmenty <- NA
  melted.data.frame$segmenty2 <- NA
  ymax <- 0.2
  for (i in seq_along(groups)) {
    melted.data.frame$segmenty <- ifelse(melted.data.frame$variable == paste0(group.by, ".", groups[i]), -(i*ymax-ymax), melted.data.frame$segmenty)
    melted.data.frame$segmenty2 <- ifelse(melted.data.frame$variable == paste0(group.by, ".", groups[i]), -(i*ymax), melted.data.frame$segmenty2)   
  }
  plot2 <- ggplot(subset(melted.data.frame, gene.set.query == "yes")) + 
    geom_segment(aes(x = value,y=segmenty,yend=segmenty2,xend=value, color = variable),
                 lwd = 1) + 
    guides(color = "none") +
    xlab("Mean Rank Order") + 
    scale_color_manual(values = col) + 
    theme(axis.title.y = element_blank(),
          axis.ticks.y = element_blank(), 
          axis.text.y = element_blank(), 
          panel.background = element_rect(fill = NA, colour = "black"))
  EnPlot <- plot1 + plot2 + plot_layout(ncol=1, heights = c(3, 1))
  return(EnPlot)
}
  
# Internal function from GSVA
compute_rank_score.mod <- function(sort_idx_vec, p){
  tmp <- rep(0, p)
  tmp[sort_idx_vec] <- abs(seq(from=p,to=1) - p/2)
  return (tmp)
}

# Modified from GSVA
#' @importFrom MatrixGenerics rowSds
.filterFeatures <- function(expr) {
  sdGenes <- rowSds(expr)
  sdGenes[sdGenes < 1e-10] <- 0
  if (any(sdGenes == 0) || any(is.na(sdGenes))) {
      expr <- expr[sdGenes > 0 & !is.na(sdGenes), ]
  }
  
  if (nrow(expr) < 2)
    stop("Less than two genes in the input assay object\n")
  
  if(is.null(rownames(expr)))
    stop("The input assay object doesn't have rownames\n")
  expr
}
