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
#'  
#'  seurat_ex <- suppressWarnings(SeuratObject::pbmc_small)
#'  ES <- enrichIt(obj = seurat_ex, gene.sets = GS)
#'  
#'  enrichmentPlot(seurat_ex, gene.set = "Cancer",
#'                 gene.sets = GS, group = "Type")
#'  }
#' @import patchwork
#' @import GSVA
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
#'
#' @return ggplot2 object mean rank gene density
enrichmentPlot <- function(obj, 
                           gene.set, 
                           gene.sets, 
                           group, 
                           colors = c("#0D0887FF","#7E03A8FF","#CC4678FF","#F89441FF","#F0F921FF")) {
  cnts <- cntEval(obj)
  cnts.filter <- as.matrix(.filterFeatures.mod(cnts, "ssgsea"))
  meta <- grabMeta(obj)
  group <- meta[, group]
  uniq.grp <- sort(as.character(unique(group)))
  colpal<- colorRampPalette(colors)(length(uniq.grp))
  gene.sets <- GS.check(gene.sets)
  gene.set <- unlist(gene.sets[[gene.set]])
  mapped.gset.idx.list <- na.omit(match(gene.set, rownames(cnts.filter)))
  gene.density <- compute.gene.density.mod(cnts.filter, seq_len(ncol(cnts)), TRUE)
  rank.scores <- rep(0, nrow(cnts.filter))
  sort.sgn.idxs <- apply(gene.density, 2, order, decreasing=TRUE)
  gsva_rnk2 <- apply(sort.sgn.idxs, 2, compute_rank_score.mod, nrow(cnts))
  output <- NULL
  for (i in seq_along(uniq.grp)) {
    means <- rowMeans(gsva_rnk2[,group == uniq.grp[i]])
    rank <- round(order(means)/2)
    output <- cbind(output, rank)
  }
  output <- data.frame(output)
  colnames(output) <- paste0("Group.", uniq.grp)
  output$GS <- NA
  output$GS[mapped.gset.idx.list] <- "yes"
  mlt <- suppressMessages(melt(output))
  plot1 <- ggplot(mlt, aes(x = value, color = variable)) + 
    geom_density(data = subset(mlt, GS == "yes"), linetype="dashed") + 
    theme_classic() + 
    scale_color_manual(values = colpal) + 
    labs(color = "Group") + 
    ylab("Rank Density") + 
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(), 
          axis.text.x = element_blank())
  mlt$segmenty <- NA
  mlt$segmenty2 <- NA
  ymax <- 0.2
  for (i in seq_along(uniq.grp)) {
    mlt$segmenty <- ifelse(mlt$variable == paste0("Group.", uniq.grp[i]), -(i*ymax-ymax), mlt$segmenty)
    mlt$segmenty2 <- ifelse(mlt$variable == paste0("Group.", uniq.grp[i]), -(i*ymax), mlt$segmenty2)   
  }
  plot2 <- ggplot(subset(mlt, GS == "yes")) + 
    geom_segment(aes(x = value,y=segmenty,yend=segmenty2,xend=value, color = variable)) + 
    guides(color = "none") +
    xlab("Mean Rank Order") + 
    scale_color_manual(values = colpal) + 
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

# Internal function from GSVA
#' @importFrom stats ecdf
compute.gene.density.mod <- function(expr, sample.idxs, rnaseq=FALSE){
  n.test.samples <- ncol(expr)
  n.genes <- nrow(expr)
  n.density.samples <- length(sample.idxs)
  
  gene.density <- NA
    gene.density <- t(apply(expr, 1, function(x, sample.idxs) {
      f <- ecdf(x[sample.idxs])
      f(x)
    }, sample.idxs))
    gene.density <- log(gene.density / (1-gene.density)) 
  
  return(gene.density)	
}

# Internal function from GSVA
#' @importFrom MatrixGenerics rowSds
.filterFeatures.mod <- function(expr, method) {
  ## filter out genes with constant expression values
  ## DelayedMatrixStats::rowSds() works for both base and 
  ## DelayedArray matrices
  sdGenes <- rowSds(expr)
  ## the following fixes this bug, see issues
  ## https://github.com/rcastelo/GSVA/issues/54
  ## https://github.com/HenrikBengtsson/matrixStats/issues/204
  sdGenes[sdGenes < 1e-10] <- 0
  if (any(sdGenes == 0) || any(is.na(sdGenes))) {
    warning(sum(sdGenes == 0 | is.na(sdGenes)),
            " genes with constant expression values throuhgout the samples.")
    if (method != "ssgsea") {
      warning("Since argument method!=\"ssgsea\", genes with constant expression values are discarded.")
      expr <- expr[sdGenes > 0 & !is.na(sdGenes), ]
    }
  } 
  
  if (nrow(expr) < 2)
    stop("Less than two genes in the input assay object\n")
  
  if(is.null(rownames(expr)))
    stop("The input assay object doesn't have rownames\n")
  expr
}