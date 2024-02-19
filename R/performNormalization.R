#' Perform Normalization on Enrichment Data
#' 
#' This function allows users to normalize the enrichment calculations 
#' by accounting for single-cell dropout and producing positive 
#' values for downstream differential enrichment analyses. A positive range
#' values is useful for several downstream analyses, like differential 
#' evaluation for log2-fold change, but will alter the original 
#' enrichment values.
#' 
#' @param input.data Enrichment output from \code{\link{escape.matrix}} or
#' \code{\link{runEscape}}.
#' @param assay Name of the assay to plot if data is a single-cell object.
#' @param gene.sets The gene set library to use to extract 
#' the individual gene set information from.
#' @param make.positive Shift enrichment values to a positive range \strong{TRUE}
#' for downstream analysis or not \strong{TRUE} (default).
#' 
#' @examples
#' GS <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
#'            Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))
#' pbmc_small <- SeuratObject::pbmc_small
#' pbmc_small <- runEscape(pbmc_small, 
#'                         gene.sets = GS, 
#'                         min.size = NULL)
#'                         
#' pbmc_small <- performNormalization(pbmc_small, 
#'                                    assay = "escape", 
#'                                    gene.sets = GS)
#'
#' @export
#' 
#' @return Single-cell object or matrix of normalized enrichment scores

performNormalization <- function(input.data,
                                 assay = NULL,
                                 gene.sets = NULL,
                                 make.positive = FALSE) {
  
  
  if(is_seurat_or_se_object(input.data)) {
    enriched <- .pull.Enrich(input.data, assay)
  } else {
    enriched <- input.data
  }
  
  #Getting the gene sets that passed filters
  egc <- .GS.check(gene.sets)
  names(egc) <- str_replace_all(names(egc), "_", "-")
  egc <- egc[names(egc) %in% colnames(enriched)]
  
  #Isolating the number of genes per cell expressed
  cnts <- .cntEval(input.data, assay = "RNA", type = "counts")
  egc.size <- lapply(egc, function(x)  {
    lapply(seq_len(ncol(cnts)), function(y) {
      values <-length(which(names(cnts[which(cnts[,y] != 0),y]) %in% x))
      values
    })
  })
  
  #Dividing the enrichment score by number of genes expressed
  lapply(seq_len(ncol(enriched)), function(x) {
    gene.set <- unlist(egc.size[colnames(enriched)[x]])
    if(any(gene.set == 0)) {
      gene.set[which(gene.set == 0)] <- 1
    }
    enriched[,x] <- enriched[,x]/gene.set
    if(any(enriched[,x] < 0) & make.positive) {
      enriched[,x] <- enriched[,x] + abs(min(enriched[,x]))
    }
    enriched[,x]
  }) -> normalized.values
  
  normalized.enriched <- do.call(cbind, normalized.values)
  colnames(normalized.enriched) <- colnames(enriched)
  
  if(is_seurat_or_se_object(input.data)) {
    input.data <- .adding.Enrich(input.data, normalized.enriched, paste0(assay, "_normalized"))
    return(input.data)
  } else {
    return(normalized.enriched)
  }
  
}
