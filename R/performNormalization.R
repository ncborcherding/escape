#' Perform Normalization on Enrichment Data
#' 
#' This function allows users to normalize the enrichment
#' calculations by accounting for single-cell dropout and
#' producing positive values for downstream differential 
#' enrichment analyses. 
#' 
#' @param input.data Enrichment output from \code{\link{escape.matrix}} or
#' \code{\link{runEscape}}.
#' @param assay Name of the assay to plot if data is a single-cell object.
#' @param gene.set.reference The gene set library to use to extract 
#' the individual gene set information from.
#' @make.positve Shift enrichment values to a positive range \strong{TRUE}
#' @param reduction.key Name of the key to use with the components.

performNormalization <- function(input.data,
                                 assay = NULL,
                                 gene.set.reference = NULL,
                                 make.positive = TRUE) {
  
  
  if(is_seurat_or_se_object(input.data)) {
    enriched <- .pull.Enrich(input.data, assay)
  } else {
    enriched <- input.data
  }
  
  #Getting the gene sets that passed filters
  egc <- .GS.check(gene.set.reference)
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
    enriched[,x]/gene.set
    if(any(enriched[,x] < 0) & make.positive) {
      enriched[,x] <- enriched[,x] + abs(min(enriched[,x]))
    }
    enriched[,x]
  }) -> normalized.values
  
  normalized.enriched <- do.call(cbind, normalized.values)
  colnames(normalized.enriched) <- colnames(enriched)
  
  #TODO add ammended enrichment values to single-cell object or return matrix
  
  if(is_seurat_or_se_object(input.data)) {
    input.data <- .adding.Enrich(input.data, normalized.enriched, assay)
    return(input.data)
  }
  return(normalized.enriched)
  
}
