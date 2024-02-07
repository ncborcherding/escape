#' Perform Normalization on Enrichment Data
#' 
#' This function allows users to normalize the enrichment
#' calculations by accounting for single-cell dropout and
#' producing positive values for downstream differential 
#' enrichment analyses. 

performNormalization <- function(input.data,
                                 assay = NULL,
                                 gene.set.reference = NULL,
                                 make.positive = TRUE) {
  
  
  egc <- .GS.check(gene.set.reference)
  names(egc) <- str_replace_all(names(egc), "_", "-")
  enriched <- .pull.Enrich(input.data, assay)
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
  
  #TODO add ammended enrichment values to single-cell object or return matrix
  
  
}
