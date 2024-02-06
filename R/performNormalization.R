

performNormalization <- function(input.data,
                                 assay = NULL
                                 gene.set.reference = NULL) {
  
  
  egc <- .GS.check(gene.set.reference)
  enriched <- .pull.Enrich(input.data, assay)
  egc <- egc[egc %in% colnames(enriched)]
  
  cnts <- .cntEval(input.data, assay = "RNA", type = "counts")
  egc.size <- lapply(egc, function(x) length(which(rownames(cnts) %in% x)))
  
  #TODO normalize by gene set size and transform to positive integers/float
  
  
}