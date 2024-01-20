#' Calculate gene set enrichment scores 
#'
#' This function allows users to input both the single-cell RNA-sequencing 
#' counts and output the enrichment scores as a matrix.  
#'
#'
#' @param input.data The count matrix, Seurat, or Single-Cell Experiment object.
#' @param gene.sets Gene sets can be a list, output from 
#' \code{\link{getGeneSets}}, or the built-in gene sets 
#' in the escape package \code{\link{escape.gene.sets}}.
#' @param method select the method to calculate enrichment, 
#' \strong{GSVA}, \strong{ssGSEA} or \strong{UCell}
#' @param groups The number of cells to separate the enrichment calculation.
#' @param min.size Minimum number of gene necessary to perform the enrichment
#' calculation
#' @param normalize Whether to divide the enrichment score by the number 
#' of genes \strong{TRUE} or report unnormalized \strong{FALSE}
#' @param ... pass arguments to GSVA, ssGSEA or UCell call
#'
#' @importFrom GSVA gsva gsvaParam ssgseaParam
#' @importFrom GSEABase GeneSetCollection 
#' @importFrom UCell ScoreSignatures_UCell
#' @importFrom AUCell AUCell_run
#' @importFrom SummarizedExperiments assay
#'
#' @examples 
#' GS <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
#'            Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))
#' pbmc_small <- SeuratObject::pbmc_small
#' ES <- escape.matrix(pbmc_small, 
#'                     gene.sets = GS, 
#'                     min.size = NULL)
#' 
#' @export
#' @author Nick Borcherding, Jared Andrews
#'
#' @seealso \code{\link{getGeneSets}} to collect gene sets.
#' @return Data frame of normalized enrichmenet scores (NES)
escape.matrix <- function(input.data, 
                          gene.sets = NULL, 
                          method = "ssGSEA", 
                          groups = 1000, 
                          min.size = 5,
                          normalize = FALSE,
                          ...) {
  
    egc <- .GS.check(gene.sets)
    cnts <- .cntEval(input.data)
    egc.size <- lapply(egc, function(x) length(which(rownames(cnts) %in% x)))
    if (!is.null(min.size)){
      remove <- unname(which(egc.size < min.size))
      egc <- egc[-remove]
      egc.size <- egc.size[-remove]
      
    }
    
    scores <- list()
    splits <- seq(1, ncol(cnts), by=groups)
    print(paste('Using sets of', groups, 'cells. Running', 
                length(splits), 'times.'))
    split.data <- .split_data.matrix(matrix=cnts, chunk.size=groups)
    
    
    for (i in seq_along(splits)) {
          last <- min(ncol(cnts), i+groups-1)
          if(method == "GSVA") {
              parameters <- .gsva.setup(split.data[[i]], egc)
          } else if (method == "ssGSEA") {
              parameters <- .ssGSEA.setup(split.data[[i]], egc)
          }
          if(method %in% c("ssGSEA", "GSVA")) {
              a <- gsva(param = parameters, 
                        verbose = FALSE,
                        ...)
          } else if(method == "UCell") {
              a <- t(suppressWarnings(
                ScoreSignatures_UCell(matrix = split.data[[i]], 
                                      features=egc,
                                      name = NULL,
                                      ...)))
          } else if (method == "AUCell") {
             a <- t(assay(suppressWarnings(
                           AUCell_run(exprMat = split.data[[i]], 
                                      geneSets = egc,
                                      normAUC = FALSE,
                                      ...))))
             
          }
          scores[[i]] <- a
    }
    scores <- do.call(cbind, scores)
    output <- t(as.matrix(scores))
    #Normalizing ssGSEA by number of genes in gene sets present
    if(normalize) {
        for(i in seq_len(ncol(output))) {
          output[,i] <- output[,i]/egc.size[[i]]
        }
    }
    return(output)
}

#' Ibex single cell calculation
#'
#' Run Ibex algorithm with Seurat or SingleCellExperiment pipelines

#' Enrichment calculation for single-cell workflows
#'
#' Run the escape-based gene-set enrichment calculation with 
#' Seurat or SingleCellExperiment pipelines
#' 
#' @examples 
#' GS <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
#'            Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))
#' pbmc_small <- SeuratObject::pbmc_small
#' pbmc_small <- runEscape(pbmc_small, 
#'                         gene.sets = GS, 
#'                         min.size = NULL)
#'
#'
#' @param input.data The count matrix, Seurat, or Single-Cell Experiment object.
#' @param gene.sets Gene sets can be a list, output from 
#' \code{\link{getGeneSets}}, or the built-in gene sets 
#' in the escape package \code{\link{escape.gene.sets}}.
#' @param method select the method to calculate enrichment, 
#' \strong{GSVA}, \strong{ssGSEA} or \strong{UCell}
#' @param groups The number of cells to separate the enrichment calculation.
#' @param min.size Minimum number of gene necessary to perform the enrichment
#' calculation
#' @param normalize Whether to divide the enrichment score by the number 
#' of genes \strong{TRUE} or report unnormalized \strong{FALSE}
#' @param new.assay.nam The new name of the assay to append to 
#' the single-cell object containing the enrichment scores.
#' @param ... pass arguments to GSVA, ssGSEA or UCell call
#' @export
#' @return Seurat or SingleCellExperiment object with escape enrichment scores
#' in the assay slot. 

runEscape <- function(input.data, 
                      gene.sets = NULL, 
                      method = "ssGSEA", 
                      groups = 1000, 
                      min.size = 5,
                      normalize = FALSE,
                      new.assay.name = "escape",
                      ...) {
  .checkSingleObject(input.data)
  enrichment <- escape.matrix(input.data = input.data,
                              gene.sets = gene.sets,
                              method = method,
                              groups = groups,
                              min.size = min.size)
  
  input.data <- .adding.Enrich(input.data, enrichment, new.assay.name)
  return(input.data)
}

.gsva.setup <- function(data, egc) {
  params.to.use <- gsvaParam(exprData = data,
                             geneSets = egc,
                             kcdf = "Poisson")
  return(params.to.use)
}

.ssGSEA.setup <- function(data, egc) {
  params.to.use <- ssgseaParam(exprData = data,
                               geneSets = egc,
                               normalize = FALSE)
  return(params.to.use)
}
