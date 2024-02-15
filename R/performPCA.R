#' Perform Principal Component Analysis on Enrichment Data
#' 
#' This function allows users to calculate the principal components 
#' for the gene set enrichment values. For single-cell data, the PCA
#' will be stored with the dimensional reductions. If a matrix is used
#' as input, the output is a list for further plotting. Alternatively,
#' users can use functions for PCA calculations based on their desired
#' workflow in lieu of using \code{\link{performPCA}}.
#'
#' @param input.data Enrichment output from \code{\link{escape.matrix}} or
#' \code{\link{runEscape}}.
#' @param assay Name of the assay to plot if data is a single-cell object.
#' @param scale Standardize the enrichment value (\strong{TRUE}) or 
#' not (\strong{FALSE})
#' @param n.dim The number of components to calculate.
#' @param reduction.name Name of the reduced dimensions object to add if 
#' data is a single-cell object.
#' @param reduction.key Name of the key to use with the components.
#'
#' @importFrom stats prcomp
#' @importFrom SeuratObject CreateDimReducObject
#' @importFrom SingleCellExperiment reducedDim reducedDim<-
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
#' @export
#' 
#' @return single-cell object or list with PCA components to plot.
performPCA <- function(input.data,
                       assay = NULL,
                       scale = TRUE,
                       n.dim = 1:10,
                       reduction.name = "escape.PCA",
                       reduction.key = "PCA") {
  
  if(is_seurat_or_se_object(input.data)) {
    enriched <- .pull.Enrich(input.data, assay)
  } else {
    enriched <- input.data
  }
  
  PCA <- prcomp(enriched, 
                scale. = scale,
                rank. = max(n.dim))
  rotation <- PCA$rotation
  eigen.values <- PCA$sdev^2
  percent.contribution <- round((eigen.values/sum(eigen.values))*100,1)
  PCA <- PCA$x
  colnames(PCA) <- paste0(reduction.key, "_", seq_len(ncol(PCA)))
  
  additional.data <- list(eigen_values = eigen.values,
                          contribution = percent.contribution, 
                          rotation = rotation)
  if(is_seurat_or_se_object(input.data)) {
    if (inherits(input.data, "Seurat")) {
      DR <- suppressWarnings(CreateDimReducObject(
                            embeddings = PCA,
                            stdev = rep(0, ncol(PCA)),
                            key = reduction.key,
                            jackstraw = NULL,
                            misc = additional.data))
      input.data[[reduction.name]] <- DR
    } else if (inherits(input.data, "SingleCellExperiment")) {
      reducedDim(input.data, reduction.name) <- PCA
      if(length(input.data@metadata) == 0) {
        input.data@metadata <- additional.data
      } else {
        input.data@metadata <- c(input.data@metadata, additional.data)
      }
    
    } 
    return(input.data)
  } else {
    PCA.results <- list(PCA = PCA,
                        eigen_values = eigen.values,
                        contribution = percent.contribution,
                        rotation = rotation)
    return(PCA.results)
  }
  
}
