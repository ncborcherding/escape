#' Calculate gene set enrichment scores for single-cell data
#'
#' This function allows users to input both the single-cell RNA-sequencing 
#' counts and any gene set pathways either from the stored data or from 
#' other sources. The enrichment calculation itself 
#' uses the two methods 1) gsva R package and the poisson distribution for RNA
#' or the \href{https://github.com/carmonalab/UCell}{UCell package}. 
#' If using the method "UCell", the core option will be ignored due to openMP
#' support issues on Macs (will be changed in the future). 
#'
#' @param obj The count matrix, Seurat, or SingleCellExperiment object.
#' @param gene.sets Gene sets from \code{\link{getGeneSets}} to use 
#' for the enrichment analysis. Alternatively a simple base R list where
#' the names of the list elements correspond to the name of the gene set
#' and the elements themselves are simple vectors of gene names representing
#' the gene set. 
#' @param method select the method to calculate enrichment, either ssGSEA or UCell
#' @param groups The number of cells to separate the enrichment calculation.
#' @param cores The number of cores to use for parallelization.
#'
#' @importFrom GSVA gsva
#' @importFrom GSEABase GeneSetCollection geneIds
#' @importFrom SingleCellExperiment counts
#' @importFrom UCell ScoreSignatures_UCell
#' @importFrom BiocParallel SnowParam
#' @importFrom Matrix summary
#'
#' 
#' @examples 
#' # download HALLMARK gene set collection
#' GS <- list(Housekeeping = c("ACTA1", "ACTN1", "GAPDH"),
#'   Cancer = c("TP53","BRCA2","ERBB2","MYC"))
#' 
#' seurat_ex <- suppressWarnings(SeuratObject::pbmc_small)
#' ES <- enrichIt(obj = seurat_ex, gene.sets = GS)
#' 
#' @export
#'
#' @author Nick Borcherding, Jared Andrews
#'
#' @seealso \code{\link{getGeneSets}} to collect gene sets.
#' @return Data frame of normalized enrichmenet scores (NES)
enrichIt <- function(obj, gene.sets = NULL, 
                     method = "ssGSEA", groups = 1000, cores = 2) {
    
    if(is.null(gene.sets)) {
        stop("Please provide the gene.sets you would like to use for 
            the enrichment analysis")
    } else {
        egc <- gene.sets
    }
    
    if (inherits(x = obj, what = "Seurat")) {
        cnts <- obj@assays[["RNA"]]@counts
    } else if (inherits(x = obj, what = "SingleCellExperiment")) {
        cnts <- counts(obj)
    } else {
        cnts <- obj
    }
    if (!inherits(cnts, what = "dgCMatrix")) {
        cnts <- Matrix::Matrix(as.matrix(cnts),sparse = T)
    }
    cnts <- cnts[tabulate(summary(cnts)$i) != 0, , drop = FALSE]
    
    if(inherits(egc, what = "GeneSetCollection")){
        egc <- GSEABase::geneIds(egc) # will return a simple list, 
        #which will work if a matrix is supplied to GSVA
    }
    scores <- list()
    wind <- seq(1, ncol(cnts), by=groups)
    print(paste('Using sets of', groups, 'cells. Running', 
                length(wind), 'times.'))
    if (method == "ssGSEA") {
        # break to groups of cells
        split.data <- split_data.matrix(matrix=cnts, chunk.size=groups)
        for (i in seq_along(wind)) {
            last <- min(ncol(cnts), i+groups-1)
            a <- suppressWarnings(gsva(split.data[[i]], egc, method = 'ssgsea', 
                ssgsea.norm = TRUE, kcdf = "Poisson", parallel.sz = cores, 
                BPPARAM = SnowParam()))
            scores[[i]] <- a
        }
    } else if (method == "UCell") {
        scores[[1]] <- suppressWarnings(ScoreSignatures_UCell(cnts, features=egc, 
                                        chunk.size = groups, ncores = 1))
    }
    scores <- do.call(cbind, scores)
    output <- data.frame(t(as.matrix(scores)))
    return(output)
}
