#' Calculate gene set enrichment scores for single-cell data
#'
#' This function allows users to input both the single-cell RNA-sequencing 
#' counts and any gene set pathways either from the stored data or from 
#' other sources. The enrichment calculation itself 
#' uses the two methods 1) gsva R package and the poisson distribution for RNA
#' or 2) the \href{https://github.com/carmonalab/UCell}{UCell package}. 
#'
#' @param obj The count matrix, Seurat, or SingleCellExperiment object.
#' @param gene.sets Gene sets from \code{\link{getGeneSets}} to use 
#' for the enrichment analysis. Alternatively a simple base R list where
#' the names of the list elements correspond to the name of the gene set
#' and the elements themselves are simple vectors of gene names representing
#' the gene set. 
#' @param method select the method to calculate enrichment, either "ssGSEA" or "UCell" 
#' @param groups The number of cells to separate the enrichment calculation.
#' @param cores The number of cores to use for parallelization.
#' @param min.size Minimum number of gene necessary to perform the enrichment
#' calculation
#' @param ssGSEA.norm normalized the enrichment score based on the range of the
#' individual gene set. If TRUE, the returned enrichment score is based may change
#' with cell composition.
#' @param ... pass arguments to ssGSEA or UCell call
#'
#' @importFrom GSVA gsva
#' @importFrom GSEABase GeneSetCollection 
#' @importFrom BiocParallel SnowParam
#'
#' 
#' @examples 
#' GS <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
#'   Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))
#' pbmc_small <- suppressWarnings(SeuratObject::pbmc_small)
#' ES <- enrichIt(obj = pbmc_small, gene.sets = GS, min.size = NULL)
#' 
#' @export
#' @author Nick Borcherding, Jared Andrews
#'
#' @seealso \code{\link{getGeneSets}} to collect gene sets.
#' @return Data frame of normalized enrichmenet scores (NES)
enrichIt <- function(obj, gene.sets = NULL, 
                     method = "ssGSEA", 
                     groups = 1000, 
                     cores = 2,
                     min.size = 5,
                     ssGSEA.norm = FALSE,
                     ...) {
    egc <- GS.check(gene.sets)
    cnts <- cntEval(obj)
    if (!is.null(min.size)){
        GS.size <- lapply(egc, function(x) length(which(rownames(cnts) %in% x)))
        remove <- unname(which(GS.size < min.size))
        if (length(remove) != 0) {
            egc <- egc[-remove]
        }
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
                ssgsea.norm = FALSE,
                kcdf = "Poisson", parallel.sz = cores, 
                BPPARAM = SnowParam()),
                ...)
            scores[[i]] <- a
        }
    } else if (method == "UCell") {
        scores[[1]] <- t(suppressWarnings(ScoreSignatures_UCell(cnts, features=egc, 
                                        chunk.size = groups, ncores = cores,
                                        ...)))
    }
    scores <- do.call(cbind, scores)
    output <- t(as.matrix(scores))
    if(method == "ssGSEA" & ssGSEA.norm) {
        output <- apply(output, 2, normalize)
    }
    output <- data.frame(output)
    return(output)
}
