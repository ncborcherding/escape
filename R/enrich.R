#' Calculate gene set enrichment scores for single-cell data
#'
#' This function allows users to input both the single-cell RNA-sequencing 
#' counts and any gene set pathways either from the stored data or from 
#' other sources. The enrichment calculation itself 
#' uses the two methods 1) gsva R package and the poisson distribution for RNA
#' or the \href{https://github.com/carmonalab/UCell}{UCell package}. 
#'
#' @param obj The count matrix, Seurat, or SingleCellExperiment object.
#' @param gene.sets Gene sets from \code{\link{getGeneSets}} to use 
#' for the enrichment analysis. Alternatively a simple base R list where
#' the names of the list elements correspond to the name of the gene set
#' and the elements themselves are simple vectors of gene names representing
#' the gene set. 
#' @param method select the method to calculate enrichment, either "ssGSEA", "UCell" or
#' "singscore"
#' @param groups The number of cells to separate the enrichment calculation.
#' @param cores The number of cores to use for parallelization.
#' @param weight.by.nFeatures Weight the normalized enrichment score by the inverse
#' of the nFeature by cell (recommended if using ssGSEA).
#'
#' @importFrom GSVA gsva
#' @importFrom GSEABase GeneSetCollection 
#' @importFrom UCell ScoreSignatures_UCell
#' @importFrom singscore rankGenes simpleScore
#' @importFrom BiocParallel SnowParam
#'
#' 
#' @examples 
#' GS <- list(Bcells = c("MS4A1", "CD79B", "CD79A),
#'   Tcells = c("CD3E","CD7","CD8A",))
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
                     method = "ssGSEA", groups = 1000, cores = 2, 
                     weight.by.nFeatures = TRUE) {
    egc <- GS.check(gene.sets)
    cnts <- cntEval(obj)
    nFeature = apply(cnts,2,function(x) sum(x > 0))
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
                                        chunk.size = groups, ncores = cores))
    } else if (method == "singscore") {
        tmp.list <- list()
        rankData <- rankGenes(as.matrix(cnts))
        for (i in seq_along(egc)) {
            tmp <- simpleScore(rankData, egc[[i]])
            colnames(tmp)[1] <- names(egc)[i]
            tmp.list[[i]] <- t(tmp[,1])
        }
        tmp.list <- do.call(rbind, tmp.list)
        rownames(tmp.list) <- names(egc)
        colnames(tmp.list) <- colnames(cnts)
        scores[[1]] <- tmp.list
    }
    scores <- do.call(cbind, scores)
    output <- t(as.matrix(scores))
    if (weight.by.nFeatures) {
        nFeature <- 0.5+normalize(nFeature)
        output <- apply(output,2 , function(x) x*nFeature)
    }
    output <- data.frame(output)
    return(output)
}
