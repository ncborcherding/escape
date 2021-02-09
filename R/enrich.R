#' Calculate gene set enrichment scores for single-cell data
#'
#' This function allows users to input both the single-cell RNA-sequencing 
#' counts and any gene set pathways either from the stored data or from 
#' other sources. The enrichment calculation itself 
#' uses the gsva R package and the poisson distribution for RNA.
#'
#' @param obj The count matrix, Seurat, or SingleCellExperiment object.
#' @param gene.sets Gene sets from \code{\link{getGeneSets}} to use 
#' for the enrichment analysis. Alternatively a simple base R list where
#' the names of the list elements correspond to the name of the gene set
#' and the elements themselves are simple vectors of gene names representing
#' the gene set. 
#' @param groups The number of cells to separate the enrichment calculation.
#' @param cores The number of cores to use for parallelization.
#'
#' @importFrom GSVA gsva
#' @importFrom GSEABase GeneSetCollection
#' @importFrom SingleCellExperiment counts
#' @importFrom BiocParallel SnowParam
#'
#' 
#' @examples 
#' GS <- getGeneSets(library = "H")
#' GS <- GS[[1]] #Reduce list size for example
#' seurat_ex <- suppressWarnings(Seurat::pbmc_small)
#' ES <- enrichIt(obj = seurat_ex, gene.sets = GS)
#'
#' @export
#'
#' @author Nick Borcherding, Jared Andrews
#'
#' @seealso \code{\link{getGeneSets}} to collect gene sets.
#' @return Data frame of normalized enrichmenet scores (NES)
enrichIt <- function(obj, gene.sets = NULL, groups = 1000, cores = 2) {
    
    if(is.null(gene.sets)) {
        stop("Please provide the gene.sets you would like to use for 
            the enrichment analysis")
    } else {
        egc <- gene.sets
    }
    
    if (inherits(x = obj, what = "Seurat")) {
        cnts <- as.matrix(as.matrix(obj@assays[["RNA"]]@counts))
    } else if (inherits(x = obj, what = "SingleCellExperiment")) {
        cnts <- as.matrix(counts(obj))
    } else {
        cnts <- obj
    }
    
    # egc <- GeneSetCollection(egc) ## maybe it's a version thing, but supplying a matrix with a GeneSet did not work for me
    # names <- NULL
    if( attr(class(egc), "package") == "GSEABase"){
        egc <- GSEABase::geneIds(egc) # will return a simple list, which will work if a matrix is supplied to GSVA
    }
    
    #for (x in seq_along(egc)) {
    #    setName <- egc[[x]]@setName
    #    names <- c(names, setName)
    #} ## this doesn't seem to serve a purpose; otherwise unlist(lapply(egc, function(x) x@setName)) does the job, too
    ## I wouldn't overwrite the inbuilt names() function, though
    
    cnts <- cnts[rowSums(cnts > 0) != 0, ] 
    # break to groups of cells
    scores <- list()
    wind <- seq(1, ncol(cnts), by=groups)
    print(paste('Using sets of', groups, 'cells. Running', 
            length(wind), 'times.'))
    
    for (i in wind) {
        last <- min(ncol(cnts), i+groups-1)
        a <- suppressWarnings(gsva(cnts[,i:last], egc, method = 'ssgsea', 
            ssgsea.norm = TRUE, kcdf = "Poisson", parallel.sz = cores, 
            BPPARAM = SnowParam()))
        scores[[i]] <- a
    }
    
    scores <- do.call(cbind, scores)
    output <- data.frame(t(scores))
    
    return(output)
}
