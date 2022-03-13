#Ucell is available at https://github.com/carmonalab/UCell 
#below is the code taken from the github repo in order
#to be compatible with Bioconductor

#' Calculate module enrichment scores from single-cell data
#' 
#' @param matrix Input matrix, 
#' @param features A list of signatures
#' @param precalc.ranks If you have pre-calculated ranks
#' @param maxRank Maximum number of genes to rank per cell; above this rank, a given gene is considered as not expressed.
#' @param assay The sce object assay where the data is to be found
#' @param chunk.size Number of cells to be processed simultaneously (lower size requires slightly more computation but reduces memory demands)
#' @param w_neg Weight on negative genes in signature. e.g. `w_neg=1` weighs equally up- and down-regulated genes, `
#' @param ncores Number of processors to parallelize computation.
#' @param ties.method How ranking ties should be resolved (passed on to [data.table::frank])
#' @param name Name suffix appended to signature names
#' @param force.gc Explicitly call garbage collector to reduce memory footprint
#' @param seed Integer seed
#' @importFrom methods is 
#' @importFrom  SingleCellExperiment altExp
#' @importFrom SummarizedExperiment assay SummarizedExperiment
#' @importFrom stringr str_remove_all
#' @return Returns input SingleCellExperiment object with UCell scores added to altExp
ScoreSignatures_UCell <- function(matrix=NULL, features, precalc.ranks=NULL, 
                                  maxRank=1500, w_neg=1, name="_UCell",
                                  assay="counts", chunk.size=1000, ncores=1, 
                                  ties.method="average", force.gc=FALSE, seed=123) {
    
    features <- check_signature_names(features)
    
    #Check type of input
    if (is(matrix, "SingleCellExperiment")) { # sce object
        if (!assay %in% names(matrix@assays)) {
            stop(sprintf("Assay %s not found in sce object.", assay))
        }
        m <- assay(matrix, assay) 
    } else if (is(matrix, "matrix") | #matrix or DF
               is(matrix, "dgCMatrix") |
               is(matrix, "data.frame")) { 
        m <- matrix
    } else {
        m <- NULL
    }
    
    if (is.null(m) & is.null(precalc.ranks)) {
        stop("Unrecognized input format.")
    }
    
    #Run on pre-calculated ranks ('m' can be NULL)
    if (!is.null(precalc.ranks)) {
        u.list <- rankings2Uscore(precalc.ranks, features=features, chunk.size=chunk.size,
                                  w_neg=w_neg, ncores=ncores, force.gc=force.gc, name=name)
    } else {
        u.list <- calculate_Uscore(m, features=features, maxRank=maxRank, chunk.size=chunk.size,
                                   w_neg=w_neg, ties.method=ties.method, ncores=ncores,
                                   force.gc=force.gc, name=name)
    }
    u.merge <- lapply(u.list,function(x) rbind(x[["cells_AUC"]]))
    u.merge <- Reduce(rbind, u.merge)
    
    if (is(matrix, "SingleCellExperiment")) {
        altExp(matrix, "UCell") <- SummarizedExperiment(assays = list("UCell" = t(u.merge)))
        return(matrix)
    } else {
        colnames(u.merge) <- str_remove_all(colnames(u.merge), "_UCell")
        return(u.merge)
    }
    
}

#' Calculate and store gene rankings for a single-cell dataset
#' @param matrix Input matrix, either stored in a [SingleCellExperiment] object or as a raw matrix. \code{dgCMatrix} format supported.
#' @param maxRank Maximum number of genes to rank per cell; above this rank, a given gene is considered as not expressed
#' @param assay Assay where the data is to be found (for input in 'sce' format)
#' @param chunk.size Number of cells to be processed simultaneously (lower size requires slightly more computation but reduces memory demands)
#' @param ncores Number of processors to parallelize computation
#' @param ties.method How ranking ties should be resolved (passed on to [data.table::frank])
#' @param force.gc Explicitly call garbage collector to reduce memory footprint
#' @param seed Integer seed
#' @return Returns a sparse matrix of pre-calculated ranks that can be used multiple times to evaluate different signatures
#' @importFrom methods is 
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment assay
#' @return Returns a sparse matrix of pre-calculated ranks that can be used multiple times to evaluate different signatures
StoreRankings_UCell <- function(matrix, maxRank=1500, chunk.size=1000,
                                ncores=1, assay='counts', ties.method="average",
                                force.gc=FALSE, seed=123) {
    
    #Check type of input
    if (is(matrix, "SingleCellExperiment")) { # sce object
        if (!assay %in% names(matrix@assays)) {
            stop(sprintf("Assay %s not found in sce object.", assay))
        }
        m <- assay(matrix, assay)
    } else if (is(matrix, "matrix") | #matrix or DF
               is(matrix, "dgCMatrix") |
               is(matrix, "data.frame")) { 
        m <- matrix
    } else {
        stop("Unrecognized input format.")
    }
    
    features <- rownames(m)[1]  #placeholder signature
    meta.list <- calculate_Uscore(m, features=features, maxRank=maxRank,
                                  chunk.size=chunk.size, ncores=ncores,
                                  ties.method=ties.method, storeRanks=TRUE,
                                  force.gc=force.gc)
    
    ranks.all <- lapply(meta.list,function(x) rbind(x[["cells_rankings"]]))
    ranks.all <- Reduce(cbind, ranks.all)
    
    return(ranks.all)
    
}

#' Calculate Mann Whitney U from a vector of ranks
#' @param rank_value A vector of ranks
#' @param maxRank Max number of features to include in ranking
#' @param sparse Whether the vector of ranks is in sparse format
#' @return Normalized AUC (as U statistic) for the vector
u_stat <- function(rank_value, maxRank=1000, sparse=FALSE){
    
    if (sparse==TRUE){
        rank_value[rank_value==0] <- maxRank+1
    }
    
    insig <- rank_value > maxRank
    if (all(insig)) {
        return(0L)
    } else {
        rank_value[insig] <- maxRank+1
        rank_sum <- sum(rank_value)
        len_sig <- length(rank_value)
        
        u_value <- rank_sum - (len_sig * (len_sig + 1))/2
        auc <- 1 - u_value/(len_sig * maxRank)
        return(auc)
    }
}

#' Calculate U scores for a list of signatures, given a rank matrix
#' @param   sig_list      A list of signatures
#' @param   ranks_matrix  Matrix of pre-computed ranks
#' @param   maxRank       Max number of features to include in ranking, for u_stat function
#' @param   sparse        Whether the vector of ranks is in sparse format
#' @param   w_neg         Weight on negative signatures
#' @return A matrix of U scores
u_stat_signature_list <- function(sig_list, ranks_matrix, maxRank=1000,
                                  sparse=FALSE, w_neg=1) {
    
    dim <- ncol(ranks_matrix)-1
    u_matrix <- vapply(sig_list, FUN.VALUE = numeric(dim), FUN=function(sig) {
        sig_neg <- grep('-$', unlist(sig), perl=TRUE, value=TRUE)
        sig_pos <- setdiff(unlist(sig), sig_neg)
        
        if (length(sig_pos)>0) {
            sig_pos <- gsub('\\+$','',sig_pos,perl=TRUE)
            u_p <- as.numeric(ranks_matrix[sig_pos,
                                           lapply(.SD, function(x) u_stat(x,maxRank = maxRank,sparse=sparse)),
                                           .SDcols=-1, on="rn"])
        } else {
            u_p <- rep(0, dim(ranks_matrix)[2]-1)
        }
        if (length(sig_neg)>0) {
            sig_neg <- gsub('-$','',sig_neg,perl=TRUE)
            u_n <- as.numeric(ranks_matrix[sig_neg,
                                           lapply(.SD, function(x) u_stat(x,maxRank = maxRank,sparse=sparse)),
                                           .SDcols=-1, on="rn"])
        } else {
            u_n <- rep(0, dim(ranks_matrix)[2]-1)
        }
        
        diff <- u_p - w_neg*u_n   #Subtract negative sets, if any
        diff[diff<0] <- 0
        return(diff)
    })
    
    rownames(u_matrix) <- colnames(ranks_matrix)[-1]
    return (u_matrix)
}

#' Calculate rankings and scores for query data and given signature set
#' @param   matrix        Input data matrix 
#' @param   features      List of signatures
#' @param   maxRank       Rank cutoff (1500) 
#' @param   chunk.size    Cells per sub-matrix (1000) 
#' @param   ncores        Number of cores to use for parallelization (1) 
#' @param   w_neg         Weight on negative signatures
#' @param   ties.method   How to break ties, for data.table::frankv method ("average")
#' @param   storeRanks    Store ranks? (FALSE) 
#' @param   force.gc      Force garbage collection? (FALSE) 
#' @param   name          Suffix for metadata columns ("_UCell") 
#' @importFrom methods is 
#' @importFrom  Matrix Matrix
#' @importFrom BiocParallel MulticoreParam bplapply
#' @return  A list of signature scores
calculate_Uscore <- function(matrix, features,  maxRank=1500, chunk.size=1000,
                             ncores=1, w_neg=1, ties.method="average", 
                             storeRanks=FALSE, force.gc=FALSE, name="") {
    
    #Make sure we have a sparse matrix
    if (!is(matrix, "dgCMatrix")) {
        matrix <- Matrix(as.matrix(matrix),sparse = TRUE)
    }
    
    #Check if all genes in signatures are present in the data matrix
    matrix <- check_genes(matrix, features)
    
    #Do not evaluate more genes than there are
    if (maxRank > nrow(matrix)) {
        maxRank <- nrow(matrix)
    }
    
    #Weight on neg signatures must be >=0
    if (is.null(w_neg)) {w_neg <- 1}
    if (w_neg<0) {stop("Weight on negative signatures (w_neg) must be >=0")}
    
    #Signatures cannot be larger than maxRank parameter
    sign.lgt <- lapply(features, length)
    if (any(sign.lgt > maxRank)) {
        stop("One or more signatures contain more genes than maxRank parameter. Increase maxRank parameter or make shorter signatures")
    }
    
    #Split into manageable chunks
    split.data <- split_data.matrix(matrix=matrix, chunk.size=chunk.size)
    
    #Parallelize?
    param <- MulticoreParam(workers=ncores)
    meta.list <- bplapply(
        X = split.data, 
        BPPARAM =  param,
        FUN = function(x) {
            cells_rankings <- data_to_ranks_data_table(x, ties.method = ties.method)
            cells_AUC <- u_stat_signature_list(features, cells_rankings, 
                                               maxRank=maxRank, sparse=FALSE, w_neg=w_neg)
            colnames(cells_AUC) <- paste0(colnames(cells_AUC),name)
            if (storeRanks==TRUE){
                gene.names <- as.character(as.matrix(cells_rankings[,1]))
                #make sparse
                cells_rankings[cells_rankings>maxRank] <- 0
                ranks.sparse <- Matrix(as.matrix(cells_rankings[,-1]),sparse = TRUE)
                dimnames(ranks.sparse)[[1]] <- gene.names
                if (force.gc) {
                    cells_rankings <- NULL
                    gc()
                }
                return(list(cells_rankings=ranks.sparse, cells_AUC=cells_AUC))
            } else {
                if (force.gc) {
                    cells_rankings <- NULL
                    gc()
                }
                return(list(cells_AUC=cells_AUC))
            }
            
        })
    return(meta.list)
}

#' Get signature scores from pre-computed rank matrix
#' @param     ranks_matrix    A rank matrix
#' @param     features        List of signatures
#' @param     chunk.size      How many cells per matrix chunk
#' @param     w_neg           Weight on negative signatures
#' @param     ncores          How many cores to use for parallelization
#' @param     force.gc        Force garbage collection to recover RAM? (FALSE)
#' @param     name            Name suffix for metadata columns ("_UCell")
#' @importFrom data.table as.data.table
#' @importFrom BiocParallel MulticoreParam bplapply
#' @return                    A list of signature scores
rankings2Uscore <- function(ranks_matrix, features, chunk.size=1000, w_neg=1,
                            ncores=1, force.gc=FALSE, name="_UCell") {
    
    #Check if all genes in signatures are present in the stored signatures
    ranks_matrix <- check_genes(ranks_matrix, features)
    
    #Weight on neg signatures must be >=0
    if (is.null(w_neg)) {w_neg <- 1}
    if (w_neg<0) {stop("Weight on negative signatures (w_neg) must be >=0")}
    
    maxRank <- max(ranks_matrix)
    split.data <- split_data.matrix(matrix=ranks_matrix, chunk.size=chunk.size)
    rm(ranks_matrix)
    
    #Parallelize?
    param <- MulticoreParam(workers=ncores)
    
    meta.list <- bplapply(
        X = split.data, 
        BPPARAM =  param,
        FUN = function(x) {
            
            dense <- as.matrix(x)
            dense <- as.data.table(dense, keep.rownames=TRUE)
            setkey(dense, "rn", physical=FALSE)
            
            cells_AUC <- u_stat_signature_list(features, dense, maxRank=maxRank, sparse=TRUE, w_neg=w_neg)
            colnames(cells_AUC) <- paste0(colnames(cells_AUC),name)
            
            if (force.gc) {
                dense <- NULL
                gc()
            }
            return(list(cells_AUC=cells_AUC))
        }
    )
    return(meta.list)
}

#' Check if all genes in signatures are found in data matrix - otherwise add zero counts in data-matrix to complete it
#' @param matrix Input data matrix
#' @param features List of genes that must be present (otherwise they are added)
#' @importFrom Matrix sparseMatrix
#' @return Same input matrix, extended to comprise any missing genes
check_genes <- function(matrix, features) {
    features <- unlist(features)
    features <- gsub("[-+]$","",features,perl=TRUE)
    missing <- setdiff(features, rownames(matrix))
    ll <- length(missing)
    
    if (ll/length(features) > 0.5) {
        warning(sprintf("Over half of genes (%s%%) in specified signatures are missing from data. Check the integrity of your dataset\n", 
                        round(100*ll/length(features))))
    }
    
    if (ll>0) {
        add.mat <- sparseMatrix(length(missing), ncol(matrix))
        rownames(add.mat) <- missing
        matrix <- rbind(matrix, add.mat)
        
        missing.concatenate <- paste(missing, collapse=",")
        warning(sprintf("The following genes were not found and will be imputed to exp=0:\n* %s",missing.concatenate))
    }
    return(matrix)
}

#' Check signature names and add standard names is missing
#' @param features List of signatures for scoring
#' @return The input list of signatures, with standard names if provided un-named
check_signature_names <- function(features) {
    defaultSigName <- paste0(rep("signature_",length(features)),seq_along(features))
    if(is.null(names(features))){
        names(features) <- defaultSigName
    } else {
        invalidNames <- names(features) == "" | duplicated(names(features))
        names(features)[invalidNames] <- defaultSigName[invalidNames]
    }
    return(features)
}


#' Calculate per-cell feature rankings
#' @param data              Expression data matrix 
#' @param ties.method       How to break ties (passed on to data.table::frankv)
#' @importFrom data.table frankv setkey
#' @return                  A data.table of ranks 
data_to_ranks_data_table <- function(data, ties.method="average") {
    dt <- as.data.table(as.matrix(data))
    rnaDT.ranks.dt <- dt[, lapply(.SD, function(x) frankv(x,ties.method=ties.method,order=c(-1L)))]
    rnaDT.ranks.rownames <- rownames(data)
    rnaDT.ranks.dt.rn <- cbind(rn=rnaDT.ranks.rownames, rnaDT.ranks.dt)
    setkey(rnaDT.ranks.dt.rn, "rn", physical = FALSE)
    return(rnaDT.ranks.dt.rn)
}

split_data.matrix <- function(matrix, chunk.size=1000) {
    ncols <- dim(matrix)[2]
    nchunks <- (ncols-1) %/% chunk.size + 1
    
    split.data <- list()
    min <- 1
    for (i in seq_len(nchunks)) {
        if (i == nchunks-1) {  #make last two chunks of equal size
            left <- ncols-(i-1)*chunk.size
            max <- min+round(left/2)-1
        } else {
            max <- min(i*chunk.size, ncols)
        }
        split.data[[i]] <- matrix[,min:max]
        min <- max+1    #for next chunk
    }
    return(split.data)
}
