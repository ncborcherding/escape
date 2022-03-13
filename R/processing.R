"%!in%" <- Negate("%in%")

#' Calculate Principal Components for the Enrichment Scores
#'
#' Using all or selected enrichment scores of individual 
#' single-cells, this function will calculate 
#' principal components using scaled values and attach 
#' to the output columns to use to graph later.
#'
#' @param enriched The output of \code{\link{enrichIt}}.
#' @param gene.sets Names of gene sets to include in the PCA
#' @param groups The column headers to use in future graphing functions.
#'
#' @importFrom dplyr select_if
#' @importFrom stats prcomp
#' 
#' @examples 
#' ES2 <- readRDS(url(
#' "https://ncborcherding.github.io/vignettes/escape_enrichment_results.rds"))
#' 
#' PCA <- performPCA(enriched = ES2, groups = c("Type", "Cluster"), 
#' gene.sets = colnames(ES2))
#'
#' @export
#' @return Data frame of principal components
#'
#' @author Nick Borcherding
#'
performPCA <- function(enriched, gene.sets = NULL, groups) {
    groups <- enriched[,colnames(enriched) %in% c(groups)]
    input <- select_if(enriched, is.numeric)
    if (!is.null(gene.sets)) {
      input <- input[,colnames(input) %in% gene.sets]
    }
    PCA <- prcomp(input, scale. = TRUE)
    merged <- merge(PCA$x, groups, by = "row.names")
    rownames(merged) <- merged[,1]
    merged <- merged[,-1]

    return(merged)
}

#split data matrix into cell chunks
#stole this from https://github.com/carmonalab/UCell
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

#' Get a collection of gene sets to perform enrichment on
#'
#' This function allows users to select libraries and specific 
#' gene.sets to form a GeneSetCollection that is a list of gene sets.
#
#' @param species The scientific name of the species of interest in 
#' order to get correcent gene nomenclature
#' @param library Individual collection(s) of gene sets, e.g. c("H", "C5").
#' See \href{https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp}{msigdbr}for
#' all MSigDB collections.
#' @param subcategory MSigDB sub-collection abbreviation, such as CGP or BP.
#' @param gene.sets Select gene sets or pathways, using specific names, 
#' example: pathways = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB"). Will only be
#' honored if library is set, too.
#'
#' @examples 
#' GS <- getGeneSets(library = "H")
#' 
#' @export
#' 
#' @importFrom GSEABase GeneSet GeneSetCollection
#' @importFrom msigdbr msigdbr msigdbr_species
#' 
#' @author Nick Borcherding, Jared Andrews
#' @return A \code{GeneSetCollection} object containing the requested \code{GeneSet} objects.
getGeneSets <- function(species = "Homo sapiens", 
                        library = NULL, 
                        subcategory = NULL,
                        gene.sets = NULL) {
    spec <- msigdbr_species()
    spec_check <- unlist(spec[spec$species_name %in% species,][,1])
    if (length(spec_check) == 0) {
        message(paste0("Please select a compatible species: ", 
                    paste(spec, collapse = ", ")))
    }
    if(!is.null(library)) {
        if (length(library) == 1) {
          if (is.null(subcategory)) {
            m_df = msigdbr(species = spec_check, category = library)
          } else {
            m_df = msigdbr(species = spec_check, category = library, subcategory = subcategory)
          }
        }
        m_df <- NULL
        for (x in seq_along(library)) {
          if (is.null(subcategory)) {
            tmp2 = msigdbr(species = spec_check, category = library[x])
          } else {
            tmp2 = msigdbr(species = spec_check, category = library, subcategory = subcategory)
          }
          m_df <- rbind(m_df, tmp2)
        }
        if(!is.null(gene.sets)) {
        m_df <- m_df[m_df$gs_name %in% gene.sets,]
        }    
    }
    gs <- unique(m_df$gs_name)
    ls <- list()
    for (i in seq_along(gs)) {
        tmp <- m_df[m_df$gs_name == gs[i],]
        tmp <- tmp$gene_symbol
        tmp <- unique(tmp)
        tmp <- GeneSet(tmp, setName=paste(gs[i]))
        ls[[i]] <- tmp
    }
    gsc <- GeneSetCollection(ls)
    return(gsc)
}

#Function for normalizing value
normalize <- function(x)
{
  (x- min(x)) /(max(x)-min(x))
}

#' @importFrom SingleCellExperiment counts
#' @importFrom Matrix summary Matrix
cntEval <- function(obj) {
  if (inherits(x = obj, what = "Seurat")) {
    cnts <- obj@assays[["RNA"]]@counts
  } else if (inherits(x = obj, what = "SingleCellExperiment")) {
    cnts <- counts(obj)
  } else {
    cnts <- obj
  }
  if (!inherits(cnts, what = "dgCMatrix")) {
    cnts <- Matrix(as.matrix(cnts),sparse = TRUE)
  }
  cnts <- cnts[tabulate(summary(cnts)$i) != 0, , drop = FALSE]
  return(cnts)
}

#' @importFrom GSEABase geneIds
GS.check <- function(gene.sets) {
  if(is.null(gene.sets)) {
    stop("Please provide the gene.sets you would like to use for 
            the enrichment analysis")
  } else {
    egc <- gene.sets
  }
  if(inherits(egc, what = "GeneSetCollection")){
    egc <- GSEABase::geneIds(egc) # will return a simple list, 
    #which will work if a matrix is supplied to GSVA
  }
  return(egc)
}

#This is to grab the meta data from a seurat or SCE object
#' @importFrom SingleCellExperiment colData 
grabMeta <- function(sc) {
  if (inherits(x=sc, what ="Seurat")) {
    meta <- data.frame(sc[[]], slot(sc, "active.ident"))
    if ("cluster" %in% colnames(meta)) {
      colnames(meta)[length(meta)] <- "cluster.active.ident"
    } else {
      colnames(meta)[length(meta)] <- "cluster"
    }
  }
  else if (inherits(x=sc, what ="SingleCellExperiment")){
    meta <- data.frame(colData(sc))
    rownames(meta) <- sc@colData@rownames
    clu <- which(colnames(meta) == "ident")
    if ("cluster" %in% colnames(meta)) {
      colnames(meta)[clu] <- "cluster.active.idents"
    } else {
      colnames(meta)[clu] <- "cluster"
    }
  }
  return(meta)
}

# Add to meta data some of the metrics calculated
#' @importFrom rlang %||%
#' @importFrom SummarizedExperiment colData colData<-
add.meta.data <- function(sc, meta, header) {
  if (inherits(x=sc, what ="Seurat")) { 
    col.name <- names(meta) %||% colnames(meta)
    sc[[col.name]] <- meta
  } else {
    rownames <- rownames(colData(sc))
    colData(sc) <- cbind(colData(sc), 
                         meta[rownames,])[, union(colnames(colData(sc)),  colnames(meta))]
    rownames(colData(sc)) <- rownames  
  }
  return(sc)
}
