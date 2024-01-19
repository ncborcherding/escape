"%!in%" <- Negate("%in%")

.checkSingleObject <- function(sc) {
  if (!inherits(x=sc, what ="Seurat") & 
      !inherits(x=sc, what ="SummarizedExperiment")){
    stop("Object indicated is not of class 'Seurat' or 
            'SummarizedExperiment', make sure you are using
            the correct data.") }
}


#Pulling a color palette for visualizations
#' @importFrom grDevices hcl.colors
#' @keywords internal
.colorizer <- function(palette = "inferno", 
                       n= NULL) {
  colors <- hcl.colors(n=n, palette = palette, fixup = TRUE)
  return(colors)
}

#split data matrix into cell chunks
#stole this from https://github.com/carmonalab/UCell
.split_data.matrix <- function(matrix, chunk.size=1000) {
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

#' @importFrom SingleCellExperiment counts
#' @importFrom Matrix summary Matrix
.cntEval <- function(obj) {
  if (inherits(x = obj, what = "Seurat")) {
    cnts <- obj@assays[["RNA"]]@counts
  } else if (inherits(x = obj, what = "SingleCellExperiment")) {
    cnts <- counts(obj)
  } else {
    cnts <- obj
  }
  cnts <- cnts[rowSums2(cnts) != 0,]
  return(cnts)
}

#Add the values to single cell object
#' @importFrom SeuratObject CreateAssayObject
#' @importFrom SingleCellExperiment reducedDim
.adding.Enrich <- function(sc, enrichment, enrichment.name) {
  if (inherits(sc, "Seurat")) {
    new.assay <- suppressWarnings(CreateAssayObject(
                                  data = as.matrix(t(enrichment))))
    sc[[enrichment.name]] <- new.assay
  } else if (inherits(sc, "SingleCellExperiment")) {
    assays(sc, enrichment.name) <- enrichment
  }
  return(sc)
}

#' @importFrom GSEABase geneIds
.GS.check <- function(gene.sets) {
  if(is.null(gene.sets)) {
    stop("Please provide the gene.sets you would like to use for 
            the enrichment analysis")
  }
  egc <- gene.sets
  if(inherits(egc, what = "GeneSetCollection")){
    egc <- GSEABase::geneIds(egc) # will return a simple list, 
    #which will work if a matrix is supplied to GSVA
  }
  return(egc)
}

#This is to grab the meta data from a seurat or SCE object
#' @importFrom SingleCellExperiment colData 
#' @importFrom SeuratObject Idents
#' @importFrom methods slot
#' @keywords internal
.grabMeta <- function(sc) {
  if (is_seurat_object(sc)) {
    meta <- data.frame(sc[[]], slot(sc, "active.ident"))
    colnames(meta)[length(meta)] <- "ident"
    
  } else if (is_se_object(sc)){
    meta <- data.frame(colData(sc))
    rownames(meta) <- sc@colData@rownames
    clu <- which(colnames(meta) == "ident")
    colnames(meta)[clu] <- "ident"
  } else {
    stop("Object indicated is not of class 'Seurat' or 
            'SummarizedExperiment', make sure you are using
            the correct data.")
  }
  return(meta)
}

# Add to meta data some of the metrics calculated
#' @importFrom rlang %||%
#' @importFrom SummarizedExperiment colData colData<-
.add.meta.data <- function(sc, meta, header) {
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