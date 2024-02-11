"%!in%" <- Negate("%in%")

is_seurat_object <- function(obj) inherits(obj, "Seurat")
is_se_object <- function(obj) inherits(obj, "SummarizedExperiment")
is_seurat_or_se_object <- function(obj) {
  is_seurat_object(obj) || is_se_object(obj)
}

.checkSingleObject <- function(sc) {
  if (!inherits(x=sc, what ="Seurat") & 
      !inherits(x=sc, what ="SummarizedExperiment")){
    stop("Object indicated is not of class 'Seurat' or 
            'SummarizedExperiment', make sure you are using
            the correct data.") }
}

#' @importFrom dplyr group_by summarise_at
#' @importFrom stringr str_sort
.orderFunction <- function(dat, order.by, group.by){
  if(order.by %!in% c("mean", "group.by")) {
    stop(paste0("Please select either 'mean' or 'group.by' for ordering."))
  }
  if(order.by == "mean") {
    summary <- dat %>%
                  group_by(dat[,group.by]) %>%
                  summarise_at(.vars = colnames(.)[1], mean) %>%
                  as.data.frame()
    summary <- summary[order(summary[,2], decreasing = TRUE),]
    dat[,group.by] <- factor(dat[,group.by], levels = summary[,1])
  }
  else if (order.by == "group.by") {
    dat[,group.by] <- factor(dat[,group.by], str_sort(unique(dat[,group.by]), numeric = TRUE))
  }
  return(dat)
}

.makeDFfromSCO <- function(input.data, 
                           assay = "escape", 
                           gene.set = NULL,
                           group.by = NULL, 
                           split.by = NULL, 
                           facet.by = NULL) {
  if(is.null(assay)){
    stop("Please add the assay name in which to plot from")
  }
  columns <- unique(c(group.by, split.by, facet.by))
  cnts <- .cntEval(input.data, 
                   assay = assay, 
                   type = "data")
  if(length(gene.set) == 1 && gene.set == "all") {
    gene.set <- rownames(cnts)
  }
  meta <- .grabMeta(input.data)
  if(length(gene.set) == 1) {
    enriched <- data.frame(cnts[gene.set,], meta[,columns])
  } else {
    enriched <- data.frame(t(cnts[gene.set,]), meta[,columns])
  }
  colnames(enriched) <- c(gene.set, columns)
  return(enriched)
}

#Prepare Data
.prepData <- function(input.data, assay, gene.set, group.by, split.by, facet.by) {
  
  if (inherits(x=input.data, what ="Seurat") || 
      inherits(x=input.data, what ="SummarizedExperiment")) {
    enriched <- .makeDFfromSCO(input.data, assay, gene.set, group.by, split.by, facet.by)
    if(length(gene.set) == 1 && gene.set == "all") {
      gene.set <- colnames(enriched)[colnames(enriched) %!in% c(group.by, split.by, facet.by)]
      gene.set <- gene.set[!grepl("meta", gene.set)]
    }
  } else if (!is_seurat_or_se_object(input.data)) {
    if(length(gene.set) == 1 && gene.set == "all") {
      gene.set <- colnames(input.data)
      gene.set <- gene.set[gene.set %!in% c(group.by, split.by, facet.by)]
    } 
      enriched <- data.frame(input.data[,c(gene.set,group.by, split.by, facet.by)])
    }
    
  colnames(enriched) <- c(gene.set, group.by, split.by, facet.by)
  return(enriched)
}

.colorby <- function(enriched,
                     plot, 
                     color.by,
                     palette, 
                     type = "fill") { 
  if(type == "fill") {
    if(inherits(enriched[,color.by], "numeric")) {
      plot <- plot +
              scale_fill_gradientn(colors = .colorizer(palette, 11)) + 
              labs(fill = color.by) 
    } else {
      col <- length(unique(enriched[,color.by]))
      plot <- plot + 
              scale_fill_manual(values = .colorizer(palette, col)) +
              labs(fill = color.by) 
    }
  } else if (type == "color") {
    if(inherits(enriched[,color.by], "numeric")) {
      plot <- plot +
              scale_color_gradientn(colors = .colorizer(palette, 11)) + 
              labs(color = color.by) 
    } else {
      col <- length(unique(enriched[,color.by]))
      plot <- plot + 
              scale_color_manual(values = .colorizer(palette, col)) +
              labs(color = color.by) 
    }
  }
  return(plot)
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

#' @importFrom SummarizedExperiment assays
#' @importFrom MatrixGenerics rowSums2
.cntEval <- function(obj, 
                     assay = "RNA", 
                     type = "counts") {
  if (inherits(x = obj, what = "Seurat")) {
    cnts <- obj@assays[[assay]][type]
  } else if (inherits(x = obj, what = "SingleCellExperiment")) {
    pos <- ifelse(assay == "RNA", "counts", assay) 
    if(assay == "RNA") {
      cnts <- assay(obj,pos)
    } else {
      cnts <- assay(altExp(obj), pos)
    }
  } else {
    cnts <- obj
  }
  cnts <- cnts[rowSums2(cnts) != 0,]
  return(cnts)
}

#Add the values to single cell object
#' @importFrom SeuratObject CreateAssayObject
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SingleCellExperiment altExps
.adding.Enrich <- function(sc, enrichment, enrichment.name) {
  if (inherits(sc, "Seurat")) {
    new.assay <- suppressWarnings(CreateAssayObject(
                                  data = as.matrix(t(enrichment))))
    
    sc[[enrichment.name]] <- new.assay
  } else if (inherits(sc, "SingleCellExperiment")) {
    altExp(sc, enrichment.name) <- SummarizedExperiment(assay = t(enrichment))
    names(assays(altExp(sc, enrichment.name))) <- enrichment.name
  }
  return(sc)
}

#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment altExps
.pull.Enrich <- function(sc, enrichment.name) {
  if (inherits(sc, "Seurat")) {
    values <- t(sc[[enrichment.name]]["data"])
  } else if (inherits(sc, "SingleCellExperiment")) {
    values <- t(assay(altExps(sc)[[enrichment.name]], "data"))
  }
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