% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrich.R
\name{enrichIt}
\alias{enrichIt}
\title{Calculate gene set enrichment scores for single-cell data}
\usage{
enrichIt(obj, gene.sets = NULL, groups = 1000, cores = 2)
}
\arguments{
\item{obj}{The count matrix, Seurat, or SingleCellExperiment object.}

\item{gene.sets}{Gene sets from \code{\link{getGeneSets}} to use 
for the enrichment analysis. Alternatively a simple base R list where
the names of the list elements correspond to the name of the gene set
and the elements themselves are simple vectors of gene names representing
the gene set.}

\item{groups}{The number of cells to separate the enrichment calculation.}

\item{cores}{The number of cores to use for parallelization.}
}
\value{
Data frame of normalized enrichmenet scores (NES)
}
\description{
This function allows users to input both the single-cell RNA-sequencing 
counts and any gene set pathways either from the stored data or from 
other sources. The enrichment calculation itself 
uses the gsva R package and the poisson distribution for RNA.
}
\examples{
# download HALLMARK gene set collection
GS <- getGeneSets(library = "H") 
GS <- GS[c(1:2)] #Reduce list size for example
seurat_ex <- suppressWarnings(SeuratObject::pbmc_small)
ES <- enrichIt(obj = seurat_ex, gene.sets = GS)

# alternatively, construct your own list of gene sets
myGS <- list(Housekeeping = c("ACTA1", "ACTN1", "GAPDH"),
  Cancer = c("TP53","BRCA2","ERBB2","MYC"))
}
\seealso{
\code{\link{getGeneSets}} to collect gene sets.
}
\author{
Nick Borcherding, Jared Andrews
}
