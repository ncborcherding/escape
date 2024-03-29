% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splitEnrichment.R
\name{splitEnrichment}
\alias{splitEnrichment}
\title{Visualize enrichment results with a split violin plot}
\usage{
splitEnrichment(
  input.data,
  assay = NULL,
  split.by = NULL,
  group.by = NULL,
  gene.set = NULL,
  order.by = NULL,
  facet.by = NULL,
  scale = TRUE,
  palette = "inferno"
)
}
\arguments{
\item{input.data}{Enrichment output from \code{\link{escape.matrix}} or
\code{\link{runEscape}}.}

\item{assay}{Name of the assay to plot if data is a single-cell object.}

\item{split.by}{Variable to form the split violin, must have 2 levels.}

\item{group.by}{Categorical parameter to plot along the x.axis. If input is
a single-cell object the default will be cluster.}

\item{gene.set}{Gene set to plot (on y-axis).}

\item{order.by}{Method to organize the x-axis - \strong{"mean"} will arrange
the x-axis by the mean of the gene.set, while \strong{"group"} will arrange
the x-axis by in alphanumerical order. Using \strong{NULL} will not reorder
the x-axis.}

\item{facet.by}{Variable to facet the plot into n distinct graphs.}

\item{scale}{Visualize raw values \strong{FALSE} or Z-transform 
enrichment values \strong{TRUE}.}

\item{palette}{Colors to use in visualization - input any 
\link[grDevices]{hcl.pals}.}
}
\value{
ggplot2 object violin-based distributions of selected gene.set
}
\description{
This function allows to the user to examine the distribution of 
enrichment across groups by generating a split violin plot.
}
\examples{
GS <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
           Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))
pbmc_small <- SeuratObject::pbmc_small
pbmc_small <- runEscape(pbmc_small, 
                        gene.sets = GS, 
                        min.size = NULL)
                        
splitEnrichment(pbmc_small, 
                assay = "escape",
                split.by = "groups",
                gene.set = "Tcells")

}
