% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualizations.R
\name{splitEnrichment}
\alias{splitEnrichment}
\title{Generate a split violin plot examine enrichment distributions}
\usage{
splitEnrichment(
  enriched,
  x.axis = NULL,
  scale.bracket = NULL,
  split = NULL,
  gene.set = NULL,
  colors = c("#0348A6", "#7AC5FF", "#C6FDEC", "#FFB433", "#FF4B20")
)
}
\arguments{
\item{enriched}{The output of \code{\link{enrichIt}}}

\item{x.axis}{Optional parameter for seperation.}

\item{scale.bracket}{This will filter the enrichment scores to remove 
extreme outliers. Values entered (1 or 2 numbers) will be the filtering 
parameter using z-scores of the selected gene.set. If only 1 value is given, 
a seocndary bracket is autommatically selected as the inverse of the number.}

\item{split}{The parameter to split, must be binary.}

\item{gene.set}{The gene set to graph on the y-axis.}

\item{colors}{The color palette for the ridge plot.}
}
\value{
ggplot2 object violin-based distributions of selected gene.set
}
\description{
This function allows to the user to examine the distribution of 
enrichment across groups by generating a split violin plot.
}
\examples{
ES2 <- readRDS(url(
"https://ncborcherding.github.io/vignettes/escape_enrichment_results.rds"))
splitEnrichment(ES2, x.axis = "cluster", split = "Type", 
gene.set = "HALLMARK_DNA_REPAIR")

}
\seealso{
\code{\link{enrichIt}} for generating enrichment scores.
}
