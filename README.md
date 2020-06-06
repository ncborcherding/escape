# escape
#### Easy single cell analysis platform for enrichment

<img align="right" src="https://github.com/ncborcherding/ncborcherding.github.io/blob/master/images/escape_hex_sticker.png" width="305" height="352">

### Introduction
Single-cell sequencing (SCS) is an emerging technology in the across the diverse array of biological fields. Part of the struggle with the high-resolution approach of SCS, is distilling the data down to meaningful scientific hypotheses. escape was created to bridge SCS results, either from raw counts or from the popular Seurat R package, with gene set enrichment analyses (GSEA), allowing users to simply and easily graph outputs. The package accesses the entire [Molecular Signature Database v7.0](https://www.gsea-msigdb.org/gsea/msigdb/search.jsp) and enables users to select single, multiple gene sets, and even libraries to perform enrichment analysis on. 

### R Packages Imported
+  dplyr
+  [facoextra](http://www.sthda.com/english/wiki/fviz-pca-quick-principal-component-analysis-data-visualization-r-software-and-data-mining) ```devtools::install_github("kassambara/factoextra")```
+  ggplots2
+  ggrepel
+  ggridges
+  grDevices
+  GSEABase
+  GSVA
+  msigdbr
+  limma
+  pheatmap
+  rlang
+  Seurat
+  SingleCellExperiment

### Installation

```devtools::install_github("ncborcherding/escape")```

### Learning To Use escape:

Vignette available [here](https://ncborcherding.github.io/vignettes/escape+vignette.html), includes 2,000 malignant and nonmalignant peripheral blood T cells from a patient with cutaenous T cell lymphoma.

### Contact
Questions, comments, suggestions, please feel free to contact Nick Borcherding via this repository, [email](mailto:ncborch@gmail.com), or using [twitter](https://twitter.com/theHumanBorch). 

