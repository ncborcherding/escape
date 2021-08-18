# escape
#### Easy single cell analysis platform for enrichment

<img align="right" src="https://github.com/ncborcherding/ncborcherding.github.io/blob/master/images/escape_hex_sticker.png" width="305" height="352">

### Introduction
Single-cell sequencing (SCS) is an emerging technology in the across the diverse array of biological fields. Part of the struggle with the high-resolution approach of SCS, is distilling the data down to meaningful scientific hypotheses. escape was created to bridge SCS results, either from raw counts or from the popular Seurat R package, with gene set enrichment analyses (GSEA), allowing users to simply and easily graph outputs. The package accesses the entire [Molecular Signature Database v7.0](https://www.gsea-msigdb.org/gsea/msigdb/search.jsp) and enables users to select single, multiple gene sets, and even libraries to perform enrichment analysis on. 

### Installation

```devtools::install_github("ncborcherding/escape")```

#### Developmental Version

```devtools::install_github("ncborcherding/escape@dev")```

A newer version of escape is in the works to address previous issues and add new functions. For the most up-to-date version, check out the dev branch. These changes include adding support of [UCell](https://github.com/carmonalab/UCell) and [singscore](https://github.com/DavisLaboratory/singscore) to the ```enrichIt()``` function. 

### Learning To Use escape:

Vignette available [here](https://ncborcherding.github.io/vignettes/escape_vignette.html), includes 2,000 malignant and nonmalignant peripheral blood T cells from a patient with cutaenous T cell lymphoma.

### Citation 
If using escape, please cite the [article](https://www.nature.com/articles/s42003-020-01625-6): Borcherding, N., Vishwakarma, A., Voigt, A.P. et al. Mapping the immune environment in clear cell renal carcinoma by single-cell genomics. Commun Biol 4, 122 (2021). https://doi.org/10.1038/s42003-020-01625-6. 


### Contact
Questions, comments, suggestions, please feel free to contact Nick Borcherding via this repository, [email](mailto:ncborch@gmail.com), or using [twitter](https://twitter.com/theHumanBorch). 

