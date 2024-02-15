# test script for pcaEnrichment.R - testcases are NOT comprehensive!

test_that("pcaEnrichment works", {
  
  seuratObj <- getdata("performPCA", "pbmc_hallmarks")
  seuratObj <- performPCA(seuratObj,
                          assay = "escape.ssGSEA",
                          n.dim = 1:10)
  expect_doppelganger(
    "pcaEnrichment_plot",
    pcaEnrichment(seuratObj, 
                  dimRed = "escape.PCA",
                  x.axis = "PC1",
                  y.axis = "PC2")
  )
  
  expect_doppelganger(
    "pcaEnrichment_hex_plot",
    pcaEnrichment(seuratObj, 
                  dimRed = "escape.PCA",
                  x.axis = "PC1",
                  y.axis = "PC2",
                  style = "hex")
  )
  
  expect_doppelganger(
    "pcaEnrichment_addFactors_plot",
    pcaEnrichment(seuratObj,
                  dimRed = "escape.PCA",
                  x.axis = "PC2",
                  y.axis = "PC3",
                  display.factors = TRUE,
                  number.of.factors = 10)
  )
  
  expect_doppelganger(
    "pcaEnrichment_facetby_plot",
    pcaEnrichment(seuratObj, 
                  dimRed = "escape.PCA",
                  x.axis = "PC2",
                  y.axis = "PC1",
                  facet.by = "groups")
  )
  
  expect_doppelganger(
    "pcaEnrichment_facetby_addFactors_plot",
    pcaEnrichment(seuratObj, 
                  dimRed = "escape.PCA",
                  x.axis = "PC1",
                  y.axis = "PC2", 
                  facet.by = "groups", 
                  display.factors = TRUE,
                  number.of.factors = 10)
  )
  
})
