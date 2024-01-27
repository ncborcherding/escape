# test script for ridgeEnrichment.R - testcases are NOT comprehensive!

test_that("ridgeEnrichment works", {
  
  seuratObj <- getdata("runEscape", "pbmc_small_ssGSEA")
  
  expect_doppelganger(
    "ridgeEnrichment_default_plot",
    ridgeEnrichment(
      seuratObj, 
      assay = "escape",
      gene.set = "Bcells"
    )
  )
  
  expect_doppelganger(
    "ridgeEnrichment_default_plot",
    ridgeEnrichment(
      seuratObj, 
      assay = "escape",
      gene.set = "Bcells",
      add.rug = TRUE
    )
  )

  expect_doppelganger(
    "splitEnrichment_mean_plot",
    splitEnrichment(
      seuratObj, 
      order.by = "mean",
      split.by = "groups",
      assay = "escape",
      gene.set = "Tcells"
    )
  )
  
  
})
