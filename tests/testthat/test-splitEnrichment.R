# test script for splitEnrichment.R - testcases are NOT comprehensive!

test_that("splitEnrichment works", {
  
  seuratObj <- getdata("runEscape", "pbmc_small_ssGSEA")
  
  expect_doppelganger(
    "splitEnrichment_default_plot",
    splitEnrichment(
      seuratObj, 
      split.by = "groups",
      assay = "escape",
      gene.set = "Tcells"
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
