# test script for scatterEnrichment.R - testcases are NOT comprehensive!

test_that("scatterEnrichment works", {
  
  seuratObj <- getdata("runEscape", "pbmc_small_ssGSEA")
  
  expect_doppelganger(
    "scatterEnrichment_default_plot",
    scatterEnrichment(
      seuratObj, 
      assay = "escape",
      x.axis = "Tcells",
      y.axis = "Bcells"
    )
  )

  expect_doppelganger(
    "scatterEnrichment_scale_plot",
    scatterEnrichment(
      seuratObj, 
      assay = "escape",
      x.axis = "Tcells",
      y.axis = "Bcells",
      scale = TRUE
    )
  )
  
  expect_doppelganger(
    "scatterEnrichment_facet_plot",
    scatterEnrichment(
      seuratObj, 
      assay = "escape",
      x.axis = "Tcells",
      y.axis = "Bcells",
      facet.by = "groups"
    )
  )
  
  
})
