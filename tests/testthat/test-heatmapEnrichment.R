# test script for heatmapEnrichment.R - testcases are NOT comprehensive!

test_that("heatmapEnrichment works", {
  
  seuratObj <- getdata("runEscape", "pbmc_small_ssGSEA")
  
  expect_doppelganger(
    "heatmapEnrichment_default_plot",
    heatmapEnrichment(
      seuratObj, 
      assay = "escape")
  )
  
  expect_doppelganger(
    "heatmapEnrichment_scale_plot",
    heatmapEnrichment(
      seuratObj, 
      assay = "escape",
      scale = TRUE
    )
  )
  
  
  expect_doppelganger(
    "heatmapEnrichment_facet_plot",
    heatmapEnrichment(
      seuratObj, 
      assay = "escape",
      facet.by = "groups"
    )
  )

  expect_doppelganger(
    "heatmapEnrichment_clusterRows_plot",
    heatmapEnrichment(
      seuratObj, 
      cluster.rows = TRUE,
      assay = "escape",
    )
  )
  
  expect_doppelganger(
    "heatmapEnrichment_clusterColumns_plot",
    heatmapEnrichment(
      seuratObj, 
      cluster.columns = TRUE,
      assay = "escape",
    )
  )
  
  expect_doppelganger(
    "geyserEnrichment_gradient_plot",
    geyserEnrichment(
      seuratObj, 
      assay = "escape",
      gene.set = "Tcells",
      color.by = "Tcells"
    )
  )
  
})
