# test script for ridgeEnrichment.R - testcases are NOT comprehensive!

test_that("ridgeEnrichment works", {
  
  seuratObj <- getdata("runEscape", "pbmc_small_ssGSEA")
  
  set.seed(42)
  expect_doppelganger(
    "ridgeEnrichment_default_plot",
    ridgeEnrichment(
      seuratObj, 
      assay = "escape",
      gene.set = "Bcells"
    )
  )
  
  set.seed(42)
  expect_doppelganger(
    "ridgeEnrichment_rugadded_plot",
    ridgeEnrichment(
      seuratObj, 
      assay = "escape",
      gene.set = "Bcells",
      add.rug = TRUE
    )
  )
  
  set.seed(42)
  expect_doppelganger(
    "ridgeEnrichment_facet_plot",
    ridgeEnrichment(
      seuratObj, 
      assay = "escape",
      gene.set = "Bcells",
      facet.by = "groups"
    )
  )

  set.seed(42)
  expect_doppelganger(
    "ridgeEnrichment_order_plot",
    ridgeEnrichment(
      seuratObj, 
      order.by = "mean",
      assay = "escape",
      gene.set = "Bcells"
    )
  )
  
  set.seed(42)
  expect_doppelganger(
    "ridgeEnrichment_gradient_plot",
    ridgeEnrichment(
      seuratObj, 
      assay = "escape",
      gene.set = "Bcells",
      color.by = "Bcells"
    )
  )
  
  set.seed(42)
  expect_doppelganger(
    "ridgeEnrichment_gradient_reorder_plot",
    ridgeEnrichment(
      seuratObj, 
      assay = "escape",
      order.by = "mean",
      gene.set = "Bcells",
      color.by = "Bcells"
    )
  )
  
  set.seed(42)
  expect_doppelganger(
    "ridgeEnrichment_gradient_facet_plot",
    ridgeEnrichment(
      seuratObj, 
      assay = "escape",
      gene.set = "Bcells",
      color.by = "Bcells",
      facet.by = "groups"
    )
  )
  
})
