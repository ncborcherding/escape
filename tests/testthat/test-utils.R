# test script for runEscape.R - testcases are NOT comprehensive!

test_that(".orderFunction works", {
  
  enrichment <- as.data.frame(getdata("runEscape", "escape.matrix_ssGSEA"))
  enrichment$grouping <- c(rep("g2", 40), rep("g1", 40))
  enrichment <- enrichment[,c(1,3)]
  
  enrichment.order1 <- .orderFunction(enrichment, order.by = "mean", group.by = "grouping")
  
  enrichment.order2 <- .orderFunction(enrichment, order.by = "group.by", group.by = "grouping")
  
  expect_equal(enrichment.order1, 
               getdata("utils", "orderFunction_mean"))
  
  expect_equal(enrichment.order2, 
               getdata("utils", "orderFunction_group"))
})


test_that(".makeDFfromSCO works", {
  seuratObj <- getdata("runEscape", "pbmc_small_ssGSEA")
  
  enriched <- .makeDFfromSCO(seuratObj, 
                             assay = "escape", 
                             group.by = NULL,
                             split.by = "groups",
                             gene.set = "Tcells")
  
  expect_equal(enriched,
               getdata("utils", "makeDFfromSCO_data.frame"))
})