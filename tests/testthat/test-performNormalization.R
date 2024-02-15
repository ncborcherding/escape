# test script for performNormalization.R - testcases are NOT comprehensive!

test_that("performNormalization works", {
  
  seuratObj <- getdata("performPCA", "pbmc_hallmarks")
  GS.hallmark <- getdata("performNormalization", "GS.Hallmark")
  
  
  seuratObj.p <- performNormalization(seuratObj, 
                                      assay = "escape.ssGSEA", 
                                      gene.sets = GS.hallmark,
                                      make.positive = TRUE)
  
  expect_equal(seuratObj.p@assays$escape.ssGSEA, 
               getdata("performNormalization", "performNormalization_positve"))
  
  seuratObj.n <- performNormalization(seuratObj, 
                                      assay = "escape.ssGSEA", 
                                      gene.sets = GS.hallmark,
                                      make.positive = FALSE)
  
  expect_equal(seuratObj.n@assays$escape.ssGSEA, 
               getdata("performNormalization", "performNormalization_nonpositive"))
  
})
