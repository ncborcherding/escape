# test script for performPCA.R - testcases are NOT comprehensive!

test_that("performPCA works", {
  
  seuratObj <- getdata("performPCA", "pbmc_hallmarks")
  
  
  output <- performPCA(seuratObj@assays$escape.ssGSEA@data)
  
  expect_equal(names(output),
               c("PCA", "eigen_values", "contribution","rotation"))
  
  expect_equal(output$PCA, 
               getdata("performPCA", "performPCA_PCAvalues"))
  
})
