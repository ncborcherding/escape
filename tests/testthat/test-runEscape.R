# test script for runEscape.R - testcases are NOT comprehensive!

test_that("runEscape works", {
  GS <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
             Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))
             
  pbmc_small <- SeuratObject::pbmc_small
  pbmc_sce <- Seurat::as.SingleCellExperiment(pbmc_small)
  
  ####################
  #Testing the methods
  ####################
  trial.ssGSEA <- escape.matrix(pbmc_small, 
                                method = "ssGSEA", 
                                gene.sets = GS, 
                                min.size = NULL)
  
  trial.GSVA <- escape.matrix(pbmc_small, 
                              method = "GSVA", 
                              gene.sets = GS, 
                              min.size = NULL)
  
  trial.UCell <- escape.matrix(pbmc_small, 
                               method = "UCell", 
                               gene.sets = GS, 
                               min.size = NULL)
  
  set.seed(123)
  trial.AUCell <- escape.matrix(pbmc_small, 
                                method = "AUCell", 
                                gene.sets = GS, 
                                min.size = NULL)
  
  expect_equal(trial.ssGSEA, 
               getdata("runEscape", "escape.matrix_ssGSEA"))
  expect_equal(trial.GSVA, 
               getdata("runEscape", "escape.matrix_GSVA"))
  expect_equal(trial.UCell, 
               getdata("runEscape", "escape.matrix_UCell"))
  expect_equal(trial.AUCell, 
               getdata("runEscape", "escape.matrix_AUCell"), 
               tolerance=1e-4)
  
  pbmc_small <- runEscape(pbmc_small,
                          method = "ssGSEA", 
                          gene.sets = GS, 
                          min.size = NULL)
  
  expect_equal(names(pbmc_small@assays),
               c("RNA", "escape"))
  
  expect_equal(t(pbmc_small@assays$escape@data),
               getdata("runEscape", "escape.matrix_ssGSEA"))
  
  
})
