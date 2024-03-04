# test script for getGeneSets.R - testcases are NOT comprehensive!

test_that("getGeneSets works", {
  
  hallmark.default <- getGeneSets(library = "H")
  
  mouse.hallmark.default <- getGeneSets(library = "H", species = "Mus musculus")
  
  C5.GO.default <- getGeneSets(library = "C5", subcategory = "BP")
  
  
  expect_equal(hallmark.default , 
               getdata("getGeneSets", "getGeneSets_default"))
  
  expect_equal(mouse.hallmark.default, 
               getdata("getGeneSets", "getGeneSets_mouse"))
  
  expect_equal(C5.GO.default, 
               getdata("getGeneSets", "getGeneSets_C5"))
})
