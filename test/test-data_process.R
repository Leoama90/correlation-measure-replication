library(testthat)
library(here)


# -- test file existence in the folder ------
test_that("file not found, where did I put it...?",{
  # give the file position, here() searches it, then store in vector file_01
  file_01 <- here("data", "raw", "meta_HMP2.csv")
  # logic, if file.exists = TRUE, test pass
  expect_true(file.exists(file_01))
  
  # give the file position, here() searches it, then store in vector file_02
  file_02 <- here("data", "raw", "otu_HMP2_16S.csv")
  # logic, if file.exists = TRUE, test pass
  expect_true(file.exists(file_02))
  
  # give the file position, here() searches it, then store in vector file_03
  file_03 <- here("data", "raw", "taxonomy_HMP2_16S.csv")
  # logic, if file.exists = TRUE, test pass
  expect_true(file.exists(file_03))
})

# -- test that files are not empty ------
test_that("file size = 0, maybe it's empty and that's a problem!",{
  # test that the file is not empty (size > 0) 
  file_01 <- here("data", "raw", "meta_HMP2.csv")
  expect_true(file.size(file_01) > 0)
  file_02 <- here("data", "raw", "otu_HMP2_16S.csv")
  expect_true(file.size(file_02) > 0)
  file_03 <- here("data", "raw", "taxonomy_HMP2_16S.csv")
  expect_true(file.size(file_03) > 0)
})

# -- test that number of columns and rows are > 0
test_that("this number of column/row is greater than zero!", {
  tab_meta <- read.csv(here("data", "raw", "meta_HMP2.csv"), row.names = 1)
  expect_true(nrow(tab_meta) > 0)
  expect_true(ncol(tab_meta) > 0)
  
  tab_otu_ <- read.csv(here("data", "raw", "otu_HMP2_16S.csv"), row.names = 1)
  expect_true(nrow(tab_otu_) > 0)
  expect_true(ncol(tab_otu_) > 0)
  
  tab_taxo <- read.csv(here("data", "raw", "taxonomy_HMP2_16S.csv"), row.names = 1)
  expect_true(nrow(tab_taxo) > 0)
  expect_true(ncol(tab_taxo) > 0)
})

# -- test that the .csv files found match the expected ones ------
test_that(".csv files found are NOT the expected ones", {
  
  # creates a vector with expected file names
  expected_files <- c("meta_HMP2.csv", "otu_HMP2_16S.csv", "taxonomy_HMP2_16S.csv")
  
  # here() returns the root directory of the project
  # list.files() searches for .csv files recursively from the root
  # basename() strips folder path, keeping only the file names
  found_files <- basename(list.files(
    path = here(),
    pattern = "\\.csv$",
    recursive = TRUE
  ))
  
  # checks that both vectors contain the same elements, regardless of order
  expect_setequal(found_files, expected_files)
})