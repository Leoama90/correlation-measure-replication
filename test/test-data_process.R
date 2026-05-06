library(testthat)
library(here)

# test that the .csv files found match the expected ones
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