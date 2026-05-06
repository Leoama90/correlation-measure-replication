library(testthat)
library(here)


test_that("I think I've found an error, sir!", {
  
  # creates a vector with the expected file name
  expected_data_frame <- c("compositional_effects_02.rds")
  
  # here() returns the root directory of the project
  # list.files() searches for .rds files recursively from the root
  # basename() strips the folder path, keeping only the file names
  found_data_frame <- basename(list.files(
    path = here(),
    pattern = "^compositional_effects_02\\.rds$",
    recursive = TRUE
  ))
  
  # checks that both vectors contain the same elements, regardless of order
  expect_setequal(found_data_frame, expected_data_frame)
})









