library(testthat)
library(here)
library(tidyverse)

source(here("scripts", "csv_analysis.R"))


dum_tibble <- tibble(       # generate a dummy tibble for this test
  `letters` = letters[1:5], # first column letters only
  `numbers` = sample(1:5))  # second column integers from 1 to 5

test_that("There is an error in the correct number of rows and columns:", {

  expect_equal(nrow(dum_tibble), 5) # test correct number of rows
  expect_equal(ncol(dum_tibble), 2) # test correct number of columns
})

test_that("The row or column names doesn't match",{
  
  expected_names <- c("letters", "numbers")
  expect_identical(names(dum_tibble), expected_names)
})



