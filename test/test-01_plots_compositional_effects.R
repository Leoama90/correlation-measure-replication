library(testthat)
library(here)

# -- test the that the preprocessed file is in the correct position ------

test_that("The file is not in the folder", {
  
  # expected full path of the preprocessed file
  expected_data_frame <- here("compositional_effects", "compositional_effects_02.rds")
  
  # search recursively and return full paths
  found_data_frame <- list.files(
    path       = here(),
    pattern    = "^compositional_effects_02\\.rds$",
    recursive  = TRUE,
    full.names = TRUE
  )
  
  # check that found path match the expected one
  expect_setequal(found_data_frame, expected_data_frame)
})

# -- test correct column name and size --------

test_that("there is a problem with the file columns: either one is empty or doesn't exist!", {
  # read the file
  df <- readRDS(here("compositional_effects", "compositional_effects_02.rds"))
  
  # required columns must be present
  expect_true(all(c("d", "pielou", "ERR_L1", "ERR_CLR") %in% names(df)))
  
  # must not be empty
  expect_gt(nrow(df), 0)
})

# -- check the data frame after the pipelines --------

test_that("there is a problem with the pipelined dataframe!",{
  # read the file
  df <- readRDS(here("compositional_effects", "compositional_effects_02.rds"))
  
  expect_true(nrow(df) == 1560)
  
})
