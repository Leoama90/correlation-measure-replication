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
  
  # check that found paths match the expected one
  expect_setequal(found_data_frame, expected_data_frame)
})

# --

