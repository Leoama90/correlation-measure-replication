# test-datasum.R

library(testthat)
library(tidyverse)
library(here)

# load the script to test
source(list.files(path       = here(), 
                  pattern    = "^zero_summary\\.R$", 
                  all.files  = FALSE,
                  full.names = TRUE,
                  recursive  = TRUE
))


# suppress cat() output so tests stay clean, keep the return value
quiet_datasum <- function(x) {
  capture.output(result <- datasum(x))
  result
}

# -------- test 1: row/column counts --------
test_that("n_col and n_row are correct", {
  # 3 columns, 5 rows
  df <- tibble(a = 1:5, b = 6:10, c = 11:15)
  result <- quiet_datasum(df)
  
  # column count should match the tibble
  expect_equal(result$n_col, 3)
  # row count should match the tibble
  expect_equal(result$n_row, 5)
})

# -------- test 2: zero count --------
test_that("total_zeros counts zeros correctly", {
  # 4 zeros spread across two columns
  df <- tibble(a = c(0, 1, 0), b = c(2, 0, 0))
  result <- quiet_datasum(df)
  
  # total zero count should be 4
  expect_equal(result$total_zeros, 4)
})

# -------- test 3: percentage calculation --------
test_that("zero_rate is calculated correctly", {
  # 2 zeros out of 4 total cells = 50%
  df <- tibble(a = c(0, 1), b = c(0, 1))
  result <- quiet_datasum(df)
  
  # check the percentage matches the expected ratio
  expect_equal(result$zero_rate, 50)
})

# -------- test 4: data.frame input --------
test_that("datasum works with a data.frame, not just a tibble", {
  # plain data.frame instead of a tibble
  df <- data.frame(a = c(0, 1, 2), b = c(3, 0, 5))
  result <- quiet_datasum(df)
  
  # conversion should not affect the zero count
  expect_equal(result$total_zeros, 2)
})