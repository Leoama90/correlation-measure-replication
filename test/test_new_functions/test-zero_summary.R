# testthat: unit testing framework for R
# https://testthat.r-lib.org/
library(testthat)
# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)
# tidyverse: data science toolkit (dplyr, ggplot2, tidyr, etc.)
# https://tidyverse.org
library(tidyverse)

# bring datasum() into scope by sourcing the file where it is defined
# (adjust the path below if datasum() lives somewhere else in your project)
source(list.files(path       = here(), 
                  pattern    = "^zero_summary\\.R$", 
                  all.files  = FALSE,
                  full.names = TRUE,
                  recursive  = TRUE
))

# -------- square matrix tests --------

test_that("datasum returns correct stats for a simple square matrix", {
  
  # build a fixed 2x2 matrix so the expected values are known in advance
  m <- matrix(c(1, 0, 3, 4), nrow = 2, ncol = 2)
  
  # run datasum(), sending its cat() output to capture.output so it
  # does not clutter the test report, and keep the invisible return value
  out <- capture.output(res <- datasum(m))
  
  # check the reported number of columns matches the matrix
  expect_equal(res$n_col, 2)
  # check the reported number of rows matches the matrix
  expect_equal(res$n_row, 2)
  # check the minimum value found (the matrix contains 0, 1, 3, 4)
  expect_equal(res$min_val, 0)
  # check the maximum value found (the matrix contains 0, 1, 3, 4)
  expect_equal(res$max_val, 4)
  # check the determinant: det(matrix(c(1,0,3,4), 2, 2)) = 1*4 - 3*0 = 4
  expect_equal(res$det_val, 4)
  # check the total number of zero entries (only one zero in the matrix)
  expect_equal(res$total_zeros, 1)
  # check the zero percentage: 1 zero out of 4 cells is 25%
  expect_equal(res$zero_rate, 25)
})

test_that("datasum prints the expected summary lines for a square matrix", {
  
  # reuse the same small matrix as in the previous test
  m <- matrix(c(1, 0, 3, 4), nrow = 2, ncol = 2)
  
  # capture everything the function prints with cat() during the call
  printed <- capture.output(datasum(m))
  
  # the header line should contain the "---" markers around the dataset name
  expect_true(any(grepl("---", printed)))
  # the column count line should mention "number of columns:"
  expect_true(any(grepl("number of columns:", printed)))
  # the zero percentage line should report "25 %"
  expect_true(any(grepl("25 %", printed)))
})

test_that("datasum's return value is invisible", {
  
  # build the smallest valid square numeric matrix, a single cell
  m <- matrix(5, nrow = 1, ncol = 1)
  
  # expect_invisible checks that the call does not auto-print at top level
  expect_invisible(datasum(m))
})

test_that("datasum truncates the variable name to its first four letters", {
  
  # use a variable name longer than four characters on purpose
  mymatrix <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)
  
  # capture the printed header, which embeds the truncated dataset_name
  printed <- capture.output(datasum(mymatrix))
  
  # the header should contain "myma", the first four letters of "mymatrix"
  expect_true(any(grepl("myma", printed)))
})

test_that("datasum ignores NA values when counting zeros but keeps them in the denominator", {
  
  # build a square matrix that contains one NA and one real zero
  m <- matrix(c(0, NA, 2, 3), nrow = 2, ncol = 2)
  
  # run datasum() while discarding its printed output
  out <- capture.output(res <- datasum(m))
  
  # only the real zero should be counted, na.rm = TRUE excludes the NA
  expect_equal(res$total_zeros, 1)
  # the denominator still uses the full cell count (nrow * ncol = 4)
  expect_equal(res$zero_rate, 25)
})

test_that("datasum works on the documentation example without error", {
  
  # fix the random seed so the random matrix is reproducible
  set.seed(1)
  # build the same kind of input used in the function's @examples block
  m <- matrix(runif(16), nrow = 4, ncol = 4)
  
  # the call should complete without raising an error
  expect_no_error(capture.output(res <- datasum(m)))
  
  # the returned list should always contain these seven named elements
  expect_named(
    res,
    c("n_col", "n_row", "min_val", "max_val", "det_val", "total_zeros", "zero_rate")
  )
})

# -------- non-square matrix tests --------
#
# min_val, max_val and det_val are only assigned inside the
# is_square_matrix branch, so datasum() is expected to error out later
# when it tries to cat() them for anything that is not a square matrix.
# these tests document that current behaviour of the function.

test_that("datasum errors on a non-square matrix because max_val is undefined", {
  
  # build a 2x3 matrix, which fails the is_square_matrix check
  m <- matrix(1:6, nrow = 2, ncol = 3)
  
  # the function should stop with an "object not found" style error
  expect_error(capture.output(datasum(m)), regexp = "max_val")
})

# -------- data.frame and tibble tests --------
#
# is.matrix() is FALSE for data.frames and tibbles, so the matrix-only
# branch never runs and min_val, max_val, det_val stay undefined for them.

test_that("datasum errors on a plain data.frame input", {
  
  # build a simple square (2x2) data.frame, still not a "matrix" for is.matrix()
  df <- data.frame(a = 1:2, b = 3:4)
  
  # expect an error since is_square_matrix is FALSE for data.frames
  expect_error(capture.output(datasum(df)), regexp = "max_val")
})

test_that("datasum prints the tibble notice and then errors on a tibble input", {
  
  # build a tibble input, which should trigger the "already a tibble" message
  tib <- tibble::tibble(a = 1:2, b = 3:4)
  
  # wrap the call in tryCatch so the error message can be inspected below
  error_message <- character(0)
  tryCatch({
    capture.output(datasum(tib))
  }, error = function(e) {
    # store the error message for the assertion outside tryCatch
    error_message <<- conditionMessage(e)
  })
  
  # the error message should reference the missing max_val object
  expect_true(grepl("max_val", error_message))
})

# -------- edge case: an all-zero square matrix --------

test_that("datasum reports 100% zero rate for an all-zero square matrix", {
  
  # build a 3x3 matrix filled entirely with zeros
  m <- matrix(0, nrow = 3, ncol = 3)
  
  # run the function, discarding its printed output
  out <- capture.output(res <- datasum(m))
  
  # every one of the 9 cells is zero
  expect_equal(res$total_zeros, 9)
  # the zero rate should therefore be exactly 100%
  expect_equal(res$zero_rate, 100)
  # the minimum and maximum of an all-zero matrix are both zero
  expect_equal(res$min_val, 0)
  expect_equal(res$max_val, 0)
  # the determinant of an all-zero matrix is zero
  expect_equal(res$det_val, 0)
})