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
source(
  list.files(
    path = here(),
    pattern = "^zero_summary\\.R$",
    all.files = FALSE,
    full.names = TRUE,
    recursive = TRUE
  )
)

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
  # check the determinant: det(matrix(c(1, 0, 3, 4), 2, 2)) = 1*4 - 3*0 = 4
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
  # min and max should also ignore the NA thanks to na.rm = TRUE
  expect_equal(res$min_val, 0)
  expect_equal(res$max_val, 3)
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
# min_val and max_val are now computed regardless of shape, so datasum()
# no longer errors on a non-square matrix; only det_val stays NA, since
# the determinant is undefined for a non-square matrix.

test_that("datasum returns min/max but NA determinant for a non-square matrix", {
  
  # build a 2x3 matrix, which fails the is_square_matrix check
  m <- matrix(1:6, nrow = 2, ncol = 3)
  
  # the call should complete without raising an error
  out <- capture.output(res <- datasum(m))
  
  # dimensions should be reported as-is
  expect_equal(res$n_col, 3)
  expect_equal(res$n_row, 2)
  # min and max are still computable on a non-square matrix
  expect_equal(res$min_val, 1)
  expect_equal(res$max_val, 6)
  # the determinant is undefined here, so it should stay NA
  expect_true(is.na(res$det_val))
  # a 2x3 matrix of 1:6 contains no zeros
  expect_equal(res$total_zeros, 0)
  expect_equal(res$zero_rate, 0)
})

# -------- data.frame and tibble tests --------
#
# is.matrix() is FALSE for data.frames and tibbles, so det_val stays NA
# for them, but min_val/max_val are now computed on the tibble version
# of x regardless.

test_that("datasum returns min/max but NA determinant for a plain data.frame", {
  
  # build a simple square (2x2) data.frame, still not a "matrix" for is.matrix()
  df <- data.frame(a = 1:2, b = 3:4)
  
  # the call should complete without raising an error
  out <- capture.output(res <- datasum(df))
  
  expect_equal(res$n_col, 2)
  expect_equal(res$n_row, 2)
  expect_equal(res$min_val, 1)
  expect_equal(res$max_val, 4)
  expect_true(is.na(res$det_val))
  expect_equal(res$total_zeros, 0)
  expect_equal(res$zero_rate, 0)
})

test_that("datasum prints the tibble notice and returns min/max for a tibble input", {
  
  # build a tibble input, which should trigger the "already a tibble" message
  tib <- tibble::tibble(a = 1:2, b = 3:4)
  
  # capture the printed output and keep the invisible return value
  printed <- capture.output(res <- datasum(tib))
  
  # the "already a tibble" notice should have been printed
  expect_true(any(grepl("already a tibble", printed)))
  # min and max should still be correctly computed
  expect_equal(res$min_val, 1)
  expect_equal(res$max_val, 4)
  # the determinant is not defined for a tibble, so it stays NA
  expect_true(is.na(res$det_val))
})

# -------- list tests --------
#
# a plain list is coerced with as.data.frame() before being turned into
# a tibble, so datasum() should work on it as long as its elements are
# equal-length numeric vectors.

test_that("datasum works on a list of equal-length numeric vectors", {
  
  # build a list with two numeric vectors of the same length
  lst <- list(a = 1:5, b = 6:10)
  
  # the call should complete without raising an error
  out <- capture.output(res <- datasum(lst))
  
  # the list should be coerced into a 5-row, 2-column tibble
  expect_equal(res$n_col, 2)
  expect_equal(res$n_row, 5)
  # min and max across both vectors (1:5 and 6:10)
  expect_equal(res$min_val, 1)
  expect_equal(res$max_val, 10)
  # a list is never square in the is.matrix() sense, so det_val stays NA
  expect_true(is.na(res$det_val))
  # 1:5 and 6:10 contain no zeros
  expect_equal(res$total_zeros, 0)
  expect_equal(res$zero_rate, 0)
})

test_that("datasum counts zeros correctly on a list containing zero values", {
  
  # build a list where one vector contains a zero
  lst <- list(a = c(0, 1, 2), b = c(3, 4, 5))
  
  out <- capture.output(res <- datasum(lst))
  
  # only one zero across the whole coerced dataset
  expect_equal(res$total_zeros, 1)
  # 1 zero out of 6 cells is ~16.67%
  expect_equal(res$zero_rate, round(1 / 6 * 100, 2))
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