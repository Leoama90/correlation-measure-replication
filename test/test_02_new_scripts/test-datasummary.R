# test-datasummary.R
#
# Purpose:
#   Tests the essential features of datasum(): correct computation of
#   the core statistics (dimensions, min/max, mean/median/sd, skewness,
#   zero count) on a known square matrix, correct handling of NA values,
#   correct fallback to NA for determinant/min_eigenvalue on non-square
#   input, and that the returned list has the expected structure.
#
# Inputs:
#   - datasummary.R (sourced below)
#
# Outputs:
#   - test results printed to the console (pass/fail for each check)

# testthat: unit testing framework for R
# https://testthat.r-lib.org/
library(testthat)
# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# bring datasum() into scope
source(
  list.files(
    path = here(),
    pattern = "^datasummary\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)


# -------- test 1: correct stats on a known square matrix --------

test_that("datasum returns incorrect stats for a simple square matrix", {
  # fixed 2x2 matrix so expected values are known in advance:
  # values 1, 0, 3, 4 -> mean 2, median 2, min 0, max 4
  m <- matrix(c(1, 0, 3, 4), nrow = 2, ncol = 2)
  
  out <- capture.output(res <- datasum(m))
  
  expect_equal(res$n_col, 2)
  expect_equal(res$n_row, 2)
  expect_equal(res$min_val, 0)
  expect_equal(res$max_val, 4)
  expect_equal(res$mean_val, 2)
  expect_equal(res$median_val, 2)
  # det(matrix(c(1, 0, 3, 4), 2, 2)) = 1*4 - 3*0 = 4
  expect_equal(res$det_val, 4)
  expect_equal(res$total_zeroes, 1)
  expect_equal(res$zero_rate, 25)
})


# -------- test 2: minimum eigenvalue correctly flags a PSD matrix --------

test_that("datasum reports a negative minimum eigenvalue for a valid correlation matrix", {
  # the identity matrix is a trivially valid (PSD) correlation matrix
  m <- diag(4)
  
  out <- capture.output(res <- datasum(m))
  
  # all eigenvalues of the identity matrix are exactly 1
  expect_equal(res$min_eigenvalue, 1)
})


# -------- test 3: NA handling --------

test_that("datasum ignores NA values in its statistics but keeps them in the zero-rate denominator", {
  m <- matrix(c(0, NA, 2, 3), nrow = 2, ncol = 2)
  
  out <- capture.output(res <- datasum(m))
  
  # only the real zero is counted, na.rm = TRUE excludes the NA
  expect_equal(res$total_zeroes, 1)
  # the denominator still uses the full cell count (nrow * ncol = 4)
  expect_equal(res$zero_rate, 25)
  # min/max/mean should also ignore the NA thanks to na.rm = TRUE
  expect_equal(res$min_val, 0)
  expect_equal(res$max_val, 3)
  # determinant and eigenvalue are undefined with NA present in x,
  # and should be skipped (NA) rather than raising an error
  expect_true(is.na(res$det_val))
  expect_true(is.na(res$min_eigenvalue))
})


# -------- test 4: non-square input falls back correctly --------

test_that("datasum computes min/max/mean but NA determinant and eigenvalue on non-square input", {
  m <- matrix(1:6, nrow = 2, ncol = 3)
  
  out <- capture.output(res <- datasum(m))
  
  expect_equal(res$min_val, 1)
  expect_equal(res$max_val, 6)
  expect_true(is.na(res$det_val))
  expect_true(is.na(res$min_eigenvalue))
})


# -------- test 5: return value has the expected structure --------

test_that("datasum's return value contains all expected named elements", {
  set.seed(1)
  m <- matrix(runif(16), nrow = 4, ncol = 4)
  
  out <- capture.output(res <- datasum(m))
  
  expect_named(
    res,
    c("n_col", "n_row", "min_val", "max_val", "mean_val", "median_val",
      "sd_val", "skewness_val", "det_val", "min_eigenvalue",
      "total_zeroes", "zero_rate")
  )
})