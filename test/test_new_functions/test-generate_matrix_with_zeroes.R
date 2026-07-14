# test_generate_matrix.R
#
# Purpose:
#   This script tests the fundamental features of generate_matrix_with_zeros.R,
#   checking matrix dimensions, correct number of zeros, diagonal handling,
#   and error handling for invalid inputs.
#
# Inputs:
#   - generate_matrix_with_zeros.R (sourced below)
#
# Outputs:
#   - test results printed to the console (pass/fail for each check)
#
# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)
# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# load the functions defined in the main script
source(here("generate_matrix_with_zeros.R"))

# -------- test 1: check that the output matrix has the requested dimensions --------
test_that("matrix has correct dimensions", {
  # generate a small matrix with a fixed seed for reproducibility
  m <- generate_matrix_with_zeros(n = 10, n_zeros = 5, seed = 1)
  # check that the number of rows matches n
  expect_equal(nrow(m), 10)
  # check that the number of columns matches n
  expect_equal(ncol(m), 10)
})

# -------- test 2: check that the matrix contains exactly the requested number of zeros --------
test_that("matrix contains the correct number of zeros", {
  # generate a matrix asking for 7 zeros
  m <- generate_matrix_with_zeros(n = 10, n_zeros = 7, seed = 2)
  # count how many entries are exactly zero
  n_zeros_found <- sum(m == 0)
  # check that the count matches what was requested
  expect_equal(n_zeros_found, 7)
})

# -------- test 3: check that the diagonal stays equal to 1 when keep_diag_one is TRUE --------
test_that("diagonal is preserved when keep_diag_one is TRUE", {
  # generate a matrix with diagonal protection enabled
  m <- generate_matrix_with_zeros(n = 10, n_zeros = 5, keep_diag_one = TRUE, seed = 3)
  # check that every diagonal value equals 1
  expect_true(all(diag(m) == 1))
})

# -------- test 4: check that values are normalized between 0 and 1 --------
test_that("matrix values are normalized between 0 and 1", {
  # generate a matrix without forcing any zero, to check normalization only
  m <- generate_matrix_with_zeros(n = 10, n_zeros = 0, keep_diag_one = FALSE, seed = 4)
  # check that no value exceeds 1
  expect_true(max(m) <= 1)
  # check that no value is below 0
  expect_true(min(m) >= 0)
})

# -------- test 5: check that an error is raised when n_zeros exceeds available positions --------
test_that("function errors when n_zeros is too large", {
  # request more zeros than positions available in a 3x3 matrix with protected diagonal
  expect_error(
    generate_matrix_with_zeros(n = 3, n_zeros = 100, keep_diag_one = TRUE)
  )
})