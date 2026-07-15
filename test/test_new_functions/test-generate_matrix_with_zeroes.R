# test_generate_matrix.R
#
# Purpose:
#   This script tests the fundamental features of generate_matrix_with_zeros.R,
#   checking matrix dimensions, symmetry, correct number of zeros, diagonal
#   handling, value range, and error handling for invalid inputs.
#
# Inputs:
#   - generate_matrix.R (sourced below)
#
# Outputs:
#   - test results printed to the console (pass/fail for each check)

# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)
# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# load the functions defined in the main script
source(list.files(
       path       = here(),
       pattern    = "^generate_matrix_with_zeroes\\.R$",
       full.names = TRUE,
       recursive  = TRUE
       ))

# -------- test 1: check that the output matrix has the requested dimensions --------
test_that("matrix has correct dimensions", {
  # generate a small matrix with a fixed seed for reproducibility
  m <- generate_matrix_with_zeros(n = 10, n_zeros = 4, seed = 1)
  # check that the number of rows matches n
  expect_equal(nrow(m), 10)
  # check that the number of columns matches n
  expect_equal(ncol(m), 10)
})

# -------- test 2: check that the matrix contains exactly the requested number of zeros --------
test_that("matrix contains the correct number of zeros", {
  # generate a matrix asking for 6 zeros (must be even, since zeros are symmetric)
  m <- generate_matrix_with_zeros(n = 10, n_zeros = 6, seed = 2)
  # count how many entries are exactly zero
  n_zeros_found <- sum(m == 0)
  # check that the count matches what was requested
  expect_equal(n_zeros_found, 6)
})

# -------- test 3: check that the diagonal is always equal to 1 --------
test_that("diagonal is always set to 1", {
  # generate a matrix with the default settings
  m <- generate_matrix_with_zeros(n = 10, n_zeros = 4, seed = 3)
  # check that every diagonal value equals 1
  expect_true(all(diag(m) == 1))
})

# -------- test 4: check that values are in the [-1, 1] range --------
test_that("matrix values are within [-1, 1]", {
  # generate a matrix without forcing any zero, to check the value range only
  m <- generate_matrix_with_zeros(n = 10, n_zeros = 0, seed = 4)
  # check that no value exceeds 1
  expect_true(max(m) <= 1)
  # check that no value is below -1
  expect_true(min(m) >= -1)
})

# -------- test 5: check that the matrix is symmetric --------
test_that("matrix is symmetric", {
  # generate a matrix with a fixed seed for reproducibility
  m <- generate_matrix_with_zeros(n = 10, n_zeros = 4, seed = 5)
  # check that the matrix equals its own transpose
  expect_equal(m, t(m))
})

# -------- test 6: check that an error is raised when n_zeros exceeds available positions --------
test_that("function errors when n_zeros is too large", {
  # request more zeros than off-diagonal positions available in a 3x3 matrix
  expect_error(
    generate_matrix_with_zeros(n = 3, n_zeros = 100)
  )
})

# -------- test 7: check that an error is raised when n_zeros is odd --------
test_that("function errors when n_zeros is odd", {
  # request an odd number of zeros, which cannot be placed symmetrically
  expect_error(
    generate_matrix_with_zeros(n = 10, n_zeros = 5)
  )
})