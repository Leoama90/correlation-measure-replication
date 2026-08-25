# test-biases_parameters_examples.R
#
# Purpose:
#   This script tests the fundamental features of biases_parameters_examples.R, 
#   checking:
#       - mean Pielou returns a valor between 0 and 1
#       - that a matrix full of ones gives mean_Pielou = 1
#       - if the dimension of the correlation matrixes are correct
#       - the length of ToyModel outputs
#
# Inputs:
#   - biases_parameters_examples.R
#
# Outputs:
#   - test results printed to the console (pass/fail for each check)

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)
# testthat: allows to write tests for the other scripts
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)
# vegan: community ecology analysis
# https://cran.r-project.org/package=vegan
library(vegan)


# -------- loads the main script to make its variables available in the test environment --------
source(
  list.files(
    path = here(),
    pattern = "^biases_parameters_examples\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)


# -------- generate a dummy matrix for testing purposes --------

# generates a 3x3 matrix with random numbers to test the function
set.seed(42)
test_matrix <- matrix(runif(9), nrow = 3, ncol = 3)


# -------- testing that the function gives a result between 0 and 1 --------

test_that("mean_Pielou returns a result strictly between 0 and 1", {
  expect_lt(mean_Pielou(test_matrix), 1)
  expect_gt(mean_Pielou(test_matrix), 0)
})


# -------- testing that a matrix full of 1s gives mean_Pielou = 1 --------

test_that("mean_Pielou returns 1 for a matrix full of 1s (maximum evenness)", {
  unitary_matrix <- matrix(rep(1, 9), ncol = 3, nrow = 3)
  expect_equal(mean_Pielou(unitary_matrix), 1)
})


# -------- test the dimensions of the correlation matrices --------

test_that("correlation matrices have the expected dimensions for each D", {
  d5_row <- nrow(cor_D5)
  d5_col <- ncol(cor_D5)
  d30_row <- nrow(cor_D30)
  d30_col <- ncol(cor_D30)
  d100_row <- nrow(cor_D100)
  d100_col <- ncol(cor_D100)
  
  expect_equal(d5_row, 5)
  expect_equal(d5_col, 5)
  expect_equal(d30_row, 30)
  expect_equal(d30_col, 30)
  expect_equal(d100_row, 100)
  expect_equal(d100_col, 100)
})


# -------- test the length of the ToyModel outputs --------

test_that("ToyModel outputs have the expected length", {
  expect_length(toy_D5_P100, 9)
  expect_length(toy_D5_P50,  9)
  expect_length(toy_D30_P50, 9)
})