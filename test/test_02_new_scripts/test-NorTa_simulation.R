# test-NorTa_simulation.R
#
# Purpose:
#   This script tests the fundamental features of norta_simulation(),
#   checking the dimensions of the generated correlation matrix and the
#   simulated dataset, the validity of the correlation matrix (PSD by
#   construction, unit diagonal), and the reproducibility of the output
#   when a seed is provided.
#
# Inputs:
#   - NorTa_simulation.R (sourced below)
#
# Outputs:
#   - test results printed to the console (pass/fail for each check)

# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)
# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# load norta_simulation() and its dependency (generate_matrix_factors)
source(
  list.files(
    path = here(),
    pattern = "^NorTa_simulation\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)


# -------- test 1: R_true has the requested dimensions --------

test_that("R_true is an n x n matrix", {
  result <- norta_simulation(n = 10, n_groups = 3, N = 50, seed = 1)
  expect_equal(dim(result$R_true), c(10, 10))
})


# -------- test 2: R_true is a valid PSD matrix --------

test_that("R_true has no negative eigenvalues", {
  result <- norta_simulation(n = 10, n_groups = 3, N = 50, seed = 2)
  eigenvalues <- eigen(result$R_true, only.values = TRUE)$values
  expect_true(all(eigenvalues >= -1e-8))
})


# -------- test 3: R_true has a unit diagonal --------

test_that("R_true has a diagonal of 1s", {
  result <- norta_simulation(n = 10, n_groups = 3, N = 50, seed = 3)
  expect_equal(diag(result$R_true), rep(1, 10), tolerance = 1e-8)
})


# -------- test 4: sim_data has the requested dimensions --------

test_that("sim_data is an N x n matrix", {
  result <- norta_simulation(n = 10, n_groups = 3, N = 50, seed = 4)
  expect_equal(dim(result$sim_data), c(50, 10))
})


# -------- test 5: the same seed produces the same R_true --------

test_that("using the same seed gives reproducible results", {
  result_a <- norta_simulation(n = 10, n_groups = 3, N = 50, seed = 42)
  result_b <- norta_simulation(n = 10, n_groups = 3, N = 50, seed = 42)
  expect_equal(result_a$R_true, result_b$R_true)
})