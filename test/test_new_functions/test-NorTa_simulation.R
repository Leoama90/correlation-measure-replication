# test-NorTa_simulation.R
#
# Purpose:
#   This script tests the fundamental features of norta_simulation(),
#   checking the dimensions of the generated correlation matrix and the
#   simulated dataset, the validity of the PSD-corrected correlation
#   matrix (non-negative eigenvalues, unit diagonal), and the
#   reproducibility of the output when a seed is provided.
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

# load NorTa_simulation() and its dependency (generate_matrix_with_zeroes)
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
  result <- norta_simulation(n = 10, n_zeroes = 4, N = 50, seed = 1)
  # the correlation matrix should be square, with side equal to n
  expect_equal(dim(result$R_true), c(10, 10))
})

# -------- test 2: R_psd is a valid PSD matrix --------
test_that("R_psd has no negative eigenvalues", {
  result <- norta_simulation(n = 10, n_zeroes = 4, N = 50, seed = 2)
  # PSD requires all eigenvalues to be >= 0 (allowing for tiny numerical
  # noise around zero)
  eigenvalues <- eigen(as.matrix(result$R_psd), only.values = TRUE)$values
  expect_true(all(eigenvalues >= -1e-8))
})

# -------- test 3: R_psd has a unit diagonal --------
test_that("R_psd has a diagonal of 1s", {
  result <- norta_simulation(n = 10, n_zeroes = 4, N = 50, seed = 3)
  # nearPD() was called with corr = TRUE, so the diagonal must stay at 1
  expect_equal(diag(as.matrix(result$R_psd)), rep(1, 10), tolerance = 1e-8)
})

# -------- test 4: sim_data has the requested dimensions --------
test_that("sim_data is an N x n matrix", {
  result <- norta_simulation(n = 10, n_zeroes = 4, N = 50, seed = 4)
  # N samples (rows) x n simulated taxa (columns)
  expect_equal(dim(result$sim_data), c(50, 10))
})

# -------- test 5: the same seed produces the same R_true --------
test_that("using the same seed gives reproducible results", {
  result_a <- norta_simulation(n = 10, n_zeroes = 4, N = 50, seed = 42)
  result_b <- norta_simulation(n = 10, n_zeroes = 4, N = 50, seed = 42)
  # identical seed and parameters should generate an identical correlation matrix
  expect_equal(result_a$R_true, result_b$R_true)
})
