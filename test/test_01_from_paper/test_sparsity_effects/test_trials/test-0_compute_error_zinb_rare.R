# test-0_compute_error_zinb_rare.R
#
# Purpose:
#   This script tests important features from the 0_compute_error_zinb_rare.R script,
#   checking:
#       - ToyModel() (from ToyModel library) returns a correlation and a NorTa Matrixes
#         with correct dimensions
#
# Inputs:
#   - data generated inside this script for purpose test only
#
# Outputs:
#   - test results printed to the console (pass/fail for each check)
# Note: 
#     This test does not source the original script, because it was too much time consuming.
#     Instead, it tests the individual building blocks the original script relies on
#     (toy_model()) on small, self-contained data.
#     The trade-off is that this duplicates part of the original script's logic here;
#     if that logic changes in 0_compute_error_zinb_rare.R, this script must be updated to
#     match, or the tests will silently keep checking the old behaviour.


# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)

# ToyModel: lightweight package with necessary function to generate metagenomics data
# https://github.com/Fuschi/ToyModel
library(ToyModel)

# VGAM: vector generalized linear models
# https://cran.r-project.org/package=VGAM
library(VGAM)


# -------- Fixture: minimal params used across multiple tests --------

# Small ZINB parameter set representing a single OTU pair
test_params <- data.frame(
  munb  = c(10, 15),
  size  = c(2, 3),
  pstr0 = c(0.2, 0.4)
)


# -------- toy_model returns the expected structure --------

test_that("toy_model returns NorTA matrix and cor_normal with correct dimensions", {
  result <- toy_model(
    n = 50,
    cor = diag(2),
    M = 1,
    qdist = VGAM::qzinegbin,
    param = test_params
  )

  # NorTA should be a 50 x 2 non-negative integer matrix
  expect_equal(dim(result$NorTA), c(50, 2))
  expect_true(all(result$NorTA >= 0))

  # cor_normal should be a 2x2 symmetric matrix with ones on the diagonal
  expect_equal(dim(result$cor_normal), c(2, 2))
  expect_equal(diag(result$cor_normal), c(1, 1))
})
