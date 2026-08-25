# test-CLR.R
#
# Purpose:
#   This script tests the fundamental features of CLR.R, checking:
#       - CLR() errors on non-matrix input
#       - CLR() runs without error on a valid positive numeric matrix
#       - CLR() preserves the input matrix dimensions
#       - CLR() output is numeric (double)
#       - each row of the CLR-transformed matrix sums to (approximately)
#         zero, as required by the CLR definition
#
# Inputs:
#   - CLR.R (sourced below)
#
# Outputs:
#   - test results printed to the console (pass/fail for each check)

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)
# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)


# -------- load the function to test --------

source(
  list.files(
    path = here(),
    pattern = "^CLR\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)


# -------- generate dummy variables for testing purposes --------

x <- 5
z <- matrix(1:25,
  nrow = 5,
  ncol = 5
)


# -------- test CLR() input validation --------

test_that("CLR() errors on non-matrix input", {
  expect_error(CLR(x))
})

test_that("CLR() runs without error on a valid positive numeric matrix", {
  expect_no_error(CLR(z))
})


# -------- test CLR() output structure --------

test_that("CLR() preserves the input matrix dimensions", {
  expect_equal(dim(CLR(z)), dim(z))
})

test_that("CLR() output is numeric (double)", {
  expect_true(is.double(CLR(z)))
})


# -------- test CLR() mathematical property --------

test_that("each row of the CLR-transformed matrix sums to (approximately) zero", {
  # by definition, CLR centers each row on its own mean log-value, so
  # rows sum to 0 up to floating point error
  row_sums <- rowSums(CLR(z))
  expect_true(all(abs(row_sums) < 1e-8))
})
