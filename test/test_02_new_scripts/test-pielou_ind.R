# test-pielou_ind.R
#
# Purpose:
#   This script tests the fundamental features of pielou_ind(): the
#   known-value cases (perfectly even vs. maximally uneven samples),
#   the 0*log(0) = 0 convention on zero-containing rows, the
#   per_sample = TRUE / FALSE behavior, and the three explicit input
#   validation errors (negative values, D < 2, zero library size).
#
# Inputs:
#   - pielou_ind.R (sourced below)
#
# Outputs:
#   - test results printed to the console (pass/fail for each check)

# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)
# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# load the script to test
source(
  list.files(
    path = here(),
    pattern = "^pielou_ind\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)


# -------- test 1: a perfectly even sample gives a Pielou index of 1 --------

test_that("a perfectly even row returns a Pielou index of 1", {
  # every taxon has the same abundance, so relative abundances are
  # uniform and Shannon entropy equals ln(D), giving Pielou = 1
  even_row <- matrix(c(5, 5, 5, 5), nrow = 1, ncol = 4)
  expect_equal(pielou_ind(even_row), 1, tolerance = 1e-8)
})


# -------- test 2: a maximally uneven sample gives a Pielou index of 0 --------

test_that("a row concentrated in a single taxon returns a Pielou index of 0", {
  # all abundance is in one taxon, so p = (1, 0, 0, 0); the 0*log(0) = 0
  # convention keeps entropy well-defined and equal to 0
  uneven_row <- matrix(c(10, 0, 0, 0), nrow = 1, ncol = 4)
  expect_equal(pielou_ind(uneven_row), 0, tolerance = 1e-8)
})


# -------- test 3: zeroes are handled without producing NA/NaN/Inf --------

test_that("rows containing some (but not all) zeroes produce a finite result", {
  # a realistic OTU-like row with a mix of zero and non-zero counts;
  # this should not trigger log(0) = -Inf issues thanks to the
  # explicit 0*log(0) = 0 convention inside pielou_ind()
  mixed_row <- matrix(c(5, 0, 3, 0, 2), nrow = 1, ncol = 5)
  result <- pielou_ind(mixed_row)
  expect_true(is.finite(result))
})


# -------- test 4: per_sample = TRUE returns one value per row --------

test_that("per_sample = TRUE returns a vector of length nrow(x)", {
  # three samples with different distributions; per_sample should give
  # one Pielou value per row, and the default (mean) should equal the
  # average of that vector
  x <- matrix(
    c(5, 5, 5, 5,
      10, 0, 0, 0,
      1, 2, 3, 4),
    nrow = 3, byrow = TRUE
  )
  per_sample_result <- pielou_ind(x, per_sample = TRUE)
  mean_result <- pielou_ind(x, per_sample = FALSE)
  expect_length(per_sample_result, nrow(x))
  expect_equal(mean(per_sample_result), mean_result, tolerance = 1e-8)
})


# -------- test 5: negative values are rejected --------

test_that("pielou_ind errors on negative values", {
  # abundances/counts cannot be negative; this must fail explicitly
  # rather than returning a nonsensical value
  x <- matrix(c(5, -1, 3, 2), nrow = 1, ncol = 4)
  expect_error(pielou_ind(x))
})


# -------- test 6: a single-column input is rejected --------

test_that("pielou_ind errors when there are fewer than 2 columns", {
  # with D = 1, ln(D) = 0, which would make the normalization divide
  # by zero; this must fail explicitly instead of returning NaN/Inf
  x <- matrix(c(5, 3, 2), nrow = 3, ncol = 1)
  expect_error(pielou_ind(x))
})


# -------- test 7: an all-zero row (library size 0) is rejected --------

test_that("pielou_ind errors on a row with library size 0", {
  # a row summing to zero has undefined relative abundances (division
  # by zero); this must fail explicitly rather than returning NaN
  x <- matrix(c(0, 0, 0, 0), nrow = 1, ncol = 4)
  expect_error(pielou_ind(x))
})