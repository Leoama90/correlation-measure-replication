# test-pseudocount.R
#
# Purpose:
#   This script tests the fundamental features of pseudocount(), checking
#   that library sizes of 0 trigger an error, that no zeroes remain in
#   the output, that non-zero proportions are left unchanged, that the
#   replaced pseudocount values match the expected formula, and that
#   output dimensions match the input.
#
# Inputs:
#   - pseudocount.R (sourced below)
#
# Outputs:
#   - test results printed to the console (pass/fail for each check)
#
# Note:
#   pseudocount() now accepts an explicit threshold_pct argument (see
#   pseudocount.R): when it is provided, the readline() prompt is skipped
#   entirely, so the real function can be called directly and
#   non-interactively in these tests. The previous version of this test
#   duplicated the function's body to avoid an infinite loop at the
#   readline() prompt when sourced; that duplicate is no longer needed
#   and has been removed, so these tests now exercise the actual
#   pseudocount() function, not a parallel copy of it.

# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)
# tidyverse: data science toolkit (dplyr, ggplot2, tidyr, etc.)
# https://tidyverse.org
library(tidyverse)
# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# load pseudocount.R: this brings the real pseudocount() function into
# scope; safe to source here, since sourcing only defines the function
# and does not call it (the readline() prompt only fires when
# pseudocount() is called without threshold_pct)
source(
  list.files(
    path = here(),
    pattern = "^pseudocount\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)

# -------- shared test data --------
# 2 samples x 3 OTUs, with one zero per row, in a different column each
# time, and library sizes of 15 and 5 respectively. Explicit column
# names avoid the as_tibble() "must have unique column names" warning
# that would otherwise appear when converting an unnamed matrix.
otu_test <- matrix(
  c(10, 0, 5,
    2, 3, 0),
  nrow = 2, byrow = TRUE,
  dimnames = list(NULL, c("OTU1", "OTU2", "OTU3"))
)

# -------- test 1: a library size of 0 triggers an error --------
test_that("pseudocount() errors when a row has a library size of 0", {
  # a row of all zeroes has a library size of 0
  otu_empty_row <- matrix(
    c(0, 0, 0,
      2, 3, 4),
    nrow = 2, byrow = TRUE,
    dimnames = list(NULL, c("OTU1", "OTU2", "OTU3"))
  )
  expect_error(pseudocount(otu_empty_row, threshold_pct = 0.5))
})

# -------- test 2: an out-of-range threshold_pct triggers an error --------
test_that("pseudocount() errors when threshold_pct is outside [0, 1]", {
  # passing threshold_pct explicitly must still be validated, exactly as
  # the interactive prompt would validate a typed value
  expect_error(pseudocount(otu_test, threshold_pct = 1.5))
  expect_error(pseudocount(otu_test, threshold_pct = -0.2))
})

# -------- test 3: no zeroes remain in the output --------
test_that("no zero values remain after pseudocount replacement", {
  result <- pseudocount(otu_test, threshold_pct = 0.5)
  expect_true(all(result != 0))
})

# -------- test 4: output dimensions match the input --------
test_that("output has the same dimensions as the input", {
  result <- pseudocount(otu_test, threshold_pct = 0.5)
  expect_equal(dim(result), dim(otu_test))
})

# -------- test 5: non-zero entries are left unchanged as proportions --------
test_that("non-zero entries are correctly converted to proportions", {
  result <- pseudocount(otu_test, threshold_pct = 0.5)
  # row 1: library size 15, so the first entry should be 10 / 15
  expect_equal(result[[1, 1]], 10 / 15, tolerance = 1e-8)
  # row 2: library size 5, so the second entry should be 3 / 5
  expect_equal(result[[2, 2]], 3 / 5, tolerance = 1e-8)
})

# -------- test 6: pseudocount values match the expected formula --------
test_that("replaced zeroes match threshold_pct * (1 / lib_size)", {
  result <- pseudocount(otu_test, threshold_pct = 0.5)
  # row 1 (lib_size = 15): pseudocount = 0.5 * (1 / 15)
  expect_equal(result[[1, 2]], 0.5 * (1 / 15), tolerance = 1e-8)
  # row 2 (lib_size = 5): pseudocount = 0.5 * (1 / 5)
  expect_equal(result[[2, 3]], 0.5 * (1 / 5), tolerance = 1e-8)
})