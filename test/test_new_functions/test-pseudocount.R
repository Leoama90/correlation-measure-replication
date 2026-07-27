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
#   pseudocount() uses a readline() command which goes into an infinite loop
#   once the script is sourced in the test. To avoid this, I duplicated the
#   body of the function in the test without the problematic command.
#   the trade-off is that this test must be changed accordignly to the original
#   script, if the script changes in the future.

# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)
# tidyverse: data science toolkit (dplyr, ggplot2, tidyr, etc.)
# https://tidyverse.org
library(tidyverse)
# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# load pseudocount.R only for reference/documentation purposes; the
# actual function used in the tests below is the non-interactive
# duplicate defined right after
source(
  list.files(
    path = here(),
    pattern = "^pseudocount\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)

# -------- non-interactive duplicate of pseudocount(), for testing --------
# same logic as pseudocount.R, but threshold_pct is a direct argument
# instead of being requested via readline()
pseudocount_for_test <- function(x, threshold_pct) {
  # Convert the input to a tibble.
  x <- as_tibble(x)
  # Compute the total counts for each row.
  lib_size <- rowSums(x)
  # Stop early if any sample has a library size of 0.
  if (any(lib_size == 0)) {
    stop("At least one row has a library size of 0; remove empty samples before computing pseudocounts.")
  }
  # Compute the detection limit for each row.
  detection_limit <- 1 / lib_size
  # Compute the row-specific pseudocount values.
  pseudo <- threshold_pct * detection_limit
  # Convert raw counts to row-wise proportions.
  y_prop <- sweep(as.matrix(x), 1, lib_size, "/")
  # Replace zeroes with the corresponding row-specific pseudocount.
  zero_mask <- y_prop == 0
  y_prop[zero_mask] <- rep(pseudo, times = ncol(y_prop))[zero_mask]
  # Return the transformed data as a tibble.
  as_tibble(y_prop)
}

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
test_that("pseudocount_for_test() errors when a row has a library size of 0", {
  # a row of all zeroes has a library size of 0
  otu_empty_row <- matrix(
    c(0, 0, 0,
      2, 3, 4),
    nrow = 2, byrow = TRUE,
    dimnames = list(NULL, c("OTU1", "OTU2", "OTU3"))
  )
  expect_error(pseudocount_for_test(otu_empty_row, threshold_pct = 0.5))
})

# -------- test 2: no zeroes remain in the output --------
test_that("no zero values remain after pseudocount replacement", {
  result <- pseudocount_for_test(otu_test, threshold_pct = 0.5)
  expect_true(all(result != 0))
})

# -------- test 3: output dimensions match the input --------
test_that("output has the same dimensions as the input", {
  result <- pseudocount_for_test(otu_test, threshold_pct = 0.5)
  expect_equal(dim(result), dim(otu_test))
})

# -------- test 4: non-zero entries are left unchanged as proportions --------
test_that("non-zero entries are correctly converted to proportions", {
  result <- pseudocount_for_test(otu_test, threshold_pct = 0.5)
  # row 1: library size 15, so the first entry should be 10 / 15
  expect_equal(result[[1, 1]], 10 / 15, tolerance = 1e-8)
  # row 2: library size 5, so the second entry should be 3 / 5
  expect_equal(result[[2, 2]], 3 / 5, tolerance = 1e-8)
})

# -------- test 5: pseudocount values match the expected formula --------
test_that("replaced zeroes match threshold_pct * (1 / lib_size)", {
  result <- pseudocount_for_test(otu_test, threshold_pct = 0.5)
  # row 1 (lib_size = 15): pseudocount = 0.5 * (1 / 15)
  expect_equal(result[[1, 2]], 0.5 * (1 / 15), tolerance = 1e-8)
  # row 2 (lib_size = 5): pseudocount = 0.5 * (1 / 5)
  expect_equal(result[[2, 3]], 0.5 * (1 / 5), tolerance = 1e-8)
})