# test-clr_pearson.R
#
# Purpose:
#   This script tests the fundamental features of clr_on_data(), checking
#   the structure of the returned list, the absence of zeroes/non-finite
#   values after the CLR transformation, the row-wise zero-sum property
#   of the CLR transform, and the basic properties of the resulting
#   Pearson correlation matrix (dimensions, symmetry, unit diagonal).
#
# Inputs:
#   - clr_pearson.R (sourced below)
#
# Outputs:
#   - test results printed to the console (pass/fail for each check)
#
# Note:
#   clr_on_data() now accepts explicit prevalence_threshold and
#   threshold_pct arguments (see clr_pearson.R): when both are provided,
#   the readline() prompts inside filt_data() and pseudocount() are
#   skipped entirely, so the real function can be called directly and
#   non-interactively in these tests. The previous version of this test
#   mocked readline() globally to avoid the prompts hanging; that mock
#   is no longer needed and has been removed, so these tests now
#   exercise clr_on_data() exactly as it would be called in a
#   reproducible notebook.
#
#
# Note_02: running this suite prints three harmless warnings.
#
# One appears on every test that reaches pseudocount(): otu_test has no
# column names (it is an unnamed matrix), so tibble::as_tibble() warns
# that it is falling back to its compatibility name-repair behaviour
# ("must have unique column names if .name_repair is omitted"). This
# does not affect the numeric result and is unrelated to the specific
# thresholds used.
#
# The other two only appear in the out-of-range prevalence_threshold
# test above: they come from datasum(), which computes min()/max() on
# the now-empty (0-column) OTU table left after every OTU fails the
# (out-of-range) prevalence filter. They occur before the actual error
# (raised later, when pseudocount() finds every sample's library size
# equal to 0) and do not affect the outcome of expect_error().

# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)
# tidyverse: data science toolkit (dplyr, ggplot2, tidyr, etc.)
# https://tidyverse.org
library(tidyverse)
# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# load clr_on_data() and its dependencies (datasum, filt_data, pseudocount)
source(
  list.files(
    path = here(),
    pattern = "^clr_pearson\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)

# suppress cat() output so tests stay clean, keep the return value;
# thresholds are passed explicitly (0.1 for both), so no prompt ever
# fires and there is nothing to capture from user input, only from
# filt_data()'s/datasum()'s printed summaries
quiet_clr_on_data <- function(x, prevalence_threshold, threshold_pct) {
  capture.output(
    result <- clr_on_data(
      x,
      prevalence_threshold = prevalence_threshold,
      threshold_pct = threshold_pct
    )
  )
  result
}


# -------- shared test data --------

# a small OTU table (5 samples x 6 OTUs), with a single zero count, chosen
# so that every OTU survives filt_data()'s prevalence and median filters
# at the 0.1 threshold supplied explicitly below
otu_test <- matrix(
  c(10, 20, 15,  5,  8, 30,
    12, 18,  0,  6,  9, 25,
    11, 22, 14,  5, 10, 28,
    9, 19, 16,  6,  8, 27,
    13, 21, 15,  5,  9, 29),
  nrow = 5, byrow = TRUE
)

# -------- test 1: the returned list has the expected structure --------

test_that("clr_on_data returns a list with the expected elements", {
  result <- quiet_clr_on_data(otu_test, prevalence_threshold = 0.1, threshold_pct = 0.1)
  # the result should contain exactly these three named elements
  expect_named(result, c("samp_filt", "y_clr", "cor_matrix"))
})


# -------- test 2: the CLR-transformed data has no zeroes or non-finite values --------

test_that("y_clr contains no NA or infinite values", {
  result <- quiet_clr_on_data(otu_test, prevalence_threshold = 0.1, threshold_pct = 0.1)
  # pseudocount() removes every zero before the log step, so log(0)/-Inf
  # cannot occur; this checks no NA/Inf slipped through regardless
  expect_true(all(is.finite(result$y_clr)))
})


# -------- test 3: each row of y_clr sums to (approximately) zero --------

test_that("y_clr rows satisfy the CLR zero-sum property", {
  result <- quiet_clr_on_data(otu_test, prevalence_threshold = 0.1, threshold_pct = 0.1)
  # centering by the row-wise mean log value should make every row sum to 0
  row_sums <- rowSums(result$y_clr)
  expect_equal(row_sums, rep(0, nrow(result$y_clr)), tolerance = 1e-8)
})


# -------- test 4: the correlation matrix has the right dimensions --------

test_that("cor_matrix is square with one row/column per surviving OTU", {
  result <- quiet_clr_on_data(otu_test, prevalence_threshold = 0.1, threshold_pct = 0.1)
  n_otu <- ncol(result$y_clr)
  # cor() on columns should produce an OTU x OTU square matrix
  expect_equal(dim(result$cor_matrix), c(n_otu, n_otu))
})


# -------- test 5: the correlation matrix is symmetric with a unit diagonal --------

test_that("cor_matrix is symmetric and has 1s on the diagonal", {
  result <- quiet_clr_on_data(otu_test, prevalence_threshold = 0.1, threshold_pct = 0.1)
  # a Pearson correlation matrix must equal its own transpose
  expect_equal(result$cor_matrix, t(result$cor_matrix))
  # every variable is perfectly correlated with itself
  expect_true(all(diag(result$cor_matrix) == 1))
})


# -------- test 6: threshold validation is only enforced for threshold_pct --------

test_that("clr_on_data propagates pseudocount()'s threshold validation error", {
  # threshold_pct outside [0, 1] should error inside pseudocount(), and
  # that error should propagate up through clr_on_data()
  expect_error(
    quiet_clr_on_data(otu_test, prevalence_threshold = 0.1, threshold_pct = -0.2)
  )
})

# -------- test 7: an out-of-range prevalence_threshold still errors (indirectly) --------

test_that("an out-of-range prevalence_threshold still errors, but not via explicit validation", {
  # filt_data() has no explicit validation branch for a directly-supplied
  # prevalence_threshold: a value above 1 simply makes the prevalence
  # filter condition (colSums(x > 0) / nrow(x) >= prevalence_threshold)
  # false for every OTU, so every OTU is silently dropped (no error at
  # that point). clr_on_data() still errors overall, but only because
  # the now-empty OTU table gives every sample a library size of 0 in
  # the subsequent pseudocount() step. This test documents that chain
  # of events rather than asserting it is the ideal behaviour.
  expect_error(
    quiet_clr_on_data(otu_test, prevalence_threshold = 1.5, threshold_pct = 0.1)
  )
})

