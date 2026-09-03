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

# override readline() so the interactive prompts inside filt_data() and
# pseudocount() are answered automatically with a valid threshold (0.1),
# instead of hanging (or looping forever) when tests run non-interactively
assign("readline", function(prompt = "") "0.1", envir = .GlobalEnv)

# suppress cat() output so tests stay clean, keep the return value
quiet_clr_on_data <- function(x) {
  capture.output(result <- clr_on_data(x))
  result
}


# -------- shared test data --------

# a small OTU table (5 samples x 6 OTUs), with a single zero count, chosen
# so that every OTU survives filt_data()'s prevalence and median filters
# at the 0.1 threshold supplied by the mocked readline() above
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
  result <- quiet_clr_on_data(otu_test)
  # the result should contain exactly these three named elements
  expect_named(result, c("samp_filt", "y_clr", "cor_matrix"))
})


# -------- test 2: the CLR-transformed data has no zeroes or non-finite values --------

test_that("y_clr contains no NA or infinite values", {
  result <- quiet_clr_on_data(otu_test)
  # pseudocount() removes every zero before the log step, so log(0)/-Inf
  # cannot occur; this checks no NA/Inf slipped through regardless
  expect_true(all(is.finite(result$y_clr)))
})


# -------- test 3: each row of y_clr sums to (approximately) zero --------

test_that("y_clr rows satisfy the CLR zero-sum property", {
  result <- quiet_clr_on_data(otu_test)
  # centering by the row-wise mean log value should make every row sum to 0
  row_sums <- rowSums(result$y_clr)
  expect_equal(row_sums, rep(0, nrow(result$y_clr)), tolerance = 1e-8)
})


# -------- test 4: the correlation matrix has the right dimensions --------

test_that("cor_matrix is square with one row/column per surviving OTU", {
  result <- quiet_clr_on_data(otu_test)
  n_otu <- ncol(result$y_clr)
  # cor() on columns should produce an OTU x OTU square matrix
  expect_equal(dim(result$cor_matrix), c(n_otu, n_otu))
})


# -------- test 5: the correlation matrix is symmetric with a unit diagonal --------

test_that("cor_matrix is symmetric and has 1s on the diagonal", {
  result <- quiet_clr_on_data(otu_test)
  # a Pearson correlation matrix must equal its own transpose
  expect_equal(result$cor_matrix, t(result$cor_matrix))
  # every variable is perfectly correlated with itself
  expect_true(all(diag(result$cor_matrix) == 1))
})

# remove the no longer needed readline created in this test 
rm(readline, envir = .GlobalEnv)