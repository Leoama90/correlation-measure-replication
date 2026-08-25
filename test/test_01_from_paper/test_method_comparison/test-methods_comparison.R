# test-methods_comparison.R
#
# Purpose:
#   This script tests the fundamental features of methods_comparison.R,
#   checking:
#       - OTU filtering reduces both rows (subject/health subsetting)
#         and columns (prevalence/abundance filtering) as expected
#       - beta_mb (raw MB coefficients) is not symmetric before
#         averaging, and adj_mb is symmetric afterwards, with non-zero
#         entries binarized to exactly 1 or -1
#       - res_cc, res_rho, and res_clr (SparCC, Rho, CLR correlation
#         matrices) are square
#       - SparCC row/column names match each other and the filtered
#         OTU table
#       - res_clr has a zeroed diagonal (self-correlations removed)
#       - adj_gl (GLASSO adjacency matrix) is symmetric, with non-zero
#         entries binarized to exactly 1 or -1, and row/column names
#         matching the filtered OTU table
#       - res_rho (Rho/propr) is symmetric, with a diagonal of all 1s
#
# Inputs:
#   - methods_comparison.R (sourced below)
#   - otu_HMP2.rds, meta_HMP2.rds, taxonomy.rds, required by
#     methods_comparison.R itself
#
# Outputs:
#   - test results printed to the console (pass/fail for each check)
#
# Note:
#   Unlike some of the other test scripts in this project, this one
#   DOES source methods_comparison.R directly, which re-runs the full
#   analysis (SPIEC-EASI GLASSO/MB estimation, SparCC, propr, and all
#   plotting/export steps) on every test run. This is slow, writes
#   several PNG files to outputs/methods_comparison_outputs/, and calls
#   browseURL() at the end, which will open an image viewer window.
#   A faster alternative (following the pattern used elsewhere in this
#   project) would load pre-computed objects instead of sourcing the
#   whole script; not done here yet.

# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)
# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)


# -------- load the script to test --------

# sourcing methods_comparison.R also loads its own dependencies
# (SpiecEasi, propr, tidyverse, etc.), so they are not re-declared here
source(
  list.files(
    path = here(),
    pattern = "^methods_comparison\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)


# -------- checks for OTU filtering --------

test_that("otu_69001_H has fewer rows than otu, after filtering by subject and health status", {
  otu_rown <- nrow(otu)
  otu_69_rown <- nrow(otu_69001_H)
  expect_true(otu_rown > otu_69_rown)
})

test_that("otu_filt has fewer columns than otu_69001_H, after prevalence and median filtering", {
  expect_true(ncol(otu_filt) < ncol(otu_69001_H))
})


# -------- checks matrices for SPIEC-EASI MB --------

test_that("beta_mb (raw MB coefficients) is not symmetric before averaging", {
  expect_false(isSymmetric.matrix(beta_mb))
})

test_that("adj_mb is symmetric after averaging beta and its transpose", {
  expect_true(isSymmetric.matrix(adj_mb))
})

test_that("non-zero entries of adj_mb are exactly 1 or -1", {
  expect_true(all(adj_mb[adj_mb != 0] %in% c(-1, 1)))
})


# -------- checks correlation matrices --------

test_that("res_cc, res_rho, and res_clr are square matrices", {
  expect_equal(nrow(res_cc), ncol(res_cc))
  expect_equal(nrow(res_rho), ncol(res_rho))
  expect_equal(nrow(res_clr), ncol(res_clr))
})


# -------- checks for SparCC --------

test_that("SparCC row and column names match filtered OTU names", {
  expect_true(identical(colnames(res_cc), rownames(res_cc)))
  expect_true(identical(colnames(res_cc), colnames(otu_filt)))
})


# -------- checks for CLR + Pearson --------

test_that("res_clr diagonal is set to zero (self-correlations removed)", {
  diag_values <- c(diag(res_clr))
  expect_true(all(diag_values == 0))
})


# -------- checks for GLASSO --------

test_that("GLASSO adjacency matrix (adj_gl) is symmetric", {
  expect_true(isSymmetric.matrix(adj_gl))
})

test_that("non-zero entries of adj_gl are exactly 1 or -1", {
  expect_true(all(adj_gl[adj_gl != 0] %in% c(-1, 1)))
})

test_that("GLASSO row and column names match filtered OTU names", {
  expect_true(identical(colnames(adj_gl), colnames(otu_filt)))
  expect_true(identical(rownames(adj_gl), colnames(otu_filt)))
})


# -------- checks for Rho --------

test_that("Rho matrix (res_rho) is symmetric", {
  expect_true(isSymmetric.matrix(res_rho))
})

test_that("Rho diagonal is all 1s (perfect self-proportionality)", {
  diag_values <- c(diag(res_rho))
  expect_true(all(diag_values == 1))
})