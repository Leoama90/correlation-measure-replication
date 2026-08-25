# test-methods_comparison_reviews.R
#
# Purpose:
#   this script tests the main outputs of the methods_comparison_review.R script,
#   checking:
#       - OTU filtering is correctly applied based on subject, health status,
#         prevalence, and median abundance
#       - phylum-level aggregation produces a valid abundance matrix
#       - OTU-level correlation matrices have the expected dimensions,
#         names, symmetry, and diagonal values
#       - phylum-level correlation matrices have the expected dimensions,
#         names, symmetry, and diagonal values
#
# Inputs:
#   - methods_comparison_review.R script (sourced below)
#   - OTU abundance and taxonomy data loaded by the script
#
# Outputs:
#   - test results printed to the console (pass/fail for each check)

# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)

# tidyverse: useful to manage data (dplyr) and make nice plots (ggplot2)
# https://cran.r-project.org/web/packages/tidyverse/index.html
library(tidyverse)

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)


# -------- load script that needs to be tested --------

source(
  list.files(
    path = here(),
    pattern = "^methods_comparison_review\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)


# -------- checks for OTU filtering --------

test_that("Filtering by subject and health status reduces the number of samples", {
  # filtering by subject and health status must reduce the number of samples
  expect_true(nrow(otu) > nrow(otu_69001_H))
})

test_that("OTU filtering reduces the number of OTUs", {
  # prevalence and median filtering must reduce the number of OTUs
  expect_true(ncol(otu_filt) < ncol(otu_69001_H))
})

test_that("All filtered OTUs meet the minimum prevalence threshold", {
  # every remaining OTU must be present in at least 33% of samples
  prevalences <- colSums(otu_filt > 0) / nrow(otu_filt)
  expect_true(all(prevalences >= 0.33))
})

test_that("All filtered OTUs meet the minimum median abundance threshold", {
  # every remaining OTU must have median non-zero abundance >= 5
  medians <- apply(otu_filt, 2, function(x) median(x[x > 0]))
  expect_true(all(medians >= 5))
})


# -------- checks for phylum aggregation --------

test_that("Phylum matrix has the same number of samples as the filtered OTU matrix", {
  # phylum matrix must have the same number of samples as otu_filt
  expect_equal(nrow(phy), nrow(otu_filt))
})

test_that("Phylum matrix contains only non-negative abundances", {
  # aggregated abundances must be non-negative
  expect_true(all(phy >= 0))
})

test_that("Phylum matrix contains valid phylum names", {
  # column names must be valid phylum names from the taxonomy table
  expect_true(all(colnames(phy) %in% taxa_filt[, "phylum"]))
})


# -------- checks for OTU-level correlation matrices --------

test_that("OTU-level correlation matrices are square", {
  # all correlation matrices must be square
  expect_equal(nrow(res_clr), ncol(res_clr))
  expect_equal(nrow(res_cc), ncol(res_cc))
  expect_equal(nrow(res_rho), ncol(res_rho))
  expect_equal(nrow(res_L1), ncol(res_L1))
})

test_that("OTU-level correlation matrices match the number of filtered OTUs", {
  # number of rows/cols must equal number of filtered OTUs
  expect_equal(nrow(res_clr), ncol(otu_filt))
  expect_equal(nrow(res_cc), ncol(otu_filt))
  expect_equal(nrow(res_rho), ncol(otu_filt))
  expect_equal(nrow(res_L1), ncol(otu_filt))
})

test_that("SparCC OTU-level matrix names match the filtered OTU names", {
  # row and column names must be consistent with each other
  expect_true(identical(colnames(res_cc), rownames(res_cc)))
  
  # names must match the filtered OTU matrix columns
  expect_true(identical(colnames(res_cc), colnames(otu_filt)))
})

test_that("Rho OTU-level matrix is symmetric", {
  # proportionality is a symmetric measure
  expect_true(isSymmetric.matrix(res_rho))
})

test_that("Rho OTU-level diagonal contains only one values", {
  # diagonal of rho must be 1 (perfect self-proportionality)
  expect_true(all(diag(res_rho) == 1))
})


# -------- checks for phylum-level correlation matrices --------

test_that("Phylum-level correlation matrices are square", {
  # all phylum-level correlation matrices must be square
  expect_equal(nrow(res_clr_phy), ncol(res_clr_phy))
  expect_equal(nrow(res_cc_phy), ncol(res_cc_phy))
  expect_equal(nrow(res_rho_phy), ncol(res_rho_phy))
  expect_equal(nrow(res_L1_phy), ncol(res_L1_phy))
})

test_that("Phylum-level correlation matrices match the number of phyla", {
  # number of rows/cols must equal number of phyla
  expect_equal(nrow(res_clr_phy), ncol(phy))
  expect_equal(nrow(res_cc_phy), ncol(phy))
  expect_equal(nrow(res_rho_phy), ncol(phy))
  expect_equal(nrow(res_L1_phy), ncol(phy))
})

test_that("SparCC phylum-level matrix names match the phylum names", {
  # row and column names must be consistent with each other
  expect_true(identical(colnames(res_cc_phy), rownames(res_cc_phy)))
  
  # names must match the phylum matrix columns
  expect_true(identical(colnames(res_cc_phy), colnames(phy)))
})

test_that("Rho phylum-level matrix is symmetric", {
  # proportionality is a symmetric measure
  expect_true(isSymmetric.matrix(res_rho_phy))
})

test_that("Rho phylum-level diagonal contains only one values", {
  # diagonal of rho must be 1 (perfect self-proportionality)
  expect_true(all(diag(res_rho_phy) == 1))
})