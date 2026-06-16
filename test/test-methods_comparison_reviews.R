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
source(here("script", "method_comparison", "methods_comparison_review.R"))


# -------- checks for OTU filtering --------

test_that("Row number after subject/health filter is not correct", {
  # filtering by subject and health status must reduce the number of samples
  expect_true(nrow(otu) > nrow(otu_69001_H))
})

test_that("Column number after OTU filtering is not correct", {
  # prevalence and median filtering must reduce the number of OTUs
  expect_true(ncol(otu_filt) < ncol(otu_69001_H))
})

test_that("OTU prevalence filter was not applied correctly", {
  # every remaining OTU must be present in at least 33% of samples
  prevalences <- colSums(otu_filt > 0) / nrow(otu_filt)
  expect_true(all(prevalences >= 0.33))
})

test_that("OTU median filter was not applied correctly", {
  # every remaining OTU must have median non-zero abundance >= 5
  medians <- apply(otu_filt, 2, function(x) median(x[x > 0]))
  expect_true(all(medians >= 5))
})


# -------- checks for phylum aggregation --------

test_that("Phylum matrix has wrong number of rows", {
  # phylum matrix must have the same number of samples as otu_filt
  expect_equal(nrow(phy), nrow(otu_filt))
})

test_that("Phylum matrix contains negative values", {
  # aggregated abundances must be non-negative
  expect_true(all(phy >= 0))
})

test_that("Phylum matrix column names do not match taxonomy phylum names", {
  # column names must be valid phylum names from the taxonomy table
  expect_true(all(colnames(phy) %in% taxa_filt[, "phylum"]))
})


# -------- checks for OTU-level correlation matrices --------

test_that("OTU-level correlation matrices are not square", {
  # all correlation matrices must be square
  expect_equal(nrow(res_clr), ncol(res_clr))
  expect_equal(nrow(res_cc) , ncol(res_cc ))
  expect_equal(nrow(res_rho), ncol(res_rho))
  expect_equal(nrow(res_L1) , ncol(res_L1 ))
})

test_that("OTU-level correlation matrices dimensions do not match otu_filt", {
  # number of rows/cols must equal number of filtered OTUs
  expect_equal(nrow(res_clr), ncol(otu_filt))
  expect_equal(nrow(res_cc) , ncol(otu_filt))
  expect_equal(nrow(res_rho), ncol(otu_filt))
  expect_equal(nrow(res_L1) , ncol(otu_filt))
})

test_that("SparCC OTU-level names do not match filtered OTU names", {
  # row and column names must be consistent with each other
  expect_true(identical(colnames(res_cc), rownames(res_cc)))
  # names must match the filtered OTU matrix columns
  expect_true(identical(colnames(res_cc), colnames(otu_filt)))
})

test_that("Rho OTU-level matrix is not symmetric", {
  # proportionality is a symmetric measure
  expect_true(isSymmetric.matrix(res_rho))
})

test_that("Rho OTU-level diagonal does not have only one values", {
  # diagonal of rho must be 1 (perfect self-proportionality)
  expect_true(all(diag(res_rho) == 1))
})


# -------- checks for phylum-level correlation matrices --------

test_that("Phylum-level correlation matrices are not square", {
  # all phylum-level correlation matrices must be square
  expect_equal(nrow(res_clr.phy), ncol(res_clr.phy))
  expect_equal(nrow(res_cc_phy ), ncol(res_cc_phy ))
  expect_equal(nrow(res_rho_phy), ncol(res_rho_phy))
  expect_equal(nrow(res_L1_phy ), ncol(res_L1_phy ))
})

test_that("Phylum-level correlation matrices dimensions do not match phy", {
  # number of rows/cols must equal number of phyla
  expect_equal(nrow(res_clr.phy), ncol(phy))
  expect_equal(nrow(res_cc_phy ), ncol(phy))
  expect_equal(nrow(res_rho_phy), ncol(phy))
  expect_equal(nrow(res_L1_phy ), ncol(phy))
})

test_that("SparCC phylum-level names do not match phylum names", {
  # row and column names must be consistent with each other
  expect_true(identical(colnames(res_cc_phy), rownames(res_cc_phy)))
  # names must match the phylum matrix columns
  expect_true(identical(colnames(res_cc_phy), colnames(phy)))
})

test_that("Rho phylum-level matrix is not symmetric", {
  # proportionality is a symmetric measure
  expect_true(isSymmetric.matrix(res_rho_phy))
})

test_that("Rho phylum-level diagonal does not have only one values", {
  # diagonal of rho must be 1 (perfect self-proportionality)
  expect_true(all(diag(res_rho_phy) == 1))
})