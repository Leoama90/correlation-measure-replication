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
  expect_true(nrow(otu) > nrow(otu.69001.H))
})

test_that("Column number after OTU filtering is not correct", {
  # prevalence and median filtering must reduce the number of OTUs
  expect_true(ncol(otu.filt) < ncol(otu.69001.H))
})

test_that("OTU prevalence filter was not applied correctly", {
  # every remaining OTU must be present in at least 33% of samples
  prevalences <- colSums(otu.filt > 0) / nrow(otu.filt)
  expect_true(all(prevalences >= 0.33))
})

test_that("OTU median filter was not applied correctly", {
  # every remaining OTU must have median non-zero abundance >= 5
  medians <- apply(otu.filt, 2, function(x) median(x[x > 0]))
  expect_true(all(medians >= 5))
})


# -------- checks for phylum aggregation --------

test_that("Phylum matrix has wrong number of rows", {
  # phylum matrix must have the same number of samples as otu.filt
  expect_equal(nrow(phy), nrow(otu.filt))
})

test_that("Phylum matrix contains negative values", {
  # aggregated abundances must be non-negative
  expect_true(all(phy >= 0))
})

test_that("Phylum matrix column names do not match taxonomy phylum names", {
  # column names must be valid phylum names from the taxonomy table
  expect_true(all(colnames(phy) %in% taxa.filt[, "phylum"]))
})


# -------- checks for OTU-level correlation matrices --------

test_that("OTU-level correlation matrices are not square", {
  # all correlation matrices must be square
  expect_equal(nrow(res.clr), ncol(res.clr))
  expect_equal(nrow(res.cc) , ncol(res.cc ))
  expect_equal(nrow(res.rho), ncol(res.rho))
  expect_equal(nrow(res.L1) , ncol(res.L1 ))
})

test_that("OTU-level correlation matrices dimensions do not match otu.filt", {
  # number of rows/cols must equal number of filtered OTUs
  expect_equal(nrow(res.clr), ncol(otu.filt))
  expect_equal(nrow(res.cc) , ncol(otu.filt))
  expect_equal(nrow(res.rho), ncol(otu.filt))
  expect_equal(nrow(res.L1) , ncol(otu.filt))
})

test_that("SparCC OTU-level names do not match filtered OTU names", {
  # row and column names must be consistent with each other
  expect_true(identical(colnames(res.cc), rownames(res.cc)))
  # names must match the filtered OTU matrix columns
  expect_true(identical(colnames(res.cc), colnames(otu.filt)))
})

test_that("Rho OTU-level matrix is not symmetric", {
  # proportionality is a symmetric measure
  expect_true(isSymmetric.matrix(res.rho))
})

test_that("Rho OTU-level diagonal does not have only one values", {
  # diagonal of rho must be 1 (perfect self-proportionality)
  expect_true(all(diag(res.rho) == 1))
})


# -------- checks for phylum-level correlation matrices --------

test_that("Phylum-level correlation matrices are not square", {
  # all phylum-level correlation matrices must be square
  expect_equal(nrow(res.clr.phy), ncol(res.clr.phy))
  expect_equal(nrow(res.cc.phy ), ncol(res.cc.phy ))
  expect_equal(nrow(res.rho.phy), ncol(res.rho.phy))
  expect_equal(nrow(res.L1.phy ), ncol(res.L1.phy ))
})

test_that("Phylum-level correlation matrices dimensions do not match phy", {
  # number of rows/cols must equal number of phyla
  expect_equal(nrow(res.clr.phy), ncol(phy))
  expect_equal(nrow(res.cc.phy ), ncol(phy))
  expect_equal(nrow(res.rho.phy), ncol(phy))
  expect_equal(nrow(res.L1.phy ), ncol(phy))
})

test_that("SparCC phylum-level names do not match phylum names", {
  # row and column names must be consistent with each other
  expect_true(identical(colnames(res.cc.phy), rownames(res.cc.phy)))
  # names must match the phylum matrix columns
  expect_true(identical(colnames(res.cc.phy), colnames(phy)))
})

test_that("Rho phylum-level matrix is not symmetric", {
  # proportionality is a symmetric measure
  expect_true(isSymmetric.matrix(res.rho.phy))
})

test_that("Rho phylum-level diagonal does not have only one values", {
  # diagonal of rho must be 1 (perfect self-proportionality)
  expect_true(all(diag(res.rho.phy) == 1))
})