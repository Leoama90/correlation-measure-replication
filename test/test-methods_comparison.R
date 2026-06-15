# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)
# tidyverse: useful to manage data (dplyr) and make nice plots (ggplot2)
# https://cran.r-project.org/web/packages/tidyverse/index.html
library(tidyverse)
# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)
# SpiecEasi: package with spiec.easi and sparCC methods
# https://github.com/zdk123/SpiecEasi
library(SpiecEasi)
# propr: package with proportionality rho method
# https://github.com/tpq/propr
library(propr)
# ggpubr: nice plots based on ggplot2
# https://cran.r-project.org/web/packages/ggpubr/index.html
library(ggpubr)
# qualpalr: generate distinct qualitative color palette
# https://cran.r-project.org/web/packages/qualpalr/index.html
library(qualpalr)
# cowplot: draw ggplot2 in new figures
# https://cran.r-project.org/web/packages/cowplot/index.html
library(cowplot)
# vegan: used to elaborates shannon entropy
# https://cran.r-project.org/web/packages/vegan/index.html
library(vegan)
# lightweight package with necessary function to generate metagenomics data
# https://github.com/Fuschi/ToyModel
library(ToyModel)
# Imputation of zeros, left-censored and missing values in compositional data
# https://cran.r-project.org/web/packages/zCompositions/index.html
library(zCompositions)

# load function to test
source(here("script", "method_comparison", "methods_comparison.R"))


# -------- checks for OTU filtering --------

test_that("The row number is not correct", {
  # get otu row number
  otu_rown    <- nrow(otu)
  # get otu.69001.H row number
  otu.69_rown <- nrow(otu.69001.H)
  # expect less rows after filtering for subject and health status
  expect_true(otu_rown > otu.69_rown)
})

test_that("The column number is not correct", {
  # expect less OTUs after prevalence and median filtering
  expect_true(ncol(otu.filt) < ncol(otu.69001.H))
})


# -------- checks matrixes for SPIEC-EASI MB --------

test_that("Matrix has not yet been symmetrized", {
  # raw MB coefficients should not be symmetric before averaging
  expect_false(isSymmetric.matrix(beta.mb))
})

test_that("Matrix is not symmetric", {
  # adjacency matrix should be symmetric after averaging beta and its transpose
  expect_true(isSymmetric.matrix(adj.mb))
})

test_that("The selected matrixes entries are not equal 1 or -1", {
  # non-zero entries should be binarized to +1 or -1
  expect_true(all(adj.mb[adj.mb != 0] %in% c(-1, 1)))
})


# -------- checks correlation matrixes --------

test_that("The matrix is not squared", {
  # all correlation matrixes must be square
  expect_equal(nrow(res.cc ), ncol(res.cc ))
  expect_equal(nrow(res.rho), ncol(res.rho))
  expect_equal(nrow(res.clr), ncol(res.clr))
})


# -------- checks for SparCC --------

test_that("SparCC row and column names match filtered OTU names", {
  # row and column names must be consistent with each other
  expect_true(identical(colnames(res.cc), rownames(res.cc)))
  # names must match the filtered OTU matrix columns
  expect_true(identical(colnames(res.cc), colnames(otu.filt)))
})


# -------- checks for CLR + Pearson --------

test_that("Diagonal does not have only zero values!", {
  # diagonal is set to zero to remove self-correlations
  diag_values <- c(diag(res.clr))
  expect_true(all(diag_values == 0))
})


# -------- checks for GLASSO --------

test_that("GLASSO adjacency matrix is not symmetric", {
  # adjacency matrix must be symmetric for an undirected network
  expect_true(isSymmetric.matrix(adj.gl))
})

test_that("GLASSO adjacency matrix entries are not 1 or -1", {
  # non-zero entries should be binarized to +1 or -1
  expect_true(all(adj.gl[adj.gl != 0] %in% c(-1, 1)))
})

test_that("GLASSO row and column names do not match filtered OTU names", {
  # column names must match the filtered OTU matrix columns
  expect_true(identical(colnames(adj.gl), colnames(otu.filt)))
  # row names must match the filtered OTU matrix columns
  expect_true(identical(rownames(adj.gl), colnames(otu.filt)))
})


# -------- checks for Rho --------

test_that("Rho matrix is not symmetric", {
  # proportionality is a symmetric measure
  expect_true(isSymmetric.matrix(res.rho))
})

test_that("Rho diagonal does not have only one values!", {
  # diagonal of rho should be 1 (perfect self-proportionality)
  diag_values <- c(diag(res.rho))
  expect_true(all(diag_values == 1))
})