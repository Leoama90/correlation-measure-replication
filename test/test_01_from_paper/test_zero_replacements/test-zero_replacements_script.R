# test-zero_replacements_script.R
#
# Purpose:
#   This script tests the fundamental features of
#   zero_replacements_script.R, checking:
#       - OTU/taxa filtering produces consistent dimensions (column
#         count changes after filtering, and matches taxa_filt rows)
#       - otu_filt column names and taxa_filt row names are identical
#       - all correlation matrices (cor_czm, cor_gbm, cor_bl, cor_65)
#         are square, symmetric, and contain no NA values
#       - otu_filt_65 no longer contains any zeros
#       - cor_df has one row per unique OTU pair (upper triangle), with
#         no duplicated pairs
#       - cor_long has exactly 4 times the rows of cor_df (one row per
#         zero-replacement method)
#       - max_abs_diff is always non-negative, and the higher flag is
#         coherent with the 0.1 threshold
#       - the output PNG file was created
#
# Inputs:
#   - zero_replacements_script.R (sourced below)
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


# -------- recall function to test --------

source(
  list.files(
    path = here(),
    pattern = "^zero_replacements_script\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)


# -------- compare column and row size after filtering --------

test_that("column count changes after filtering, and matches taxa_filt rows", {
  # count number of columns of otu_69001_H
  before_filt_col_num <- ncol(otu_69001_H)
  # count column number after filtering
  after_filt_col_num <- ncol(otu_filt)
  # compare the two numbers
  expect_true(before_filt_col_num != after_filt_col_num)
  # get taxa_filt number of rows
  taxa_row <- nrow(taxa_filt)
  # taxa_filt row number compared with otu_filt column number
  expect_true(after_filt_col_num == taxa_row)
})


# -------- compare the names --------

test_that("otu_filt column names match taxa_filt row names", {
  # get otu_filt column names
  otu_filt_colname <- colnames(otu_filt)
  # get taxa_filt row names
  taxa_filt_rowname <- rownames(taxa_filt)
  # compare taxa_filt row names with otu_filt column names
  expect_identical(taxa_filt_rowname, otu_filt_colname)
})


# -------- correlation matrices are square --------

test_that("all correlation matrices are square", {
  # get row and col num for cor_czm matrix
  czm_row_num <- nrow(cor_czm)
  czm_col_num <- ncol(cor_czm)
  # test the ncol equal nrow of cor_czm, to check if it is a square matrix
  expect_identical(czm_row_num, czm_col_num)
  
  # get row and col num for cor_gbm matrix
  gbm_row_num <- nrow(cor_gbm)
  gbm_col_num <- ncol(cor_gbm)
  # test the ncol equal nrow of cor_gbm, to check if it is a square matrix
  expect_identical(gbm_row_num, gbm_col_num)
  
  # get row and col num for cor_bl matrix
  bl_row_num <- nrow(cor_bl)
  bl_col_num <- ncol(cor_bl)
  # test the ncol equal nrow of cor_bl, to check if it is a square matrix
  expect_identical(bl_row_num, bl_col_num)
  
  # get row and col num for cor_65 matrix
  cor65_row_num <- nrow(cor_65)
  cor65_col_num <- ncol(cor_65)
  # test the ncol equal nrow of cor_65, to check if it is a square matrix
  expect_identical(cor65_row_num, cor65_col_num)
})


# -------- symmetry test --------

test_that("correlation matrices are symmetric", {
  expect_true(isSymmetric(cor_czm))
  expect_true(isSymmetric(cor_gbm))
  expect_true(isSymmetric(cor_bl))
  expect_true(isSymmetric(cor_65))
})


# -------- NA test --------

test_that("correlation matrices contain no NA values", {
  expect_false(anyNA(cor_czm))
  expect_false(anyNA(cor_gbm))
  expect_false(anyNA(cor_bl))
  expect_false(anyNA(cor_65))
})


# -------- zero replacement test --------

test_that("otu_filt_65 no longer contains zeros", {
  expect_false(any(otu_filt_65 == 0))
})


# -------- cor_df dimensions test --------

test_that("cor_df has the expected number of rows (one per unique OTU pair)", {
  # expected number of pairs in the upper triangle
  expected_pairs <- ncol(otu_filt) * (ncol(otu_filt) - 1) / 2
  expect_identical(nrow(cor_df), as.integer(expected_pairs))
})


# -------- no duplicate pairs in cor_df --------

test_that("cor_df contains no duplicate OTU pairs", {
  # concatenate var1 and var2 to create a unique pair identifier
  pairs <- paste(cor_df$var1, cor_df$var2, sep = "_")
  expect_false(any(duplicated(pairs)))
})


# -------- cor_long dimensions test --------

test_that("cor_long has 4 times the rows of cor_df (one per method)", {
  # cor_long should have 4x the rows of cor_df (one per method)
  expect_identical(nrow(cor_long), nrow(cor_df) * 4L)
})


# -------- max_abs_diff non-negativity test --------

test_that("max_abs_diff is always non-negative", {
  expect_true(all(tbl$max_abs_diff >= 0))
})


# -------- higher flag coherence test --------

test_that("higher flag is coherent with the 0.1 threshold", {
  expect_identical(tbl$higher, tbl$max_abs_diff > .1)
})


# -------- PNG output test --------

test_that("PNG file was created", {
  expect_true(file.exists(here("Plots", "correlation_differences_between_zeroRepl.png")))
})