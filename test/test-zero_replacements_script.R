# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)
# tidyverse: useful to manage data (dplyr) and make nice plots (ggplot2)
# https://cran.r-project.org/web/packages/tidyverse/index.html
library(tidyverse)
# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# recall function to test
source(here("script", "zero_replacements", "zero_replacements_script.R"))


# -------- Compare column and row size after filtering --------
test_that("the column number or row is not the expected one", {
  
  # count number of columns of otu_69001_H
  before_filt_col_num <- ncol(otu_69001_H)
  # count column number after filtering
  after_filt_col_num  <- ncol(otu_filt)    
  # compare the two numbers
  expect_true(before_filt_col_num != after_filt_col_num)
  # get taxa_filt number of rows
  taxa_row <- nrow(taxa_filt)
  # taxa_filt row number compared with otu_filt column number
  expect_true(after_filt_col_num == taxa_row)
  
})

# -------- Compare the names --------
test_that("the names doesn't match!", {
  
  # get otu_filt column names
  otu_filt_colname  <- colnames(otu_filt)
  # get taxa_filt row names
  taxa_filt_rowname <- rownames(taxa_filt)
  # compare taxa_filt row names with otu_filt column names
  expect_identical(taxa_filt_rowname, otu_filt_colname)
  
})

# -------- Correlation matrixes test --------
test_that("There is a problem with the matrixes", {
  
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
  # test the ncol equal nrow of cor_bl, to check if it is a square matrix
  expect_identical(cor65_row_num, cor65_col_num)
  
})

# -------- cor_65 square matrix test --------
test_that("cor_65 is not a square matrix", {
  
  # test the ncol equal nrow of cor_65, to check if it is a square matrix
  expect_identical(nrow(cor_65), ncol(cor_65))
  
})

# -------- Symmetry test --------
test_that("correlation matrices are not symmetric", {
  
  expect_true(isSymmetric(cor_czm))
  expect_true(isSymmetric(cor_gbm))
  expect_true(isSymmetric(cor_bl))
  expect_true(isSymmetric(cor_65))
  
})

# -------- NA test --------
test_that("correlation matrices contain NA values", {
  
  expect_false(anyNA(cor_czm))
  expect_false(anyNA(cor_gbm))
  expect_false(anyNA(cor_bl))
  expect_false(anyNA(cor_65))
  
})

# -------- Zero replacement test --------
test_that("otu_filt_65 still contains zeros", {
  
  expect_false(any(otu_filt_65 == 0))
  
})

# -------- cor_df dimensions test --------
test_that("cor_df has unexpected number of rows", {
  
  # expected number of pairs in the upper triangle
  expected_pairs <- ncol(otu_filt) * (ncol(otu_filt) - 1) / 2
  expect_identical(nrow(cor_df), as.integer(expected_pairs))
  
})

# -------- No duplicate pairs in cor_df --------
test_that("cor_df contains duplicate OTU pairs", {
  
  # concatenate var1 and var2 to create a unique pair identifier
  pairs <- paste(cor_df$var1, cor_df$var2, sep = "_")
  expect_false(any(duplicated(pairs)))
  
})

# -------- cor_long dimensions test --------
test_that("cor_long has unexpected number of rows", {
  
  # cor_long should have 4x the rows of cor_df (one per method)
  expect_identical(nrow(cor_long), nrow(cor_df) * 4L)
  
})

# -------- max_abs_diff non-negativity test --------
test_that("max_abs_diff contains negative values", {
  
  expect_true(all(tbl$max_abs_diff >= 0))
  
})

# -------- higher flag coherence test --------
test_that("higher flag is not coherent with threshold 0.1", {
  
  expect_identical(tbl$higher, tbl$max_abs_diff > .1)
  
})

# -------- PNG output test --------
test_that("PNG file was not created", {
  
  expect_true(file.exists(here("Plots", "correlation_differences_between_zeroRepl.png")))
  
})

