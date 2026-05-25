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
  
  # count number of columns of otu.69001.H
  before_filt_col_num <- ncol(otu.69001.H)
  # count column number after filtering
  after_filt_col_num  <- ncol(otu.filt)    
  # compare the two numbers
  expect_true(before_filt_col_num != after_filt_col_num)
  # get taxa.filt number of rows
  taxa_row <- nrow(taxa.filt)
  # taxa.filt row number compared with otu.filt column number
  expect_true(after_filt_col_num == taxa_row)
  
})

# -------- Compare the names --------
test_that("the names doesn't match!", {
  
  # get otu.filt column names
  otu.filt_colname <- colnames(otu.filt)
  # get taxa.filt row names
  taxa.filt_rowname <- rownames(taxa.filt)
  # compare taxa.filt row names with otu.filt column names
  expect_identical(taxa.filt_rowname, otu.filt_colname)
  
})

# -------- Correlation matrixes test --------
test_that("There is a problem with the matrixes", {
  
  # get row and col num for cor.czm matrix
  czm_row_num <- nrow(cor.czm)
  czm_col_num <- ncol(cor.czm)
  # test the ncol equal nrow of cor.czm, to check if it is a square matrix
  expect_identical(czm_row_num, czm_col_num)
  
  # get row and col num for cor.gbm matrix
  gbm_row_num <- nrow(cor.gbm)
  gbm_col_num <- ncol(cor.gbm)
  # test the ncol equal nrow of cor.gbm, to check if it is a square matrix
  expect_identical(gbm_row_num, gbm_col_num)
  
  # get row and col num for cor.bl matrix
  bl_row_num <- nrow(cor.bl)
  bl_col_num <- ncol(cor.bl)
  # test the ncol equal nrow of cor.bl, to check if it is a square matrix
  expect_identical(bl_row_num, bl_col_num)
  
  # get row and col num for cor.65 matrix
  cor65_row_num <- nrow(cor.65)
  cor65_col_num <- ncol(cor.65)
  # test the ncol equal nrow of cor.bl, to check if it is a square matrix
  expect_identical(cor65_row_num, cor65_col_num)
  
})

# -------- cor.65 square matrix test --------
test_that("cor.65 is not a square matrix", {
  
  # test the ncol equal nrow of cor.65, to check if it is a square matrix
  expect_identical(nrow(cor.65), ncol(cor.65))
  
})

# -------- Symmetry test --------
test_that("correlation matrices are not symmetric", {
  
  expect_true(isSymmetric(cor.czm))
  expect_true(isSymmetric(cor.gbm))
  expect_true(isSymmetric(cor.bl))
  expect_true(isSymmetric(cor.65))
  
})

# -------- NA test --------
test_that("correlation matrices contain NA values", {
  
  expect_false(anyNA(cor.czm))
  expect_false(anyNA(cor.gbm))
  expect_false(anyNA(cor.bl))
  expect_false(anyNA(cor.65))
  
})

# -------- Zero replacement test --------
test_that("otu.filt.65 still contains zeros", {
  
  expect_false(any(otu.filt.65 == 0))
  
})

# -------- cor_df dimensions test --------
test_that("cor_df has unexpected number of rows", {
  
  # expected number of pairs in the upper triangle
  expected_pairs <- ncol(otu.filt) * (ncol(otu.filt) - 1) / 2
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

