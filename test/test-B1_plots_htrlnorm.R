# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)

# tidyverse: useful to manage data (dplyr) and make nice plots (ggplot2)
# https://cran.r-project.org/web/packages/tidyverse/index.html
library(tidyverse)

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# load script to test
source(here("script", "sparsity_effects", "B1_plots_htrlnorm.R"))

# -------- test on the dataframe df --------

# check that the "df" data frame does not have the phi_2 column
test_that("Data frame still has the phi_2 column", {
  expect_false("phi_2" %in% colnames(df))
})

# check that the "df" data frame does not have the phi_1 column
test_that("Data frame still has the phi_1 column", {
  expect_false("phi_1" %in% colnames(df))
})

# check that the renaming phi_1 -> phi worked
test_that("Data frame has the phi column after renaming", {
  expect_true("phi" %in% colnames(df))
})

# check that ERR_CLR column was created
test_that("Data frame has the ERR_CLR column", {
  expect_true("ERR_CLR" %in% colnames(df))
})

# check that ERR_CLR is computed correctly as abs(cor_normal - cor_NorTA_PCLR)
test_that("ERR_CLR is correctly computed as absolute error", {
  expected <- abs(df$cor_normal - df$cor_NorTA_PCLR)
  expect_equal(df$ERR_CLR, expected)
})

# check that no ERR_CLR value exceeds 1 (mirrors the stop() guard in the script)
test_that("No ERR_CLR value exceeds 1", {
  expect_true(all(df$ERR_CLR <= 1))
})

# check that ERR_CLR is always non-negative (it's an absolute value)
test_that("ERR_CLR values are all non-negative", {
  expect_true(all(df$ERR_CLR >= 0))
})

# -------- test on the plot p --------

# check that the plot object p exists and is a ggplot
test_that("p is a ggplot object", {
  expect_true(inherits(p, "gg"))
})

# check that the plot has a geom_tile layer
test_that("Plot contains a geom_tile layer", {
  geom_classes <- sapply(p$layers, function(l) class(l$geom)[1])
  expect_true("GeomTile" %in% geom_classes)
})

# -------- test the output file --------

# check that the output file error.rds was created
test_that("error.rds file was saved", {
  expect_true(file.exists("error.rds"))
})

# check that error.rds contains a ggplot object
test_that("error.rds contains a ggplot object", {
  saved_plot <- readRDS("error.rds")
  expect_true(inherits(saved_plot, "gg"))
})