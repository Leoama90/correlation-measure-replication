# test-B1_plots_htrlnorm.R
#
# Purpose:
#   This script tests the fundamental features of B1_plots_htrlnorm.R,
#   checking:
#       - df no longer has the phi_2 column, and phi_1 was renamed to phi
#       - ERR_CLR is computed correctly as abs(cor_normal - cor_NorTA_PCLR),
#         is always non-negative, and never exceeds 1
#       - p is a valid ggplot object containing a geom_tile layer
#       - the output file error_htlrlnorm.rds was saved at the expected
#         path and contains a ggplot object
#
# Inputs:
#   - B1_plots_htrlnorm.R (sourced below)
#   - Sparsity_Effects_htrlnorm_2.rds, required by the sourced script
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


# -------- load script to test --------

source(
  list.files(
    path = here(),
    pattern = "^B1_plots_htrlnorm\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)


# -------- test on the dataframe df --------

test_that("df no longer has the phi_2 column", {
  expect_false("phi_2" %in% colnames(df))
})

test_that("df no longer has the phi_1 column (renamed to phi)", {
  expect_false("phi_1" %in% colnames(df))
})

test_that("df has the phi column after renaming", {
  expect_true("phi" %in% colnames(df))
})

test_that("df has the ERR_CLR column", {
  expect_true("ERR_CLR" %in% colnames(df))
})

test_that("ERR_CLR is correctly computed as absolute error", {
  expected <- abs(df$cor_normal - df$cor_NorTA_PCLR)
  expect_equal(df$ERR_CLR, expected)
})

test_that("no ERR_CLR value exceeds 1", {
  expect_true(all(df$ERR_CLR <= 1))
})

test_that("ERR_CLR values are all non-negative", {
  expect_true(all(df$ERR_CLR >= 0))
})


# -------- test on the plot p --------

test_that("p is a ggplot object", {
  expect_true(inherits(p, "gg"))
})

test_that("plot contains a geom_tile layer", {
  geom_classes <- sapply(p$layers, function(l) class(l$geom)[1])
  expect_true("GeomTile" %in% geom_classes)
})


# -------- test the output file --------

test_that("error_htlrlnorm.rds file was saved at the expected path", {
  expect_true(file.exists(here("01_from_paper", "sparsity_effects", "error_htlrlnorm.rds")))
})

test_that("error_htlrlnorm.rds contains a ggplot object", {
  saved_plot <- readRDS(here("01_from_paper", "sparsity_effects", "error_htlrlnorm.rds"))
  expect_true(inherits(saved_plot, "gg"))
})