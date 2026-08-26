# test-B1_plots_zinb_ecdf.R
#
# Purpose:
#   This script tests the fundamental features of B1_plots_zinb_ecdf.R,
#   checking:
#       - df is read correctly as a list-like object (data frame)
#       - pstr0_2 was dropped, and pstr0_1 was renamed to pstr0
#       - df has the ERR_CLR column, with all values below 1
#       - p is a valid ggplot object with a geom_tile layer
#       - p's underlying data has a MEAN_ERR_CLR column, bounded in [0, 1]
#       - the plot's x and y axis labels match the expected "r" and phi
#       - the output file error_zinb.rds was saved at the expected path
#         and contains a ggplot object
#
# Inputs:
#   - B1_plots_zinb_ecdf.R (sourced below)
#   - Sparsity_Effects_zinbin_rare_ecdf_2.rds, required by the sourced script
#
# Outputs:
#   - test results printed to the console (pass/fail for each check)

# ggpubr: nice plots based on ggplot2
# https://cran.r-project.org/web/packages/ggpubr/index.html
library(ggpubr)
# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)
# tidyverse: useful to manage data (dplyr) and make nice plots (ggplot2)
# https://cran.r-project.org/web/packages/tidyverse/index.html
library(tidyverse)

# load script to test
source(
  list.files(
    path = here(),
    pattern = "^B1_plots_zinb_ecdf\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)


# -------- test that the dataset was correctly read --------

test_that("df was correctly read as a list-like object", {
  expect_true(is.list(df))
})


# -------- control the mutate pipeline --------

test_that("the column pstr0_2 was dropped", {
  expect_false(hasName(df, "pstr0_2"))
})

test_that("the column pstr0_1 was renamed to pstr0", {
  expect_false(hasName(df, "pstr0_1"))
  expect_true(hasName(df, "pstr0"))
})

test_that("the ERR_CLR column is present", {
  expect_true(hasName(df, "ERR_CLR"))
})

test_that("all values in the ERR_CLR column are below 1", {
  expect_true(all(df$ERR_CLR < 1))
})


# -------- control the aggregated data and plot --------

test_that("p is a ggplot object", {
  expect_true(inherits(p, "ggplot"))
})

test_that("MEAN_ERR_CLR column is present in plot data", {
  expect_true(hasName(p$data, "MEAN_ERR_CLR"))
})

test_that("MEAN_ERR_CLR values are all within [0, 1]", {
  expect_true(all(p$data$MEAN_ERR_CLR >= 0))
  expect_true(all(p$data$MEAN_ERR_CLR <= 1))
})

test_that("the plot uses geom_tile", {
  geom_classes <- sapply(p$layers, function(l) class(l$geom)[1])
  expect_true("GeomTile" %in% geom_classes)
})

test_that("x-axis label is 'r'", {
  expect_equal(p$labels$x, "r")
})

test_that("y-axis label is phi", {
  expect_equal(as.character(p$labels$y), "phi")
})


# -------- control the saved RDS file --------

test_that("the RDS file was saved at the expected path", {
  expect_true(file.exists(here("01_from_paper", "sparsity_effects", "error_zinb.rds")))
})

test_that("the saved RDS file contains a ggplot object", {
  saved_p <- readRDS(here("01_from_paper", "sparsity_effects", "error_zinb.rds"))
  expect_true(inherits(saved_p, "ggplot"))
})