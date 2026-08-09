# ggpubr: nice plots based on ggplot2
# https://cran.r-project.org/web/packages/ggpubr/index.html
library(ggpubr)

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# tidyverse: useful to manage data (dplyr) and make nice plots (ggplot2)
# https://cran.r-project.org/web/packages/tidyverse/index.html
library(tidyverse)

# load function to test
source(
  list.files(
    path = here(),
    pattern = "^B1_plots_zinb_ecdf\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)

# -------- test that the database was correctly read --------
test_that("The file was not correctly read", {
  expect_true(is.list(df))
})

# -------- control the mutate pipeline --------
# check that pstr0_2 was dropped and pstr0_1 was renamed
test_that("The column pstr0_2 was not dropped", {
  expect_false(hasName(df, "pstr0_2"))
})

test_that("The column pstr0_1 was not renamed to pstr0", {
  expect_false(hasName(df, "pstr0_1"))
  expect_true(hasName(df, "pstr0"))
})

# check for the presence of the ERR_CLR column
test_that("The ERR_CLR column is not present", {
  expect_true(hasName(df, "ERR_CLR"))
})

# checks that all values in ERR_CLR are less than 1 (safety check mirrors the stop() in script)
test_that("There is a value in the ERR_CLR column that is higher than 1", {
  expect_true(all(df$ERR_CLR < 1))
})

# -------- control the aggregated data and plot --------
# check that p is a ggplot object
test_that("p is not a ggplot object", {
  expect_true(inherits(p, "ggplot"))
})

# check that the aggregated data inside p contains MEAN_ERR_CLR
test_that("MEAN_ERR_CLR column is missing from plot data", {
  expect_true(hasName(p$data, "MEAN_ERR_CLR"))
})

# check that MEAN_ERR_CLR values are all within [0, 1]
test_that("MEAN_ERR_CLR contains values outside [0, 1]", {
  expect_true(all(p$data$MEAN_ERR_CLR >= 0))
  expect_true(all(p$data$MEAN_ERR_CLR <= 1))
})

# check that the plot uses geom_tile
test_that("The plot does not use geom_tile", {
  geom_classes <- sapply(p$layers, function(l) class(l$geom)[1])
  expect_true("GeomTile" %in% geom_classes)
})

# check axis labels
test_that("x-axis label is not 'r'", {
  expect_equal(p$labels$x, "r")
})

test_that("y-axis label is missing or incorrect", {
  expect_false(is.null(p$labels$y))
})

# -------- control the saved RDS file --------
test_that("The RDS file was not saved", {
  expect_true(file.exists(here("script", "sparsity_effects", "error_zinb.rds")))
})

test_that("The saved RDS file does not contain a ggplot object", {
  saved_p <- readRDS(here("script", "sparsity_effects", "error_zinb.rds"))
  expect_true(inherits(saved_p, "ggplot"))
})
