# test-Z_final_plot.R
#
# Purpose:
#   This script tests the fundamental features of Z_final_plot.R,
#   checking:
#       - exactly one file is found for each expected .rds pattern
#         (error_zinb.rds, example_couple.rds), with no missing or
#         duplicated matches
#       - the legend.position theme override to "top" is correctly
#         stored on the loaded heatmap plot
#       - combining the two panels with ggarrange() produces a
#         renderable object, and that it actually draws without error
#         when printed to a real graphics device
#
# Inputs:
#   - error_zinb.rds, example_couple.rds (located via list.files())
#
# Outputs:
#   - test results printed to the console (pass/fail for each check)
#
# Note:
#   This script does not source Z_final_plot.R directly. Instead, it
#   loads the same two saved plot objects the script depends on, and
#   re-applies the same theme override and ggarrange() call locally, so
#   the composition logic can be tested without also re-running the
#   PNG export and browseURL() calls at the end of the original script
#   (the latter would open an image viewer window on every test run).

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)
# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)
# tidyverse: useful to manage data (dplyr) and make nice plots (ggplot2)
# https://cran.r-project.org/web/packages/tidyverse/index.html
library(tidyverse)
# ggpubr: tools for publication-ready ggplot2 figures (ggarrange)
# https://cran.r-project.org/web/packages/ggpubr/index.html
library(ggpubr)

# resolve file paths once so every test block can reuse them without
# repeating the list.files() call
path_err <- list.files(
  path = here(),
  pattern = "^error_zinb\\.rds$",
  full.names = TRUE,
  recursive = TRUE
)
path_couple <- list.files(
  path = here(),
  pattern = "example_couple.rds",
  full.names = TRUE,
  recursive = TRUE
)


# -------- input files are found --------

test_that("exactly one file is found for each RDS pattern", {
  # fail if the file is missing or if duplicates exist under different subdirectories
  expect_length(path_err, 1)
  expect_length(path_couple, 1)
})


# -------- legend position override --------

test_that("legend.position is set to 'top' after theme override", {
  p_err_zinb <- readRDS(path_err) +
    theme(legend.position = "top")

  # verify that the override was actually stored in the plot's theme layer
  expect_equal(p_err_zinb$theme$legend.position, "top")
})


# -------- combined plot is a valid ggarrange object --------

test_that("ggarrange returns a renderable object without errors", {
  p_err_zinb <- readRDS(path_err) +
    theme(legend.position = "top")
  p_couple <- readRDS(path_couple)

  # combine the two panels with a shared legend, mirroring the main script call
  pall <- ggpubr::ggarrange(p_err_zinb, p_couple,
    common.legend = TRUE,
    legend        = "top",
    labels        = c("A", "")
  )

  # the returned object must belong to a renderable class
  expect_true(inherits(pall, c("ggarrange", "gg")))

  # render to a temporary PNG to confirm the full print pipeline works end-to-end;
  # this catches errors that only surface during actual drawing (e.g. missing aesthetics)
  expect_no_error({
    grDevices::png(tempfile(fileext = ".png"), width = 480, height = 240)
    print(pall)
    grDevices::dev.off()
  })
})
