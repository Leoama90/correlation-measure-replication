# test-demo_pseudocount.R
#
# Purpose:
#   This script tests the fundamental features of demo_pseudocount.R,
#   checking:
#       - demo_otu has the expected known structure (dimensions, values,
#         one zero per row in a different column)
#       - demo_result has no zeroes left, matching the number of rows
#         and columns of demo_otu
#       - non-zero entries of demo_result are correctly converted to
#         proportions (raw count / row library size)
#       - the pseudocount replacing each zero matches the expected
#         formula (threshold_pct * detection_limit), for the threshold
#         used in this test
#
# Inputs:
#   - demo_pseudocount.R (sourced below)
#
# Outputs:
#   - test results printed to the console (pass/fail for each check)
#
# Note:
#   pseudocount() always asks for its threshold interactively via
#   readline(), with no non-interactive override available (unlike
#   filt_data()). To run demo_pseudocount.R (and this test) without
#   hanging, readline() is overridden to always return "0.5", assigned
#   directly into .GlobalEnv (a plain top-level assignment in a testthat
#   file may instead land in a local test environment, which
#   pseudocount() - defined in .GlobalEnv via source()'s default envir -
#   would never see), then restored afterwards so the override does not
#   leak into other test files run later in the same session.

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)
# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)


# -------- load the demo script, with readline() mocked to answer "0.5" --------

assign("readline", function(prompt = "") "0.5", envir = .GlobalEnv)

capture.output(
  source(
    list.files(
      path = here(),
      pattern = "^demo_pseudocounts\\.R$",
      full.names = TRUE,
      recursive = TRUE
    )
  )
)




# -------- test 1: demo_otu has the expected known structure --------

test_that("demo_otu has the expected dimensions (3x3), values, and one zero per row", {
  expect_equal(dim(demo_otu), c(3, 3))
  expect_equal(colnames(demo_otu), c("OTU1", "OTU2", "OTU3"))
  expect_equal(unname(demo_otu[1, ]), c(10, 0, 5))
  expect_equal(unname(demo_otu[2, ]), c(2, 3, 0))
  expect_equal(unname(demo_otu[3, ]), c(0, 1, 20))
})


# -------- test 2: demo_result has no zeroes and matching dimensions --------

test_that("demo_result has no zeroes left and matches demo_otu's dimensions", {
  expect_equal(dim(demo_result), dim(demo_otu))
  expect_true(all(demo_result != 0))
})


# -------- test 3: non-zero entries are correctly converted to proportions --------

test_that("non-zero entries are correctly converted to proportions", {
  # row 1: library size 15, so the first entry should be 10 / 15
  expect_equal(demo_result[[1, 1]], 10 / 15, tolerance = 1e-8)
  # row 2: library size 5, so the second entry should be 3 / 5
  expect_equal(demo_result[[2, 2]], 3 / 5, tolerance = 1e-8)
})


# -------- test 4: pseudocount values match the expected formula --------

test_that("replaced zeroes match threshold_pct * (1 / lib_size), for threshold 0.5", {
  # row 1 (lib_size = 15): pseudocount = 0.5 * (1 / 15)
  expect_equal(demo_result[[1, 2]], 0.5 * (1 / 15), tolerance = 1e-8)
  # row 2 (lib_size = 5): pseudocount = 0.5 * (1 / 5)
  expect_equal(demo_result[[2, 3]], 0.5 * (1 / 5), tolerance = 1e-8)
})

rm(readline, envir = .GlobalEnv)