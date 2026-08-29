# test-demo_pielou_ind.R
#
# Purpose:
#   This script tests the fundamental features of demo_pielou_ind.R,
#   checking:
#       - mean_pielou matches pielou_ind() called on the same data
#         without per_sample
#       - per_sample_pielou has one value per sample (length N = 25)
#         and matches pielou_ind() called with per_sample = TRUE
#       - all per-sample Pielou values are bounded within [0, 1]
#       - the demo dataset (piel_demo_data) actually contains zeroes,
#         so the zero-handling behaviour described by the demo is
#         genuinely being exercised
#
# Inputs:
#   - demo_pielou_ind.R (sourced below)
#
# Outputs:
#   - test results printed to the console (pass/fail for each check)
#
# Note:
#   This script DOES source demo_pielou_ind.R directly, since the demo
#   itself is lightweight (n = 10, N = 25) and does not re-run any slow
#   simulation. capture.output() is used to suppress its cat()/print()
#   messages so the test output stays clean.

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)
# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)


# -------- load the demo script to test --------

capture.output(
  source(
    list.files(
      path = here(),
      pattern = "^demo_pielou_ind\\.R$",
      full.names = TRUE,
      recursive = TRUE
    )
  )
)


# -------- test 1: mean_pielou matches a direct call to pielou_ind() --------

test_that("mean_pielou matches pielou_ind(piel_demo_data)", {
  expect_equal(mean_pielou, pielou_ind(piel_demo_data))
})


# -------- test 2: per_sample_pielou has one value per sample --------

test_that("per_sample_pielou has length N = 25 and matches a direct call", {
  expect_length(per_sample_pielou, 25)
  expect_equal(per_sample_pielou, pielou_ind(piel_demo_data, per_sample = TRUE))
})


# -------- test 3: per-sample values are bounded within [0, 1] --------

test_that("all per-sample Pielou values are bounded within [0, 1]", {
  expect_true(all(per_sample_pielou >= 0))
  expect_true(all(per_sample_pielou <= 1))
})


# -------- test 4: the demo dataset genuinely contains zeroes --------

test_that("piel_demo_data contains at least one zero", {
  expect_true(any(piel_demo_data == 0))
})