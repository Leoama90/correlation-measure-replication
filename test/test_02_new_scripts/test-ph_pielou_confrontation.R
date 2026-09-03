# test-ph_pielou_confrontation.R
#
# Purpose:
#   This script tests the essential features of the first section of
#   ph_pielou_confrontation.R (the sigma = 0.1-0.2 scenario), checking:
#       - ph_values stays within the expected 3.5-10.5 range
#       - pielou_list is a list of data.frame elements
#       - pielou_list has the expected number of elements (one per pH value)
#       - the ph column in every element of pielou_list is within 3.5-10.5
#
# Inputs:
#   - ph_pielou_confrontation.R (sourced below)
#
# Outputs:
#   - test results printed to the console (pass/fail for each check)
#
# Note:
#   The original script repeats the same logic five times (once per
#   tolerance scenario: 0.1-0.2, 0.3-0.4, 0.5-0.6, 0.7-0.8, 0.9-1.0),
#   with identical structure and only sigma_min/sigma_max changing. This
#   test file only checks the first (0.1-0.2) scenario's objects
#   (ph_values, pielou_list), since the other four are built the exact
#   same way and would only duplicate this coverage.
#   data_sim_ph_driven.R already has its own test script, so its
#   internal behaviour is not re-tested here.

# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)
# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)


# -------- load the script to test --------

source(
  list.files(
    path = here(),
    pattern = "^ph_pielou_confrontation\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)


# -------- test 1: ph_values stays within the expected range --------

test_that("ph_values is within 3.5 and 10.5", {
  expect_true(all(ph_values > 3.5))
  expect_true(all(ph_values < 10.5))
})


# -------- test 2: pielou_list is a list of data.frames --------

test_that("pielou_list is a list and its elements are data.frames", {
  expect_true(is.list(pielou_list))
  expect_true(all(sapply(pielou_list, is.data.frame)))
})


# -------- test 3: pielou_list has the expected number of elements --------

test_that("pielou_list has one element per value in ph_values", {
  expect_equal(length(pielou_list), length(ph_values))
})


# -------- test 4: ph values inside pielou_list are within range --------

test_that("the ph column in every element of pielou_list is between 3.5 and 10.5", {
  for (i in seq_along(pielou_list)) {
    expect_true(all(pielou_list[[i]]$ph > 3.5))
    expect_true(all(pielou_list[[i]]$ph < 10.5))
  }
})