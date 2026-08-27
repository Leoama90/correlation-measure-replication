# test-demo_data_sim_ph_driven.R
#
# Purpose:
#   This script tests the fundamental features of
#   demo_data_sim_ph_driven.R, checking:
#       - demo_data has the expected list structure with all 7 elements
#       - sim_counts and sim_data have the dimensions implied by the
#         demo's parameters (n = 20, N = 50)
#       - total_zero_rate matches the actual fraction of zeroes in
#         sim_counts
#       - best_taxon has the lowest phi_per_taxon and worst_taxon has
#         the highest, consistent with what the script prints
#       - sim_counts columns are zero-padded OTU names matching n
#
# Inputs:
#   - demo_data_sim_ph_driven.R (sourced below)
#
# Outputs:
#   - test results printed to the console (pass/fail for each check)
#
# Note:
#   This script DOES source demo_data_sim_ph_driven.R directly, since
#   the demo itself is lightweight (n = 20, N = 50) and does not re-run
#   any slow simulation, unlike some of the other test scripts in this
#   project. capture.output() is used to suppress its cat()/print()
#   messages so the test output stays clean.

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)
# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)


# -------- load the demo script --------

capture.output(
  source(
    list.files(
      path = here(),
      pattern = "^demo_data_sim_ph_driven\\.R$",
      full.names = TRUE,
      recursive = TRUE
    )
  )
)


# -------- test 1: demo_data has the expected structure --------

test_that("demo_data has all 7 expected elements", {
  expect_named(
    demo_data,
    c("mat", "groups", "ph_optima", "phi_per_taxon", "sim_data", "sim_counts", "total_zero_rate")
  )
})


# -------- test 2: dimensions match the demo's parameters (n = 20, N = 50) --------

test_that("sim_counts and sim_data have dimensions matching n = 20, N = 50", {
  expect_equal(dim(demo_data$sim_counts), c(50, 20))
  expect_equal(dim(demo_data$sim_data), c(50, 20))
})


# -------- test 3: total_zero_rate matches sim_counts --------

test_that("total_zero_rate matches the actual fraction of zeroes in sim_counts", {
  expect_equal(demo_data$total_zero_rate, mean(demo_data$sim_counts == 0))
})


# -------- test 4: best/worst taxon match phi_per_taxon extremes --------

test_that("best_taxon and worst_taxon correspond to the min/max of phi_per_taxon", {
  expect_equal(best_taxon, which.min(demo_data$phi_per_taxon))
  expect_equal(worst_taxon, which.max(demo_data$phi_per_taxon))
  # the best match must have a lower (or equal) zero-inflation probability
  # than the worst match
  expect_true(demo_data$phi_per_taxon[best_taxon] <= demo_data$phi_per_taxon[worst_taxon])
})


# -------- test 5: OTU column names are zero-padded and match n --------

test_that("sim_counts has zero-padded OTU column names matching n = 20", {
  expect_equal(colnames(demo_data$sim_counts), sprintf("OTU_%02d", 1:20))
})