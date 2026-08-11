# test-demo_datasummary.R
#
# Purpose:
#   This script tests the fundamental features of demo_datasummary.R,
#   checking that the dummy tibble it builds has the expected known
#   structure, and that datasum() correctly reports its statistics both
#   when run on the tibble directly (determinant/eigenvalue undefined)
#   and on its matrix version (determinant/eigenvalue computed, and
#   both exactly 0 due to the all-zero column).
#
# Inputs:
#   - demo_datasummary.R (sourced below)
#
# Outputs:
#   - test results printed to the console (pass/fail for each check)

# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)
# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# source demo_datasummary.R: this builds dum_tibble and runs datasum()
# on it (twice), leaving dum_tibble available in this environment.
# capture.output() suppresses its cat()/print() messages so the test
# output stays clean.
capture.output(
  source(
    list.files(
      path = here(),
      pattern = "^demo_datasummary\\.R$",
      full.names = TRUE,
      recursive = TRUE
    )
  )
)


# -------- test 1: dum_tibble has the expected known structure --------

test_that("dum_tibble does not have the expected columns and values", {
  expect_named(dum_tibble, c("a", "b", "c", "d", "e"))
  expect_equal(dum_tibble$a, 1:5)
  expect_equal(dum_tibble$b, 2:6)
  expect_equal(dum_tibble$c, (1:5)^2 - 2)
  expect_equal(dum_tibble$d, rep(0, 5))
  expect_equal(dum_tibble$e, dum_tibble$b - dum_tibble$a)
})


# -------- test 2: datasum() on the tibble itself has no determinant --------

test_that("datasum() on dum_tibble directly does not return NA determinant and eigenvalue", {
  out <- capture.output(res <- datasum(dum_tibble))
  expect_true(is.na(res$det_val))
  expect_true(is.na(res$min_eigenvalue))
  # 5 of the 25 cells belong to the all-zero column d
  expect_equal(res$total_zeroes, 5)
})


# -------- test 3: datasum() on the matrix version computes the determinant --------

test_that("datasum() on as.matrix(dum_tibble) should return a zero determinant and eigenvalue", {
  out <- capture.output(res <- datasum(as.matrix(dum_tibble)))
  # the all-zero column d makes the matrix singular by construction
  expect_equal(res$det_val, 0)
})
