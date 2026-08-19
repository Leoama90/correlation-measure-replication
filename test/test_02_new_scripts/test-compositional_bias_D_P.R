# test-compositional_bias_D_P.R
#
# Purpose:
#   Tests the essential helper functions defined in
#   compositional_bias_D_P.R: l1_normalize(), clr_transform(),
#   mae_matrix(), and tune_diversity(). Does not test the full grid
#   simulation itself (too slow to run on every test invocation).
#
# Inputs:
#   - a small test matrix created inside this script
#
# Outputs:
#   - test results printed to the console (pass/fail for each check)
#
# Note:
#   compositional_bias_D_P.R runs its entire grid simulation as a
#   top-level side effect when sourced (not wrapped in a guarded
#   main()/if block), so a plain source() here would re-run the full
#   15x15 grid before any test executes. To avoid this, the helper
#   functions are redefined locally, identical to the originals. This
#   duplicates their bodies (a known trade-off, flagged here rather
#   than hidden) until compositional_bias_D_P.R is refactored to guard
#   its top-level execution (e.g. behind a run_analysis() wrapper or
#   an interactive() check).

# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)
# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# bring pielou_ind() into scope, needed by tune_diversity() below
source(
  list.files(
    path = here(),
    pattern = "^pielou_ind\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)


# -------- local, non-executing copies of the helper functions --------

# identical to l1_normalize() in compositional_bias_D_P.R
l1_normalize <- function(x) {
  sweep(x, 1, rowSums(x), "/")
}

# identical to clr_transform() in compositional_bias_D_P.R
clr_transform <- function(x) {
  log_x <- log(x)
  log_x - rowMeans(log_x)
}

# identical to mae_matrix() in compositional_bias_D_P.R
mae_matrix <- function(R_est, R_true) {
  mean(abs(R_est - R_true))
}

# identical to tune_diversity() in compositional_bias_D_P.R
tune_diversity <- function(x, target_p, tol = 0.005, max_log_factor = 8) {
  baseline_p <- pielou_ind(x)
  
  if (target_p > baseline_p) {
    warning("target_p (", round(target_p, 3), ") exceeds the dataset's ",
            "baseline diversity (", round(baseline_p, 3), "); skipping.")
    return(list(data = NULL, achieved_p = NA))
  }
  
  diversity_gap <- function(log_factor) {
    x_scaled <- x
    x_scaled[, 1] <- x_scaled[, 1] * (10^log_factor)
    pielou_ind(x_scaled) - target_p
  }
  
  root <- uniroot(diversity_gap, interval = c(0, max_log_factor),
                  tol = tol)$root
  
  x_tuned <- x
  x_tuned[, 1] <- x_tuned[, 1] * (10^root)
  
  list(data = x_tuned, achieved_p = pielou_ind(x_tuned))
}


# -------- shared test data --------

test_matrix <- matrix(
  c(
    3, 4, 3,
    1, 2, 9,
    7, 5, 1
  ),
  nrow = 3, ncol = 3, byrow = TRUE
)


# -------- tests for l1_normalize() --------

test_that("l1_normalize returns the expected matrix", {
  l1_test_matrix <- l1_normalize(test_matrix)
  
  # precalculated matrix that is expected after applying the l1_normalization
  expected_matrix <- matrix(
    c(
      3 / 10, 4 / 10, 3 / 10,
      1 / 12, 2 / 12, 9 / 12,
      7 / 13, 5 / 13, 1 / 13
    ),
    nrow = 3, ncol = 3, byrow = TRUE
  )
  
  expect_equal(l1_test_matrix, expected_matrix, tolerance = 1e-12)
})

test_that("l1_normalize preserves the matrix dimensions", {
  l1_test_matrix <- l1_normalize(test_matrix)
  expect_equal(dim(l1_test_matrix), dim(test_matrix))
})

test_that("each row of the L1-normalized matrix sums to one", {
  l1_test_matrix <- l1_normalize(test_matrix)
  expect_equal(rowSums(l1_test_matrix), rep(1, nrow(test_matrix)), tolerance = 1e-12)
})


# -------- tests for clr_transform() --------

test_that("clr_transform preserves the matrix dimensions", {
  clr_test_matrix <- clr_transform(test_matrix)
  expect_equal(dim(clr_test_matrix), dim(test_matrix))
})

test_that("each row of the CLR-transformed matrix sums to (approximately) zero", {
  # CLR subtracts the row-wise mean of the logs, so each row's values
  # are centered around 0 by construction; small floating point error
  # is expected, hence a tolerance rather than exact equality to 0
  clr_test_matrix <- clr_transform(test_matrix)
  row_sums <- rowSums(clr_test_matrix)
  expect_true(all(abs(row_sums) < 1e-10))
})


# -------- tests for mae_matrix() --------

test_that("mae_matrix returns 0 for two identical matrices", {
  m <- diag(4)
  expect_equal(mae_matrix(m, m), 0)
})

test_that("mae_matrix returns the correct value for a known difference", {
  # every entry differs by exactly 0.5, so the mean absolute error is 0.5
  R_est <- matrix(0.5, nrow = 3, ncol = 3)
  R_true <- matrix(0, nrow = 3, ncol = 3)
  expect_equal(mae_matrix(R_est, R_true), 0.5)
})


# -------- tests for tune_diversity() --------

test_that("tune_diversity achieves a target diversity within tolerance", {
  set.seed(1)
  x <- matrix(runif(50, min = 1, max = 10), nrow = 10, ncol = 5)
  
  target_p <- 0.7
  tuned <- tune_diversity(x, target_p, tol = 0.005)
  
  expect_false(is.null(tuned$data))
  expect_equal(tuned$achieved_p, target_p, tolerance = 0.01)
})

test_that("tune_diversity returns NULL data when target_p exceeds baseline diversity", {
  set.seed(1)
  # a perfectly even dataset has baseline diversity very close to 1;
  # requesting a target above 1 is always unachievable
  x <- matrix(1, nrow = 10, ncol = 5)
  
  expect_warning(
    tuned <- tune_diversity(x, target_p = 1.5),
    "exceeds the dataset's baseline diversity"
  )
  expect_true(is.null(tuned$data))
  expect_true(is.na(tuned$achieved_p))
})