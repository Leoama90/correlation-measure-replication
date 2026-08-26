# test-generate_matrix_factors.R
#
# Purpose:
#   Tests the fundamental properties of generate_matrix_factors():
#   that the output is a valid correlation matrix (symmetric, unit
#   diagonal, positive semi-definite), that taxa in different groups
#   are exactly uncorrelated, that n_zeroes matches the actual zero
#   count in the matrix, and that invalid n_groups values are rejected.
#
# Inputs:
#   - generate_matrix_factors.R (sourced below)
#
# Outputs:
#   - test results printed to the console (pass/fail for each check)

# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)
# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# bring generate_matrix_factors() into scope
source(here("02_new_scripts", "generate_matrix_factors.R"))


# -------- test 1: output is a valid correlation matrix --------

test_that("mat is symmetric, has a unit diagonal, and is positive semi-definite", {
  result <- generate_matrix_factors(n = 20, n_groups = 4, seed = 42)
  mat <- result$mat
  
  # symmetry: mat and its transpose must be (numerically) identical
  expect_equal(mat, t(mat), tolerance = 1e-8)
  
  # unit diagonal, as required for a correlation matrix
  expect_equal(diag(mat), rep(1, 20), tolerance = 1e-8)
  
  # PSD: all eigenvalues must be non-negative (allowing tiny negative
  # noise from floating point arithmetic)
  eigenvalues <- eigen(mat, only.values = TRUE)$values
  expect_true(all(eigenvalues > -1e-8))
})


# -------- test 2: taxa in different groups are exactly uncorrelated --------

test_that("taxa assigned to different groups have exactly zero correlation", {
  result <- generate_matrix_factors(n = 20, n_groups = 4, seed = 42)
  mat <- result$mat
  groups <- result$groups
  
  # pick one pair of taxa known to be in different groups
  different_group_pair <- which(groups != groups[1])[1]
  expect_equal(mat[1, different_group_pair], 0, tolerance = 1e-8)
})


# -------- test 3: n_zeroes matches the actual zero count in mat --------

test_that("n_zeroes matches the number of (near-)zero off-diagonal entries in mat", {
  result <- generate_matrix_factors(n = 20, n_groups = 4, seed = 42)
  mat <- result$mat
  
  tol <- 1e-8
  actual_zeroes <- sum(abs(mat[upper.tri(mat)]) < tol) * 2
  expect_equal(result$n_zeroes, actual_zeroes)
})


# -------- test 4: n_groups = 1 produces no zeroes --------

test_that("n_groups = 1 leaves every taxon in the same group, with no zero entries", {
  result <- generate_matrix_factors(n = 10, n_groups = 1, seed = 42)
  expect_equal(result$n_zeroes, 0)
})


# -------- test 5: n_groups = n produces the identity matrix --------

test_that("n_groups = n makes every taxon its own group, giving the identity matrix", {
  result <- generate_matrix_factors(n = 10, n_groups = 10, seed = 42)
  expect_equal(result$mat, diag(10), tolerance = 1e-8)
})


# -------- test 6: invalid n_groups is rejected --------

test_that("n_groups outside [1, n] raises an error", {
  expect_error(generate_matrix_factors(n = 10, n_groups = 0))
  expect_error(generate_matrix_factors(n = 10, n_groups = 11))
})


# -------- test 7: seed makes the output reproducible --------

test_that("using the same seed produces identical output", {
  result_a <- generate_matrix_factors(n = 15, n_groups = 3, seed = 123)
  result_b <- generate_matrix_factors(n = 15, n_groups = 3, seed = 123)
  expect_equal(result_a$mat, result_b$mat)
  expect_equal(result_a$groups, result_b$groups)
})

# -------- test 8: check that one of the listed output is a matrix --------

test_that("mat is one of the element of the final output list and is a matrix",{
  # generate a dummy list to check that the <list_name>$mat is a matrix
  dum_mtrx <- generate_matrix_factors(20, 5, min_loading = 0.3, max_loading = 1.0, seed = 42)
  
  expect_true(is.matrix(dum_mtrx$mat))
})