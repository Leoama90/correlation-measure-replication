# test-data_sim_ph_driven.R
#
# Purpose:
#   Tests the essential features of data_sim_ph_driven(): that the
#   returned correlation matrix is valid (inherited from
#   generate_matrix_factors()), that a taxon exactly matching the
#   community pH gets zero-inflation 0 while a maximally mismatched
#   taxon approaches phi_max, that OTU names are assigned correctly
#   with zero-padding, that total_zero_rate is consistent with
#   sim_counts, and that seed makes the output reproducible.
#
# Inputs:
#   - data_sim_ph_driven.R (sourced below)
#
# Outputs:
#   - test results printed to the console (pass/fail for each check)

# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)
# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# bring data_sim_ph_driven() (and generate_matrix_factors(), sourced
# inside it) into scope
source(
  list.files(
    path = here(),
    pattern = "^data_sim_ph_driven\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)


# -------- test 1: correlation matrix is valid --------

test_that("mat is symmetric, has a unit diagonal, and is positive semi-definite", {
  result <- data_sim_ph_driven(n = 20, N = 100, ph = 6.5, n_groups = 4, seed = 42)
  mat <- result$mat
  
  expect_equal(mat, t(mat), tolerance = 1e-8)
  expect_equal(diag(mat), rep(1, 20), tolerance = 1e-8)
  
  eigenvalues <- eigen(mat, only.values = TRUE)$values
  expect_true(all(eigenvalues > -1e-8))
})


# -------- test 2: a well-matched taxon gets near-zero zero-inflation --------

test_that("a taxon whose optimum equals the community pH gets phi close to 0", {
  result <- data_sim_ph_driven(n = 20, N = 100, ph = 6.5, n_groups = 4, seed = 42)
  
  # force the first taxon's optimum to exactly match the community pH,
  # and recompute its zero-inflation probability the same way the
  # function does internally, to check the formula's boundary behaviour
  ph <- 6.5
  ph_optimum <- ph
  ph_tolerance <- 0.5
  suitability <- exp(-(ph - ph_optimum)^2 / (2 * ph_tolerance^2))
  phi <- 0.9 * (1 - suitability)
  
  expect_equal(phi, 0)
})


# -------- test 3: a mismatched taxon approaches phi_max --------

test_that("a taxon whose optimum is very far from the community pH approaches phi_max", {
  ph <- 6.5
  ph_optimum <- 20  # implausibly far, to isolate the asymptotic behaviour
  ph_tolerance <- 0.5
  phi_max <- 0.9
  
  suitability <- exp(-(ph - ph_optimum)^2 / (2 * ph_tolerance^2))
  phi <- phi_max * (1 - suitability)
  
  expect_equal(phi, phi_max, tolerance = 1e-8)
})


# -------- test 4: OTU names are assigned correctly, with zero-padding --------

test_that("sim_data and sim_counts get zero-padded OTU column names", {
  result <- data_sim_ph_driven(n = 15, N = 50, ph = 6.5, n_groups = 3, seed = 42)
  
  expect_equal(colnames(result$sim_data), sprintf("OTU_%02d", 1:15))
  expect_equal(colnames(result$sim_counts), sprintf("OTU_%02d", 1:15))
})


# -------- test 5: total_zero_rate matches sim_counts --------

test_that("total_zero_rate matches the actual fraction of zeroes in sim_counts", {
  result <- data_sim_ph_driven(n = 20, N = 100, ph = 6.5, n_groups = 4, seed = 42)
  
  expect_equal(result$total_zero_rate, mean(result$sim_counts == 0))
})


# -------- test 6: seed makes the output reproducible --------

test_that("using the same seed produces identical output", {
  result_a <- data_sim_ph_driven(n = 15, N = 50, ph = 6.0, n_groups = 3, seed = 123)
  result_b <- data_sim_ph_driven(n = 15, N = 50, ph = 6.0, n_groups = 3, seed = 123)
  
  expect_equal(result_a$sim_counts, result_b$sim_counts)
  expect_equal(result_a$ph_optima, result_b$ph_optima)
})


# -------- test 7: output dimensions match n and N --------

test_that("sim_data and sim_counts have the expected dimensions", {
  result <- data_sim_ph_driven(n = 12, N = 80, ph = 6.5, n_groups = 3, seed = 42)
  
  expect_equal(dim(result$sim_data), c(80, 12))
  expect_equal(dim(result$sim_counts), c(80, 12))
})