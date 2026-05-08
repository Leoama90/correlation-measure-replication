library(testthat)
library(here)
library(vegan)

# -- load the simulation results to test --------
results <- readRDS(here("compositional_effects", "compositional_effects_01.rds"))

# -- test output structure --------
test_that("results is a data frame with the expected columns", {
  expect_s3_class(results, "data.frame")
  expect_named(results, c("d", "m", "err_l1", "err_clr", "pielou"))
})

test_that("results has the expected number of rows", {
  # dimensions has 5 values (5, 55, 105, 155, 205 -> seq(5,200,by=50) = 5 values)
  # magnification has 2 values -> 5 * 2 = 10 rows
  expected_rows <- length(seq(5, 200, by = 50)) * length(c(1, 2))
  expect_equal(nrow(results), expected_rows)
})

# -- test parameter coverage --------
test_that("all dimension values are present in results", {
  expect_setequal(unique(results$d), seq(5, 200, by = 50))
})

test_that("all magnification values are present in results", {
  expect_setequal(unique(results$m), c(1, 2))
})

test_that("each (d, m) combination appears exactly once", {
  # duplicated rows would indicate a parallelization merge error
  combos <- paste(results$d, results$m)
  expect_false(any(duplicated(combos)))
})

# -- test metric validity --------
test_that("ERR_L1 and ERR_CLR are non-negative", {
  # errors are absolute values, so they must never be negative
  expect_true(all(results$ERR_L1 >= 0))
  expect_true(all(results$ERR_CLR >= 0))
})

test_that("Pielou evenness is bounded between 0 and 1", {
  # Pielou index is a normalised entropy, strictly in [0, 1]
  expect_true(all(results$Pielou >= 0))
  expect_true(all(results$Pielou <= 1))
})

test_that("no NA values exist in the results", {
  # any NA signals a failed simulation run
  expect_false(anyNA(results))
})

# -- test output file --------
test_that("output .rds file exists in the expected location", {
  expect_true(file.exists(here("compositional_effects", "compositional_effects_01.rds")))
})