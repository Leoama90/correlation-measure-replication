# testthat: allows to write tests for the other scripts
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)
# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)
# vegan: tools for descriptive community ecology (diversity indices, ordination)
# https://cran.r-project.org/web/packages/vegan/index.html
library(vegan)


# -------- load the simulation results to test --------

# search recursively for the target .rds file and load it directly
results <- readRDS(list.files(path       = here(),
                              pattern    = "compositional_effects_01.rds",
                              full.names = TRUE,
                              recursive  = TRUE))


# -------- test output structure --------

test_that("results is not a data frame with the expected columns", {
  # confirm the object inherits the data.frame class
  expect_s3_class(results, "data.frame")
  # confirm column names match the expected output schema exactly
  expect_named(results, c("d", "m", "ERR_L1", "ERR_CLR", "Pielou"))
})
test_that("results has not the expected number of rows", {
  # dimensions has 5 values (5, 55, 105, 155, 205 -> seq(5,200,by=50) = 5 values)
  # magnification has 2 values -> 5 * 2 = 10 rows
  expected_rows <- length(seq(5, 200, by = 50)) * length(c(1, 2))
  expect_equal(nrow(results), expected_rows)
})


# -------- test output file --------

test_that("output .rds file does not exists in the expected location", {
  # search for the dataset path then keeps only the file .rds file name
  exp_file_name <- basename(list.files(path       = here(),
                                       pattern    = "compositional_effects_01.rds",
                                       full.names = TRUE,
                                       recursive  = TRUE))
  file_name <- "compositional_effects_01.rds"
  
  expect_equal(exp_file_name, file_name)
})


# -------- test parameter coverage --------

test_that("all dimension values are not present in results", {
  # every value in seq(5, 200, by = 50) must appear at least once
  expect_setequal(unique(results$d), seq(5, 200, by = 50))
})
test_that("all magnification values are not present in results", {
  # both magnification levels must be represented in the output
  expect_setequal(unique(results$m), c(1, 2))
})
test_that("each (d, m) combination appears more than once", {
  # duplicated rows would indicate a parallelization merge error
  duos <- paste(results$d, results$m)
  expect_false(any(duplicated(duos)))
})


# -------- test metric validity --------

test_that("ERR_L1 and ERR_CLR are negative", {
  # errors are absolute values, so they must never be negative
  expect_true(all(results$ERR_L1  >= 0))
  expect_true(all(results$ERR_CLR >= 0))
})
test_that("Pielou evenness is not bounded between 0 and 1!", {
  # Pielou index is a normalised entropy, strictly in [0, 1]
  expect_true(all(results$Pielou >= 0))
  expect_true(all(results$Pielou <= 1))
})
test_that("there is an unexpected NA value in the results!", {
  # any NA signals a failed simulation run
  expect_false(anyNA(results))
})