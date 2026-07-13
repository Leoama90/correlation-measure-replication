# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)
# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)


# ---- test TRIU function for errors ------------------------------------------

test_that("TRIU does not raise an error on invalid input", {
  
  # load function from script
  source(here("script", "method_comparison", "TRIU.R"))
  
  # generate non-matrix dummy variable
  x <- 5
  # generate non-numeric matrix
  z_char <- matrix(c("a","b","c","d"), nrow = 2, ncol = 2)
  # generate non-symmetric matrix
  z_nonsym <- matrix(1:9, nrow = 3, ncol = 3)
  
  # check non-matrix input raises an error
  expect_error(TRIU(x))
  # check non-numeric matrix raises an error
  expect_error(TRIU(z_char))
  # check non-symmetric matrix raises an error
  expect_error(TRIU(z_nonsym))
})


# ---- test TRIU function output ----------------------------------------------

test_that("TRIU does not return the expected upper triangular vector", {
  
  # load function from script
  source(here("script", "method_comparison", "TRIU.R"))

  # generate symmetric dummy matrix
  z <- matrix(c(1, 2, 3,
                2, 4, 5,
                3, 5, 6),
              nrow = 3, ncol = 3)
  
  # apply function to dummy matrix
  result <- TRIU(z)
  
  # check output is a vector
  expect_true(is.vector(result))
  # check output is numeric
  expect_true(is.numeric(result))
  # check output length equals number of upper triangular elements (excluding diagonal)
  expect_equal(length(result), (nrow(z) * (nrow(z) - 1)) / 2)
  # check output values match expected upper triangular values
  expect_equal(result, c(2, 3, 5))
})