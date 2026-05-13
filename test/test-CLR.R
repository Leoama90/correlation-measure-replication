# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)


# -- test CLR function for errors --------

test_that("there is a problem with the script",{
  
  source(here("script", "method_comparison", "CLR.R"))  
  
  # generate dummy variables for testing purpose
  x <- 5
  z <- matrix(1:25, 
              nrow = 5,
              ncol = 5)
  # check with non-matrix variable
  expect_error(CLR(x))
  # check with matrix variable
  expect_no_error(CLR(z))
  # verify transformation of matrix entries
  expect_true(all(CLR(z)) > 0)
  expect_true(is.double(CLR(z)))
  
})