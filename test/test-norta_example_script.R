# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# Collection of packages for data wrangling and visualization
# https://cran.r-project.org/web/packages/tidyverse/index.html
library(tidyverse)

# NorTA (Normal To Anything): simulates correlated count data via Gaussian copula
# https://github.com/nome-repo/ToyModel
library(ToyModel)

# Quantile functions for zero-inflated and extended distributions (qzinegbin)
# https://cran.r-project.org/web/packages/VGAM/index.html
library(VGAM)

# load script to test
source(here("script", "NorTa example", "norta_example_script.R"))


# -------- test for the R matrix --------

# test that R matrix is symmetric
test_that("R is not symmetric",{
  
  expect_true(isSymmetric.matrix(R))
})

# test that the matrix entry [2, 4] has the proper value
test_that("Matrix entry [2, 4] is different from -0.9",{
  expect_true(R[2, 4] == -0.9)
})

# testing that [4, 2] entry is equal to [2, 4] entry
test_that("Matrix entries in positions [4, 2] and [2, 4] are not identical",{
  expect_equal(R[4, 2], R[2, 4])
})




