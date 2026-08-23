# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)

# tidyverse: useful to manage data (dplyr) and make nice plots (ggplot2)
# https://cran.r-project.org/web/packages/tidyverse/index.html
library(tidyverse)

# ToyModel: lightweight package with necessary function to generate metagenomics data
# https://github.com/Fuschi/ToyModel
library(ToyModel)

# VGAM: vector generalized linear models
# https://cran.r-project.org/package=VGAM
library(VGAM)


# -------- Fixture: minimal params used across multiple tests --------

# Small ZINB parameter set representing a single OTU pair
test_params <- data.frame(
  munb  = c(10, 15),
  size  = c(2, 3),
  pstr0 = c(0.2, 0.4)
)


# -------- toy_model returns the expected structure --------

test_that("toy_model returns NorTA matrix and cor_normal with correct dimensions", {
  result <- toy_model(
    n = 50,
    cor = diag(2),
    M = 1,
    qdist = VGAM::qzinegbin,
    param = test_params
  )

  # NorTA should be a 50 x 2 non-negative integer matrix
  expect_equal(dim(result$NorTA), c(50, 2))
  expect_true(all(result$NorTA >= 0))

  # cor_normal should be a 2x2 symmetric matrix with ones on the diagonal
  expect_equal(dim(result$cor_normal), c(2, 2))
  expect_equal(diag(result$cor_normal), c(1, 1))
})
