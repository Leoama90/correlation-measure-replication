# SpiecEasi: package with spiec.easi and sparCC methods
# https://github.com/zdk123/SpiecEasi
library(SpiecEasi)
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


# ---- Fixture: minimal params used across multiple tests ----

# Small ZINB parameter set representing a single OTU pair
test_params <- data.frame(
  munb  = c(10 , 15 ),
  size  = c(2  , 3  ),
  pstr0 = c(0.2, 0.4)
)


# -------- toy_model returns the expected structure --------

test_that("toy_model returns NorTA matrix and cor_normal with correct dimensions", {
  
  result <- toy_model(n     = 50,
                      cor   = diag(2),
                      M     = 1,
                      qdist = VGAM::qzinegbin,
                      param = test_params)
  
  # NorTA should be a 50 x 2 non-negative integer matrix
  expect_equal(dim(result$NorTA), c(50, 2))
  expect_true(all(result$NorTA >= 0))
  
  # cor_normal should be a 2x2 symmetric matrix with ones on the diagonal
  expect_equal(dim(result$cor_normal),  c(2, 2))
  expect_equal(diag(result$cor_normal), c(1, 1))
})


# -------- fitdistr extracts the three expected ZINB parameters --------

test_that("fitdistr returns munb, size, and pstr0 for a ZINB-distributed vector", {
  
  set.seed(42)
  # Generate a small ZINB sample to fit
  x <- VGAM::rzinegbin(n     = 200, 
                       munb  = 10, 
                       size  = 2, 
                       pstr0 = 0.3)
  
  fitted <- SpiecEasi::fitdistr(as.numeric(x), "zinegbin")$par
  
  # All three parameter names must be present
  expect_true(all(c("munb", "size", "pstr0") %in% names(fitted)))
  
  # Fitted values should be positive and pstr0 in [0, 1)
  expect_true(fitted["munb"]  > 0)
  expect_true(fitted["size"]  > 0)
  expect_true(fitted["pstr0"] >= 0 && fitted["pstr0"] < 1)
})


# -------- CLR transformation produces zero-sum rows --------

test_that("clr() produces rows that sum to zero (up to floating-point tolerance)", {
  
  set.seed(42)
  # Simulate a small count matrix (rows = samples, cols = OTUs)
  counts <- toy_model(n     = 30,
                      cor   = diag(5),
                      M     = 1,
                      qdist = VGAM::qzinegbin,
                      param = data.frame(
                      munb  = rep(10, 5),
                      size  = rep(2,  5),
                      pstr0 = rep(0.1, 5)
                      ))$NorTA
  
  clr_matrix <- ToyModel::clr(counts)
  
  # Each row of the CLR-transformed matrix must sum to ~0
  row_sums <- rowSums(clr_matrix)
  expect_true(all(abs(row_sums) < 1e-10))
})


# -------- injection of the OTU pair into background does not alter other columns --------

test_that("replacing columns 25 and 125 leaves all other columns unchanged", {
  
  set.seed(42)
  bg   <- toy_model(n     = 100,
                    cor   = diag(200),
                    M     = 1,
                    qdist = VGAM::qzinegbin,
                    param = data.frame(
                    munb  = rep(10,  200),
                    size  = rep(2,   200),
                    pstr0 = rep(0.2, 200)
                  ))
  
  pair <- toy_model(n     = 100,
                    cor   = 0.6,
                    M     = 1,
                    qdist = VGAM::qzinegbin,
                    param = test_params)
  
  modified        <- bg$NorTA
  modified[, 25]  <- pair$NorTA[, 1]
  modified[, 125] <- pair$NorTA[, 2]
  
  # Columns other than 25 and 125 must be identical to the original background
  other_cols <- setdiff(1:200, c(25, 125))
  expect_equal(modified[, other_cols], bg$NorTA[, other_cols])
})

