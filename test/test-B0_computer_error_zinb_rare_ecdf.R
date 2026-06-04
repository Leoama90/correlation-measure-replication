# doSNOW: parallel backend for foreach, supports progress bars via snow clusters
# https://cran.r-project.org/web/packages/doSNOW/index.html
library(doSNOW)

# foreach: parallel foreach loops
# https://cran.r-project.org/package=foreach
library(foreach)

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# progress: displays text progress bars for long-running loops
# https://cran.r-project.org/web/packages/progress/index.html
library(progress)

# SpiecEasi: package with spiec.easi and sparCC methods
# https://github.com/zdk123/SpiecEasi
library(SpiecEasi)

# tidyverse: useful to manage data (dplyr) and make nice plots (ggplot2)
# https://cran.r-project.org/web/packages/tidyverse/index.html
library(tidyverse)

# ToyModel: lightweight package with necessary function to generate metagenomics data
# https://github.com/Fuschi/ToyModel
library(ToyModel)

# VGAM: vector generalized linear models
# https://cran.r-project.org/package=VGAM
library(VGAM)

# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)

# load function to test
source(here("script", "sparsity_effects", "B0_computer_error_zinb_rare_ecdf.R"))

#test_that("..."{
  
  
#})

test_that("Column number after OTU filtering is not correct", {
  # prevalence and median filtering must reduce the number of OTUs
  expect_true(ncol(otu.filt) < ncol(otu.69001.H))
})

test_that("OTU prevalence filter was not applied correctly", {
  # every remaining OTU must be present in at least 33% of samples
  prevalences <- colSums(otu.filt > 0) / nrow(otu.filt)
  expect_true(all(prevalences >= 0.33))
})

test_that("OTU median filter was not applied correctly", {
  # every remaining OTU must have median non-zero abundance >= 5
  medians <- apply(otu.filt, 2, function(x) median(x[x > 0]))
  expect_true(all(medians >= 5))
})
