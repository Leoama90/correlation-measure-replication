# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)

# load function to test
source(here("script", "sparsity_effects", "B0_computer_error_zinb_rare_ecdf.R"))

#test_that("..."{
  
  
#})

test_that("Column number after OTU filtering is not correct", {
  # prevalence and median filtering must reduce the number of OTUs
  expect_true(ncol(otu_filt) < ncol(otu_69001_H))
})

test_that("OTU prevalence filter was not applied correctly", {
  # every remaining OTU must be present in at least 33% of samples
  prevalences <- colSums(otu_filt > 0) / nrow(otu_filt)
  expect_true(all(prevalences >= 0.33))
})

test_that("OTU median filter was not applied correctly", {
  # every remaining OTU must have median non-zero abundance >= 5
  medians <- apply(otu_filt, 2, function(x) median(x[x > 0]))
  expect_true(all(medians >= 5))
})
