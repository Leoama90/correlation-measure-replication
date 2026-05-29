# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)

# tidyverse: useful to manage data (dplyr) and make nice plots (ggplot2)
# https://cran.r-project.org/web/packages/tidyverse/index.html
library(tidyverse)

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# SpiecEasi: package with spiec.easi and sparCC methods
# https://github.com/zdk123/SpiecEasi
library(SpiecEasi)

# propr: package with proportionality rho method
# https://github.com/tpq/propr
library(propr)

# ggpubr: nice plots based on ggplot2
# https://cran.r-project.org/web/packages/ggpubr/index.html
library(ggpubr)

# qualpalr: generate distinct qualitative color palette
# https://cran.r-project.org/web/packages/qualpalr/index.html
library(qualpalr)

# cowplot: draw ggplot2 in new figures
# https://cran.r-project.org/web/packages/cowplot/index.html
library(cowplot)

# vegan: used to elaborates shannon entropy
# https://cran.r-project.org/web/packages/vegan/index.html
library(vegan)

# lightweight package with necessary function to generate metagenomics data
# https://github.com/Fuschi/ToyModel
library(ToyModel)

# Imputation of zeros, left-censored and missing values in compositional data
# https://cran.r-project.org/web/packages/zCompositions/index.html
library(zCompositions)

# load function to test
source(here("script", "method_comparison", "methods_comparison.R"))

test_that("The row number is not correct",{
  
  # get otu row number
  otu_rown    <- nrow(otu)
  # get otu.69001.H row number
  otu.69_rown <- nrow(otu.69001.H)
  # expect less row number after first filter
  expect_true(otu_rown > otu.69_rown)
  
})

# -------- checks matrixes for Spiec-easi MB --------

test_that("Matrix has not yet been symmetrized",{
  
  # testing the non-symmetrized matrix is not symmetric
  expect_false(isSymmetric.matrix(beta.mb))
  
})

test_that("Matrix is not symmetric",{
  
  # testing the symmetry of adj.mb matrix
  expect_true(isSymmetric.matrix(adj.mb))
  
})

test_that("The selected matrixes entries are not equal 1 or -1",{
  
  expect_equal(adj.mb[adj.mb >  0] <-   1,  1)
  expect_equal(adj.mb[adj.mb <  0] <-  -1, -1)
  
})

# -------- checks for SparCC --------

test_that("SparCC row and column names match filtered OTU names", {
  # Use identical() to compare two vectors, not == which returns a vector
  expect_true(identical(colnames(res.cc), rownames(res.cc)))
  expect_true(identical(colnames(res.cc), colnames(otu.filt)))
})








