# -------- the begin of all the tests! --------
# this is script has the collection of all the libraries needed in this project
# its aim is to simplify the writing of tests by having all the libraries here 
# commented and ready to use

# testthat: unit testing framework for R
# https://cran.r-project.org/web/packages/testthat/index.html
library(testthat)
# tidyverse: useful to manage data (dplyr) and make nice plots (ggplot2)
# https://cran.r-project.org/web/packages/tidyverse/index.html
library(tidyverse)
# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# recall function to test
source(here("script", "zero_replacements", "zero_replacements_script.R"))

test_that("the column number is not the expected one", {
  
    
  # count number of columns of otu.69001.H
  before_filt_col_num <- ncol(otu.69001.H)
  # count column number after filtering
  after_filt_col_num  <- ncol(otu.filt)    
  # compare the two numbers
  expect_true(before_filt_col_num != after_filt_col_num)
  
})

