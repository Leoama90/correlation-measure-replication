# data_analysis_and_statics
#
# This function takes as input a data set (tibble, data.frame, or matrix)
# and prints out the number of columns, number of rows, and the count
# and percentage of zero values it contains
# it is intended to gather some data useful to build further code
#
# Input:
#   - a tibble, data.frame, or matrix with numeric columns
#
# Output:
#   - printed column count, row count, number of zeroes and their
#     percentage; invisibly returns these same values as a list
#
# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)
# tidyverse: data science toolkit (dplyr, ggplot2, tidyr, etc.)
# [https://tidyverse.org](https://tidyverse.org)
library(tidyverse)


# -------- body of the function --------

datasum <- function(x) {
  
  # give the dataset the first four letters of its original name
  dataset_name <- substr(deparse(substitute(x)), 1, 4)
  
  # transform the input in a tibble
  if (!is_tibble(x)) {
    x <- as_tibble(x)
  } else {
    cat("\n")
    cat("this is already a tibble, no need to transform it")
    cat("\n")
  }
  
  # gather the number of column, number of rows and the zeroes percentage
  cat("\n")
  cat("the column number of", dataset_name, "is", ncol(x), "\n")
  cat("the row number of", dataset_name, "is", nrow(x), "\n")
  cat("\n")
  
  
  # -------- zero counts --------
  
  # total count of zeros across the whole dataset
  total_zeros <- sum(x == 0, na.rm = TRUE)
  cat("the total number of zeroes in", dataset_name, "is", total_zeros, "\n")
  
  zero_rate <- round(total_zeros / (ncol(x) * nrow(x)) * 100, 2)
  cat("the percentage of zeroes in", dataset_name, "is:", zero_rate, "%\n")
  
  
  # -------- return the invisible values --------
  
  invisible(list(
    n_col = ncol(x),
    n_row = nrow(x),
    total_zeros = total_zeros,
    zero_rate = zero_rate
  ))
}