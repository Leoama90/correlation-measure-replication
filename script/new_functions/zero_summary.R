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
#
#' Summarize a dataset's dimensions and zero content
#'
#' @param x a data.frame, tibble, or matrix with numeric columns.
#'
#' @return An invisible list containing: n_col (number of columns),
#'   n_row (number of rows), total_zeros (count of zero values), min_val 
#'   (minimum value of the matrix), max_val (maximum value of the matrix),
#'   det_val (determinant value of the matrix)
#'   and zero_rate (percentage of zero values in the dataset).
#'
#' @examples
#' datasum(matrix(runif(16), nrow = 4, ncol = 4))


# -------- body of the function --------
datasum <- function(x) {
  
  # give the dataset the first four letters of its original name
  dataset_name <- substr(deparse(substitute(x)), 1, 4)
  
  # check, before any conversion, whether x is a square numeric matrix
  is_square_matrix <- is.matrix(x) && is.numeric(x) && nrow(x) == ncol(x)
  
  # if x qualifies, compute matrix-specific properties while it is still a matrix
  if (is_square_matrix) {
    min_val <- min(x)
    max_val <- max(x)
    det_val <- det(x)
  }
  
  # transform the input in a tibble
  if (!is_tibble(x)) {
    x <- as_tibble(x)
  } else {
    cat("\n")
    cat("this is already a tibble, no need to transform it")
    cat("\n")
  }
  
  # print a header line with the given dataset_name
  cat("\n---", dataset_name, "---\n")
  # print the number of columns
  cat("number of columns:", ncol(x), "\n")
  # print the number of rows
  cat("number of rows:   ", nrow(x), "\n")
  # print the maximum value in the matrix
  cat("maximum value:    ", max_val, "\n")
  # print the minimum value in the matrix
  cat("minimum value:    ", min_val, "\n")
  # print the determinant of the matrix
  cat("determinant value ", det_val, "\n")
  
  
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
    min_val = min_val,
    max_val = max_val,
    det_val = det_val,
    total_zeros = total_zeros,
    zero_rate = zero_rate
  ))
}