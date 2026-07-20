# datasummary.R
#
# This function takes as input a data set (tibble, data.frame, list, or matrix)
# and prints out the number of columns, number of rows, and the count
# and percentage of zero values it contains
# it is intended to gather some data useful to build further code
#
# Input:
#   - a tibble, data.frame, list, or matrix with numeric columns
#     (does not need to be square; determinant is only computed
#     when the input is a square numeric matrix)
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
#' @param x a data.frame, tibble, list, or matrix with numeric columns.
#'   Does not need to be square.
#'
#' @return An invisible list containing: n_col (number of columns),
#'   n_row (number of rows), total_zeroes (count of zero values), min_val
#'   (minimum value of the dataset), max_val (maximum value of the dataset),
#'   det_val (determinant value, only if x is a square numeric matrix)
#'   and zero_rate (percentage of zero values in the dataset).
#'
#' @examples
#' datasum(matrix(runif(16), nrow = 4, ncol = 4))
#' datasum(matrix(runif(12), nrow = 3, ncol = 4))
#' datasum(list(a = 1:5, b = 6:10))
#' @export


# -------- body of the function --------

datasum <- function(x) {
  # give the dataset the first four letters of its original name
  dataset_name <- substr(deparse(substitute(x)), 1, 4)

  # check, before any conversion, whether x is a square numeric matrix
  # (needed for the determinant, which is only defined for square matrices)
  is_square_matrix <- is.matrix(x) && is.numeric(x) && nrow(x) == ncol(x)

  # default value for the determinant: only defined for square matrices
  det_val <- NA

  # if x qualifies, compute the determinant while it is still a matrix
  if (is_square_matrix) {
    det_val <- det(x)
  } else {
    cat("the dataset does not have the same number of columns and rows.
The determinant can't be computed\n")
  }

  # print on screen the dataset type before converting it in a tibble
  cat("the type of dataset before becoming a tibble was", typeof(x))

  # if x is a list (but not already a data.frame/tibble), try to coerce it
  # into a data.frame so the rest of the function can work on it
  if (is.list(x) && !is.data.frame(x)) {
    x <- as.data.frame(x)
  }

  # transform the input in a tibble
  if (!is_tibble(x)) {
    x <- as_tibble(x)
  } else {
    cat("\n")
    cat("this is already a tibble, no need to transform it")
    cat("\n")
  }

  # min and max can be computed on any numeric dataset, square or not,
  # so they are calculated here on the tibble version of x
  min_val <- min(x, na.rm = TRUE)
  max_val <- max(x, na.rm = TRUE)

  # print a header line with the given dataset_name
  cat("\n---", dataset_name, "---\n")
  # print the number of columns
  cat("number of columns:", ncol(x), "\n")
  # print the number of rows
  cat("number of rows:   ", nrow(x), "\n")
  # print the maximum value in the dataset
  cat("maximum value:    ", max_val, "\n")
  # print the minimum value in the dataset
  cat("minimum value:    ", min_val, "\n")
  # print the determinant of the matrix (NA if x is not a square matrix)
  cat("determinant value ", det_val, "\n")


  # -------- zero counts --------

  # total count of zeros across the whole dataset
  total_zeroes <- sum(x == 0, na.rm = TRUE)
  cat("the total number of zeroes in", dataset_name, "is", total_zeroes, "\n")

  zero_rate <- round(total_zeroes / (ncol(x) * nrow(x)) * 100, 2)
  cat("the percentage of zeroes in", dataset_name, "is:", zero_rate, "%\n")


  # -------- return the invisible values --------

  invisible(list(
    n_col = ncol(x),
    n_row = nrow(x),
    min_val = min_val,
    max_val = max_val,
    det_val = det_val,
    total_zeroes = total_zeroes,
    zero_rate = zero_rate
  ))
}
