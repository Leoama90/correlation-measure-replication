# datasummary.R
#
# This function takes as input a data set (tibble, data.frame, list, or matrix)
# and prints out the number of columns, number of rows, and the count
# and percentage of zero values it contains, along with basic central
# tendency, dispersion, and shape statistics (mean, median, standard
# deviation, skewness). For square numeric matrices, it also computes
# the determinant and the minimum eigenvalue, the latter being a more
# direct diagnostic of positive semi-definiteness than the determinant
# alone (relevant when summarizing correlation matrices such as those
# produced by generate_matrix_factors()).
# it is intended to gather some data useful to build further code
#
# Input:
#   - a tibble, data.frame, list, or matrix with numeric columns
#     (does not need to be square; determinant and minimum eigenvalue
#     are only computed when the input is a square numeric matrix)
#
# Output:
#   - printed column count, row count, mean, median, standard deviation,
#     skewness, number of zeroes and their percentage, determinant, and
#     minimum eigenvalue (the latter two only for square matrices);
#     invisibly returns these same values as a list
#

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)
# tidyverse: data science toolkit (dplyr, ggplot2, tidyr, etc.)
# [https://tidyverse.org](https://tidyverse.org)
library(tidyverse)

#
#' Summarize a dataset's dimensions, central tendency, shape, and zero content
#'
#' @param x a data.frame, tibble, list, or matrix with numeric columns.
#'   Does not need to be square.
#'
#' @return An invisible list containing: n_col (number of columns),
#'   n_row (number of rows), min_val (minimum value), max_val (maximum
#'   value), mean_val (mean value), median_val (median value), sd_val
#'   (standard deviation), skewness_val (skewness), det_val (determinant,
#'   only if x is a square numeric matrix), min_eigenvalue (smallest
#'   eigenvalue, only if x is a square numeric matrix), total_zeroes
#'   (count of zero values), and zero_rate (percentage of zero values).
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
  # (needed for the determinant and the minimum eigenvalue, both only
  # defined for square matrices)
  is_square_matrix <- is.matrix(x) && is.numeric(x) && nrow(x) == ncol(x)

  # default values: only defined for square matrices
  det_val <- NA
  min_eigenvalue <- NA

  # if x qualifies, compute the determinant and minimum eigenvalue while
  # it is still a matrix; both det() and eigen() error out on matrices
  # containing NA/Inf, so skip them (leaving NA) if any such value is
  # present, instead of letting the function crash
  if (is_square_matrix) {
    if (anyNA(x) || any(is.infinite(x))) {
      cat("the matrix contains NA or infinite values.\n
  The determinant and minimum eigenvalue can't be computed\n")
    } else {
      det_val <- det(x)
      eigenvalues <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
      min_eigenvalue <- min(eigenvalues)
    }
  } else {
    cat("The dataset does not have the same number of columns and rows.\nThe determinant and minimum eigenvalue can't be computed\n")
  }

  # print on screen the dataset type before converting it in a tibble
  cat("The type of dataset before becoming a tibble was", typeof(x))

  # if x is a list (but not already a data.frame/tibble), try to coerce it
  # into a data.frame so the rest of the function can work on it
  if (is.list(x) && !is.data.frame(x)) {
    x <- as.data.frame(x)
  }

  # transform the input in a tibble; .name_repair = "unique" silences the
  # tibble warning that would otherwise appear for matrices without column
  # names (assigning default names like "V1", "V2", ... instead)
  if (!is_tibble(x)) {
    x <- as_tibble(x, .name_repair = "unique")
  } else {
    cat("\n")
    cat("this is already a tibble, no need to transform it")
    cat("\n")
  }

  # min, max, mean, median, sd, and skewness can be computed on any
  # numeric dataset, square or not, so they are calculated here on the
  # tibble version of x
  min_val <- min(x, na.rm = TRUE)
  max_val <- max(x, na.rm = TRUE)
  mean_val <- mean(as.matrix(x), na.rm = TRUE)
  median_val <- median(as.matrix(x), na.rm = TRUE)
  sd_val <- sd(as.matrix(x), na.rm = TRUE)

  # skewness: third standardized moment, computed manually (no external
  # package needed); positive values indicate a right-skewed (heavy
  # right tail) distribution, as typically found in raw metagenomic
  # counts before CLR transformation
  skewness_val <- mean((as.matrix(x) - mean_val)^3, na.rm = TRUE) / sd_val^3

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
  # print the mean value in the dataset
  cat("mean value:       ", mean_val, "\n")
  # print the median value in the dataset
  cat("median value:     ", median_val, "\n")
  # print the standard deviation of the dataset
  cat("standard deviation:", sd_val, "\n")
  # print the skewness of the dataset
  cat("skewness:         ", skewness_val, "\n")
  # print the determinant of the matrix (NA if x is not a square matrix)
  cat("determinant value ", det_val, "\n")
  # print the minimum eigenvalue (NA if x is not a square matrix);
  # a non-negative value indicates the matrix is (numerically) PSD
  cat("minimum eigenvalue", min_eigenvalue, "\n")


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
    mean_val = mean_val,
    median_val = median_val,
    sd_val = sd_val,
    skewness_val = skewness_val,
    det_val = det_val,
    min_eigenvalue = min_eigenvalue,
    total_zeroes = total_zeroes,
    zero_rate = zero_rate
  ))
}
