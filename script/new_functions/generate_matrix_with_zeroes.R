# generate_matrix.R
#
# Purpose:
#   This script generates a normalized square matrix filled with random
#   values, allowing the user to control exactly how many entries are
#   forced to zero, while optionally keeping the diagonal equal to 1.
#
# Inputs:
#   - No external files required. All inputs are defined inline:
#     * n <- 40
#     * n_zeros <- 10
#
# Outputs:
#   - the generated matrix printed to the console
#   - a summary of the matrix properties (dimensions, min, max, determinant,
#     number of zeros) printed to the console

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)
# tidyverse: data science toolkit (dplyr, ggplot2, tidyr, etc.)
# https://tidyverse.org
library(tidyverse)

#' Generate a normalized matrix with a specific number of zeros
#'
#' @param n integer. Dimension of the matrix (n x n).
#' @param n_zeros integer. Number of entries to force to zero.
#' @param keep_diag_one logical. If TRUE, keeps the diagonal equal to 1
#'   and never places zeros on it. Default is TRUE.
#' @param seed integer or NULL. Optional seed for reproducibility of the
#'   random number generation. Default is NULL.
#'
#' @return A numeric n x n matrix, normalized in [0, 1], with exactly
#'   n_zeros entries set to zero.
#'
#' @examples
#' generate_matrix_with_zeros(n = 5, n_zeros = 3, seed = 42)
#
# generate_matrix_with_zeros(): builds a normalized n x n matrix with a
# chosen number of zero entries, optionally protecting the diagonal
generate_matrix_with_zeros <- function(n, n_zeros, keep_diag_one = TRUE, seed = NULL) {
  
  # if a seed is provided, fix it for reproducible random numbers
  if (!is.null(seed)) set.seed(seed)
  
  # generate n*n random numbers and arrange them into a square matrix
  mat <- matrix(runif(n^2, min = 0, max = 100), nrow = n, ncol = n)
  
  # normalize the matrix by dividing every value by the matrix's maximum
  norm_mat <- mat / max(mat)
  
  # if requested, force the diagonal values to be 1
  if (keep_diag_one) diag(norm_mat) <- 1
  
  # build a data frame with every possible (row, col) position in the matrix
  all_positions <- expand.grid(row = 1:n, col = 1:n)
  
  # if the diagonal must stay untouched, remove diagonal positions from the pool
  if (keep_diag_one) {
    all_positions <- all_positions[all_positions$row != all_positions$col, ]
  }
  
  # stop with a clear error if more zeros are requested than available positions
  if (n_zeros > nrow(all_positions)) {
    stop("requested n_zeros (", n_zeros, ") exceeds available positions (",
         nrow(all_positions), ")")
  }
  
  # randomly pick n_zeros distinct positions from the available pool
  chosen <- all_positions[sample(nrow(all_positions), n_zeros), ]
  
  # loop over the chosen positions and set each one to zero
  for (i in seq_len(nrow(chosen))) {
    norm_mat[chosen$row[i], chosen$col[i]] <- 0
  }
  
  # return the resulting matrix
  return(norm_mat)
}

#' Print a summary of matrix properties
#'
#' @param mat numeric matrix. The matrix to summarize.
#' @param label character. A short label to identify this summary in the
#'   console output. Default is an empty string.
#'
#' @return Invisibly returns NULL. Called for its side effect of printing
#'   a summary to the console.
#'
#' @examples
#' m <- generate_matrix_with_zeros(n = 5, n_zeros = 3, seed = 42)
#' print_matrix_summary(m, label = "example matrix")
#
# print_matrix_summary(): prints dimensions, min/max, determinant and zero
# count of a given matrix, labeled for readability
print_matrix_summary <- function(mat, label = "") {
  # print a header line with the given label
  cat("\n---", label, "---\n")
  # print the number of columns
  cat("number of columns:", ncol(mat), "\n")
  # print the number of rows
  cat("number of rows:   ", nrow(mat), "\n")
  # print the maximum value in the matrix
  cat("maximum value:    ", max(mat), "\n")
  # print the minimum value in the matrix
  cat("minimum value:    ", min(mat), "\n")
  # print the determinant of the matrix
  cat("determinant:      ", det(mat), "\n")
  # print how many entries are exactly zero
  cat("number of zeros:  ", sum(mat == 0), "\n")
}

# set the matrix dimension
n <- 40
# set how many zeros should be inserted
n_zeros <- 10
# generate the matrix using the function defined above
norm_mat <- generate_matrix_with_zeros(n = n, n_zeros = n_zeros,
                                       keep_diag_one = TRUE, seed = 42)
# print the matrix to the console
print(norm_mat)
# print the summary of the generated matrix
print_matrix_summary(norm_mat, label = "generated matrix")