# generate_matrix.R
#
# Purpose:
#   This script generates a symmetric matrix that mimics a correlation
#   matrix: values in [-1, 1], diagonal equal to 1, and a controlled
#   number of off-diagonal zero entries placed symmetrically.
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

#' Generate a symmetric correlation-like matrix with a specific number of zeros
#'
#' @param n integer. Dimension of the matrix (n x n).
#' @param n_zeros integer. Number of off-diagonal entries to force to zero.
#'   Must be an even number, since zeros are placed symmetrically in pairs
#'   (position [i,j] and its mirror [j,i]).
#' @param seed integer or NULL. Optional seed for reproducibility of the
#'   random number generation. Default is NULL.
#'
#' @return A symmetric numeric n x n matrix, with values in [-1, 1], a
#'   diagonal of 1s, and exactly n_zeros off-diagonal entries set to zero.
#'
#' @examples
#' generate_matrix_with_zeros(n = 5, n_zeros = 4, seed = 42)
#
# generate_matrix_with_zeros(): builds a symmetric n x n matrix with values
# in [-1, 1], a diagonal of 1s, and a chosen number of symmetric zero pairs
generate_matrix_with_zeros <- function(n, n_zeros, seed = NULL) {
  
  # if a seed is provided, fix it for reproducible random numbers
  if (!is.null(seed)) set.seed(seed)
  
  # stop if n_zeros is not even, since zeros must be placed in symmetric pairs
  if (n_zeros %% 2 != 0) {
    stop("n_zeros must be an even number, since zeros are placed symmetrically")
  }
  
  # count how many entries are available in the upper triangle (off-diagonal)
  n_upper <- n * (n - 1) / 2
  
  # stop with a clear error if more zero-pairs are requested than available
  if (n_zeros / 2 > n_upper) {
    stop("requested n_zeros (", n_zeros, ") exceeds available off-diagonal ",
         "positions (", n_upper * 2, ")")
  }
  
  # generate random values in [-1, 1] only for the upper triangle
  upper_vals <- runif(n_upper, min = -1, max = 1)
  
  # start from a matrix of zeros
  corr_mat <- matrix(0, nrow = n, ncol = n)
  
  # fill the upper triangle with the random values
  corr_mat[upper.tri(corr_mat)] <- upper_vals
  
  # mirror the upper triangle onto the lower triangle to enforce symmetry
  corr_mat <- corr_mat + t(corr_mat)
  
  # build a data frame with every (row, col) position in the upper triangle
  upper_positions <- which(upper.tri(corr_mat), arr.ind = TRUE)
  
  # randomly pick n_zeros/2 distinct upper-triangle positions to zero out
  chosen <- upper_positions[sample(nrow(upper_positions), n_zeros / 2), , drop = FALSE]
  
  # loop over the chosen positions and zero out both [i,j] and its mirror [j,i]
  for (i in seq_len(nrow(chosen))) {
    corr_mat[chosen[i, "row"], chosen[i, "col"]] <- 0
    corr_mat[chosen[i, "col"], chosen[i, "row"]] <- 0
  }
  
  # set the diagonal to 1, as in a correlation matrix
  diag(corr_mat) <- 1
  
  # return the resulting matrix
  return(corr_mat)
}
