# generate_matrix_factors.R
#
# Purpose:
#   Generates a valid correlation matrix (symmetric, PSD, unit diagonal)
#   with an exact, controlled zero pattern, using a latent-factor (Gram
#   matrix) construction instead of generating random values and
#   correcting them afterwards (as generate_matrix_with_zeroes.R +
#   nearPD() does). Zeroes arise here by construction, not by projection,
#   so no PSD-correction step is needed and no zeroes are lost.
#
# Design note:
#   Exact control over the *number* of zeroes (as in
#   generate_matrix_with_zeroes()) is not generally achievable with a
#   block/factor structure, since the count of zero pairs is fixed by
#   the combinatorics of the group sizes. This function instead takes
#   the number of groups as input, and returns the resulting number of
#   zero pairs alongside the matrix, so the actual sparsity achieved is
#   always known exactly (not just approximately, as with nearPD).
#
# Inputs:
#   - No external files required. All inputs are defined inline:
#     * n <- 40
#     * n_groups <- 5
#
# Outputs:
#   - a list with three elements: mat (the resulting n x n symmetric,
#     PSD, unit-diagonal correlation matrix), groups (the group
#     assignment vector of length n, one group id per taxon), and
#     n_zeroes (the exact number of off-diagonal zero entries in mat,
#     counting both symmetric halves)

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)
# tidyverse: data science toolkit (dplyr, ggplot2, tidyr, etc.)
# https://tidyverse.org
library(tidyverse)

#' Generate a valid correlation matrix with an exact block-zero pattern
#'
#' Builds a correlation matrix as R = D^(-1/2) (Lambda Lambda') D^(-1/2),
#' where Lambda is an n x n_groups matrix of latent-factor loadings with
#' each taxon (row) loading only on its assigned group's column. Since
#' Sigma = Lambda %*% t(Lambda) is a Gram matrix, it is guaranteed to be
#' positive semi-definite for any Lambda, with no PSD-correction step
#' required. Off-diagonal entries between taxa in different groups are
#' exactly zero, since their loading vectors do not share any non-zero
#' column.
#'
#' @param n integer. Dimension of the matrix (n x n), i.e. number of taxa.
#' @param n_groups integer. Number of latent groups (factors) to split
#'   the n taxa into. Taxa within the same group are correlated with
#'   each other; taxa in different groups have exactly zero correlation.
#'   Must be between 1 (fully correlated, no zeroes) and n (fully
#'   independent, diagonal matrix).
#' @param min_loading numeric. Minimum absolute value for a taxon's
#'   random loading on its group. Default 0.3.
#' @param max_loading numeric. Maximum absolute value for a taxon's
#'   random loading on its group. Default 1.
#' @param seed integer or NULL. Optional seed for reproducibility of the
#'   random number generation. Default is NULL.
#'
#' @return A list with mat (the resulting n x n symmetric, PSD,
#'   unit-diagonal correlation matrix), groups (the group assignment
#'   vector of length n, one group id per taxon), and n_zeroes (the
#'   exact number of off-diagonal zero entries in mat, counting both
#'   symmetric halves).
#'
#' @examples
#' # 40 taxa split into 5 groups: taxa in different groups are
#' # exactly uncorrelated, no PSD correction is ever needed
#' generate_matrix_factors(n = 40, n_groups = 5, seed = 42)
#'
#' # n_groups = 1: every taxon shares the same group, so no zeroes
#' generate_matrix_factors(n = 10, n_groups = 1, seed = 42)
#'
#' # n_groups = n: every taxon is its own group, so the matrix is
#' # the identity (no correlation at all between distinct taxa)
#' generate_matrix_factors(n = 10, n_groups = 10, seed = 42)
#'
#' @export


# -------- body of the function --------

generate_matrix_factors <- function(n, n_groups, min_loading = 0.3,
                                    max_loading = 1, seed = NULL) {
  
  # if a seed is provided, fix it for reproducible random numbers
  if (!is.null(seed)) set.seed(seed)
  
  # stop if n_groups is out of the valid range [1, n]
  if (n_groups < 1 || n_groups > n) {
    stop("n_groups must be between 1 and n (got n_groups = ", n_groups,
         ", n = ", n, ")")
  }
  
  # randomly assign each of the n taxa to one of the n_groups groups;
  # this assignment determines the zero pattern of the final matrix
  groups <- sample(rep_len(1:n_groups, n))
  
  # start from a matrix of zeroes: rows = taxa, columns = latent groups
  Lambda <- matrix(0, nrow = n, ncol = n_groups)
  
  # fill in each taxon's loading only in the column of its own group;
  # sign is randomized so correlations within a group can be negative
  loadings <- runif(n, min = min_loading, max = max_loading) * sample(c(-1, 1), n, replace = TRUE)
  for (i in seq_len(n)) {
    Lambda[i, groups[i]] <- loadings[i]
  }
  
  # Gram matrix: guaranteed positive semi-definite for any Lambda, with
  # no PSD-correction step needed, unlike generate_matrix_with_zeroes()
  Sigma <- Lambda %*% t(Lambda)
  
  # rescale to a correlation matrix (unit diagonal), dividing each
  # entry by the geometric mean of its row/column standard deviations
  d <- 1 / sqrt(diag(Sigma))
  corr_mat <- Sigma * outer(d, d)
  
  # clean up numerical noise on the diagonal (should already be ~1)
  diag(corr_mat) <- 1
  
  # count the exact number of off-diagonal zero entries (both symmetric
  # halves), i.e. the pairs of taxa placed in different groups
  tol <- 1e-8
  n_zeroes <- sum(abs(corr_mat[upper.tri(corr_mat)]) < tol) * 2
  
  # return the matrix, the group assignment, and the exact zero count
  list(
    mat = corr_mat,
    groups = groups,
    n_zeroes = n_zeroes
  )
}