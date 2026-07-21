# NorTa_simulation.R
#
# Purpose:
#   This script generates a synthetic correlation matrix with
#   generate_matrix_with_zeroes.R, projects it onto the nearest valid
#   positive semi-definite (PSD) correlation matrix, and uses the "Normal
#   To Anything" (NorTA) approach to simulate a metagenomic dataset with
#   that known correlation structure.
#
# Input:
#   - no external files required; the correlation matrix is generated
#     inline via generate_matrix_with_zeroes()
#
# Output:
#   - a simulated metagenomic dataset with a known correlation structure
#
# Scripts used:
#   - generate_matrix_with_zeroes.R

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)
# Matrix: tools for working with dense and sparse matrices, including nearPD()
# https://cran.r-project.org/web/packages/Matrix/index.html
library(Matrix)
# mvtnorm: multivariate normal distribution functions, including rmvnorm()
# https://cran.r-project.org/web/packages/mvtnorm/index.html
library(mvtnorm)
# tidyverse: data science toolkit (dplyr, ggplot2, tidyr, etc.)
# https://tidyverse.org
library(tidyverse)

# load the script that allows to generate matrices with a specified number
# of zeroes
source(
  list.files(
    path = here(),
    pattern = "^generate_matrix_with_zeroes\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)

#' Simulate a metagenomic dataset with a known correlation matrix
#'
#' Generates a synthetic correlation matrix with a controlled number of
#' zeroes, projects it onto the nearest valid positive semi-definite (PSD)
#' correlation matrix, and simulates multivariate normal data with that
#' exact correlation structure, following the NorTA approach.
#'
#' @param n integer. Dimension of the correlation matrix (n x n), i.e. the
#'   number of simulated taxa.
#' @param n_zeroes integer. Number of off-diagonal entries to force to
#'   zero in the correlation matrix. Must be an even number (see
#'   generate_matrix_with_zeroes()).
#' @param N integer. Number of samples to simulate.
#' @param seed integer or NULL. Optional seed for reproducibility.
#'   Default is NULL.
#'
#' @return An invisible list with R_true (the original correlation matrix
#'   before PSD correction), R_psd (the corrected, valid PSD correlation
#'   matrix actually used for simulation), and sim_data (the simulated
#'   N x n matrix of multivariate normal data).
#'
#' @examples
#' norta_simulation(n = 10, n_zeroes = 4, N = 100, seed = 42)
#'
#' @export
norta_simulation <- function(n, n_zeroes, N, seed = NULL) {
  
  # generate a synthetic n x n correlation matrix with the requested
  # number of symmetric zero pairs
  R_true <- generate_matrix_with_zeroes(n = n, n_zeroes = n_zeroes, seed = seed)
  
  # project the generated matrix onto the nearest valid positive
  # semi-definite (PSD) correlation matrix (eigenvalues >= 0), since
  # R_true is not guaranteed to be PSD by construction
  R_psd <- nearPD(R_true, corr = TRUE)$mat
  
  # simulate N samples from a multivariate normal distribution with the
  # corrected correlation structure R_psd
  sim_data <- rmvnorm(n = N, mean = rep(0, nrow(R_psd)), sigma = as.matrix(R_psd))
  
  # return the invisible values
  invisible(list(
    R_true = R_true,
    R_psd = R_psd,
    sim_data = sim_data
  ))
}