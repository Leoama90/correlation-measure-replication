# NorTa_simulation.R
#
# Purpose:
#   This script generates a synthetic correlation matrix with
#   generate_matrix_with_zeroes.R, projects it onto the nearest valid
#   positive semi-definite (PSD) correlation matrix, and uses the "Normal
#   To Anything" (NorTA) approach to simulate a sparse metagenomic count
#   dataset with that known correlation structure. Sparsity is introduced
#   by mapping the correlated normal data through a zero-inflated
#   negative binomial (ZINB) quantile function.
#
#   Since the PSD projection is not guaranteed to preserve the exact
#   zero pattern requested via n_zeroes, this script also measures how
#   many zero pairs survive the correction, instead of assuming they do.
#
# Input:
#   - no external files required; the correlation matrix is generated
#     inline via generate_matrix_with_zeroes()
#
# Output:
#   - an invisible list with six elements: R_true (the correlation
#     matrix before PSD correction), R_psd (the valid PSD correlation
#     matrix used for simulation), sim_data (the intermediate,
#     continuous multivariate normal data with that correlation
#     structure), sim_counts (the final sparse count dataset, obtained
#     by mapping sim_data through a ZINB marginal distribution),
#     n_zeroes_true (zero pairs actually found in R_true) and
#     n_zeroes_psd (zero pairs surviving in R_psd after PSD correction)
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
# VGAM: provides qzinegbin(), the quantile function of the zero-inflated
# negative binomial distribution, used here to introduce sparsity
# https://cran.r-project.org/web/packages/VGAM/index.html
library(VGAM)

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

#' Simulate a sparse metagenomic dataset with a known correlation matrix
#'
#' Generates a synthetic correlation matrix with a controlled number of
#' zeroes, projects it onto the nearest valid positive semi-definite (PSD)
#' correlation matrix, simulates multivariate normal data with that exact
#' correlation structure, and finally maps it through a zero-inflated
#' negative binomial (ZINB) quantile function to obtain a sparse count
#' dataset, following the NorTA approach. Since the PSD projection can
#' alter the requested zero pattern, the number of zero pairs actually
#' surviving in R_psd is also measured and returned.
#'
#' @param n integer. Dimension of the correlation matrix (n x n), i.e. the
#'   number of simulated taxa. To match a real dataset, pass
#'   ncol(real_data).
#' @param n_zeroes integer. Number of off-diagonal entries to force to
#'   zero in the correlation matrix. Must be an even number (see
#'   generate_matrix_with_zeroes()).
#' @param N integer. Number of samples to simulate. To match a real
#'   dataset, pass nrow(real_data).
#' @param mu numeric. Mean parameter of the negative binomial component
#'   of the ZINB distribution, shared by all simulated taxa. Default 20.
#' @param size numeric. Dispersion (size) parameter of the negative
#'   binomial component; smaller values mean higher variance/overdispersion.
#'   Default 30.
#' @param phi numeric. Zero-inflation probability of the ZINB
#'   distribution (0 to 1), controlling the extra sparsity added on top
#'   of the negative binomial's own zeroes. Default 0.3.
#' @param seed integer or NULL. Optional seed for reproducibility.
#'   Default is NULL.
#'
#' @return An invisible list with R_true (the original correlation matrix
#'   before PSD correction), R_psd (the corrected, valid PSD correlation
#'   matrix actually used for simulation), sim_data (the intermediate
#'   N x n matrix of continuous multivariate normal data), sim_counts
#'   (the final N x n matrix of sparse ZINB counts), n_zeroes_true (zero
#'   pairs actually found in R_true, should equal n_zeroes) and
#'   n_zeroes_psd (zero pairs surviving in R_psd after the PSD
#'   projection, typically fewer than n_zeroes_true).
#'
#' @examples
#' # basic usage, with default ZINB parameters
#' norta_simulation(n = 10, n_zeroes = 4, N = 100, seed = 42)
#'
#' \dontrun{
#' # matching the dimensions of a real dataset
#' norta_simulation(
#'   n = ncol(real_data), n_zeroes = 20,
#'   N = nrow(real_data), seed = 42
#' )
#' }
#'
#' @export


# -------- body of the function --------

norta_simulation <- function(n, n_zeroes, N, mu = 20, size = 30, phi = 0.3,
                             seed = NULL) {
  
  # generate a synthetic n x n correlation matrix with the requested
  # number of symmetric zero pairs
  R_true <- generate_matrix_with_zeroes(n = n, n_zeroes = n_zeroes, seed = seed)
  
  # project the generated matrix onto the nearest valid positive
  # semi-definite (PSD) correlation matrix (eigenvalues >= 0), since
  # R_true is not guaranteed to be PSD by construction
  R_psd <- nearPD(R_true, corr = TRUE)$mat
  
  
  # -------- measure how many zeroes survive the PSD projection --------
  
  # tolerance below which an entry is still considered "zero"; after
  # nearPD's correction exact 0s are unlikely to remain, since the
  # projection introduces small floating point deviations
  zero_tol <- 1e-8
  
  # count zero pairs actually present in R_true (should equal n_zeroes,
  # since generate_matrix_with_zeroes() enforces it by construction)
  n_zeroes_true <- sum(abs(R_true[upper.tri(R_true)]) < zero_tol) * 2
  
  # count how many of the off-diagonal entries are still (approximately)
  # zero in R_psd, after the PSD projection has been applied
  n_zeroes_psd <- sum(abs(as.matrix(R_psd)[upper.tri(R_psd)]) < zero_tol) * 2
  
  # report the requested vs. actual vs. surviving zero counts, so the
  # loss introduced by the PSD correction is visible at simulation time
  cat("Requested zero pairs:", n_zeroes, "\n")
  cat("Zero pairs in R_true (before PSD correction):", n_zeroes_true, "\n")
  cat("Zero pairs surviving in R_psd (after PSD correction):", n_zeroes_psd, "\n")
  
  
  # -------- NorTA step 1: simulate correlated normal data --------
  
  # simulate N samples from a multivariate normal distribution with the
  # corrected correlation structure R_psd; since R_psd is a correlation
  # matrix (unit diagonal), each column has a standard normal marginal
  sim_data <- rmvnorm(n = N, mean = rep(0, nrow(R_psd)), sigma = as.matrix(R_psd))
  
  
  # -------- NorTA step 2: map normal data to sparse ZINB counts --------
  
  # convert each standard normal value to its rank/percentile (uniform
  # on [0,1]) via the normal CDF; this step preserves the correlation
  # structure, since pnorm() is a monotonic transformation
  sim_unif <- pnorm(sim_data)
  
  # apply the ZINB quantile function (inverse CDF) to turn each uniform
  # value into a sparse count, introducing zeroes via phi on top of the
  # negative binomial's own dispersion (size, mu)
  sim_counts <- qzinegbin(sim_unif, size = size, munb = mu, pstr0 = phi)
  
  # restore the original matrix shape (qzinegbin flattens its input)
  dim(sim_counts) <- dim(sim_data)
  
  # return the invisible values
  invisible(list(
    R_true = R_true,
    R_psd = R_psd,
    sim_data = sim_data,
    sim_counts = sim_counts,
    n_zeroes_true = n_zeroes_true,
    n_zeroes_psd = n_zeroes_psd
  ))
}