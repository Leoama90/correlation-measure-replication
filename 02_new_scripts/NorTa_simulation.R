# NorTa_simulation.R
#
# Purpose:
#   This script generates a synthetic correlation matrix with
#   generate_matrix_factors.R and uses the "Normal To Anything" (NorTA)
#   approach to simulate a sparse metagenomic count dataset with that
#   known correlation structure. Sparsity is introduced by mapping the
#   correlated normal data through a zero-inflated negative binomial
#   (ZINB) quantile function.
#
#   The correlation matrix is built via a latent-factor (Gram matrix)
#   construction: Sigma = Lambda %*% t(Lambda) is guaranteed positive
#   semi-definite by construction, for any Lambda, so no PSD-correction
#   step (e.g. nearPD()) is needed. This avoids the issue found with the
#   previous approach (generate_matrix_with_zeroes.R + nearPD()), where
#   the PSD projection was found to eliminate nearly all requested
#   zeroes. Here, zero pairs correspond to taxa placed in different
#   latent groups, and are exact by construction, with no loss.
#
# Input:
#   - no external files required; the correlation matrix is generated
#     inline via generate_matrix_factors()
#
# Output:
#   - an invisible list with five elements: R_true (the correlation
#     matrix used for simulation), groups (the latent group assignment
#     for each taxon), n_zeroes (the exact number of zero pairs in
#     R_true), sim_data (the intermediate, continuous multivariate
#     normal data with that correlation structure), and sim_counts
#     (the final sparse count dataset, obtained by mapping sim_data
#     through a ZINB marginal distribution)
#
# Scripts used:
#   - generate_matrix_factors.R

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)
# mvtnorm: multivariate normal distribution functions, including rmvnorm()
# https://cran.r-project.org/web/packages/mvtnorm/index.html
library(mvtnorm)
# VGAM: provides qzinegbin(), the quantile function of the zero-inflated
# negative binomial distribution, used here to introduce sparsity
# https://cran.r-project.org/web/packages/VGAM/index.html
library(VGAM)

# load the script that allows to generate a valid correlation matrix
# with an exact, controlled zero pattern via a latent-factor structure
source(
  list.files(
    path = here(),
    pattern = "^generate_matrix_factors\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)

#' Simulate a sparse metagenomic dataset with a known correlation matrix
#'
#' Generates a valid correlation matrix with an exact, controlled zero
#' pattern via a latent-factor (Gram matrix) construction, simulates
#' multivariate normal data with that exact correlation structure, and
#' finally maps it through a zero-inflated negative binomial (ZINB)
#' quantile function to obtain a sparse count dataset, following the
#' NorTA approach.
#'
#' @param n integer. Dimension of the correlation matrix (n x n), i.e. the
#'   number of simulated taxa. To match a real dataset, pass
#'   ncol(real_data).
#' @param n_groups integer. Number of latent groups to split the n taxa
#'   into (see generate_matrix_factors()). Taxa in the same group are
#'   correlated with each other; taxa in different groups have exactly
#'   zero correlation. Must be between 1 and n.
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
#' @return An invisible list with R_true (the correlation matrix used
#'   for simulation, valid and PSD by construction), groups (the latent
#'   group assignment for each of the n taxa), n_zeroes (the exact
#'   number of off-diagonal zero pairs in R_true), sim_data (the
#'   intermediate N x n matrix of continuous multivariate normal data),
#'   and sim_counts (the final N x n matrix of sparse ZINB counts).
#'
#' @examples
#' # basic usage, with default ZINB parameters
#' norta_simulation(n = 10, n_groups = 3, N = 100, seed = 42)
#'
#' \dontrun{
#' # matching the dimensions of a real dataset
#' norta_simulation(
#'   n = ncol(real_data), n_groups = 5,
#'   N = nrow(real_data), seed = 42
#' )
#' }
#'
#' @export


# -------- body of the function --------

norta_simulation <- function(n, n_groups, N, mu = 20, size = 30, phi = 0.3,
                             seed = NULL) {
  
  # generate a valid n x n correlation matrix with an exact zero
  # pattern determined by the latent group assignment; PSD by
  # construction, so no correction step is needed here
  factors_result <- generate_matrix_factors(n = n, n_groups = n_groups, seed = seed)
  R_true <- factors_result$mat
  groups <- factors_result$groups
  n_zeroes <- factors_result$n_zeroes
  
  # report the achieved sparsity, which is exact (no loss), unlike the
  # previous nearPD-based approach
  cat("Number of latent groups:", n_groups, "\n")
  cat("Zero pairs in R_true (exact, by construction):", n_zeroes, "\n")
  
  
  # -------- NorTA step 1: simulate correlated normal data --------
  
  # simulate N samples from a multivariate normal distribution with the
  # correlation structure R_true; since R_true is a correlation matrix
  # (unit diagonal), each column has a standard normal marginal
  sim_data <- rmvnorm(n = N, mean = rep(0, nrow(R_true)), sigma = R_true)
  
  
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
    groups = groups,
    n_zeroes = n_zeroes,
    sim_data = sim_data,
    sim_counts = sim_counts
  ))
}