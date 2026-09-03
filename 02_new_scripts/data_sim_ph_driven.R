# data_sim_ph_driven.R
#
# Purpose:
#   Simulates a sparse metagenomic count dataset (following the NorTA
#   approach, as in NorTa_simulation.R) where sparsity is not imposed
#   as a single fixed zero-inflation probability shared by all taxa,
#   but derives mechanistically from a single environmental pH value
#   for the simulated community and each taxon's own pH niche. 
#   Each taxon is assigned a random optimal pH (mu_i) and tolerance (sigma_i);
#   taxa whose niche is far from the community's pH are given a higher
#   zero-inflation probability (more likely to be absent from a sample),
#   while taxa well matched to the environmental pH are given a low
#   zero-inflation probability. This is a standalone script, independent
#   of NorTa_simulation.R, so the two zero-inflation strategies (fixed
#   phi vs pH-driven per-taxon phi) can be compared directly.
#
#   For a given community, the pH is a single value (not varying
#   sample-to-sample), representing one simulated environmental
#   condition; running the script for several pH values produces a
#   family of datasets across an environmental gradient.
#
# Inputs:
#   - No external files required. All inputs are defined inline:
#     * n <- 40
#     * N <- 200
#     * ph <- 6.5
#
# Outputs:
#   - a list with: mat (the n x n correlation matrix used for the
#     underlying correlated structure, from generate_matrix_factors()),
#     groups (the latent group assignment from generate_matrix_factors()),
#     ph_optima (each taxon's optimal pH), phi_per_taxon (each taxon's
#     resulting zero-inflation probability, given the community pH),
#     sim_data (the intermediate continuous multivariate normal data),
#     sim_counts (the final sparse count matrix), and total_zero_rate
#     (the overall fraction of zero entries in sim_counts)
#
# Used scripts:
#   - generate_matrix_factors.R
# Note:
#     columns are labeled OTU_XX for consistency with the naming convention of the
#     real HMP2 dataset used elsewhere in this project,
#     though no sequence clustering is actually performed here,
#     each column is an abstract simulated taxon.

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


# -------- bring generate_matrix_factors() into scope --------

source(
  list.files(
    path = here(),
    pattern = "^generate_matrix_factors\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)


#' Simulate a sparse metagenomic dataset with pH-driven, per-taxon sparsity
#'
#' Assigns each taxon a random optimal pH (its ecological niche centre)
#' and tolerance, then derives a per-taxon zero-inflation probability
#' from the distance between the community's environmental pH and each
#' taxon's optimum: taxa poorly matched to the environment are more
#' likely to be absent (zero) from a given sample, while well-matched
#' taxa are rarely absent. The underlying correlation structure is
#' generated with generate_matrix_factors(), and correlated counts are
#' simulated via the NorTA approach (as in NorTa_simulation.R), but with
#' a per-taxon zero-inflation vector instead of a single shared phi.
#'
#' @param n integer. Number of simulated taxa.
#' @param N integer. Number of samples to simulate.
#' @param ph numeric. The environmental pH of the simulated community
#'   (a single value, shared by all samples in this dataset).
#' @param n_groups integer. Number of latent correlation groups, passed
#'   to generate_matrix_factors() (see that script for details).
#' @param ph_min numeric. Lower bound of the plausible pH range from
#'   which taxa optima are drawn. Default 5.5 (approximate lower bound
#'   of human colon pH).
#' @param ph_max numeric. Upper bound of the plausible pH range from
#'   which taxa optima are drawn. Default 7.5 (approximate upper bound
#'   of human colon pH).
#' @param sigma_min numeric. Minimum pH tolerance (niche width) a taxon
#'   can be assigned. Default 0.3.
#' @param sigma_max numeric. Maximum pH tolerance (niche width) a taxon
#'   can be assigned. Default 0.8.
#' @param phi_max numeric. Upper bound on zero-inflation probability
#'   (0 to 1): even a taxon maximally mismatched to the environmental pH
#'   is capped at this probability, so it is never certain to be absent
#'   from every sample. Default 0.9.
#' @param mu numeric. Mean parameter of the negative binomial component
#'   of the ZINB distribution, shared by all simulated taxa. Default 20.
#' @param size numeric. Dispersion (size) parameter of the negative
#'   binomial component. Default 30.
#' @param seed integer or NULL. Optional seed for reproducibility.
#'   Default is NULL.
#'
#' @return An invisible list with mat (the correlation matrix used),
#'   groups (latent group assignment), ph_optima (each taxon's optimal
#'   pH), phi_per_taxon (each taxon's resulting zero-inflation
#'   probability), sim_data (intermediate continuous normal data),
#'   sim_counts (the final sparse count matrix), and total_zero_rate
#'   (overall fraction of zero entries in sim_counts).
#'
#' @examples
#' # a community well within the typical colon pH range
#' data_sim_ph_driven(n = 40, N = 200, ph = 6.5, n_groups = 5, seed = 42)
#'
#' # a more acidic environment: taxa with high-pH optima become sparser
#' data_sim_ph_driven(n = 40, N = 200, ph = 5.6, n_groups = 5, seed = 42)
#'
#' @export


# -------- body of the function --------

data_sim_ph_driven <- function(n, N, ph, n_groups, ph_min = 5.5, ph_max = 7.5,
                               sigma_min = 0.3, sigma_max = 0.8,
                               phi_max = 0.7, mu = 20, size = 30,
                               seed = NULL) {
  # if a seed is provided, fix it for reproducible random numbers
  if (!is.null(seed)) set.seed(seed)


  # -------- correlation structure (reused from generate_matrix_factors) --------

  factors_result <- generate_matrix_factors(n = n, n_groups = n_groups, seed = seed)
  R_true <- factors_result$mat
  groups <- factors_result$groups


  # -------- per-taxon pH niche and resulting zero-inflation --------

  # each taxon gets a random optimal pH and a random tolerance (niche width)
  ph_optima <- runif(n, min = ph_min, max = ph_max)
  ph_tolerance <- runif(n, min = sigma_min, max = sigma_max)

  # zero-inflation probability per taxon: 0 when the community pH exactly
  # matches the taxon's optimum, approaching phi_max as the mismatch grows,
  # following a Gaussian suitability curve (same functional form used for
  # correlation decay in generate_matrix_ph.R, applied here to sparsity
  # instead of correlation)
  suitability <- exp(-(ph - ph_optima)^2 / (2 * ph_tolerance^2))
  phi_per_taxon <- phi_max * (1 - suitability)


  # -------- NorTA step 1: simulate correlated normal data --------

  sim_data <- mvtnorm::rmvnorm(n = N, mean = rep(0, n), sigma = R_true)


  # -------- NorTA step 2: map normal data to sparse ZINB counts --------

  # convert each standard normal value to its rank/percentile (uniform
  # on [0,1]) via the normal CDF; this step preserves the correlation
  # structure, since pnorm() is a monotonic transformation
  sim_unif <- pnorm(sim_data)

  # expand phi_per_taxon (length n, one value per taxon/column) to match
  # sim_unif's column-major flattened layout: sim_unif is N x n, so its
  # flattened form lists all N values of column 1, then all N values of
  # column 2, and so on. rep(phi_per_taxon, each = N) produces exactly
  # that pattern - each taxon's phi repeated N times in a row - so every
  # element gets the zero-inflation probability of its own taxon/column,
  # regardless of which sample/row it belongs to
  phi_expanded <- rep(phi_per_taxon, each = N)

  # apply the ZINB quantile function with a per-element (per-taxon)
  # zero-inflation probability, instead of the single shared phi used
  # in NorTa_simulation.R
  sim_counts <- VGAM::qzinegbin(sim_unif, size = size, munb = mu, pstr0 = phi_expanded)

  # restore the original matrix shape (qzinegbin flattens its input)
  dim(sim_counts) <- dim(sim_data)


  # -------- assign OTU names --------

  # zero-padded names (OTU_01, OTU_02, ...), width chosen automatically
  # based on n, so column order stays numerically correct even when
  # sorted alphabetically (e.g. "OTU_02" before "OTU_10")
  otu_names <- sprintf(paste0("OTU_%0", nchar(n), "d"), seq_len(n))
  colnames(sim_data) <- otu_names
  colnames(sim_counts) <- otu_names


  # -------- summarize achieved sparsity --------

  total_zero_rate <- mean(sim_counts == 0)


  # -------- return the invisible values --------

  invisible(list(
    mat = R_true,
    groups = groups,
    ph_optima = ph_optima,
    phi_per_taxon = phi_per_taxon,
    sim_data = sim_data,
    sim_counts = sim_counts,
    total_zero_rate = total_zero_rate
  ))
}
