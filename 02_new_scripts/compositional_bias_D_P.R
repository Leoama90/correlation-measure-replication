# compositional_bias_D_P.R
#
# Purpose:
#   Replicates Section 3.1 of Fuschi et al. (2025): quantifies the
#   compositional bias introduced by L1 and CLR normalizations on
#   Pearson correlation, as a function of dimensionality (D, number of
#   taxa) and within-sample diversity (P, Pielou index), using a fully
#   uncorrelated ground-truth correlation matrix (identity) so that any
#   non-zero estimated correlation is attributable purely to the
#   compositional constraint, not to genuine correlation structure.
#
#   For each (D, P) combination, this script:
#     1. generates N independent standard normal taxa (ground truth
#        correlation = identity matrix)
#     2. shifts the data to be strictly positive (required to compute
#        the Pielou index and to apply L1/CLR)
#     3. tunes the diversity of the dataset to the target P, by scaling
#        one taxon's column (as in the paper), via a numerical root
#        search on the scaling factor
#     4. applies L1 and CLR normalizations, computes their correlation
#        matrices, and the mean absolute error (MAE) against the
#        ground-truth identity matrix
#
#   Scale note: the paper uses a 40 x 39 grid (D from 5 to 200 step 5,
#   P from 0.025 to 0.975 step 0.025) with N = 10,000 samples per
#   dataset (1560 datasets total). This script uses a reduced 15 x 15
#   grid with N = 5,000, a deliberate, declared time-driven compromise,
#   not a methodological limitation of the approach itself.
#
# Inputs:
#   - No external files required. All inputs are defined inline:
#     * D_values: 15 dimensionality values between 5 and 200
#     * P_values: 15 target diversity values between 0.05 and 0.9
#     * N <- 5000
#
# Outputs:
#   - compositional_bias_results.rds: a tibble with columns D, P_target,
#     P_achieved, mae_L1, mae_CLR, one row per (D, P) combination
#   - compositional_bias_heatmap.png: two-panel heatmap of MAE (log10
#     scale) as a function of D and P, for L1 (left) and CLR (right),
#     analogous to Figure 2A of the paper
#
# Used scripts:
#   - pielou_ind.R
#
# Known limitation:
#   clr_transform() below duplicates the CLR logic already present
#   inside clr_on_data() (clr_pearson.R). It is kept as a small local
#   helper here, rather than extracted into a shared clr.R script, due
#   to time constraints; refactoring into a single shared function is a
#   natural follow-up if time allows.

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# mvtnorm: multivariate normal distribution functions, including rmvnorm()
# https://cran.r-project.org/web/packages/mvtnorm/index.html
library(mvtnorm)

# tidyverse: data science toolkit (dplyr, ggplot2, tidyr, etc.)
# https://tidyverse.org
library(tidyverse)


# bring pielou_ind() into scope
source(
  list.files(
    path = here(),
    pattern = "^pielou_ind\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)


# -------- local helper functions --------

# L1 normalization: convert each sample (row) to relative abundances
l1_normalize <- function(x) {
  # divide each element by the sum of its corresponding sample
  sweep(x, 1, rowSums(x), "/")
}

# CLR transformation: log of each value minus the row-wise mean of the
# logs (equivalent to dividing by the geometric mean of each sample);
# see the "Known limitation" note above regarding duplication with
# clr_on_data() in clr_pearson.R
clr_transform <- function(x) {
  # compute the natural logarithm of every element in the matrix
  log_x <- log(x)

  # subtract the mean log-abundance of each sample from every taxon
  log_x - rowMeans(log_x)
}

# mean absolute error between an estimated correlation matrix and the
# ground-truth matrix, averaged over all D^2 entries, as defined in the
# paper (Section 3.1)
mae_matrix <- function(R_est, R_true) {
  # compute the average absolute difference between the two matrices
  mean(abs(R_est - R_true))
}

# tunes the diversity of a positive dataset to a target Pielou index,
# by scaling its first column by a factor >= 1 (searched on a log10
# scale via uniroot()); increasing the factor concentrates abundance
# in that column, monotonically decreasing diversity from the dataset's
# baseline (factor = 1) towards 0. Returns NA (with a warning) if
# target_P is above the achievable baseline, since only decreasing
# diversity is supported by this scaling strategy (matching the paper's
# approach of applying a multiplicative factor to skew the data).
tune_diversity <- function(x, target_p, tol = 0.005, max_log_factor = 8) {
  # calculate the baseline Pielou diversity before scaling
  baseline_p <- pielou_ind(x)

  # check whether the requested diversity is higher than the baseline
  if (target_p > baseline_p) {
    # issue a warning because this scaling strategy cannot increase diversity
    warning(
      "target_p (", round(target_p, 3), ") exceeds the dataset's ",
      "baseline diversity (", round(baseline_p, 3), "); skipping."
    )

    # return an empty dataset and NA for the achieved diversity
    return(list(data = NULL, achieved_p = NA))
  }

  # define the function whose root corresponds to the target diversity
  diversity_gap <- function(log_factor) {
    # create a copy of the dataset that can be modified
    x_scaled <- x
    # multiply the first taxon's abundance by the scaling factor
    x_scaled[, 1] <- x_scaled[, 1] * (10^log_factor)
    # calculate the difference between achieved and target diversity
    pielou_ind(x_scaled) - target_p
  }

  # find the scaling factor that produces the target diversity
  root <- uniroot(diversity_gap,
    interval = c(0, max_log_factor),
    tol = tol
  )$root

  # create a copy of the original dataset for the final scaling
  x_tuned <- x
  # apply the scaling factor corresponding to the numerical root
  x_tuned[, 1] <- x_tuned[, 1] * (10^root)
  # return the tuned dataset and its achieved Pielou diversity
  list(data = x_tuned, achieved_p = pielou_ind(x_tuned))
}


# -------- experimental grid (reduced scale, see header note) --------

# set the random seed to make the simulation reproducible
set.seed(42)
# define the dimensionality values used in the reduced experimental grid
D_values <- round(seq(5, 200, length.out = 15))
# define the target Pielou diversity values used in the reduced grid
P_values <- seq(0.05, 0.9, length.out = 15)
# define the number of samples generated for each dataset
N <- 5000


# -------- main loop over the (D, P) grid --------

# report the size of the experimental grid and the number of samples
cat(
  "Running compositional bias analysis over a", length(D_values), "x",
  length(P_values), "grid (D x P), N =", N, "samples per dataset.\n"
)

# report that the simulation is a reduced replication of the paper
cat("This replicates a reduced version of Section 3.1 in Fuschi et al. (2025).\n\n")

# initialize the list used to store the results for each dataset
results <- list()
# initialize the index used to fill the results list
i <- 1
# iterate over all dimensionality values
for (D in D_values) {
  # create the ground-truth identity correlation matrix for the current D
  R_true <- diag(D)

  # generate independent standard normal variables for all taxa
  raw_data <- rmvnorm(n = N, mean = rep(0, D), sigma = R_true)

  # shift each taxon so that its minimum becomes zero, then add one
  shifted_data <- sweep(raw_data, 2, apply(raw_data, 2, min), "-") + 1

  # iterate over all target diversity values
  for (P_target in P_values) {
    # tune the dataset to the current target Pielou diversity
    tuned <- tune_diversity(shifted_data, P_target)
    # skip this (D, P) combination if the target diversity was not achievable
    if (is.null(tuned$data)) {
      # move to the next position in the results list
      i <- i + 1
      # skip the remaining operations for this P value
      next
    }

    # apply L1 normalization to the tuned dataset
    l1_data <- l1_normalize(tuned$data)
    # apply the CLR transformation to the tuned dataset
    clr_data <- clr_transform(tuned$data)
    # calculate the Pearson correlation matrix after L1 normalization
    R_L1 <- cor(l1_data)
    # calculate the Pearson correlation matrix after CLR transformation
    R_CLR <- cor(clr_data)
    # calculate the MAE between the L1 correlation matrix and the ground truth
    mae_L1 <- mae_matrix(R_L1, R_true)
    # calculate the MAE between the CLR correlation matrix and the ground truth
    mae_CLR <- mae_matrix(R_CLR, R_true)

    # store the results for the current (D, P) combination
    results[[i]] <- tibble(
      # store the current dimensionality
      D = D,
      # store the requested Pielou diversity
      P_target = P_target,
      # store the diversity actually achieved after tuning
      P_achieved = tuned$achieved_p,
      # store the MAE obtained after L1 normalization
      mae_L1 = mae_L1,
      # store the MAE obtained after CLR transformation
      mae_CLR = mae_CLR
    )

    # increment the result index
    i <- i + 1
  }

  # report that all P values for the current dimensionality are complete
  cat("Completed D =", D, "\n")
}

# combine all individual result rows into a single tibble
results_df <- bind_rows(results)


# -------- save the results --------

# create the output directory if it does not already exist
dir.create(here("outputs"), showWarnings = FALSE)

# save the complete results tibble as an RDS file
saveRDS(results_df, here("outputs", "compositional_bias_results.rds"))

# report the output path and number of generated result rows
cat(
  "\nResults saved to outputs/compositional_bias_results.rds (",
  nrow(results_df), "rows ).\n"
)


# -------- heatmap: MAE(D, P) for L1 and CLR, analogous to Figure 2A --------

# reshape the results from wide to long format for plotting
results_long <- results_df %>%
  # keep only the variables required for the heatmap
  select(D, P_target, mae_L1, mae_CLR) %>%
  # combine the L1 and CLR MAE columns into a single method column
  pivot_longer(
    cols = c(mae_L1, mae_CLR),
    names_to = "method",
    values_to = "mae"
  ) %>%
  # replace the original column names with readable method labels
  mutate(method = recode(method, mae_L1 = "L1", mae_CLR = "CLR"))

# construct the heatmap using dimensionality, diversity and MAE
p_heatmap <- ggplot(results_long, aes(x = D, y = P_target, fill = mae)) +
  # represent each (D, P) combination as a rectangular tile
  geom_tile() +
  # create a separate heatmap panel for L1 and CLR
  facet_wrap(~method) +
  # map the MAE values to a viridis colour scale on a logarithmic axis
  scale_fill_viridis_c(name = "MAE", trans = "log10") +
  # apply a simple black-and-white theme
  theme_bw() +
  # label the x-axis with dimensionality
  xlab("Dimensionality (D)") +
  # label the y-axis with within-dataset diversity
  ylab("Within-dataset diversity (P)") +
  # add the main title and subtitle to the heatmap
  ggtitle("Compositional bias of L1 vs CLR normalization",
    subtitle = "Reduced-scale replication of Fuschi et al. (2025), Figure 2A"
  )

# open a PNG graphics device for saving the heatmap
png(here("outputs", "compositional_bias_heatmap.png"),
  width = 2400, height = 1200, res = 300
)

# render the heatmap on the active graphics device
print(p_heatmap)

# close the graphics device and save the PNG file
dev.off()

# report the location of the generated heatmap
cat("Heatmap saved to outputs/compositional_bias_heatmap.png\n")
