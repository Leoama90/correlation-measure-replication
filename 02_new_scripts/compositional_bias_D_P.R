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
  sweep(x, 1, rowSums(x), "/")
}

# CLR transformation: log of each value minus the row-wise mean of the
# logs (equivalent to dividing by the geometric mean of each sample);
# see the "Known limitation" note above regarding duplication with
# clr_on_data() in clr_pearson.R
clr_transform <- function(x) {
  log_x <- log(x)
  log_x - rowMeans(log_x)
}

# mean absolute error between an estimated correlation matrix and the
# ground-truth matrix, averaged over all D^2 entries, as defined in the
# paper (Section 3.1)
mae_matrix <- function(R_est, R_true) {
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
  baseline_p <- pielou_ind(x)
  
  if (target_p > baseline_p) {
    warning("target_p (", round(target_p, 3), ") exceeds the dataset's ",
            "baseline diversity (", round(baseline_p, 3), "); skipping.")
    return(list(data = NULL, achieved_p = NA))
  }
  
  # function whose root (in log10 factor) gives the target diversity
  diversity_gap <- function(log_factor) {
    x_scaled <- x
    x_scaled[, 1] <- x_scaled[, 1] * (10^log_factor)
    pielou_ind(x_scaled) - target_p
  }
  
  # search for the scaling factor between 1 (log10 = 0) and a large
  # upper bound; diversity_gap is monotonically decreasing in this range
  root <- uniroot(diversity_gap, interval = c(0, max_log_factor),
                  tol = tol)$root
  
  x_tuned <- x
  x_tuned[, 1] <- x_tuned[, 1] * (10^root)
  
  list(data = x_tuned, achieved_p = pielou_ind(x_tuned))
}


# -------- experimental grid (reduced scale, see header note) --------

set.seed(42)

D_values <- round(seq(5, 200, length.out = 15))
P_values <- seq(0.05, 0.9, length.out = 15)
N <- 5000


# -------- main loop over the (D, P) grid --------

cat("Running compositional bias analysis over a", length(D_values), "x",
    length(P_values), "grid (D x P), N =", N, "samples per dataset.\n")
cat("This replicates a reduced version of Section 3.1 in Fuschi et al. (2025).\n\n")

results <- list()
i <- 1

for (D in D_values) {
  
  # ground-truth correlation matrix: identity (all true correlations
  # are exactly 0), so any non-zero estimated correlation after L1/CLR
  # is attributable purely to the compositional constraint
  R_true <- diag(D)
  
  # generate the raw, uncorrelated multivariate normal data once per D
  # (reused across all P values for this D, then re-tuned each time)
  raw_data <- rmvnorm(n = N, mean = rep(0, D), sigma = R_true)
  
  # shift to strictly positive values, per column, so the correlation
  # structure is unaffected (Pearson correlation is invariant to any
  # per-column additive shift)
  shifted_data <- sweep(raw_data, 2, apply(raw_data, 2, min), "-") + 1
  
  for (P_target in P_values) {
    
    tuned <- tune_diversity(shifted_data, P_target)
    
    # skip this (D, P) combination if the target diversity was not
    # achievable (see tune_diversity() for when this happens)
    if (is.null(tuned$data)) {
      i <- i + 1
      next
    }
    
    # apply both normalizations and compute their correlation matrices
    l1_data <- l1_normalize(tuned$data)
    clr_data <- clr_transform(tuned$data)
    
    R_L1 <- cor(l1_data)
    R_CLR <- cor(clr_data)
    
    # compare each against the ground-truth identity matrix
    mae_L1 <- mae_matrix(R_L1, R_true)
    mae_CLR <- mae_matrix(R_CLR, R_true)
    
    results[[i]] <- tibble(
      D = D,
      P_target = P_target,
      P_achieved = tuned$achieved_p,
      mae_L1 = mae_L1,
      mae_CLR = mae_CLR
    )
    i <- i + 1
  }
  
  cat("Completed D =", D, "\n")
}

# combine all rows into a single tibble, dropping any skipped (NULL) entries
results_df <- bind_rows(results)


# -------- save the results --------

dir.create(here("outputs"), showWarnings = FALSE)
saveRDS(results_df, here("outputs", "compositional_bias_results.rds"))

cat("\nResults saved to outputs/compositional_bias_results.rds (",
    nrow(results_df), "rows ).\n")


# -------- heatmap: MAE(D, P) for L1 and CLR, analogous to Figure 2A --------

# reshape to long format for faceted plotting
results_long <- results_df %>%
  select(D, P_target, mae_L1, mae_CLR) %>%
  pivot_longer(cols = c(mae_L1, mae_CLR), names_to = "method", values_to = "mae") %>%
  mutate(method = recode(method, mae_L1 = "L1", mae_CLR = "CLR"))

p_heatmap <- ggplot(results_long, aes(x = D, y = P_target, fill = mae)) +
  geom_tile() +
  facet_wrap(~method) +
  scale_fill_viridis_c(name = "MAE", trans = "log10") +
  theme_bw() +
  xlab("Dimensionality (D)") +
  ylab("Within-dataset diversity (P)") +
  ggtitle("Compositional bias of L1 vs CLR normalization",
          subtitle = "Reduced-scale replication of Fuschi et al. (2025), Figure 2A")

png(here("outputs", "compositional_bias_heatmap.png"), width = 2400, height = 1200, res = 300)
print(p_heatmap)
dev.off()

cat("Heatmap saved to outputs/compositional_bias_heatmap.png\n")