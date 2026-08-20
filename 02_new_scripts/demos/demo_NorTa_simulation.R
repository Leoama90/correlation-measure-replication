# demo_NorTa_simulation.R
#
# Purpose:
#   Demonstrates the usage of norta_simulation() (defined in
#   NorTa_simulation.R), walking step by step through the full NorTA
#   pipeline: generating a known correlation structure, simulating
#   correlated normal data from it, and mapping that data through a
#   zero-inflated negative binomial (ZINB) distribution to obtain a
#   sparse count dataset. Also verifies, empirically, the two key
#   properties the pipeline relies on: that pnorm() preserves the rank
#   correlation of the underlying normal data, and that the final count
#   data carries a correlation signal recognizably related to R_true.
#
# Input:
#   - no external files required; all parameters are defined inline
#
# Output:
#   - printed explanation of each pipeline step, intermediate objects,
#     and empirical checks of the correlation structure at each stage
#
# Used scripts:
#   - NorTa_simulation.R (which itself sources generate_matrix_factors.R)

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# -------- load the NorTa_simulation.R script --------
source(
  list.files(
    path = here(),
    pattern = "^NorTa_simulation\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)


# -------- step 0: what this demo is about --------

cat("#-------------------------------------------------------#\n")
cat("This demo simulates a sparse metagenomic count dataset with a\n")
cat("KNOWN correlation structure, using the NorTA (Normal To Anything)\n")
cat("approach. Since we control the ground-truth correlation matrix,\n")
cat("we can later check how well any correlation-estimation method\n")
cat("recovers it from the (sparse, non-normal) count data alone.\n")
cat("#-------------------------------------------------------#\n\n")


# -------- step 1: run the full pipeline --------

cat("Simulating 6 taxa, split into 2 correlation groups, over 300\n")
cat("samples, with default ZINB parameters (mu = 20, size = 30,\n")
cat("phi = 0.3, i.e. 30% extra zero-inflation on top of the negative\n")
cat("binomial's own dispersion).\n\n")

result <- norta_simulation(n = 6, n_groups = 2, N = 300, seed = 42)

cat("\nnorta_simulation() returns a list with 5 elements:\n")
cat("  - R_true:     the known n x n correlation matrix used to simulate\n")
cat("  - groups:     which of the n_groups each taxon was assigned to\n")
cat("  - n_zeroes:   exact number of zero correlation pairs in R_true\n")
cat("  - sim_data:   intermediate N x n correlated normal data\n")
cat("  - sim_counts: final N x n sparse count data (the simulated dataset)\n")


# -------- step 2: inspect the ground-truth correlation structure --------

cat("\n#-------------------------------------------------------#\n")
cat("Step 1 recap: the ground-truth correlation matrix (R_true).\n")
cat("This is what we are trying to recover later from sim_counts alone.\n")
cat("#-------------------------------------------------------#\n\n")

cat("Group assignment:\n")
print(result$groups)

cat("\nR_true (rounded):\n")
print(round(result$R_true, 2))


# -------- step 3: correlated normal data (NorTA step 1) --------

cat("\n#-------------------------------------------------------#\n")
cat("Step 2: rmvnorm() draws 300 samples from a multivariate normal\n")
cat("distribution with covariance R_true. Since R_true has a unit\n")
cat("diagonal, each simulated taxon (column) is individually a\n")
cat("standard N(0,1) variable, but the columns are correlated with\n")
cat("each other exactly as specified in R_true.\n")
cat("#-------------------------------------------------------#\n\n")

cat("Sample correlation of sim_data (should be close to R_true, since\n")
cat("it is estimated from only 300 samples, not the infinite-sample\n")
cat("limit):\n")
print(round(cor(result$sim_data), 2))

cat("\nMean absolute difference between the sample correlation of\n")
cat("sim_data and the true R_true (should be small, shrinking further\n")
cat("as N increases):\n")
cat(round(mean(abs(cor(result$sim_data) - result$R_true)), 4), "\n")


# -------- step 4: uniform ranks via pnorm() (NorTA step 2a) --------

cat("\n#-------------------------------------------------------#\n")
cat("Step 3: pnorm() converts each normal value to its percentile,\n")
cat("i.e. a value in [0, 1]. This is the 'Normal To Anything' trick:\n")
cat("pnorm() is a MONOTONIC transformation, so it changes the shape of\n")
cat("each taxon's marginal distribution, but preserves the RANK order\n")
cat("of values within each taxon, and therefore preserves the rank\n")
cat("correlation between taxa (though not necessarily the exact\n")
cat("Pearson correlation, since that is scale-sensitive).\n")
cat("#-------------------------------------------------------#\n\n")

sim_unif <- pnorm(result$sim_data)

cat("Range of sim_unif values (should be within [0, 1]):\n")
cat("min =", round(min(sim_unif), 4), ", max =", round(max(sim_unif), 4), "\n")

cat("\nSpearman (rank) correlation of sim_unif vs sim_data, for the\n")
cat("first pair of taxa (should match almost exactly, since pnorm()\n")
cat("is monotonic and rank correlation only depends on ordering):\n")
cat("  sim_data:  ", round(cor(result$sim_data[, 1], result$sim_data[, 2], method = "spearman"), 4), "\n")
cat("  sim_unif:  ", round(cor(sim_unif[, 1], sim_unif[, 2], method = "spearman"), 4), "\n")


# -------- step 5: sparse ZINB counts (NorTA step 2b) --------

cat("\n#-------------------------------------------------------#\n")
cat("Step 4: qzinegbin() maps each uniform value to a count, using the\n")
cat("zero-inflated negative binomial quantile function. Every taxon\n")
cat("shares the same distribution parameters here (mu, size, phi), so\n")
cat("sparsity is uniform across taxa in this simulation (unlike\n")
cat("data_sim_ph_driven.R, where sparsity varies per taxon).\n")
cat("#-------------------------------------------------------#\n\n")

cat("First 6 rows of the final simulated count data (sim_counts):\n")
print(head(result$sim_counts))

observed_zero_rate <- mean(result$sim_counts == 0)
cat("\nObserved overall zero rate in sim_counts:", round(observed_zero_rate, 3), "\n")
cat("(expected to be somewhat above phi = 0.3, since the negative\n")
cat("binomial component can itself also produce zero counts)\n")


# -------- step 6: does the correlation signal survive to the counts? --------

cat("\n#-------------------------------------------------------#\n")
cat("Step 5: checking whether the known correlation structure is still\n")
cat("recognizable in the final sparse count data, using plain Pearson\n")
cat("correlation directly on sim_counts (no CLR or other correction).\n")
cat("This is the 'naive' baseline the rest of the project's pipeline\n")
cat("(clr_pearson.R) is meant to improve upon.\n")
cat("#-------------------------------------------------------#\n\n")

R_naive <- cor(result$sim_counts)

cat("Naive Pearson correlation on raw sim_counts (rounded):\n")
print(round(R_naive, 2))

cat("\nMean absolute error between this naive estimate and R_true:\n")
cat(round(mean(abs(R_naive - result$R_true)), 4), "\n")

cat("\nThis error is expected to be non-trivial: the ZINB mapping is a\n")
cat("non-linear transformation, so Pearson correlation on raw counts\n")
cat("does not exactly preserve R_true, even before considering the\n")
cat("compositional bias that CLR is specifically designed to correct.\n")