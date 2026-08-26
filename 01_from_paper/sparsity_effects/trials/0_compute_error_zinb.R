# 0_computer_error_zinb.R
#
# Purpose:
#   ZINB variant of the asymmetric sparsity simulation, where pstr0_1 and
#   pstr0_2 vary independently, mirroring the structure of the htrlnorm
#   random pseudo-count trial from the previous script.
#
# Input:
#   - HMP2 OTU abundance table
#   - HMP2 sample metadata
#   - The script subsets subject 69-001 in healthy status
#   - OTUs are filtered by prevalence and abundance
#   - ZINB parameters are fitted to each filtered OTU
#   - The 10th–90th percentile range of fitted parameters is used to define
#     a realistic simulation space
#
#   Compared with earlier ZINB scripts:
#   - pstr0_1 and pstr0_2 vary independently
#   - background OTUs are sampled with uniform draws within the empirical range
#
# Output:
#   - A results file saved in:
#     here("01_from_paper", "sparsity_effects", "trials", "Sparsity_Effects_zinbin_rare.rds")
#
# doSNOW: parallel backend for foreach, supports progress bars via snow clusters
# https://cran.r-project.org/web/packages/doSNOW/index.html
library(doSNOW)

# foreach: parallel foreach loops
# https://cran.r-project.org/package=foreach
library(foreach)

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# progress: displays text progress bars for long-running loops
# https://cran.r-project.org/web/packages/progress/index.html
library(progress)

# tidyverse: useful to manage data (dplyr) and make nice plots (ggplot2)
# https://cran.r-project.org/web/packages/tidyverse/index.html
library(tidyverse)

# ToyModel: lightweight package with necessary function to generate metagenomics data
# https://github.com/Fuschi/ToyModel
library(ToyModel)

# VGAM: vector generalized linear models
# https://cran.r-project.org/package=VGAM
library(VGAM)

# SpiecEasi: package with spiec.easi and sparCC methods
# https://github.com/zdk123/SpiecEasi
library(SpiecEasi)



# -------- create new folder --------

dir.create(here("01_from_paper", "sparsity_effects", "trials"), showWarnings = FALSE, recursive = TRUE)


# -------- read data --------

# Load OTU abundance matrix (samples x OTUs)
otu <- readRDS(list.files(
  path       = here(),
  pattern    = "^otu_HMP2\\.rds$",
  full.names = TRUE,
  recursive  = TRUE
))
# Load sample metadata
meta <- readRDS(list.files(
  path       = here(),
  pattern    = "^meta_HMP2\\.rds$",
  full.names = TRUE,
  recursive  = TRUE
))


# -------- filter data --------

# Subset samples belonging to subject 69-001 in the healthy state
otu_69001_H <- otu[meta$SubjectID == "69-001" & meta$CL4_2 == "Healthy", ]

# Remove OTUs present in fewer than 33% of samples (low prevalence)
otu_filt <- otu_69001_H[, colSums(otu_69001_H > 0) / nrow(otu_69001_H) >= 0.33]
# Further remove OTUs whose median non-zero count is below 5 (too rare)
otu_filt <- otu_filt[, apply(otu_filt, 2, function(x) median(x[x > 0]) >= 5)]
# Free memory by removing objects no longer needed
rm(otu, otu_69001_H, meta)

# Fit a Zero-Inflated Negative Binomial (ZINB) distribution to each filtered OTU,
# then summarize the fitted parameters across OTUs by computing
# the 10th and 90th percentiles — these define a realistic parameter range for simulation
HMP2_quantile_params <- apply(otu_filt, 2, function(x) {
  SpiecEasi::fitdistr(as.numeric(x), "zinegbin")$par
}) %>% apply(1, function(x) quantile(x, probs = c(0.1, 0.9)))

# Initialize a progress bar to track the total number of parallel iterations
nIteration <- 100
pb <- progress_bar$new(
  format = "[:bar] :elapsed | eta: :eta",
  total  = nIteration * 7600,
  width  = 60
)
# Define the progress callback function used by doSNOW workers
progress <- function(n) {
  pb$tick()
}

set.seed(42)
result <- data.frame()


# -------- OUTER LOOP: repeat the full experiment nIteration times to assess variability --------

for (iter in 1:nIteration) {
  # Draw random ZINB parameters for 200 background OTUs,
  # sampling uniformly within the empirical [10%, 90%] range estimated from HMP2
  params_random_HMP2 <- data.frame(
    # Mean of the Negative Binomial component
    "munb" = runif(
      n = 200,
      min = HMP2_quantile_params["10%", "munb"],
      max = HMP2_quantile_params["90%", "munb"]
    ),
    # Dispersion (size) parameter of the Negative Binomial component
    "size" = runif(
      n = 200,
      min = HMP2_quantile_params["10%", "size"],
      max = HMP2_quantile_params["90%", "size"]
    ),
    # Zero-inflation probability for each OTU
    "pstr0" = runif(
      n = 200,
      min = HMP2_quantile_params["10%", "pstr0"],
      max = HMP2_quantile_params["90%", "pstr0"]
    )
  )

  # Simulate a background dataset of 10,000 samples x 200 OTUs
  # using the NorTA method (Normal-To-Anything) with an identity correlation matrix
  # (i.e., all OTUs are truly uncorrelated in the latent Gaussian space)
  random_HMP2 <- toy_model(
    n = 10^4,
    cor = diag(200),
    M = 1,
    qdist = qzinegbin,
    param = params_random_HMP2
  )
  # Store the latent (Gaussian) correlation matrix of the background dataset
  random_cor0_HMP2 <- random_HMP2$cor_normal

  # Build the full grid of parameter combinations for the supervised pair of OTUs:
  # all combinations of zero-inflation levels (pstr0_1, pstr0_2) and true correlations (cor)
  params_set <- expand_grid(
    "pstr0_1" = seq(0, 0.95, by = 0.05), # Zero-inflation of OTU 1: from 0% to 95%
    "pstr0_2" = seq(0, 0.95, by = 0.05), # Zero-inflation of OTU 2: from 0% to 95%
    "cor"     = seq(-0.9, 0.9, by = 9.1) # True latent correlation: from -0.9 to 0.9
  ) %>%
    # Randomly assign NB mean and dispersion to each combination (within HMP2 range)
    mutate(
      "munb_1" = runif(
        n = n(),
        min = HMP2_quantile_params["10%", "munb"],
        max = HMP2_quantile_params["90%", "munb"]
      ),
      "munb_2" = runif(
        n = n(),
        min = HMP2_quantile_params["10%", "munb"],
        max = HMP2_quantile_params["90%", "munb"]
      ),
      "size_1" = runif(
        n = n(),
        min = HMP2_quantile_params["10%", "size"],
        max = HMP2_quantile_params["90%", "size"]
      ),
      "size_2" = runif(
        n = n(),
        min = HMP2_quantile_params["10%", "size"],
        max = HMP2_quantile_params["90%", "size"]
      )
    ) %>%
    as.data.frame()

  # Set up a parallel cluster with 6 workers for the inner loop
  cl <- makeCluster(6)
  registerDoSNOW(cl)


  # -------- INNER LOOP: process each parameter combination in parallel --------

  df <- foreach(
    i = 1:nrow(params_set),
    .combine = "rbind",
    .packages = c("ToyModel"),
    .options.snow = list(progress = progress)
  ) %dopar% {
    # Simulate a pair of correlated OTUs with the i-th parameter combination,
    # using the specified true latent correlation and ZINB marginal distributions
    couple <- toy_model(
      n = 10^4,
      cor = params_set[i, "cor"],
      M = 1,
      qdist = VGAM::qzinegbin,
      param = data.frame(
        "munb"  = c(params_set[i, "munb_1"], params_set[i, "munb_2"]),
        "size"  = c(params_set[i, "size_1"], params_set[i, "size_2"]),
        "pstr0" = c(params_set[i, "pstr0_1"], params_set[i, "pstr0_2"])
      )
    )

    # Inject the simulated pair into the background dataset:
    # replace columns 25 and 125 of the background NorTA matrix
    # with the two OTUs of the supervised pair
    random_HMP2_NorTA_i <- random_HMP2$NorTA
    random_HMP2_NorTA_i[, 25] <- couple$NorTA[, 1]
    random_HMP2_NorTA_i[, 125] <- couple$NorTA[, 2]

    # Apply the CLR (Centered Log-Ratio) transformation to the full count matrix
    # to account for compositionality, then compute Pearson correlations
    cor_PCLR <- random_HMP2_NorTA_i %>%
      ToyModel::clr() %>%
      cor()

    # Return a one-row data frame with all input parameters and the estimated correlation
    data.frame(
      "iteration" = iter,
      # True input correlation in the latent Gaussian space
      "cor_input" = params_set[i, "cor"],
      # Realised latent correlation between the pair (may differ slightly from input)
      "cor_normal" = couple$cor_normal[1, 2],
      "munb_1" = params_set[i, "munb_1"],
      "munb_2" = params_set[i, "munb_2"],
      "size_1" = params_set[i, "size_1"],
      "size_2" = params_set[i, "size_2"],
      "pstr0_1" = params_set[i, "pstr0_1"],
      "pstr0_2" = params_set[i, "pstr0_2"],
      # CLR-based Pearson correlation between the two injected OTUs
      "cor_NorTA_PCLR" = cor_PCLR[25, 125]
    )
  }


  # -------- END INNER LOOP --------

  # Shut down the parallel cluster after the inner loop completes
  stopCluster(cl)
  # Accumulate results from this iteration into the master data frame
  result <- bind_rows(result, df)
}


# -------- END OUTER LOOP --------


# Save the final results to disk
saveRDS(result, here("script", "sparsity_effects", "trials", "Sparsity_Effects_zinbin_rare.rds"))

# UI reminder where to search the generated plot
cat("the .png file has been saved in the 'trials' folder with the name 'Sparsity_Effects_zinbin_rare.rds'")
