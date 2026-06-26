# B0_computer_error_htrlnorm.R
#
# Purpose:
#   Study how sparsity affects CLR correlation estimation, using a hurdle
#   truncated log-normal distribution as the marginal model for OTU counts.
#
# Input:
#   - HMP2 OTU abundance table
#   - HMP2 sample metadata
#   - The script subsets subject 69-001 in healthy status
#   - OTU-level hurdle truncated log-normal parameters are fitted from the
#     filtered HMP2 data, and the 10th–90th percentile range is used to
#     avoid outlier parameter values
#
#   Additional model inputs are generated internally:
#   - a regression model relating mean and maximum log-abundance across OTUs
#   - random parameter sets sampled within the observed HMP2 range
#   - a background compositional dataset of 200 uncorrelated OTUs
#
# Outputs:
#   - Sparsity_Effects_htrlnorm.rds
#
#   The script simulates OTU pairs across a grid of sparsity and input
#   correlation values, inserts them into a background compositional
#   context, applies CLR transformation, and records the resulting
#   Pearson correlation.
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

# Read the OTU counts and subject metadata
otu <- readRDS(list.files(
  path       = here(),
  pattern    = "otu_HMP2.rds",
  full.names = TRUE,
  recursive  = TRUE
))
meta <- readRDS(list.files(
  path       = here(),
  pattern    = "meta_HMP2.rds",
  full.names = TRUE,
  recursive  = TRUE
))

# select samples belonging to subject 69-001 in healthy status
otu_69001_H <- otu[meta$SubjectID == "69-001" & meta$CL4_2 == "Healthy", ]

# remove rarest OTUs: keep only those present in at least 33% of samples
otu_filt <- otu_69001_H[, colSums(otu_69001_H > 0)/nrow(otu_69001_H) >= 0.33]

# keep only OTUs whose median non-zero abundance is at least 5
otu_filt <- otu_filt[, apply(otu_filt, 2, function(x) median(x[x > 0]) >= 5)]

# remove objects no longer needed to free memory
rm(otu, otu_69001_H, meta)

# fit ZINB parameters from real data for each filtered OTU and
# save the 10th and 90th percentile of the parameter distributions
HMP2_quantile_params <- apply(otu_filt, 2, function(x){
  ToyModel::mle.htrlnorm(x)$estimate
}) %>% apply(1, function(x) quantile(x, probs = c(0.1, 0.9)))

# fit a linear model between mean and max abundance (log scale) across OTUs
df_mean_max <- tibble(
  "y" = apply(log(otu_filt + 1), 2, max),
  "x" = apply(log(otu_filt + 1), 2, mean)
)
model <- lm(y~x, data = df_mean_max)

# function to predict new maximum values from mean values,
# accounting for prediction uncertainty via leverage-based standard errors
predict_max <- function(new_xs){
  
  # predict mean values for new inputs
  new_data <- data.frame(x = new_xs)
  predicted_values   <- predict(model, new_data, interval = "none")
  
  # calculate residual variance from the fitted model
  residuals_variance <- sum(residuals(model)^2) / model$df.residual
  
  # number of observations in the original training data
  n   <- length(model$model$x)
  
  # mean of the original independent variable
  x_bar <- mean(model$model$x)
  
  # compute leverage for each new point (influence on the regression line)
  # hii = 1/n + (xi - x̄)^2 / Σ(xi - x̄)^2
  leverages <- 1/n + ((new_data$x - x_bar)^2 / sum((model$model$x - x_bar)^2))
  
  # compute standard error of prediction for each new input
  std_error_prediction <- sqrt(residuals_variance * (1 + leverages))
  
  # sample predicted maximum values adding gaussian noise
  simulated_ys <- rnorm(nrow(new_data), mean = predicted_values, sd = std_error_prediction)
  return(simulated_ys)
}

# set total number of outer iterations
nIteration <- 100

# initialize progress bar
pb <- progress_bar$new(
  format = "[:bar] :elapsed | eta: :eta",
  total  = nIteration * 7600,
  width  = 60)
progress <- function(n){pb$tick()}

# set seed for reproducibility
set.seed(42)
result <- data.frame()
for(iter in 1:nIteration){
  
  # -------- START OUTER LOOP --------

  
  # randomly sample ZINB parameters for 200 OTUs within the observed HMP2 range
  params_random_HMP2 <- data.frame(
    "meanlog" = runif(n = 200,
                    min = HMP2_quantile_params["10%", "meanlog"],
                    max = HMP2_quantile_params["90%", "meanlog"]),
    "sdlog"   = runif(n = 200,
                    min = HMP2_quantile_params["10%", "sdlog"],
                    max = HMP2_quantile_params["90%", "sdlog"]),
    "phi"     = runif(n = 200,
                    min = HMP2_quantile_params["10%", "phi"],
                    max = HMP2_quantile_params["90%", "phi"]))
  
  # predict maximum abundance values for the randomly sampled OTUs
  params_random_HMP2$b <- predict_max(params_random_HMP2$meanlog)
  
  # simulate a background dataset of 200 uncorrelated OTUs (identity correlation matrix)
  random_HMP2 <- ToyModel::toy_model(n = 10^4, cor = diag(200), M = 1,
                                     qdist = ToyModel::qhtrlnorm,
                                     param = params_random_HMP2)
  # store the underlying normal correlation matrix of the background dataset
  # random_cor0_HMP2 <- random_HMP2$cor_normal
  
  # generate all combinations of sparsity (phi) and correlation values to test
  params_set <- expand_grid(
    "phi_1" = seq( 0  , 0.95,  by = 0.05),
    "phi_2" = seq( 0  , 0.95,  by = 0.05),
    "cor"   = seq(-0.9, 0.9 ,  by = 0.1 )
  ) %>%
    # assign random log-normal parameters to each combination
    mutate("meanlog_1" = runif(n = n(),
                             min = HMP2_quantile_params["10%", "meanlog"],
                             max = HMP2_quantile_params["90%", "meanlog"]),
           "meanlog_2" = runif(n = n(),
                             min = HMP2_quantile_params["10%", "meanlog"],
                             max = HMP2_quantile_params["90%", "meanlog"]),
           "sdlog_1"   = runif(n = n(),
                             min = HMP2_quantile_params["10%", "sdlog"],
                             max = HMP2_quantile_params["90%", "sdlog"]),
           "sdlog_2"   = runif(n = n(),
                             min = HMP2_quantile_params["10%", "sdlog"],
                             max = HMP2_quantile_params["90%", "sdlog"])) %>%
    as.data.frame() %>%
    # predict maximum values for each OTU pair
    mutate(b_1 = predict_max(meanlog_1),
           b_2 = predict_max(meanlog_2))
  
  # initialize a parallel cluster with 6 workers
  cl <- makeCluster(6)
  registerDoSNOW(cl)
  
  # -------- START INNER LOOP --------

  
  # iterate over all parameter combinations in parallel
  df <- foreach(i = 1:nrow(params_set), 
                .combine      = "rbind",
                .packages     = c("ToyModel"),
                .options.snow = list(progress = progress)) %dopar% {
                  
                  # simulate a pair of OTUs with the given sparsity and correlation
                  couple <-
                    ToyModel::toy_model(n = 10^4, 
                                        cor = params_set[i, "cor"], 
                                        M = 1,
                                        qdist = ToyModel::qhtrlnorm,
                                        param = data.frame(
                                        "meanlog" = c(params_set[i, "meanlog_1"], params_set[i, "meanlog_2"]),
                                        "sdlog"   = c(params_set[i, "sdlog_1"],   params_set[i, "sdlog_2"]),
                                        "phi"     = c(params_set[i, "phi_1"],     params_set[i, "phi_2"])
                                        ))
                  
                  # inject the simulated OTU pair into the background dataset
                  # at positions 25 and 125
                  random_HMP2_NorTA_i        <- random_HMP2$NorTA
                  random_HMP2_NorTA_i[, 25]  <- couple$NorTA[, 1]
                  random_HMP2_NorTA_i[, 125] <- couple$NorTA[, 2]
                  
                  # apply CLR transformation and compute Pearson correlation matrix
                  cor_PCLR <- random_HMP2_NorTA_i %>% ToyModel::clr() %>% cor
                  
                  # store results for this parameter combination
                  data.frame("iteration" = iter,
                             "cor_input" = params_set[i, "cor"],
                             "cor_normal"= couple$cor_normal[1, 2],
                             "meanlog_1" = params_set[i, "meanlog_1"],
                             "meanlog_2" = params_set[i, "meanlog_2"],
                             "sdlog_1"   = params_set[i, "sdlog_1"],
                             "sdlog_2"   = params_set[i, "sdlog_2"],
                             "b_1"       = params_set[i, "b_1"],
                             "b_2"       = params_set[i, "b_2"],
                             "phi_1"     = params_set[i, "phi_1"],
                             "phi_2"     = params_set[i, "phi_2"],
                             "cor_NorTA_PCLR" = cor_PCLR[25, 125])
                }
  # -------- END INNER LOOP --------
  
  
  # shut down the parallel cluster
  stopCluster(cl)
  # append results of this iteration to the global result dataframe
  result <- bind_rows(result, df)
  
  # -------- END OUTER LOOP --------
  
}
# save the final results to disk
saveRDS(result, here("script", "sparsity_effects", "Sparsity_Effects_htrlnorm.rds"))