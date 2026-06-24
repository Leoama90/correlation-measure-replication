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


# -------- create new folder --------

dir.create(here("script", "sparsity_effects", "trials"), showWarnings = FALSE, recursive = TRUE)


# -------- read data --------

# Load OTU abundance matrix (samples x OTUs)
otu <- readRDS(list.files(
  path       = here(),
  pattern    = "otu_HMP2.rds",
  full.names = TRUE,
  recursive  = TRUE
))
# Load sample metadata
meta <- readRDS(list.files(
  path       = here(),
  pattern    = "meta_HMP2.rds",
  full.names = TRUE,
  recursive  = TRUE
))


# -------- filter data --------

# Subset samples belonging to subject 69-001 in the healthy state
otu_69001_H <- otu[meta$SubjectID == "69-001" & meta$CL4_2 == "Healthy", ]

# Remove OTUs present in fewer than 33% of samples (low prevalence)
otu_filt <- otu_69001_H[, colSums(otu_69001_H > 0)/nrow(otu_69001_H) >= 0.33]
# Further remove OTUs whose median non-zero count is below 5 (too rare)
otu_filt <- otu_filt[, apply(otu_filt, 2, function(x) median(x[x > 0]) >= 5)]
# Free memory by removing objects no longer needed
rm(otu, otu_69001_H, meta)


# -------- create the fits --------

# Fit Hurdle Truncated Log-Normal (htrlnorm) parameters to each filtered OTU,
# then summarise their distribution by extracting the 10th and 90th percentiles
HMP2_quantile_params <- apply(otu_filt, 2, function(x){
  ToyModel::mle.htrlnorm(x)$estimate
}) %>% apply(1, function(x) quantile(x, probs = c(0.1, 0.9)))

# Build a data frame pairing each OTU's log-mean with its log-maximum,
# used to learn a realistic mean–max relationship from the real data
df_mean_max <- tibble(
  "y" = apply(log(otu_filt + 1), 2, max ),
  "x" = apply(log(otu_filt + 1), 2, mean)
)

# Fit a simple linear model: log(max) ~ log(mean)
model <- lm(y~x, data = df_mean_max)

# Function to predict new maximum values (b) from given mean values,
# adding prediction uncertainty via random draws from the error distribution
predict_max <- function(new_xs){
  
  # Predict fitted values at the requested mean values
  new_data  <- data.frame(x = new_xs)
  predicted_values <- predict(model, new_data, interval = "none")
  
  # Estimate the residual variance of the fitted model
  residuals_variance <- sum(residuals(model)^2) / model$df.residual
  # Number of observations used to fit the original model
  n     <- length(model$model$x)
  # Mean of the original predictor variable
  x_bar <- mean(model$model$x)
  
  # Compute leverage for each new point (how far it is from the training mean)
  leverages <- 1/n + ((new_data$x - x_bar)^2 / sum((model$model$x - x_bar)^2))
  
  # Compute the full prediction standard error (model + residual uncertainty)
  std_error_prediction <- sqrt(residuals_variance * (1 + leverages))
  
  # Sample a random predicted maximum for each new mean, adding realistic noise
  simulated_ys <- rnorm(nrow(new_data), mean = predicted_values, sd = std_error_prediction)
  return(simulated_ys)
}

# Set the total number of outer iterations and initialise a progress bar
n_iteration <- 100
pb <- progress_bar$new(
  format = "[:bar] :elapsed | eta: :eta",
  total  = n_iteration * 7600,
  width  = 60)
# Callback function used to tick the progress bar from parallel workers
progress <- function(n){pb$tick()}

#set sedd for reproducibility
set.seed(42)
result <- data.frame()
for(iter in 1:n_iteration){
  
  
  # -------- OUTER LOOP --------
  
  # Sample random htrlnorm parameters for 200 background OTUs
  # uniformly within the HMP2 10th–90th percentile range
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
  
  # Predict a realistic maximum count (b) for each background OTU from its meanlog
  params_random_HMP2$b <- predict_max(params_random_HMP2$meanlog)
  
  # Simulate n = 10,000 samples for the 200 background OTUs with zero correlation
  # (identity matrix), using the htrlnorm quantile function (NorTA approach)
  random_HMP2 <- toy_model(n     = 10^4, 
                           cor   = diag(200), 
                           M     = 1,
                           qdist = ToyModel::qhtrlnorm,
                           param = params_random_HMP2)
  # Store the Pearson correlation in normal space for reference
  random_cor_zero_HMP2 <- random_HMP2$cor_normal
  
  # Build a full factorial grid of sparsity levels (phi_1, phi_2) and
  # latent correlations (cor), then attach random HMP2-calibrated meanlog and sdlog
  params_set <- expand_grid(
    # zero-inflation probability of OTU 1
    "phi_1" = seq(  0  , 0.95, by = 0.05),   
    # zero-inflation probability of OTU 2
    "phi_2" = seq(  0  , 0.95, by = 0.05),
    # target latent Pearson correlation
    "cor"   = seq(- 0.9, 0.9 , by = 0.1 )    
  ) %>%
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
    # Add a predicted maximum count for each of the two OTUs in the pair
    mutate(b_1 = predict_max(meanlog_1),
           b_2 = predict_max(meanlog_2))
  
  # Initialise a parallel cluster with 6 workers
  cl <- makeCluster(6)
  registerDoSNOW(cl)
  
  
  # -------- INNER LOOP (parallel) --------
  
  df <- foreach(i = 1:nrow(params_set), 
                .combine  = "rbind",
                .packages = c("ToyModel"),
                .options.snow = list(progress = progress)) %dopar% {
                  
                # Generate a pair of correlated OTUs with the i-th parameter set
                couple <- toy_model(n   = 10^4, 
                                    cor = params_set[i, "cor"], 
                                    M   = 1,
                                    qdist = ToyModel::qhtrlnorm,
                                    param = data.frame(
                                    "meanlog" = c(params_set[i, "meanlog_1"], params_set[i, "meanlog_2"]),
                                    "sdlog"   = c(params_set[i, "sdlog_1"],   params_set[i, "sdlog_2"  ]),
                                    "phi"     = c(params_set[i, "phi_1"],     params_set[i, "phi_2"    ])
                                    ))
                  
                  # Copy the background NorTA matrix and replace columns 25 and 125
                  # with the two simulated OTUs of interest
                  random_HMP2_NorTA_i <- random_HMP2$NorTA
                  random_HMP2_NorTA_i[, 25]  <- couple$NorTA[, 1]
                  random_HMP2_NorTA_i[, 125] <- couple$NorTA[, 2]
                  
                  # Replace exact zeros with small random pseudo-values in [0.065, 0.65]
                  # to avoid log(0) issues before CLR transformation
                  rand_pseudo <-  matrix(runif(length(random_HMP2_NorTA_i), min = 0.065, max = 0.65), 
                  nrow = nrow(otu_filt), 
                  ncol = ncol(otu_filt))
                  
                  random_HMP2_NorTA_i <- random_HMP2_NorTA_i + (random_HMP2_NorTA_i == 0) * rand_pseudo
                  
                  # Apply the Centred Log-Ratio (CLR) transform and compute the
                  # full Pearson correlation matrix across all 200 OTUs
                  cor_PCLR <- random_HMP2_NorTA_i %>% ToyModel::clr() %>% cor
                  
                  # Return a one-row summary with all input params and the estimated
                  # CLR correlation between the two target OTUs (cols 25 and 125)
                  data.frame("iteration"      = iter,
                             "cor_input"      = params_set[i, "cor"],
                             "cor_normal"     = couple$cor_normal[1, 2],
                             "meanlog_1"      = params_set[i, "meanlog_1"],
                             "meanlog_2"      = params_set[i, "meanlog_2"],
                             "sdlog_1"        = params_set[i, "sdlog_1"],
                             "sdlog_2"        = params_set[i, "sdlog_2"],
                             "b_1"            = params_set[i, "b_1"],
                             "b_2"            = params_set[i, "b_2"],
                             "phi_1"          = params_set[i, "phi_1"],
                             "phi_2"          = params_set[i, "phi_2"],
                             "cor_NorTA_PCLR" = cor_PCLR[25, 125])
                }
  
  
  # -------- END INNER LOOP --------
  
  # shut down the parallel cluster after each outer iteration
  stopCluster(cl)
  # append this iteration's results to the accumulator data frame
  result <- bind_rows(result, df)
  
  
  # -------- END OUTER LOOP --------
}


# -------- save results in a .rds file --------

# Save the complete results table to disk as an R serialised object
saveRDS(result, here("script", "sparsity_effects", "trials", "Sparsity_Effects_htrlnorm_rand_pseudo.rds"))