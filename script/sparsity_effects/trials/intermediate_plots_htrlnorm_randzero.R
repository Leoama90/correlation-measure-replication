# cowplot: draw ggplot2 in new figures
# https://cran.r-project.org/web/packages/cowplot/index.html
library(cowplot)

# ggpubr: nice plots based on ggplot2
# https://cran.r-project.org/web/packages/ggpubr/index.html
library(ggpubr)

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# tidyverse: useful to manage data (dplyr) and make nice plots (ggplot2)
# https://cran.r-project.org/web/packages/tidyverse/index.html
library(tidyverse)


# -------- DATA LOADING AND PREPARATION --------

# Recursively search the project directory for the simulation results file
path <- list.files(
  path       = here(),
  pattern    = "Sparsity_Effects_htrlnorm.rds",
  full.names = TRUE,
  recursive  = TRUE
)

# Load the RDS file and compute derived summary columns
df <- readRDS(path) %>%
  # Compute absolute error between the normal and CLR-transformed correlation
  mutate(err_clr  = abs(cor_normal - cor_NorTA_PCLR)) %>%
  # Compute the average sparsity across the two OTUs
  mutate(phi_mean = 0.5 * (phi_1 + phi_2)) %>%
  # Compute the maximum sparsity across the two OTUs
  mutate(phi_max  = pmax(phi_1, phi_2)) %>%
  # Compute the minimum sparsity across the two OTUs
  mutate(phi_min  = pmin(phi_1, phi_2))

# Sanity check: absolute error must never exceed 1
if (any(df$err_clr > 1)) stop("Error greater than 1 detected in err_clr.")



# -------- colour palette --------


# Build a diverging colour ramp from the 11-colour Spectral palette (reversed)
my_palette <- RColorBrewer::brewer.pal(n = 11, "Spectral") %>%
  rev() %>%
  grDevices::colorRampPalette()


# -------- helper: shared tile plot theme --------


# Common scale and theme layers reused across the three heatmaps
tile_scales <- list(
  # Colour gradient with custom breakpoints to emphasise the 0–0.5 range
  scale_fill_gradientn(
    name    = "Absolute Error",
    colours = my_palette(12),
    values  = c(seq(0, 0.5, by = 0.05), 1),
    limits  = c(0, 1)
  ),
  # Remove tick marks from the colour-bar guide
  guides(fill = guide_colorbar(ticks.colour = NA)),
  # Rotate x-axis labels for readability
  theme(
    legend.text  = element_text(size = 10),
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ),
  # Axis labels
  xlab("Correlation"),
  ylab("Zero %")
)


# -------- heatmaps: CLR error by sparsity summary statistic --------

# -------- Minimum sparsity across the pair --------
# Aggregate mean absolute error by input correlation and minimum phi
p_htrlnorm_min <- df %>%
  summarise(err_clr_mean = mean(err_clr), .by = c(cor_input, phi_min)) %>%
  ggplot(aes(x = cor_input, y = phi_min, fill = err_clr_mean)) +
  geom_tile() +
  theme_bw()  +
  tile_scales +
  # Set axis breaks and remove padding around the tiles
  scale_y_continuous(breaks = seq( 0.1, 0.9, 0.1), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(-0.8, 0.8, 0.2), expand = c(0, 0))

# -------- Mean sparsity across the pair --------
# Aggregate and keep only multiples of 5 % to reduce visual clutter
p_htrlnorm_mean <- df %>%
  summarise(err_clr_mean = mean(err_clr), .by = c(cor_input, phi_mean)) %>%
  filter((100 * phi_mean) %% 5 == 0) %>%
  ggplot(aes(x = cor_input, y = factor(phi_mean), fill = err_clr_mean)) +
  geom_tile() +
  theme_bw() +
  tile_scales +
  # Use discrete y-axis because phi_mean was coerced to factor
  scale_y_discrete(breaks   = seq( 0.1, 0.9, 0.1), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(-0.8, 0.8, 0.2), expand = c(0, 0))

# -------- Maximum sparsity across the pair --------
# Aggregate mean absolute error by input correlation and maximum phi
p_htrlnorm_max <- df %>%
  summarise(err_clr_mean = mean(err_clr), .by = c(cor_input, phi_max)) %>%
  ggplot(aes(x = cor_input, y = phi_max, fill = err_clr_mean)) +
  geom_tile() +
  theme_bw() +
  tile_scales +
  scale_y_continuous(breaks = seq( 0.1, 0.9, 0.1), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(-0.8, 0.8, 0.2), expand = c(0, 0))

# Save the three-panel heatmap figure
png(here("script", "sparsity_effects", "trials", "effects_htrlnorm_zerorand.png"), width = 3600, height = 1500, res = 300)
print(ggarrange(
  p_htrlnorm_min, p_htrlnorm_mean, p_htrlnorm_max,
  ncol         = 3,
  labels       = c("Phi Min", "Phi Mean", "Phi Max"),
  common.legend = TRUE,
  legend       = "bottom"
))
dev.off()


# -------- examples using HMP2 data --------

# -------- Load and filter the HMP2 OTU table --------

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

# Subset to healthy samples from subject 69-001
otu_69001_h <- otu[meta$SubjectID == "69-001" & meta$CL4_2 == "Healthy", ]

# Retain OTUs present in at least 33 % of samples
otu_filt <- otu_69001_h[, colSums(otu_69001_h > 0) / nrow(otu_69001_h) >= 0.33]

# Further restrict to OTUs whose non-zero median count is at least 5
otu_filt <- otu_filt[, apply(otu_filt, 2, function(x) median(x[x > 0]) >= 5)]

# Remove large objects that are no longer needed
rm(otu, otu_69001_h, meta)


# -------- Estimate HMP2 parameter ranges --------

# Fit htrlnorm parameters for each retained OTU and extract the 10th/90th percentiles
hmp2_quantile_params <- apply(otu_filt, 2, function(x) {
  ToyModel::mle.htrlnorm(x)$estimate
}) %>%
  apply(1, function(x) quantile(x, probs = c(0.1, 0.9)))


# -------- Linear model: log-mean → log-max --------

# Build a data frame of per-OTU log-scale mean and maximum values
df_mean_max <- tibble(
  # log-scale maximum
  y = apply(log(otu_filt + 1), 2, max),   
  # log-scale mean
  x = apply(log(otu_filt + 1), 2, mean)   
)

# Fit a simple linear model to predict max from mean (used to simulate b parameters)
model <- lm(y ~ x, data = df_mean_max)

# Predict new maximum values given a vector of log-scale means,
# drawing randomly from the predictive distribution to capture uncertainty
predict_max <- function(new_xs) {
  new_data  <- data.frame(x = new_xs)
  
  # Point predictions from the fitted model
  predicted_values <- predict(model, new_data, interval = "none")
  
  # Residual variance from the model fit
  residuals_variance <- sum(residuals(model)^2) / model$df.residual
  # number of observations used in fitting
  n     <- length(model$model$x)   
  # mean of the training predictor
  x_bar <- mean(model$model$x)     
  
  # Leverage (hat value) for each new observation
  leverages <- 1 / n +
    ((new_data$x - x_bar)^2 / sum((model$model$x - x_bar)^2))
  
  # Predictive standard error (accounts for both residual variance and leverage)
  std_error_prediction <- sqrt(residuals_variance * (1 + leverages))
  
  # Draw one simulated response per new observation
  rnorm(nrow(new_data), mean = predicted_values, sd = std_error_prediction)
}


# -------- Simulation function --------

# Simulate a pair of correlated OTUs embedded in a 200-OTU HMP2-like community,
# then return the Normal, NorTA, and CLR-transformed NorTA representations.
example <- function(n, cor, meanlog_1, meanlog_2, sdlog_1, sdlog_2, phi_1, phi_2) {
  
  # Draw random htrlnorm parameters for 200 background OTUs from the HMP2 ranges
  params_random_hmp2 <- data.frame(
    meanlog = runif(200,
                    min = hmp2_quantile_params["10%", "meanlog"],
                    max = hmp2_quantile_params["90%", "meanlog"]),
    sdlog   = runif(200,
                    min = hmp2_quantile_params["10%", "sdlog"],
                    max = hmp2_quantile_params["90%", "sdlog"]),
    phi     = runif(200,
                    min = hmp2_quantile_params["10%", "phi"],
                    max = hmp2_quantile_params["90%", "phi"])
  ) %>%
    # Simulate the b (upper-truncation) parameter from the log-mean via the linear model
    mutate(b = predict_max(meanlog))
  
  # Simulate 200 independent background OTUs using NorTA
  random_hmp2 <- ToyModel::toy_model(
    n      = n,
    cor    = diag(200),
    M      = 1,
    qdist  = ToyModel::qhtrlnorm,
    param  = params_random_hmp2
  )
  
  # Build the parameter set for the two focal OTUs
  params_set <- data.frame(
    phi     = c(phi_1, phi_2),
    meanlog = c(meanlog_1, meanlog_2),
    sdlog   = c(sdlog_1, sdlog_2)
  ) %>%
    mutate(b = predict_max(meanlog))
  
  # Simulate the correlated pair using the specified correlation structure
  couple <- ToyModel::toy_model(
    n      = n,
    cor    = cor,
    M      = 1,
    qdist  = ToyModel::qhtrlnorm,
    param  = params_set
  )
  
  # Replace the first two columns of the background community with the focal pair
  random_hmp2_norta       <- random_hmp2$NorTA
  random_hmp2_norta[, 1]  <- couple$NorTA[, 1]
  random_hmp2_norta[, 2]  <- couple$NorTA[, 2]
  
  # Add a small uniform pseudo-count to all zeros (avoids log(0) in CLR)
  rand_pseudo <- matrix(
    runif(length(random_hmp2_norta), min = 0.065, max = 0.65),
    nrow  = nrow(random_hmp2_norta),
    ncol  = ncol(random_hmp2_norta)
  )
  random_hmp2_norta <- random_hmp2_norta + (random_hmp2_norta == 0) * rand_pseudo
  
  # Apply CLR transformation and retain only the two focal OTUs
  norta_clr <- random_hmp2_norta %>%
    ToyModel::clr() %>%
    (\(x) x[, 1:2])()
  
  list(
    Normal    = couple$normal,
    NorTA     = couple$NorTA,
    NorTA_CLR = norta_clr
  )
}


# -------- Plotting function --------

# Build a three-panel scatter plot (Normal / NorTA / CLR) for one simulated example
plot_example <- function(ex) {
  
  # Helper: scatter plot with regression line and Pearson r label
  scatter <- function(data, title) {
    data %>%
      as_tibble() %>%
      ggscatter(
        x          = "V1",
        y          = "V2",
        add        = "reg.line",
        add.params = list(color = "blue", fill = "lightgray"),
        conf.int   = TRUE
      ) +
      stat_cor(aes(label = after_stat(r.label)), color = "red") +
      ggtitle(title) +
      theme(axis.title = element_blank())
  }
  
  list(
    norm      = scatter(ex$Normal,    "Normal"),
    NorTA     = scatter(ex$NorTA,     "NorTA"),
    NorTA_CLR = scatter(ex$NorTA_CLR, "NorTA_CLR")
  )
}


# -------- generate the example plots --------

# -------- Example 1: high negative correlation, high sparsity --------
e1 <- example(
  n         = 1000, 
  cor       = -0.9,
  meanlog_1 = 1, 
  meanlog_2 = 3,
  sdlog_1   = 1, 
  sdlog_2   = 1.3,
  phi_1     = 0.85, 
  phi_2     = 0.9
) %>% plot_example()

pall1 <- ggarrange(e1$norm, e1$NorTA, e1$NorTA_CLR, ncol = 3) %>%
  annotate_figure("n = 1000, cor = -0.9, meanlog = [1, 3], sdlog = [1, 1.3], phi = [0.85, 0.9]")

# -------- Example 2: high positive correlation, high sparsity --------
e2 <- example(
  n         = 1000, 
  cor       = 0.9,
  meanlog_1 = 4, 
  meanlog_2 = 2,
  sdlog_1   = 1.5, 
  sdlog_2   = 0.9,
  phi_1     = 0.9, 
  phi_2     = 0.75
) %>% plot_example()

pall2 <- ggarrange(e2$norm, e2$NorTA, e2$NorTA_CLR, ncol = 3) %>%
  annotate_figure("n = 1000, cor = 0.9, meanlog = [4, 2], sdlog = [1.5, 0.9], phi = [0.9, 0.75]")

# -------- Example 3: high negative correlation, low sparsity --------
e3 <- example(
  n         = 1000, 
  cor       = -0.9,
  meanlog_1 = 1.5, 
  meanlog_2 = 5,
  sdlog_1   = 1, 
  sdlog_2   = 1.5,
  phi_1     = 0.1, 
  phi_2     = 0.05
) %>% plot_example()

pall3 <- ggarrange(e3$norm, e3$NorTA, e3$NorTA_CLR, ncol = 3) %>%
  annotate_figure("n = 1000, cor = -0.9, meanlog = [1.5, 5], sdlog=[1, 1.5], phi=[0.1, 0.05]")

# -------- Example 4: high positive correlation, low sparsity --------
e4 <- example(
  n         = 1000, 
  cor       = 0.9,
  meanlog_1 = 2.5, 
  meanlog_2 = 5.5,
  sdlog_1   = 1.7, 
  sdlog_2   = 2,
  phi_1     = 0.1, 
  phi_2     = 0.15
) %>% plot_example()

pall4 <- ggarrange(e4$norm, e4$NorTA, e4$NorTA_CLR, ncol = 3) %>%
  annotate_figure("n = 1000, cor = 0.9, meanlog = [2.5, 5.5], sdlog = [1.7, 2], phi = [0.1, 0.15]")

# Stack all four annotated panels into a single tall figure
pall <- plot_grid(pall1, pall2, pall3, pall4, nrow = 4)


# save the final figure as a high-resolution PNG
png(here("script", "sparsity_effects", "trials", "examples_htrlnorm_zerorand.png"), width = 3600, height = 4800, res = 300)
# allows to generate the .png file
print(pall)
# close the graphics device and write the file to disk
dev.off()

# Sparsity Effects Analysis — htrlnorm
# Examines how sparsity (zero inflation) affects CLR-transformed correlations,
# and generates illustrative examples using simulated HMP2-like microbiome data.