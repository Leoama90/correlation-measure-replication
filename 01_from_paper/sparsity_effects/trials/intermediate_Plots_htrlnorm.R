# intermediate_plots_htrlnorm.R
#
# Purpose:
#   Examines how sparsity (zero inflation) affects the accuracy of
#   CLR-transformed correlations, by summarising simulation results
#   from a two-OTU htrlnorm/NorTA design as heatmaps of mean absolute
#   error (Normal vs. NorTA-CLR correlation) across input correlation
#   and zero-inflation levels (phi_min, phi_mean, phi_max). Also
#   generates illustrative scatterplot examples (Normal / NorTA /
#   NorTA-CLR) for four combinations of correlation sign and sparsity
#   level, using simulated pairs of OTUs embedded in a 200-OTU
#   background community, with per-OTU parameters drawn from empirical
#   ranges estimated on real HMP2 microbiome data.
#
# Input:
#   - Sparsity_Effects_htrlnorm.rds: pre-computed simulation results
#     with columns cor_input, cor_normal, cor_NorTA_PCLR, phi_1, phi_2
#     (located anywhere under the project root, found recursively)
#   - otu_HMP2.rds: real HMP2 OTU count table (samples x OTUs),
#     expected to be a single file found recursively under the
#     project root
#   - meta_HMP2.rds: real HMP2 sample metadata, used to subset to
#     healthy samples from subject 69-001, also expected to be a
#     single file found recursively
#
# Output:
#   - effects_htrlnorm.png: three-panel heatmap of mean absolute CLR
#     error vs. input correlation and sparsity (phi min / mean / max)
#   - examples_htrlnorm.png: four-example, three-panel scatterplot
#     figure (Normal / NorTA / NorTA-CLR) illustrating the effect of
#     correlation sign and sparsity level
#
# Note:
#   The two png() calls currently save to different folders
#   ("script/sparsity_effects/trials" for the heatmap, versus
#   "01_from_paper/sparsity_effects/trials" for the examples figure);
#   check this is intentional, or unify the output path.
#
# Packages used:
#   - ToyModel: provides toy_model(), qhtrlnorm(), mle.htrlnorm(), and
#     clr(), used to simulate and transform the htrlnorm-based data
#     (not on CRAN; assumed already installed/available in the project
#     environment)

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


# -------- data searching and filtering --------

# Search recursively in the project folder for the simulation results file
path <- list.files(
  path       = here(),
  pattern    = "^Sparsity_Effects_htrlnorm\\.rds$",
  full.names = TRUE,
  recursive  = TRUE
)

# Load the file and compute derived columns for error and sparsity summaries
df <- readRDS(path) %>%
  # Compute the absolute error between the normal and CLR-transformed correlation
  mutate("ERR_CLR" = abs(cor_normal - cor_NorTA_PCLR)) %>%
  # Compute the average sparsity across the two OTUs
  mutate("phi_mean" = 0.5 * (phi_1 + phi_2)) %>%
  # Compute the maximum sparsity across the two OTUs
  mutate("phi_max" = pmax(phi_1, phi_2)) %>%
  # Compute the minimum sparsity across the two OTUs
  mutate("phi_min" = pmin(phi_1, phi_2))

# Stop execution if any absolute error exceeds 1 (sanity check)
if (any(df$ERR_CLR > 1)) stop("Find Error greater than 1")

# Build a diverging color palette (Spectral, reversed) interpolated to 12 colors
myPalette <-
  RColorBrewer::brewer.pal(n = 11, "Spectral") %>%
  rev() %>%
  grDevices::colorRampPalette()

# -------- Heatmap of mean CLR error by correlation and minimum sparsity --------

# Average the CLR error over all replicates for each (cor_input, phi_min) combination
p.htrlnorm.min <- df %>%
  reframe(
    ERR_CLR_MEAN = mean(ERR_CLR),
    .by = c(cor_input, phi_min)
  ) %>%
  ggplot(aes(
    x = cor_input,
    y = phi_min,
    fill = ERR_CLR_MEAN
  )) +
  # Draw a heatmap tile for each (correlation, sparsity) cell
  geom_tile() +
  theme_bw() +

  # Apply the custom diverging color scale with fine-grained control near 0
  scale_fill_gradientn(
    name    = "Absolute Error",
    colours = myPalette(12),
    values  = c(seq(0, 0.5, by = 0.05), 1),
    limits  = c(0, 1)
  ) +

  # Remove tick marks from the color bar
  guides(fill = guide_colorbar(ticks.colour = NA)) +
  theme(legend.text = element_text(size = 10)) +

  # Set y-axis breaks at every 10% sparsity level
  scale_y_continuous(
    breaks = seq(0.1, 0.9, 0.1),
    expand = c(0, 0)
  ) +

  # Set x-axis breaks at every 0.2 correlation step
  scale_x_continuous(
    breaks = seq(-0.8, 0.8, 0.2),
    expand = c(0, 0)
  ) +

  # Rotate x-axis labels 90 degrees for readability
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  xlab("Correlation") +
  ylab("Zero %")

# -------- Heatmap of mean CLR error by correlation and mean sparsity --------

# Average the CLR error for each (cor_input, phi_mean) combination,
# then keep only phi_mean values that are multiples of 5% to reduce visual clutter
p.htrlnorm.mean <- df %>%
  reframe(
    ERR_CLR_MEAN = mean(ERR_CLR),
    .by = c(cor_input, phi_mean)
  ) %>%
  filter((100 * phi_mean) %% 5 == 0) %>%
  ggplot(aes(
    x = cor_input,
    # Treat phi_mean as a discrete axis since filtered values are evenly spaced
    y = factor(phi_mean),
    fill = ERR_CLR_MEAN
  )) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradientn(
    name    = "Absolute Error",
    colours = myPalette(12),
    values  = c(seq(0, 0.5, by = 0.05), 1),
    limits  = c(0, 1)
  ) +
  guides(fill = guide_colorbar(ticks.colour = NA)) +
  theme(legend.text = element_text(size = 10)) +
  scale_y_discrete(breaks = seq(0.1, 0.9, 0.1), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(-0.8, 0.8, 0.2), expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  xlab("Correlation") +
  ylab("Zero %")

# -------- Heatmap of mean CLR error by correlation and maximum sparsity --------

# Average the CLR error for each (cor_input, phi_max) combination
p.htrlnorm.max <- df %>%
  reframe(
    ERR_CLR_MEAN = mean(ERR_CLR),
    .by = c(cor_input, phi_max)
  ) %>%
  ggplot(aes(
    x = cor_input,
    y = phi_max,
    fill = ERR_CLR_MEAN
  )) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradientn(
    name    = "Absolute Error",
    colours = myPalette(12),
    values  = c(seq(0, 0.5, by = 0.05), 1),
    limits  = c(0, 1)
  ) +
  guides(fill = guide_colorbar(ticks.colour = NA)) +
  theme(legend.text = element_text(size = 10)) +
  scale_y_continuous(breaks = seq(0.1, 0.9, 0.1), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(-0.8, 0.8, 0.2), expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  xlab("Correlation") +
  ylab("Zero %")

# Open a PNG graphics device at high resolution for saving the combined figure
png(here("script", "sparsity_effects", "trials", "effects_htrlnorm.png"), width = 3600, height = 1500, res = 300)
# Arrange the three heatmaps side by side with a shared legend at the bottom
print(ggarrange(p.htrlnorm.min, p.htrlnorm.mean, p.htrlnorm.max,
  ncol          = 3,
  labels        = c("Phi Min", "Phi Mean", "Phi Max"),
  common.legend = T,
  legend        = "bottom"
))
# Close the graphics device and write the file to disk
dev.off()


# -------- Some Examples with htrlnorm --------

# Read HMP2
# this code allows R to search for the target file without specifying the folder names
f1 <- list.files(
  path       = here(),
  pattern    = "^otu_HMP2\\.rds$",
  full.names = TRUE,
  recursive  = TRUE
)
# blocks if the code finds more than one file with same name
stopifnot(length(f1) == 1)
otu <- readRDS(f1)

# this code allows R to search for the target file without specifying the folder names
f2 <- list.files(
  path       = here(),
  pattern    = "^meta_HMP2\\.rds$",
  full.names = TRUE,
  recursive  = TRUE
)
# blocks if the code finds more than one file with same name
stopifnot(length(f2) == 1)
meta <- readRDS(f2)

# filters the data
otu_69001_H <- otu[meta$SubjectID == "69-001" & meta$CL4_2 == "Healthy", ]
otu_filt <- otu_69001_H[, colSums(otu_69001_H > 0) / nrow(otu_69001_H) >= 0.33]
otu_filt <- otu_filt[, apply(otu_filt, 2, function(x) median(x[x > 0]) >= 5)]
# remove variables to free memory
rm(otu, otu_69001_H, meta)


# Quantiles
HMP2.quantile.params <- apply(otu_filt, 2, function(x) {
  ToyModel::mle.htrlnorm(x)$estimate
}) %>% apply(1, function(x) quantile(x, probs = c(0.1, 0.9)))
# Fitting a linear model between the mean and max in log scale for each OTU
df_mean_max <- tibble(
  "y" = apply(log(otu_filt + 1), 2, max),
  "x" = apply(log(otu_filt + 1), 2, mean)
)
model <- lm(y ~ x, data = df_mean_max)

# function to predict new maximum values from the means
predict_max <- function(new_xs) {
  # predict Mean Values and Standard Errors
  new_data <- data.frame(x = new_xs)
  predicted_values <- predict(model, new_data, interval = "none")

  # calculate standard error of prediction
  residuals_variance <- sum(residuals(model)^2) / model$df.residual
  # number of observations in original data
  n <- length(model$model$x)
  # mean of original independent variable
  x_bar <- mean(model$model$x)

  # leverage for each new_x (hii = 1/n + (xi - x̄)^2 / Σ(xi - x̄)^2)
  leverages <- 1 / n + ((new_data$x - x_bar)^2 / sum((model$model$x - x_bar)^2))

  # standard error of prediction for each new_x
  std_error_prediction <- sqrt(residuals_variance * (1 + leverages))


  # generate random values from normal distributions with these means and variances
  simulated_ys <- rnorm(nrow(new_data), mean = predicted_values, sd = std_error_prediction)
  return(simulated_ys)
}

example <- function(n, cor, meanlog_1, meanlog_2, sdlog_1, sdlog_2, phi_1, phi_2) {
  # generate random parameters for 200 background OTUs, sampled uniformly
  # within the 10th–90th percentile range estimated from HMP2 data
  params_random_HMP2 <- data.frame(
    "meanlog" = runif(
      n = 200,
      min = HMP2.quantile.params["10%", "meanlog"],
      max = HMP2.quantile.params["90%", "meanlog"]
    ),
    "sdlog" = runif(
      n = 200,
      min = HMP2.quantile.params["10%", "sdlog"],
      max = HMP2.quantile.params["90%", "sdlog"]
    ),
    "phi" = runif(
      n = 200,
      min = HMP2.quantile.params["10%", "phi"],
      max = HMP2.quantile.params["90%", "phi"]
    )
  )

  # predict the maximum log-value for each OTU based on its meanlog
  params_random_HMP2$b <- predict_max(params_random_HMP2$meanlog)

  # simulate 200 independent background OTUs using the htrlnorm distribution
  random_HMP2 <- ToyModel::toy_model(
    n = n,
    cor = diag(200),
    M = 1,
    qdist = ToyModel::qhtrlnorm,
    param = params_random_HMP2
  )

  # build the parameter set for the two OTUs of interest
  params_set <- data.frame(
    "phi" = c(phi_1, phi_2),
    "meanlog" = c(meanlog_1, meanlog_2),
    "sdlog" = c(sdlog_1, sdlog_2)
  ) %>%
    # predict the b parameter for each of the two OTUs
    mutate(b = predict_max(meanlog))

  # simulate the correlated pair of OTUs using the specified correlation
  couple <-
    ToyModel::toy_model(
      n = n,
      cor = cor,
      M = 1,
      qdist = ToyModel::qhtrlnorm,
      param = params_set
    )

  # replace the first two columns of the background matrix with the correlated pair
  random_HMP2_NorTA <- random_HMP2$NorTA
  random_HMP2_NorTA[, 1] <- couple$NorTA[, 1]
  random_HMP2_NorTA[, 2] <- couple$NorTA[, 2]

  # apply CLR transformation to the full 200-OTU compositional matrix
  NorTA_CLR <- ToyModel::clr(random_HMP2_NorTA)

  # retain only the two OTUs of interest after CLR transformation
  NorTA_CLR <- NorTA_CLR[, 1:2]

  # return the latent normal data, the NorTA counts, and the CLR-transformed counts
  return(list(
    "Normal" = couple$normal,
    "NorTA" = couple$NorTA,
    "NorTA_CLR" = NorTA_CLR
  ))
}

plot_example <- function(example) {
  # plot the latent normal variables with regression line and Pearson correlation
  p_norm <- example$Normal %>%
    as_tibble() %>%
    ggscatter(
      x = "V1",
      y = "V2",
      add = "reg.line",
      add.params = list(color = "blue", fill = "lightgray"),
      conf.int = TRUE
    ) +

    # annotate with the Pearson correlation coefficient
    stat_cor(aes(label = ..r.label..), color = "red") +

    ggtitle("Normal") +

    theme(axis.title = element_blank())

  # plot the NorTA-simulated counts with regression line and Pearson correlation
  p_NorTA <- example$NorTA %>%
    as_tibble() %>%
    ggscatter(
      x = "V1",
      y = "V2",
      add = "reg.line",
      add.params = list(color = "blue", fill = "lightgray"),
      conf.int = TRUE
    ) +

    stat_cor(aes(label = ..r.label..), color = "red") +

    ggtitle("NorTA") +

    theme(axis.title = element_blank())

  # plot the CLR-transformed NorTA counts with regression line and Pearson correlation
  p_NorTA_CLR <- example$NorTA_CLR %>%
    as_tibble() %>%
    ggscatter(
      x = "V1",
      y = "V2",
      add = "reg.line",
      add.params = list(color = "blue", fill = "lightgray"),
      conf.int = TRUE
    ) +

    stat_cor(aes(label = ..r.label..), color = "red") +

    ggtitle("NorTA_CLR") +

    theme(axis.title = element_blank())

  # return the three scatter plots as a named list
  return(list("norm" = p_norm, "NorTA" = p_NorTA, "NorTA_CLR" = p_NorTA_CLR))
}

# example 1: high sparsity, strong negative correlation
e1 <- example(
  n = 1000,
  cor = -0.9,
  meanlog_1 = 1,
  meanlog_2 = 3,
  sdlog_1 = 1,
  sdlog_2 = 1.3,
  phi_1 = 0.85,
  phi_2 = 0.9
) %>% plot_example()
# arrange the three panels side by side and add a descriptive title
pall1 <- ggarrange(e1$norm, e1$NorTA, e1$NorTA_CLR, ncol = 3)
pall1 <- annotate_figure(pall1, "n = 1000, cor = -0.9, meanlog = [1, 3], sdlog = [1, 1.3], phi = [0.85, 0.9]")

# example 2: high sparsity, strong positive correlation
e2 <- example(
  n = 1000,
  cor = 0.9,
  meanlog_1 = 4,
  meanlog_2 = 2,
  sdlog_1 = 1.5,
  sdlog_2 = 0.9,
  phi_1 = 0.9,
  phi_2 = 0.75
) %>% plot_example()
pall2 <- ggarrange(e2$norm, e2$NorTA, e2$NorTA_CLR, ncol = 3)
pall2 <- annotate_figure(pall2, "n = 1000, cor = -0.9, meanlog = [4, 2], sdlog = [1.5, 0.9], phi = [0.9, 0.75]")

# example 3: low sparsity, strong negative correlation
e3 <- example(
  n = 1000,
  cor = -0.9,
  meanlog_1 = 1.5,
  meanlog_2 = 5,
  sdlog_1 = 1,
  sdlog_2 = 1.5,
  phi_1 = 0.1,
  phi_2 = 0.05
) %>% plot_example()
pall3 <- ggarrange(e3$norm, e3$NorTA, e3$NorTA_CLR, ncol = 3)
pall3 <- annotate_figure(pall3, "n = 1000, cor = -0.9, meanlog = [1.5, 5], sdlog=[1, 1.5], phi = [0.1, 0.05]")

# example 4: low sparsity, strong positive correlation
e4 <- example(
  n = 1000,
  cor = 0.9,
  meanlog_1 = 2.5,
  meanlog_2 = 5.5,
  sdlog_1 = 1.7,
  sdlog_2 = 2,
  phi_1 = 0.1,
  phi_2 = 0.15
) %>% plot_example()
pall4 <- ggarrange(e4$norm, e4$NorTA, e4$NorTA_CLR, ncol = 3)
pall4 <- annotate_figure(pall4, "n = 1000, cor = 0.9, meanlog = [2.5, 5.5], sdlog = [1.7, 2], phi=[0.1,0.15]")

# combine all four examples into a single figure with 4 rows
pall <- plot_grid(pall1, pall2, pall3, pall4, nrow = 4)

# save the final figure as a high-resolution PNG
png(here("01_from_paper", "sparsity_effects", "trials", "examples_htrlnorm.png"), width = 3600, height = 4800, res = 300)
# this allows to generate the .png file
print(pall)
# close the graphics device and write the file to disk
dev.off()
