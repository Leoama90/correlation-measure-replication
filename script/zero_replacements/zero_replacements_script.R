# Collection of packages for data wrangling and visualization
# https://cran.r-project.org/web/packages/tidyverse/index.html
library(tidyverse)
# ggplot2 extensions for publication-ready graphics (scatter, stat_cor)
# https://cran.r-project.org/web/packages/ggpubr/index.html
library(ggpubr)
# Imputation of zeros, left-censored and missing values in compositional data
# https://cran.r-project.org/web/packages/zCompositions/index.html
library(zCompositions)
# Load custom function
source(here("script", "method_comparison", "CLR.R"))

# -------- read and filter data --------

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
# Read the taxonomy table
taxa <- readRDS(list.files(path       = here(),
                           pattern    = "taxonomy.rds",
                           full.names = TRUE,
                           recursive  = TRUE))

# Select samples belonging to 69-001 subject in health status
otu_69001_H <- otu[meta$SubjectID == "69-001" & meta$CL4_2 == "Healthy", ]

# Remove rarest OTUs using prevalence and median of non-zero values
otu_filt  <- otu_69001_H[, colSums(otu_69001_H > 0) / nrow(otu_69001_H) >= 0.33]
otu_filt  <- otu_filt[, apply(otu_filt, 2, function(x) median(x[x > 0]) >= 5)]
# Subset taxonomy to retain only the filtered OTUs
taxa_filt <- taxa[colnames(otu_filt), ]

# -------- different zero strategies --------

# Compute CLR-based correlation matrices under four zero-replacement strategies:
# CZM: Count Zero Multiplicative
cor_czm <- cor(CLR(as.matrix(zCompositions::cmultRepl(otu_filt, method = "CZM", label = 0)))) 
# GBM: Geometric Bayesian Multiplicative
cor_gbm <- cor(CLR(as.matrix(zCompositions::cmultRepl(otu_filt, method = "GBM", label = 0))))  
# BL: Beta-binomial Log-ratio
cor_bl  <- cor(CLR(as.matrix(zCompositions::cmultRepl(otu_filt, method = "BL",  label = 0)))) 
# Simple scalar replacement: zeros set to 0.65 (65% of detection limit)
otu_filt_65 <- otu_filt
otu_filt_65[otu_filt_65 == 0] <- 0.65
cor_65 <- cor(CLR(otu_filt_65))

# Extract upper-triangle indices to avoid duplicate pairs
otu_names <- colnames(otu_filt)
idx <- which(upper.tri(cor_czm), arr.ind = TRUE)

# Build a tidy tibble with one row per OTU pair and one column per method
cor_df <- tibble(
  var1 = otu_names[idx[, 1]],
  var2 = otu_names[idx[, 2]],
  CZM  = cor_czm[upper.tri(cor_czm)],
  GBM  = cor_gbm[upper.tri(cor_gbm)],
  BL   = cor_bl[upper.tri(cor_bl)],
  PC65 = cor_65[upper.tri(cor_65)]
)

# Reshape to long format for grouped summaries
cor_long <- cor_df %>%
  pivot_longer(cols      = c(CZM, GBM, BL, PC65),
               names_to  = "Method",
               values_to = "Correlation")

# For each OTU pair, compute the max absolute spread across methods
# and flag pairs where the spread exceeds 0.1
tbl <- cor_long %>%
  group_by(var1, var2) %>%
  summarise(max_abs_diff = max(Correlation) - min(Correlation),
            mean_corr    = mean(Correlation), .groups = "drop") %>%
  mutate(higher = ifelse(max_abs_diff > 0.1, TRUE, FALSE))

# -------- plot --------

# Save histogram of max absolute differences as PNG
png(width = 1600, height = 1200, res = 300,
    filename = here("Plots", "correlation_differences_between_zeroRepl.png"))

# create histogram plot
p <- tbl %>%
  ggplot(aes(x = max_abs_diff)) +
  # Draw histogram with 50 bins and green fill
  geom_histogram(bins = 50, fill = "lightgreen", color = "darkgreen") +
  # Apply a clean black-and-white theme
  theme_bw() +
  xlab("Maximum Absolute Difference in Correlation\nbetween Zero-Replacement Strategies\n(CZM, GBM, BL, 65% detection threshold)")

print(p)
# Close the PNG graphics device and write the file to disk
dev.off()