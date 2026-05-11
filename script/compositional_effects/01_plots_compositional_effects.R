# tidyverse: data science toolkit (dplyr, ggplot2, tidyr, etc.)
# [https://tidyverse.org](https://tidyverse.org)
library(tidyverse)

# ggpubr: ggplot2 publication-ready plots
# [https://cran.r-project.org/package=ggpubr](https://cran.r-project.org/package=ggpubr)
library(ggpubr)

# here: project-oriented file paths
# [https://cran.r-project.org/package=here](https://cran.r-project.org/package=here)
library(here)

# -- Create output directory if it doesn't exist --------
dir.create(here("Plots"), showWarnings = FALSE, recursive = TRUE)

# -- Data loading --------

# Load preprocessed compositional effects data
df <- readRDS(here("compositional_effects", "compositional_effects_02.rds"))

# -- Grid sampling: match Pielou values to target grid --------

# Iterate over all combinations of d and target Pielou evenness values
df_sort <- tibble()

for (di in seq(5, 200, by = 5)) {
  for (ei in seq(0.025, 0.975, by = 0.025)) {
    
    # Subset rows for current value of d
    sub_df <- df %>% filter(d == di)
    
    # Select the row whose Pielou index is closest to the target ei
    row_i <- sub_df[which.min(abs(sub_df$pielou - ei)), ]
    
    # Append the rounded target Pielou value as a new column
    row_i <- c(row_i, "pielou_round" = ei)
    
    df_sort <- bind_rows(df_sort, row_i)
  }
}

# Remove loop variables no longer needed
rm(di, ei, sub_df, row_i)

# -- Pipeline to elaborate data to be readable --------

df_sort <- df_sort %>%
  # Compute absolute deviation between rounded and actual Pielou values
  mutate(pielou_error = abs(pielou_round - pielou)) %>%
  # Flag rows where the deviation is small enough to be considered valid
  mutate(pielou_error_logical = pielou_error < 0.005, .after = pielou_error) %>%
  # Log-transform L1 error and floor at -2 (i.e., MAE < 0.01 treated as equal)
  mutate(log_err_l1 = log10(ERR_L1)) %>%
  mutate(log_err_l1 = if_else(log_err_l1 < -2, -2, log_err_l1)) %>%
  # Log-transform CLR error and apply the same floor
  mutate(log_err_clr = log10(ERR_CLR)) %>%
  mutate(log_err_clr = if_else(log_err_clr < -2, -2, log_err_clr))

# Subset rows that did not match closely to any target Pielou value (quality control)
df_control <- df_sort %>%
  filter(pielou_error_logical == FALSE)

# -- Colour palette --------

# Build an 11-colour diverging palette (Spectral reversed → cool-to-warm)
my_palette <- RColorBrewer::brewer.pal(11, "Spectral") %>%
  rev() %>%
  grDevices::colorRampPalette()

# -- Sanity check: verify grid dimensions --------

# Pivot to wide format and check that dimensions match the expected grid size
df_sort %>%
  dplyr::select(d, pielou_round, log_err_clr) %>%
  pivot_wider(names_from = d, values_from = log_err_clr) %>%
  column_to_rownames("pielou_round") %>%
  dim()

# -- Heatmap: L1 bias --------

p_l1 <- ggplot(df_sort, aes(x = d, y = pielou_round, fill = log_err_l1)) +
  geom_tile() +
  theme_bw() +
  # Map fill to log-scale MAE with readable break labels
  scale_fill_gradientn(
    "MAE",
    limits  = c(-2, 0),
    colours = my_palette(11),
    breaks  = c(-0.1, -1, -2),
    labels  = c(">1", "0.1", "< 0.01")
  ) +
  # Remove tick marks from the colour bar
  guides(fill = guide_colorbar(ticks.colour = NA)) +
  theme(legend.text = element_text(size = 10)) +
  scale_y_continuous(breaks = seq(0.05, 0.95, 0.05), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(10, 200, 10),      expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab(expression(bar(P))) +
  xlab("D") +
  ggtitle("L1 Bias")

# -- Heatmap: CLR bias --------

p_clr <- ggplot(df_sort, aes(x = d, y = pielou_round, fill = log_err_clr)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradientn(
    "MAE",
    limits  = c(-2, 0),
    colours = my_palette(11),
    breaks  = c(-0.1, -1, -2),
    labels  = c(">1", "0.1", "< 0.01")
  ) +
  guides(fill = guide_colorbar(ticks.colour = NA)) +
  theme(legend.text = element_text(size = 10)) +
  scale_y_continuous(breaks = seq(0.05, 0.95, 0.05), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(10, 200, 10),      expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab(expression(bar(P))) +
  xlab("D") +
  ggtitle("CLR Bias")

# -- Export: side-by-side heatmaps --------

png(filename = here("Plots", "Normalization_Bias.png"), width = 6000, height = 3000, res = 600)
print(ggarrange(                                         
  p_l1, p_clr,
  labels        = c("L1", "CLR"),
  common.legend = TRUE,
  label.y       = 0.125,
  label.x       = c(0.03, -0.02)
))
dev.off()

# -- Line plot: CLR MAE vs dimensionality (mean only) --------

p_clr_dim <- df_sort %>%
  # Compute mean CLR error for each value of d
  reframe(err_clr_mean_d = mean(ERR_CLR), .by = d) %>%
  ggplot(aes(x = d, y = err_clr_mean_d)) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_y_continuous(breaks = c(0.01, 0.02, 0.05, 0.1, 0.15, 0.2)) +
  ylab("MAE") +
  xlab("D") +
  theme(plot.margin = unit(c(2, 1, 1, 1), "cm"))

png(filename = here("Plots", "CLR_Compositional.png"), width = 1200, height = 1200, res = 300)
print(p_clr_dim)                                         
dev.off()

# -- Export: combined multi-panel figure --------

# Stack the two heatmaps vertically with a shared legend
p_l1_clr <- ggarrange(p_l1, p_clr, common.legend = TRUE, ncol = 1)

# Combine heatmap panel and line plot side by side
p_all <- ggarrange(
  p_l1_clr,
  p_clr_dim + theme(
    axis.text  = element_text(size = 14),
    axis.title = element_text(size = 16)
  ),
  ncol    = 2,
  widths  = c(0.35, 0.65),
  labels  = c("A", "B"),
  label.x = c(0.05, 0.90)
)

png(filename = here("Plots", "Normalization_Bias_all.png"), width = 6000, height = 4500, res = 600)
print(p_all)                                             
dev.off()

# -- Line plot: CLR MAE with 10th-90th percentile error bars --------

# Summarise mean and decile bounds of CLR error for each d
df_percentiles <- df_sort %>%
  group_by(d) %>%
  summarise(
    err_clr_mean = mean(ERR_CLR),
    err_clr_p10  = quantile(ERR_CLR, 0.10),
    err_clr_p90  = quantile(ERR_CLR, 0.90)
  )

# Plot mean CLR error with inter-decile range as error bars
p_clr_dim <- ggplot(df_percentiles, aes(x = d, y = err_clr_mean)) +
  geom_errorbar(
    aes(ymin = err_clr_p10, ymax = err_clr_p90),
    width = 2,
    color = "red"
  ) +
  geom_point(size = 0.5, color = "black") +
  theme_bw() +
  scale_y_continuous(breaks = c(0.01, 0.02, 0.05, 0.1, 0.15, 0.2)) +
  ylab("MAE") +
  xlab("D") +
  theme(plot.margin = unit(c(2, 1, 1, 1), "cm"))

png(filename = here("Plots", "CLR_Compositional_percentiles.png"), width = 1200, height = 1200, res = 300)
print(p_clr_dim)                                         
dev.off()