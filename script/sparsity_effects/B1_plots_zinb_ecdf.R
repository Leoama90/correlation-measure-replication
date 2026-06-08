# ggpubr: nice plots based on ggplot2
# https://cran.r-project.org/web/packages/ggpubr/index.html
library(ggpubr)

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# tidyverse: useful to manage data (dplyr) and make nice plots (ggplot2)
# https://cran.r-project.org/web/packages/tidyverse/index.html
library(tidyverse)


# Load the dataset, drop the second zero-inflation parameter column,
# rename the first one, and compute the absolute correlation error for CLR
df <- readRDS(here("script", "sparsity_effects", "Sparsity_Effects_zinbin_rare_ecdf_2.rds")) %>%
  dplyr::select(-pstr0_2) %>% rename(pstr0 = pstr0_1) %>%
  mutate("ERR_CLR" = abs(cor_normal - cor_NorTA_PCLR))

# Safety check: stop execution if any CLR error exceeds 1 (would indicate a data issue)
if(any(df$ERR_CLR > 1)) stop("Find Error greater than 1")

# Build a 12-color diverging palette (Spectral, reversed) with black appended,
# then wrap it in colorRampPalette for smooth interpolation
myPalette <- 
  c(RColorBrewer::brewer.pal(n = 11, "Spectral")) %>% rev() %>% c(., "#000000") %>%
  grDevices::colorRampPalette()

# Aggregate: compute mean CLR error for each (cor_input, pstr0) combination,
# then build a heatmap tile plot with a custom gradient fill
p <- df %>%
  reframe(MEAN_ERR_CLR = mean(ERR_CLR), .by = c(cor_input, pstr0)) %>%
  ggplot(aes(x = cor_input, y = pstr0, fill = MEAN_ERR_CLR)) +
  geom_tile() + theme_bw() +
  # Apply the custom palette with non-linear breakpoints to emphasize low-error regions
  scale_fill_gradientn("MAE", colours = myPalette(12),
                       values = c(seq(0,.5, by = .05), 1),
                       limits = c(0, 1),
                       labels = c("0", "0.25", "0.5", "0.75", "1")) +
  # Remove colorbar tick marks for a cleaner legend appearance
  guides(fill = guide_colorbar(ticks.colour = NA)) +
  theme(legend.text = element_text(size     = 10)) +
  # Set y-axis (phi) breaks every 0.1, no padding
  scale_y_continuous(breaks = seq(0, .9, .1), expand   = c(0,0))  +
  # Set x-axis (r) breaks every 0.2, no padding
  scale_x_continuous(breaks = seq(-.8, .8, .2), expand = c(0, 0)) +
  # Rotate x-axis labels 90° for readability
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  # Label axes: r for correlation input, phi (Greek letter) for zero-inflation
  xlab("r") + ylab(expression(phi))  +
  theme(axis.title = element_text(size = 12)) +
  # Place the colorbar legend at the top of the plot
  theme(legend.position = "top")

# Save the ggplot object to an RDS file for later reuse
saveRDS(p, here("script", "sparsity_effects", "error_zinb.rds"))