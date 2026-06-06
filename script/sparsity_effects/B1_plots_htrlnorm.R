# ggpubr: nice plots based on ggplot2
# https://cran.r-project.org/web/packages/ggpubr/index.html
library(ggpubr)

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# tidyverse: useful to manage data (dplyr) and make nice plots (ggplot2)
# https://cran.r-project.org/web/packages/tidyverse/index.html
library(tidyverse)

# Load the dataset, drop phi_2, rename phi_1 to phi, and compute CLR absolute error
df <- readRDS(here("script", "sparsity_effects", "Sparsity_Effects_htrlnorm_2.rds")) %>%
  dplyr::select(-phi_2) %>% rename(phi = phi_1) %>%
  mutate("ERR_CLR" = abs(cor_normal - cor_NorTA_PCLR))

# Stop execution if any error value exceeds 1
if(any(df$ERR_CLR > 1)) stop("Find Error greater than 1")

# Build a 12-color diverging palette (Spectral reversed + black), then make it interpolatable
myPalette <- 
  c(RColorBrewer::brewer.pal(n = 11, "Spectral")) %>% rev() %>% c(., "#000000") %>%
  grDevices::colorRampPalette()

p <- df %>%
  # Compute mean CLR error grouped by input correlation and sparsity level
  summarise(MEAN_ERR_CLR = mean(ERR_CLR), .by = c(cor_input, phi)) %>%
  # Map cor_input to x-axis, phi to y-axis, and mean error to fill color
  ggplot(aes(x = cor_input, y = phi, fill = MEAN_ERR_CLR)) +
  # Draw a heatmap tile for each (cor_input, phi) combination
  geom_tile() + theme_bw() +
  # Apply the custom palette with non-linear breakpoints and fixed [0,1] range
  scale_fill_gradientn("MAE", colours = myPalette(12),
                       values = c(seq(0,.5,by = .05),1),
                       limits = c(0,1),
                       labels = c("0","0.25","0.5","0.75","1")) +
  # Remove tick marks from the colorbar
  guides(fill = guide_colorbar(ticks.colour = NA)) +
  # Set legend text size
  theme(legend.text = element_text(size = 10)) +
  # Set y-axis breaks every 0.1 with no padding
  scale_y_continuous(breaks = seq(0,.9,.1), expand = c(0,0))  +
  # Set x-axis breaks every 0.2 with no padding
  scale_x_continuous(breaks = seq(-.8,.8,.2), expand = c(0,0)) +
  # Rotate x-axis labels 90° for readability
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # Label the axes
  xlab("Correlation") + ylab("Zero %") 

# Save the plot object to an RDS file for later use
saveRDS(p, here("script", "sparsity_effects", "error_htlrlnorm.rds"))