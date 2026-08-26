# Z_final_plot.R
#
# Purpose:
#   Assemble the final sparsity figure by combining two previously saved
#   plot objects into a single publication-ready panel.
#
# Input:
#   - A saved heatmap plot object:
#     here("script", "sparsity_effects", "error_zinb.rds")
#   - A saved real-data example plot object:
#     here("01_from_paper", "sparsity_effects", "example_couple.rds")
#
# Output:
#   - A high-resolution PNG figure:
#     here("Plots", "sparsity.png")
# Notes:
#   - The final figure places the simulation-based heatmap and the real-data
#     example side by side as panels A and B.
#   - The heatmap legend is moved to the top before composition.

# ggpubr: nice plots based on ggplot2
# https://cran.r-project.org/web/packages/ggpubr/index.html
library(ggpubr)

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# tidyverse: useful to manage data (dplyr) and make nice plots (ggplot2)
# https://cran.r-project.org/web/packages/tidyverse/index.html
library(tidyverse)


# Load a saved ggplot object for error visualization and move legend to top
p_err_zinb <- readRDS(list.files(
  path = here(),
  pattern = "^error_zinb\\.rds$",
  full.names = TRUE,
  recursive = TRUE
)) +
  theme(legend.position = "top")

# Load a saved ggplot object showing an example couple (pair of elements)
p_couple <- readRDS(list.files(
  path = here(),
  pattern = "^example_couple\\.rds$",
  full.names = TRUE,
  recursive = TRUE
))

# Combine the two plots side by side, sharing a single legend placed at the top
pall <- ggarrange(p_err_zinb, p_couple, common.legend = T, legend = "top", labels = c("A", ""))

# Open a PNG graphics device with high resolution (400 dpi, 4800x2400 px)
png(filename = here("Plots", "sparsity.png"), width = 4800, height = 2400, res = 400)

# UI reminder where to search the generated plot
cat("the .png file has been saved in the 'Plots' folder with the name 'sparsity.png' ")

# Render the combined plot to the PNG file
print(pall)
# Close the graphics device and save the file
dev.off()

# opens the output PNG in the default image viewer once the script completes
browseURL(list.files(
  path       = here(),
  pattern    = "^sparsity\\.png$",
  full.names = TRUE,
  recursive  = TRUE
))
