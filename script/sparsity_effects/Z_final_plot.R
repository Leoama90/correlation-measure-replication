# tidyverse: useful to manage data (dplyr) and make nice plots (ggplot2)
# https://cran.r-project.org/web/packages/tidyverse/index.html
library(tidyverse)

# igraph: used to create and analyze network graphs
# https://cran.r-project.org/web/packages/igraph/index.html
library(igraph)

# ggpubr: nice plots based on ggplot2
# https://cran.r-project.org/web/packages/ggpubr/index.html
library(ggpubr)

# magick: advanced image processing and manipulation (read, write, transform images)
# https://cran.r-project.org/web/packages/magick/index.html
library(magick)

# Set working directory to the folder containing the current R script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load a saved ggplot object for error visualization and move legend to top
p_err <- readRDS("error.rds") + theme(legend.position="top")
# Load a saved ggplot object showing an example couple (pair of elements)
p_couple <- readRDS("example_couple.rds")

# Combine the two plots side by side, sharing a single legend placed at the top
pall <- ggarrange(p_err, p_couple, common.legend=T, legend="top", labels = c("A", ""))

# Open a PNG graphics device with high resolution (400 dpi, 4800x2400 px)
png(filename = here("Plots", "sparsity.png"), width=4800, height=2400, res=400)
# Render the combined plot to the PNG file
pall
# Close the graphics device and save the file
dev.off()