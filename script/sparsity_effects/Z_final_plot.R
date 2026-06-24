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
p_err_zinb <- readRDS(list.files(path       = here(),
                                 pattern    = "error_zinb.rds",
                                 full.names = TRUE,
                                 recursive  = TRUE
                                 )) +  
  theme(legend.position = "top")

# Load a saved ggplot object showing an example couple (pair of elements)
p_couple   <- readRDS(list.files(path       = here(),
                                 pattern    = "example_couple.rds",
                                 full.names = TRUE,
                                 recursive  = TRUE
                                 ))

# Combine the two plots side by side, sharing a single legend placed at the top
pall <- ggarrange(p_err_zinb, p_couple, common.legend = T, legend = "top", labels = c("A", ""))

# Open a PNG graphics device with high resolution (400 dpi, 4800x2400 px)
png(filename = here("Plots", "sparsity.png"), width = 4800, height = 2400, res = 400)

# Render the combined plot to the PNG file
print(pall)
# Close the graphics device and save the file
dev.off()