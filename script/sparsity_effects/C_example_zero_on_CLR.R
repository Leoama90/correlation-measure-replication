# C_example_zero_on_CLR.R
#
# Purpose:
#   Illustrate how sparsity distorts CLR correlation estimates on real data
#   by visualizing two hand-picked OTU pairs with high zero rates and
#   extreme correlations.
#
# Input:
#   - HMP2 OTU abundance table
#   - HMP2 sample metadata
#   - The script keeps OTUs present in at least 25% of samples
#   - CLR transformation is applied to the filtered OTU table
#   - Pairwise Pearson correlations are computed on CLR-transformed data
#
#   The script then builds a long-format table of OTU-pair correlations and
#   zero rates, and selects two focal OTU pairs for visualization.
#
# Output:
#   - A 2x2 panel figure saved as an RDS object:
#     here("script", "sparsity_effects", "example_couple.rds")
#
# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# ggpubr: nice plots based on ggplot2
# https://cran.r-project.org/web/packages/ggpubr/index.html
library(ggpubr)

# tidyverse: useful to manage data (dplyr) and make nice plots (ggplot2)
# https://cran.r-project.org/web/packages/tidyverse/index.html
library(tidyverse)

# load CLR script
source(here("script", "method_comparison", "CLR.R"))

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

# Keep only OTUs present in at least 25% of samples (prevalence filter)
otu.filt <- otu[, colSums(otu > 0)/nrow(otu) >= 0.25]
# Further keep only OTUs whose median count (among non-zero samples) is >= 5 (abundance filter)
otu.filt <- otu.filt[, apply(otu.filt, 2, \(x)median(x[x > 0]) >= 5)]

# Apply CLR transformation to the filtered OTU table
otu.filt.CLR <- CLR(otu.filt)
# Compute the Pearson correlation matrix on CLR-transformed data (PCLR = Pearson CLR)
PCLR     <- cor(otu.filt.CLR)

# Compute prevalence (fraction of samples where count > 0) for each OTU
otu_prev <- setNames(colSums(otu.filt > 0)/nrow(otu.filt), colnames(otu.filt))

# Build a long-format data frame of all pairwise OTU correlations with their prevalences
info     <- data.frame(PCLR) %>% rownames_to_column("OTU_I") %>%
  # Pivot to long format: one row per OTU pair with their correlation value
  pivot_longer(!OTU_I, 
               names_to  = "OTU_J", 
               values_to = "cor") %>%
  # Keep only the upper triangle of the correlation matrix (avoid duplicates and self-pairs)
  filter(OTU_I > OTU_J) %>%
  # Attach prevalence values for each OTU in the pair
  mutate(prev_I = otu_prev[OTU_I],
         prev_J = otu_prev[OTU_J]) %>%
  # Convert prevalences to zero-rates (% of samples where the OTU is absent)
  mutate(zero_I = 100 * round(1 - prev_I, 2),
         zero_J = 100 * round(1 - prev_J, 2))

# Subset pairs where both OTUs are low-prevalence (<=50%) and strongly correlated (|r|>=0.4)
info_filt <- info %>%
  filter(prev_I   <= 0.5, 
         prev_J   <= 0.5, 
         abs(cor) >= 0.4)

# Define the OTU pair with the most negative (min) CLR correlation to visualize
idx.max   <- c("OTU_269", "OTU_97")
# Define the OTU pair with the most positive (max) CLR correlation to visualize
idx.min   <- c("OTU_269", "OTU_313")

# Classify each sample by detection pattern for the 'min' pair (OTU_269 vs OTU_313)
detection <- otu.filt %>% as_tibble() %>%
  select(tidyselect::all_of(idx.min)) %>%
  # Label each sample: which OTUs are zero vs. both detected
  mutate(detection = case_when(
    OTU_269 == 0 & OTU_313 == 0~"Both are 0",
    OTU_269 == 0~"OTU_269 is 0",
    OTU_313 == 0~"OTU_313 is 0",
    OTU_269 > 0 & OTU_313 > 0~"Both are > 0",
  )) %>%
  dplyr::select(detection)

# Define a named color palette for the four detection categories
palette_named <- c(
  "Both are 0"   = "#D73027",  # red
  "OTU_269 is 0" = "#1A9850",  # green
  "OTU_313 is 0" = "#4575B4",  # blue
  "Both are > 0" = "black"
)

# Apply 80% opacity to all colors
palette_named <- alpha(palette_named, 0.8)

# Scatter plot of raw counts for the min-correlation pair, colored by detection pattern
p.min <- otu.filt %>% as_tibble() %>%
  select(tidyselect::all_of(idx.min)) %>%
  cbind(detection) %>%
  ggscatter(x = "OTU_269", 
            y = "OTU_313", 
            add   = "reg.line",       # add linear regression line
            color = "detection",
            size  = 2, 
            palette  = palette_named,
            add.params = list(color = "#800080", fill = "lightgray"),
            
            conf.int = TRUE) +        
  # show confidence interval around regression
  # Annotate with Pearson r value
  stat_cor(aes(label = after_stat(r.label)), color = rgb(0.5, 0, 0), label.x.npc = 0) +
  xlab(expression("Count OTU 269 ("*phi*"~73%)")) +
  ylab(expression("Count OTU 313 ("*phi*"~55%)")) +
  theme(axis.title   = element_text(size = 12)) +
  theme(legend.title = element_blank())

# Same scatter plot but using CLR-transformed values instead of raw counts
p.min.clr <- otu.filt.CLR %>% as_tibble %>%
  dplyr::select(tidyselect::all_of(idx.min)) %>%
  cbind(detection) %>%
  ggscatter(x = "OTU_269", 
            y = "OTU_313", 
            add   = "reg.line", 
            color = "detection",
            size  = 2, 
            palette  = palette_named,
            add.params = list(color = "#800080", fill = "lightgray"),
            conf.int = TRUE) +
  
  stat_cor(aes(label = after_stat(r.label)), color = rgb(0.5, 0, 0), label.x.npc = 0) +
  xlab(expression("CLR OTU 269 ("*phi*"~73%)")) +
  ylab(expression("CLR OTU 313 ("*phi*"~55%)")) +
  theme(axis.title   = element_text(size = 12)) +
  theme(legend.title = element_blank())

# Classify each sample by detection pattern for the 'max' pair (OTU_269 vs OTU_97)
detection <- otu.filt %>% as_tibble() %>%
  select(tidyselect::all_of(idx.max)) %>%
  mutate(detection = case_when(
    OTU_269 == 0 & OTU_97 == 0~"Both are 0",
    OTU_269 == 0~"OTU_269 is 0",
    OTU_97  == 0~"OTU_97 is 0",
    OTU_269 > 0 & OTU_97 > 0~"Both are > 0", 
  )) %>%
  dplyr::select(detection)

# New named palette for the max pair (same structure, different OTU label)
palette_named <- c(
  "Both are 0"   = "#D73027",
  "OTU_269 is 0" = "#1A9850",
  "OTU_97 is 0"  = "#4575B4",
  "Both are > 0"  = "black"
)
palette_named <- alpha(palette_named, 0.8)

# Scatter plot of raw counts for the max-correlation pair
p.max <- otu.filt %>% as_tibble() %>%
  select(tidyselect::all_of(idx.max)) %>%
  cbind(detection) %>%
  ggscatter(x = "OTU_269", 
            y = "OTU_97", 
            add   = "reg.line", 
            color = "detection",
            size  = 1.5, 
            palette  = palette_named,
            add.params = list(color = "#800080", fill = "lightgray"),
            conf.int = TRUE) +
  stat_cor(aes(label = after_stat(r.label)), color = rgb(0.5, 0, 0), label.x.npc = 0) +
  xlab(expression("Count_OTU_269 ("*phi*"~73%)")) +
  ylab(expression("Count_OTU_97  ("*phi*"~56%)")) +
  theme(axis.title   = element_text(size = 12)) +
  theme(legend.title = element_blank())

# Same scatter plot but using CLR-transformed values for the max pair
p.max.clr <- otu.filt.CLR %>% as_tibble %>%
  dplyr::select(tidyselect::all_of(idx.max)) %>%
  cbind(detection) %>%
  ggscatter(x = "OTU_269", 
            y = "OTU_97", 
            add   = "reg.line",
            color = "detection",
            size  = 1.5, 
            palette  = palette_named,
            add.params = list(color = "#800080", fill = "lightgray"),
            
            conf.int = TRUE) +
  
  stat_cor(aes(label = after_stat(r.label)), color = rgb(0.5, 0, 0), label.x.npc = 0) +
  xlab(expression("CLR_OTU_269 ("*phi*"~73%)")) +
  ylab(expression("CLR_OTU_97  ("*phi*"~56%)")) +
  theme(axis.title   = element_text(size = 12)) +
  theme(legend.title = element_blank())

# Arrange all four scatter plots into a 2x2 panel figure (panels B, C, D, E)
pall <- ggarrange(
  # Top row: raw counts vs CLR for the min-correlation pair
  ggarrange(p.min, p.min.clr, common.legend = T, legend = "top", ncol = 2, labels = c("B", "C"), 
            label.x = 0.06, label.y = 1.07), 
  # Bottom row: raw counts vs CLR for the max-correlation pair
  ggarrange(p.max, p.max.clr, common.legend = T, legend = "top", ncol = 2, labels = c("D", "E"),
            label.x = 0.06, label.y = 1.07),
  nrow = 2)

# Save the combined figure as an RDS object for later use
saveRDS(pall,here("script", "sparsity_effects", "example_couple.rds"))

# UI reminder where to search the generated plot
cat("the .rds file has been saved in the 'sparsity_effects' folder with the name 'example_couple.rds' ")