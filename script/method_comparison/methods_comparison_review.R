# methods_comparison_review.R
#
# Purpose:
#   Compare four correlation methods for microbial compositional data
#   (Pearson + CLR, Pearson + L1, SparCC, and Rho proportionality) at
#   two taxonomic resolutions — OTU level and phylum level — using real
#   HMP2 data, producing a publication-ready 6-panel scatter plot figure.
#
# Input:
#   - Preprocessed HMP2 objects stored as RDS files:
#     * otu_HMP2.rds
#     * meta_HMP2.rds
#     * taxonomy.rds
#   - The script subsets samples from subject 69-001 in healthy status.
#   - OTUs are filtered by prevalence and abundance to retain only
#     well-represented taxa before correlation estimation.
#   - Phylum-level abundances are derived by aggregating filtered OTU
#     counts within each phylum per sample.
#
# Methods compared:
#   - Pearson + CLR: Pearson correlation after Centered Log-Ratio
#     transformation (reference method)
#   - Pearson + L1:  Pearson correlation on L1-normalized (relative)
#     abundances
#   - SparCC:        log-ratio correlations robust to compositionality
#   - Rho:           symmetric proportionality measure (propr package)
#
# Outputs:
#   - A single PNG figure (Methods_comparison_review.png) with 6 scatter
#     plots arranged in a 2-row x 3-column grid:
#     * Row 1 (OTU level):    panels A (CLR vs L1), B (CLR vs Rho),
#                             C (CLR vs SparCC)
#     * Row 2 (Phylum level): panels D (CLR vs L1), E (CLR vs Rho),
#                             F (CLR vs SparCC)
#   - Each panel includes a regression line with confidence interval
#     and a Pearson R label anchored to the top-left corner.
#
# Converts base R plots and grobs into ggplot2 objects (as.grob)
# https://cran.r-project.org/web/packages/ggplotify/index.html
library(ggplotify)

# ggplot2 extensions for publication-ready graphics (scatter, stat_cor)
# https://cran.r-project.org/web/packages/ggpubr/index.html
library(ggpubr)

# Low-level graphics: viewport and grob management
# https://cran.r-project.org/web/packages/grid/index.html
library(grid)

# Arranges multiple plots in grids (grid.arrange, tableGrob)
# https://cran.r-project.org/web/packages/gridExtra/index.html
library(gridExtra)

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# Proportionality analysis for compositional data (rho, phi)
# https://cran.r-project.org/web/packages/propr/index.html
library(propr)

# Sparse microbial network estimation (SPIEC-EASI: MB and GLASSO)
# https://github.com/zdk123/SpiecEasi  ← Bioconductor/GitHub
library(SpiecEasi)

# Collection of packages for data wrangling and visualization
# https://cran.r-project.org/web/packages/tidyverse/index.html
library(tidyverse)


# -------- load other scripts needed for the test (custom functions) --------

# loads Centered Log-Ratio transform function
source(here("script", "method_comparison", "CLR.R"))
# loads upper-triangle matrix extraction function
source(here("script", "method_comparison", "TRIU.R"))
# loads signed network layout function
source(here("script", "method_comparison", "layout_signed.R"))


# -------- READ AND FILTER DATA --------

# Load OTU abundance matrix (samples x OTUs)
otu <- readRDS(list.files(
  path       = here(),
  pattern    = "otu_HMP2.rds",
  full.names = TRUE,
  recursive  = TRUE
))
# Load sample metadata
meta <- readRDS(list.files(
  path       = here(),
  pattern    = "meta_HMP2.rds",
  full.names = TRUE,
  recursive  = TRUE
))
# Load OTU taxonomic annotations
taxa <- readRDS(list.files(
  path       = here(),
  pattern    = "taxonomy.rds",
  full.names = TRUE,
  recursive  = TRUE
))

# keep only samples from subject 69-001 in healthy status
otu_69001_H <- otu[meta$SubjectID == "69-001" & meta$CL4_2 == "Healthy", ]

# remove OTUs present in fewer than 33% of samples
otu_filt  <- otu_69001_H[, colSums(otu_69001_H > 0) / nrow(otu_69001_H) >= .33]
# further remove OTUs with median non-zero abundance < 5
otu_filt  <- otu_filt[, apply(otu_filt, 2, function(x) median(x[x > 0]) >= 5)]
# subset taxonomy to keep only filtered OTUs
taxa_filt <- taxa[colnames(otu_filt), ]


# -------- CORRELATION METHODS (OTU LEVEL) --------

# Pearson on L1-normalized (relative) abundances
res_L1 <- cor(otu_filt / rowSums(otu_filt), method = "pearson")

# set seed for reproducibility before SparCC, which uses internal bootstrapping (uses internal bootstrap)
set.seed(42)
# SparCC correlation (designed for compositional count data)
sparcc_res <- sparcc(otu_filt)
# access the correlation matrix; $Cor field name verified against installed SpiecEasi version
res_cc <- sparcc_res$Cor
# restore OTU names to matrix
colnames(res_cc) <- rownames(res_cc) <- colnames(otu_filt)

# Rho proportionality (compositional association metric)
res_rho <- propr(counts = otu_filt, metric = "rho")@matrix

# Pearson on CLR-transformed abundances (used as reference method)
res_clr <- cor(CLR(otu_filt), method = "pearson")


# -------- SCATTER PLOTS (OTU LEVEL) --------

# Pearson + CLR vs Pearson+L1
p0 <- ggscatter(data.frame("PearsonCLR" = TRIU(res_clr),
                           "PearsonL1"  = TRIU(res_L1)),
                        x   = "PearsonCLR", 
                        y   = "PearsonL1",
                        add = "reg.line", 
                        conf.int   = TRUE,
                        add.params = list(color = "red", fill = "lightgray")) +
  # overlay Pearson R label, anchored to the top-left corner of the panel
  ggpubr::stat_cor(aes(label = after_stat(r.label)),
                   label.x = -Inf, 
                   label.y =  Inf,
                   hjust   = -0.1, 
                   vjust   =  1.5, 
                   size    =  6.0) +
  theme_bw() +
  xlab("Pearson + CLR") + 
  ylab("Pearson+L1" ) +
  theme(plot.title = element_text(hjust = 0.5))

# Pearson + CLR vs SparCC
p1 <- ggscatter(data.frame("PearsonCLR" = TRIU(res_clr),
                           "SparCC"     = TRIU(res_cc)),
                        x   = "PearsonCLR", 
                        y   = "SparCC",
                        add = "reg.line", 
                        conf.int   = TRUE,
                        add.params = list(color = "red", fill = "lightgray")) +
  # overlay Pearson R label, anchored to the top-left corner of the panel
  stat_cor(aes(label = after_stat(r.label)),
                   label.x  = -Inf, 
                   label.y  =  Inf,
                   hjust    = -0.1,
                   vjust    =  1.5, 
                   size     =  6.0) +
  theme_bw() +
  xlab("Pearson + CLR") +
  theme(plot.title = element_text(hjust = 0.5))

# Pearson + CLR vs Rho
p2 <- ggscatter(data.frame("PearsonCLR" = TRIU(res_clr),
                           "Rho"        = TRIU(res_rho)),
                        x   = "PearsonCLR", 
                        y   = "Rho",
                        add = "reg.line", 
                        conf.int = TRUE,
                        add.params = list(color = "red", fill = "lightgray")) +
  # overlay Pearson R label, anchored to the top-left corner of the panel
  stat_cor(aes(label = after_stat(r.label)),
             label.x = -Inf, 
             label.y =  Inf,
             hjust   = -0.1, 
             vjust   =  1.5, 
             size    =  6.0) +
  theme_bw() +
  xlab("Pearson + CLR") +
  theme(plot.title = element_text(hjust = 0.5))


# -------- AGGREGATE TO PHYLUM LEVEL --------

# aggregate filtered OTU counts to phylum level by joining taxonomy,
# summing abundances within each phylum per sample, and reshaping to a wide matrix
phy <- otu_filt %>%
  # convert matrix to tibble, keep row names as column
  as_tibble(rownames = "sample_id") %>%
  # reshape from wide to long format
  pivot_longer(-sample_id, names_to = "OTU", values_to = "abundance") %>%
  # attach taxonomy info to each OTU row
  left_join(as_tibble(taxa, rownames = "OTU"), by = "OTU") %>%
  # remove OTUs with no phylum annotation
  filter(!is.na(phylum)) %>%
  # group by sample and phylum
  group_by(sample_id, phylum) %>%
  # sum OTU abundances within each phylum per sample
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  # reshape back to wide format (one column per phylum)
  pivot_wider(names_from = phylum, values_from = abundance, values_fill = 0) %>%
  # restore sample IDs as row names
  column_to_rownames("sample_id") %>%
  # convert to matrix for downstream functions
  as.matrix()


# -------- CORRELATION METHODS (PHYLUM LEVEL) --------

# Pearson on L1-normalized phylum abundances
res_L1_phy <- cor(phy / rowSums(phy), method = "pearson")

# set seed for reproducibility before SparCC, which uses internal bootstrapping
set.seed(42)
# SparCC on phylum abundances
sparcc_phy_res <- sparcc(phy)
res_cc_phy     <- sparcc_phy_res$Cor
# restore phylum names to matrix
colnames(res_cc_phy) <- rownames(res_cc_phy) <- colnames(phy)

# Rho proportionality on phylum abundances
res_rho_phy <- propr(counts = phy, metric = "rho")@matrix

# Pearson + CLR on phylum abundances (reference method)
res_clr_phy <- cor(CLR(phy), method = "pearson")


# -------- SCATTER PLOTS (PHYLUM LEVEL) --------

# Pearson + CLR vs Pearson+L1
p0_phy <- ggscatter(data.frame("PearsonCLR" = TRIU(res_clr_phy),
                               "PearsonL1"  = TRIU(res_L1_phy)),
                            x   = "PearsonCLR", 
                            y   = "PearsonL1",
                            add = "reg.line", 
                            conf.int = TRUE,
                            add.params = list(color = "red", fill = "lightgray")) +
  # overlay Pearson R label, anchored to the top-left corner of the panel
  stat_cor(aes(label = after_stat(r.label)),
             label.x = - Inf, 
             label.y =   Inf,
             hjust   = - 0.1, 
             vjust   =   1.5, 
             size    =   6.0) +
  theme_bw() +
  xlab("Pearson + CLR") + 
  ylab("Pearson+L1")

# Pearson + CLR vs SparCC
p1_phy <- ggscatter(data.frame("PearsonCLR" = TRIU(res_clr_phy),
                               "SparCC"     = TRIU(res_cc_phy)),
                            x   = "PearsonCLR", 
                            y   = "SparCC",
                            add = "reg.line", 
                            conf.int = TRUE,
                            add.params = list(color = "red", fill = "lightgray")) +
  # overlay Pearson R label, anchored to the top-left corner of the panel
  stat_cor(aes(label = after_stat(r.label)),
                   label.x = - Inf, 
                   label.y =   Inf,
                   hjust   = - 0.1, 
                   vjust   =   1.5, 
                   size    =   6.0) +
  theme_bw() +
  xlab("Pearson + CLR")

# Pearson + CLR vs Rho
p2.phy <- ggscatter(data.frame("PearsonCLR" = TRIU(res_clr_phy),
                               "Rho"        = TRIU(res_rho_phy)),
                            x   = "PearsonCLR", 
                            y   = "Rho",
                            add = "reg.line", 
                            conf.int = TRUE,
                            add.params = list(color = "red", fill = "lightgray")) +
  # overlay Pearson R label, anchored to the top-left corner of the panel
  stat_cor(aes(label = after_stat(r.label)),
                   label.x = - Inf, 
                   label.y =   Inf,
                   hjust   = - 0.1, 
                   vjust   =   1.5, 
                   size    =   6.0 ) +
  theme_bw() +
  xlab("Pearson + CLR") +
  theme(plot.title = element_text(hjust = 0.5))


# -------- COMBINE AND EXPORT PLOTS --------

# blank plot used as row label "OTU level"
label_otu <- ggplot() +
  theme_void() +
  annotate("text", x = 0.5, y = 0.5, label = "OTU level",    angle = 90, size = 6)

# blank plot used as row label "Phylum level"
label_phylum <- ggplot() +
  theme_void() +
  annotate("text", x = 0.5, y = 0.5, label = "Phylum level", angle = 90, size = 6)

# assemble all 8 panels in a 2-row x 4-column grid
combined_plots <- ggarrange(
  # row 1: label + OTU-level plots (L1, Rho, SparCC)
  label_otu, p0, p2, p1,
  # row 2: label + phylum-level plots
  label_phylum, p0_phy, p2.phy, p1_phy,
  ncol = 4,
  nrow = 2,
  # label column narrower than plot columns
  widths = c(0.15, 1, 1, 1), legend = "none",
  # panel labels (A-F, skipping label columns)
  labels = c("", "A", "B", "C", "", "D", "E", "F")
)

# build output path relative to the project root and render the plot to the PNG device
out_path <- here("Plots", "Methods_comparison_review.png")
png(filename = out_path, width = 13.8 * 600, height = 9 * 600, res = 600)
print(combined_plots)
dev.off()

message("Plot saved to: ", out_path)