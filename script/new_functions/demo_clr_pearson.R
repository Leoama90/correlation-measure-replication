# demo_clr_pearson.R
#
# Purpose:
#   Demonstrates the usage of clr_on_data(), loading a real metagenomic
#   OTU count table and running the full filtering + pseudocount + CLR +
#   Pearson correlation pipeline on it.
#
# Inputs:
#   - clr_pearson.R (sourced below)
#   - otu_HMP2.rds: a metagenomic OTU count table (samples x OTUs)
#
# Outputs:
#   - the CLR-transformed OTU table (result$y_clr)
#   - the OTU x OTU Pearson correlation matrix (result$cor_matrix)

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)


# -------- load the clr_pearson.R script --------

# load clr_on_data() and its dependencies (datasum, filt_data, pseudocount)
source(
  list.files(
    path = here(),
    pattern = "^clr_pearson\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)


# -------- read raw data and run the pipeline --------

# Load the OTU abundance matrix (samples x OTUs) from its .rds file
otu <- readRDS(
  list.files(
    path = here(),
    pattern = "^otu_HMP2\\.rds$",
    full.names = TRUE,
    recursive = TRUE
  )
)
# run the full pipeline: filtering, pseudocount, CLR, Pearson correlation
result <- clr_on_data(otu)