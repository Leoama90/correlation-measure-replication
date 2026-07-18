# clr_pearson.R
#
# Purpose:
#   Main pipeline script that ties together data summary, filtering,
#   pseudocount handling, CLR (Centered Log-Ratio) transformation, and
#   Pearson correlation estimation on a metagenomic OTU count table.
#
# Input:
#   - a metagenomic OTU count table (samples x OTUs), read and prepared
#     via the sourced scripts below
#
# Output:
#   - the CLR-transformed OTU table
#   - the OTU x OTU Pearson correlation matrix

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)
# tidyverse: data science toolkit (dplyr, ggplot2, tidyr, etc.)
# [https://tidyverse.org](https://tidyverse.org)
library(tidyverse)

# -------- recall scripts for filtering, pseudocount add and summarizing --------

# bring datasum() into scope
source(
  list.files(
    path       = here(),
    pattern    = "^datasummary\\.R$",
    full.names = TRUE,
    recursive  = TRUE
  )
)

# bring filt_data() into scope
source(
  list.files(
    path       = here(),
    pattern    = "^filtering\\.R$",
    full.names = TRUE,
    recursive  = TRUE
  )
)

# bring pseudocount() into scope
source(
  list.files(
    path       = here(),
    pattern    = "^pseudocounts\\.R$",
    full.names = TRUE,
    recursive  = TRUE
  )
)

