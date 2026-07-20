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
    path = here(),
    pattern = "^datasummary\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)

# bring filt_data() into scope
source(
  list.files(
    path = here(),
    pattern = "^filt_data\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)

# bring pseudocount() into scope
source(
  list.files(
    path = here(),
    pattern = "^pseudocounts\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)

#' Run the CLR + Pearson pipeline on a metagenomic OTU table
#'
#' Filters the input OTU table, replaces zeroes with row-specific
#' pseudocounts, applies the Centered Log-Ratio (CLR) transformation,
#' and computes the OTU x OTU Pearson correlation matrix on the
#' transformed data.
#'
#' @param y a data.frame, tibble, or matrix with samples on rows and
#'   OTUs on columns (raw counts).
#'
#' @return an invisible list with samp_filt (the filtered OTU table),
#'   y_clr (the CLR-transformed OTU table), and cor_matrix (the OTU x
#'   OTU Pearson correlation matrix).
#'
#' @export

clr_on_data <- function(y) {
  # -------- step 1: filter rare/low-quality OTUs --------
  # filt_data() prints a before/after summary and asks the user for a
  # prevalence threshold interactively
  filt_result <- filt_data(y)
  samp_filt <- filt_result$samp_filt

  # -------- step 2: replace zeroes with row-specific pseudocounts --------
  # pseudocount() converts counts to row-wise proportions and asks the
  # user for a percentage threshold of the detection limit interactively;
  # the returned tibble no longer contains any zeroes
  y_prop <- as.matrix(pseudocount(samp_filt))

  # -------- step 3: CLR transformation --------
  # log of every value (safe: y_prop has no zeroes left, thanks to
  # pseudocount()), then subtract the row-wise mean of the logs, which
  # is equivalent to dividing by the geometric mean of each sample
  log_data <- log(y_prop)
  log_geo_mean <- rowMeans(log_data)
  y_clr <- log_data - log_geo_mean

  # -------- step 4: Pearson correlation on CLR-transformed data --------
  # cor() computes correlations between columns (OTUs) by default,
  # yielding the OTU x OTU correlation matrix
  cor_matrix <- cor(y_clr, method = "pearson")

  # -------- return the invisible values --------
  invisible(list(
    samp_filt = samp_filt,
    y_clr = y_clr,
    cor_matrix = cor_matrix
  ))
}

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
