# clr_pearson.R
#
# Purpose:
#   Main pipeline function that ties together data summary, filtering,
#   pseudocount handling, CLR (Centered Log-Ratio) transformation, and
#   Pearson correlation estimation on a metagenomic OTU count table.
#
# Input:
#   - a metagenomic OTU count table (samples x OTUs)
#   - prevalence_threshold: a numeric value between 0 and 1, passed on to
#     filt_data(). If NULL (default), filt_data() asks for it interactively;
#     if provided, the interactive prompt is skipped.
#   - threshold_pct: a numeric value between 0 and 1, passed on to
#     pseudocount(). If NULL (default), pseudocount() asks for it
#     interactively; if provided, the interactive prompt is skipped.
#
# Output:
#   - the CLR-transformed OTU table
#   - the OTU x OTU Pearson correlation matrix
#
# Used scripts:
#   - datasummary.R
#   - filt_data.R
#   - pseudocount.R
#
# Note:
#   This script only defines clr_on_data(). To run the pipeline on real
#   data, see demo_clr_pearson.R. clr_on_data() can be run either
#   interactively (default, both thresholds requested via prompt) or
#   non-interactively (passing prevalence_threshold and/or threshold_pct
#   explicitly), which is useful for reproducible notebooks and
#   sensitivity analyses over a range of thresholds.

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)
# tidyverse: data science toolkit (dplyr, ggplot2, tidyr, etc.)
# https://tidyverse.org
library(tidyverse)


# -------- recall scripts for filtering, pseudocount adding and summarizing --------

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
    pattern = "^pseudocount\\.R$",
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
#' @param x a data.frame, tibble, or matrix with samples on rows and
#'   OTUs on columns (raw counts).
#' @param prevalence_threshold numeric or NULL. Forwarded to
#'   filt_data(): the minimum fraction of samples in which an OTU must
#'   be present (non-zero) to be kept. If NULL (default), filt_data()
#'   asks for it interactively; if provided, the interactive prompt is
#'   skipped, enabling non-interactive/automated use.
#' @param threshold_pct numeric or NULL. Forwarded to pseudocount(): the
#'   fraction of each row's detection limit used as the pseudocount for
#'   that row. If NULL (default), pseudocount() asks for it
#'   interactively; if provided, the interactive prompt is skipped,
#'   enabling non-interactive/automated use.
#'
#' @return an invisible list with samp_filt (the filtered OTU table),
#'   y_clr (the CLR-transformed OTU table), and cor_matrix (the OTU x
#'   OTU Pearson correlation matrix).
#'
#' @examples
#' \dontrun{
#' # interactive use: asks for both thresholds at the prompt
#' clr_on_data(otu_table)
#'
#' # non-interactive use: same pipeline, but repeatable and safe to run
#' # inside a knitted/rendered notebook
#' clr_on_data(otu_table, prevalence_threshold = 0.33, threshold_pct = 0.5)
#' }
#' @export


# -------- body of the function --------

clr_on_data <- function(x, prevalence_threshold = NULL, threshold_pct = NULL) {
  
  # -------- step 1: filter rare/low-quality OTUs --------
  
  # filt_data() prints a before/after summary; if prevalence_threshold is
  # NULL it asks the user for it interactively, otherwise it uses the
  # value passed in directly and skips the prompt
  filt_result <- filt_data(x, prevalence_threshold = prevalence_threshold)
  samp_filt <- filt_result$samp_filt
  
  # -------- step 2: replace zeroes with row-specific pseudocounts --------
  
  # pseudocount() converts counts to row-wise proportions; if
  # threshold_pct is NULL it asks the user for a percentage threshold of
  # the detection limit interactively, otherwise it uses the value passed
  # in directly and skips the prompt. Either way, the returned tibble no
  # longer contains any zeroes
  y_prop <- as.matrix(pseudocount(samp_filt, threshold_pct = threshold_pct))
  
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