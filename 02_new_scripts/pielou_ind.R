# pielou_ind.R
#
# Purpose:
#   Computes the within-dataset diversity P as defined in Fuschi et al.
#   (2025), Section 2.1: the mean Pielou evenness index across samples,
#   where each sample's index is the Shannon entropy of its relative
#   abundances, normalized by ln(D) (D = number of taxa/columns). 
#   This version accepts zero values (common in real OTU count tables), 
#   following the standard convention that 0*log(0) = 0 (the limit of x*log(x)
#   as x approaches 0 from above), so it can be applied directly to raw
#   or filtered real data, not only to artificially shifted, strictly positive
#   simulated data.
#
# Inputs:
#   - No external files required. Any non-negative numeric matrix/
#     data.frame of samples (rows) x taxa (columns) can be passed to
#     pielou_ind().
#
# Outputs:
#   - a single numeric value: the Pielou index averaged across samples
#     (or, if per_sample = TRUE, a numeric vector with one value per
#     sample)

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)


#' Compute the average Pielou evenness index (within-dataset diversity P)
#'
#' For each sample (row), computes the Shannon entropy of its relative
#' abundances, normalized by ln(D), then averages this value across all
#' samples (unless per_sample = TRUE, in which case the per-sample
#' values are returned instead). A value close to 1 indicates taxa are
#' evenly distributed within samples; a value close to 0 indicates
#' abundance concentrated in very few taxa. Rows summing to zero, and
#' datasets with a single column (D = 1, where ln(D) = 0), are rejected
#' with an explicit error rather than silently returning NaN/Inf.
#'
#' @param x a numeric matrix or data.frame of non-negative values, with
#'   samples on rows and taxa on columns. Must have at least 2 columns.
#' @param per_sample logical. If FALSE (default), returns the single
#'   value averaged across samples. If TRUE, returns the vector of
#'   per-sample Pielou indices instead, without averaging.
#'
#' @return numeric. The mean Pielou index across all samples (default),
#'   or a numeric vector with one Pielou index per sample (if
#'   per_sample = TRUE).
#'
#' @examples
#' # works on strictly positive data
#' pielou_ind(matrix(runif(50, 1, 10), nrow = 10, ncol = 5))
#'
#' # also works on data containing zeros, e.g. real OTU count tables
#' pielou_ind(matrix(c(5, 0, 3, 0, 2, 1, 4, 4, 0), nrow = 3, ncol = 3))
#'
#' # per-sample values instead of the dataset-wide average
#' pielou_ind(matrix(runif(50, 1, 10), nrow = 10, ncol = 5), per_sample = TRUE)
#'
#' @export


# -------- body of the function --------

pielou_ind <- function(x, per_sample = FALSE) {
  
  # coerce to a plain matrix, so row/column operations below are simple
  x <- as.matrix(x)
  
  # negative values are not valid abundances/counts; this function
  # only expects non-negative data (zeros are fine, negatives are not)
  if (any(x < 0)) {
    stop("pielou_ind() requires non-negative values")
  }
  
  # a single taxon gives ln(D) = 0, which would make the normalization
  # divide by zero; the Pielou index is undefined in this case
  D <- ncol(x)
  if (D < 2) {
    stop("pielou_ind() requires at least 2 columns (taxa); ",
         "got ", D, ", and ln(1) = 0 would make the index undefined")
  }
  
  # samples with a library size of 0 (all-zero row) have undefined
  # relative abundances (division by zero); reject them explicitly
  lib_size <- rowSums(x)
  if (any(lib_size == 0)) {
    stop("pielou_ind() found at least one row with a library size ",
         "of 0; remove empty samples before calling it")
  }
  
  # row-wise relative abundances (each row sums to 1)
  p <- x / lib_size
  
  # per-entry contribution to the Shannon entropy sum, with the
  # standard convention 0*log(0) = 0 (avoids NaN from log(0) on
  # zero-abundance entries, which are common in real OTU tables)
  p_logp <- ifelse(p == 0, 0, p * log(p))
  
  # per-sample Shannon entropy, normalized by ln(D) (Pielou evenness)
  H <- -rowSums(p_logp)
  pielou_per_sample <- H / log(D)
  
  # return either the per-sample vector or its average across samples,
  # depending on what the caller asked for
  if (per_sample) {
    return(pielou_per_sample)
  }
  mean(pielou_per_sample)
}