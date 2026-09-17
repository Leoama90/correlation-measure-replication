# pseudocount.R
#
# Purpose:
#   Convert count data to row-wise proportions and replace zeroes with
#   row-specific pseudocounts computed from a user-defined threshold.
#
# Input:
#   - x: a numeric matrix or data frame of non-negative counts.
#   - threshold_pct: a numeric value between 0 and 1. If NULL (default),
#     requested interactively; if provided, the interactive prompt is
#     skipped, enabling non-interactive/automated use (e.g. reproducible
#     notebooks, sensitivity analyses over a range of thresholds).
#
# Output:
#   - A tibble with counts converted to proportions and zeroes replaced by
#     row-specific pseudocounts.
#
# Description:
#   The function computes library size by row, obtains a threshold (either
#   passed in directly or requested interactively), calculates row-specific
#   pseudocounts, and substitutes zeroes in the proportional table with
#   those values.

# tidyverse: data science toolkit (dplyr, ggplot2, tidyr, etc.)
# https://tidyverse.org
library(tidyverse)

#' Replace zeroes with row-wise pseudocounts
#'
#' This function computes row-wise library sizes, obtains a threshold
#' between 0 and 1 (either passed in directly or requested interactively),
#' converts counts to proportions, and replaces zeroes with row-specific
#' pseudocount values.
#'
#' @param x A numeric matrix or data frame of non-negative counts. Rows
#'   (samples) must have a library size greater than 0.
#' @param threshold_pct numeric or NULL. Fraction (0 to 1) of the row's
#'   detection limit (1 / library size) used as the pseudocount for that
#'   row. If NULL (default), the threshold is requested interactively from
#'   the user; if provided, the interactive prompt is skipped, enabling
#'   non-interactive/automated use (e.g. reproducible notebooks).
#'
#' @return A tibble with counts converted to proportions and zeroes replaced by
#'   row-specific pseudocounts.
#'
#' @examples
#' \dontrun{
#' x <- matrix(c(10, 0, 5, 2, 3, 0), nrow = 2, byrow = TRUE)
#'
#' # interactive use: asks the user for a threshold at the prompt
#' pseudocount(x)
#' # if the user enters 0.5 at the prompt: each zero is replaced with
#' # half of that row's detection limit (1 / library size)
#' # e.g. row 1 has library size 15 -> detection limit 1/15 -> the zero
#' # becomes 0.5 * (1/15) ≈ 0.033
#'
#' # non-interactive use: same result as entering 0.5 above, but repeatable
#' # and safe to run inside a knitted/rendered notebook
#' pseudocount(x, threshold_pct = 0.5)
#' }
#' @export


# -------- body of the function --------

pseudocount <- function(x, threshold_pct = NULL) {
  # Convert the input to a tibble.
  x <- as_tibble(x)
  # Compute the total counts for each row.
  lib_size <- rowSums(x)
  # Stop early if any sample has a library size of 0, since the detection
  # limit (1 / lib_size) would be undefined (Inf) for that row.
  if (any(lib_size == 0)) {
    stop("At least one row has a library size of 0; remove empty samples before computing pseudocounts.")
  }
  # If no threshold was passed in, ask the user for one interactively
  # (keep asking until a valid number between 0 and 1 is given). If a
  # threshold was passed in, skip the prompt entirely and validate it
  # instead, so the function can be used non-interactively (e.g. in a
  # rendered notebook, or for sensitivity analyses over several thresholds).
  if (is.null(threshold_pct)) {
    repeat {
      threshold_input <- readline(
        paste0(
          "Percentage threshold for the pseudocounts: each zero is replaced ",
          "with this fraction of the sample's detection limit (1/library size).\n",
          "Enter a number between 0 and 1 (e.g. 0.5 = replace zeroes with half ",
          "the smallest detectable proportion for that sample): "
        )
      )
      threshold_pct <- suppressWarnings(as.numeric(threshold_input))
      # exit only if the input is valid.
      if (!is.na(threshold_pct) && threshold_pct >= 0 && threshold_pct <= 1) {
        break
      }
      # Warn the user if the input is invalid.
      cat("Please provide a number between 0 and 1.\n")
    }
  } else {
    # A threshold was supplied directly: validate it with the same rule
    # used in the interactive branch, instead of silently accepting
    # anything (e.g. a typo like 5 instead of 0.5).
    if (!is.numeric(threshold_pct) || is.na(threshold_pct) ||
        threshold_pct < 0 || threshold_pct > 1) {
      stop("threshold_pct must be a single number between 0 and 1 (got ",
           threshold_pct, ")")
    }
  }
  # Compute the detection limit for each row.
  detection_limit <- 1 / lib_size
  # Compute the row-specific pseudocount values.
  pseudo <- threshold_pct * detection_limit
  # Convert raw counts to row-wise proportions.
  y_prop <- sweep(as.matrix(x), 1, lib_size, "/")
  # Replace zeroes with the corresponding row-specific pseudocount.
  # This is done on a plain matrix, not a tibble: tibbles do not support
  # indexing/assignment via a logical matrix, only base matrices do.
  zero_mask <- y_prop == 0
  # rep() relies on R's column-major storage: repeating the row-wise
  # pseudo vector `ncol` times produces the same flattened order as a
  # matrix where every column equals `pseudo`, so indexing by zero_mask
  # (also flattened column-major) pairs each zero with its own row's value.
  y_prop[zero_mask] <- rep(pseudo, times = ncol(y_prop))[zero_mask]
  # Return the transformed data as a tibble.
  as_tibble(y_prop)
}