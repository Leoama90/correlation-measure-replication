# pseudocount.R
#
# Purpose:
#   Convert count data to row-wise proportions and replace zeroes with
#   row-specific pseudocounts computed from a user-defined threshold.
#
# Input:
#   - x: a numeric matrix or data frame of non-negative counts.
#   - threshold_pct: a numeric value between 0 and 1 (requested interactively).
#
# Output:
#   - A tibble with counts converted to proportions and zeroes replaced by
#     row-specific pseudocounts.
#
# Description:
#   The function computes library size by row, asks the user for a threshold,
#   calculates row-specific pseudocounts, and substitutes zeroes in the
#   proportional table with those values.

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)
# tidyverse: data science toolkit (dplyr, ggplot2, tidyr, etc.)
# https://tidyverse.org
library(tidyverse)

#' Replace zeroes with row-wise pseudocounts
#'
#' This function computes row-wise library sizes, asks the user for a threshold
#' between 0 and 1, converts counts to proportions, and replaces zeroes with
#' row-specific pseudocount values.
#'
#' @param x A numeric matrix or data frame of non-negative counts. Rows
#'   (samples) must have a library size greater than 0.
#'
#' @return A tibble with counts converted to proportions and zeroes replaced by
#'   row-specific pseudocounts.
#'
#' @examples
#' \dontrun{
#' x <- matrix(c(10, 0, 5, 2, 3, 0), nrow = 2, byrow = TRUE)
#' pseudocount(x)
#' }
#'
#' @export

# -------- body of the function --------

pseudocount <- function(x) {
  # Convert the input to a tibble.
  x <- as_tibble(x)

  # Compute the total counts for each row.
  lib_size <- rowSums(x)

  # Stop early if any sample has a library size of 0, since the detection
  # limit (1 / lib_size) would be undefined (Inf) for that row.
  if (any(lib_size == 0)) {
    stop("At least one row has a library size of 0; remove empty samples before computing pseudocounts.")
  }

  # Ask the user for a threshold between 0 and 1.
  repeat {
    threshold_input <- readline("Give me the percentage threshold for the pseudocounts (number between 0 and 1): ")
    threshold_pct <- suppressWarnings(as.numeric(threshold_input))

    # Stop only if the input is valid.
    if (!is.na(threshold_pct) && threshold_pct >= 0 && threshold_pct <= 1) {
      break
    }

    # Warn the user if the input is invalid.
    cat("Please provide a number between 0 and 1.\n")
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
  y_prop[zero_mask] <- rep(pseudo, times = ncol(y_prop))[zero_mask]

  # Return the transformed data as a tibble.
  as_tibble(y_prop)
}
